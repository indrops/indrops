import os, subprocess
import itertools
import operator
from collections import defaultdict, OrderedDict
import errno

# cPickle is a faster version of pickle that isn't installed in python3
# inserted try statement just in case
try:
   import cPickle as pickle
except:
   import pickle

from io import BytesIO

import numpy as np
import re
import shutil
import gzip

# product: product(A, B) returns the same as ((x,y) for x in A for y in B).
# combination: Return r length subsequences of elements from the input iterable.
from itertools import product, combinations
import time

import yaml

import tempfile
import string
from contextlib import contextmanager

# -----------------------
#
# Helper functions
#
# -----------------------

def string_hamming_distance(str1, str2):
    """
    Fast hamming distance over 2 strings known to be of same length.
    In information theory, the Hamming distance between two strings of equal 
    length is the number of positions at which the corresponding symbols 
    are different.

    eg "karolin" and "kathrin" is 3.
    """
    return sum(itertools.imap(operator.ne, str1, str2))

___tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def rev_comp(seq):
    return ''.join(___tbl[s] for s in seq[::-1])


def to_fastq(name, seq, qual):
    """
    Return string that can be written to fastQ file
    """
    return '@'+name+'\n'+seq+'\n+\n'+qual+'\n'

def to_fastq_lines(bc, umi, seq, qual, read_name=''):
    """
    Return string that can be written to fastQ file
    """
    reformated_name = read_name.replace(':', '_')
    name = '%s:%s:%s' % (bc, umi, reformated_name)
    return to_fastq(name, seq, qual)

def from_fastq(handle):
    while True:
        name = next(handle).rstrip()[1:] #Read name
        seq = next(handle).rstrip() #Read seq
        next(handle) #+ line
        qual = next(handle).rstrip() #Read qual
        if not name or not seq or not qual:
            break
        yield name, seq, qual

def seq_neighborhood(seq, n_subs=1):
    """
    Given a sequence, yield all sequences within n_subs substitutions of 
    that sequence by looping through each combination of base pairs within
    each combination of positions.
    """
    for positions in combinations(range(len(seq)), n_subs):
    # yields all unique combinations of indices for n_subs mutations
        for subs in product(*("ATGCN",)*n_subs):
        # yields all combinations of possible nucleotides for strings of length
        # n_subs
            seq_copy = list(seq)
            for p, s in zip(positions, subs):
                seq_copy[p] = s
            yield ''.join(seq_copy)

def build_barcode_neighborhoods(barcode_file, expect_reverse_complement=True):
    """
    Given a set of barcodes, produce sequences which can unambiguously be
    mapped to these barcodes, within 2 substitutions. If a sequence maps to 
    multiple barcodes, get rid of it. However, if a sequences maps to a bc1 with 
    1change and another with 2changes, keep the 1change mapping.
    """

    # contains all mutants that map uniquely to a barcode
    clean_mapping = dict()

    # contain single or double mutants 
    mapping1 = defaultdict(set)
    mapping2 = defaultdict(set)
    
    #Build the full neighborhood and iterate through barcodes
    with open(barcode_file, 'rU') as f:
        # iterate through each barcode (rstrip cleans string of whitespace)
        for line in f:
            barcode = line.rstrip()
            if expect_reverse_complement:
                barcode = rev_comp(line.rstrip())

            # each barcode obviously maps to itself uniquely
            clean_mapping[barcode] = barcode

            # for each possible mutated form of a given barcode, either add
            # the origin barcode into the set corresponding to that mutant or 
            # create a new entry for a mutant not already in mapping1
            # eg: barcodes CATG and CCTG would be in the set for mutant CTTG
            # but only barcode CATG could generate mutant CANG
            for n in seq_neighborhood(barcode, 1):
                mapping1[n].add(barcode)
            
            # same as above but with double mutants
            for n in seq_neighborhood(barcode, 2):
                mapping2[n].add(barcode)   
    
    # take all single-mutants and find those that could only have come from one
    # specific barcode
    for k, v in mapping1.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]
    
    for k, v in mapping2.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]
    del mapping1
    del mapping2
    return clean_mapping

def check_dir(path):
    """
    Checks if directory already exists or not and creates it if it doesn't
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def print_to_stderr(msg, newline=True):
    """
    Wrapper to eventually write to stderr
    """
    sys.stderr.write(str(msg))
    if newline:
        sys.stderr.write('\n')

def worker_filter(iterable, worker_index, total_workers):
    return (p for i,p in enumerate(iterable) if (i-worker_index)%total_workers==0)

class FIFO():
    """
    A context manager for a named pipe.
    """
    def __init__(self, filename="", suffix="", prefix="tmp_fifo_dir", dir=None):
        if filename:
            self.filename = filename
        else:
            self.tmpdir = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
            self.filename = os.path.join(self.tmpdir, 'fifo')

    def __enter__(self):
        if os.path.exists(self.filename):
            os.unlink(self.filename)
        os.mkfifo(self.filename)
        return self

    def __exit__(self, type, value, traceback):
        os.remove(self.filename)
        if hasattr(self, 'tmpdir'):
            shutil.rmtree(self.tmpdir)

# -----------------------
#
# Core objects
#
# -----------------------

class IndropsProject():

    def __init__(self, project_yaml_file_handle, read_only=False):

        self.yaml = yaml.load(project_yaml_file_handle)

        self.name = self.yaml['project_name']
        self.project_dir = self.yaml['project_dir']

        self.libraries = OrderedDict()
        self.runs = OrderedDict()

        self.read_only = read_only

        for run in self.yaml['sequencing_runs']:
            """
            After filtering, each sequencing run generates between 1 ... X files with filtered reads.
               X = (N x M)
             - N: The run is often split into several files (a typical NextSeq run is split into L001,
                  L002, L003, L004 which match different lanes, but this can also be done artificially.
             - M: The same run might contain several libraries. The demultiplexing can be handled by the script (or externally).
                  If demultiplexing is done externally, there will be a different .fastq file for each library.
            """
            version = run['version']

            filtered_filename = '{library_name}_{run_name}'
            if run['version'] == 'v3':
                filtered_filename += '_{library_index}'
            # Prepare to iterate over run split into several files
            if 'split_affixes' in run:
                filtered_filename += '_{split_affix}'
                split_affixes = run['split_affixes']
            else:
                split_affixes = ['']

            filtered_filename += '.fastq'

            # Prepare to iterate over libraries
            if 'libraries' in run:
                run_libraries = run['libraries']
            elif 'library_name' in run:
                run_libraries = [{'library_name' : run['library_name'], 'library_prefix':''}]
            else:
                raise Exception('No library name or libraries specified.')

            if run['version']=='v1' or run['version']=='v2':
                for affix in split_affixes:
                    for lib in run_libraries:
                        lib_name = lib['library_name']
                        if lib_name not in self.libraries:
                            self.libraries[lib_name] = IndropsLibrary(name=lib_name, project=self, version=run['version'])
                        else:
                            assert self.libraries[lib_name].version == run['version']

                        if version == 'v1':
                            metaread_filename = os.path.join(run['dir'],run['fastq_path'].format(split_affix=affix, read='R1', library_prefix=lib['library_prefix']))
                            bioread_filename = os.path.join(run['dir'],run['fastq_path'].format(split_affix=affix, read='R2', library_prefix=lib['library_prefix']))
                        elif version == 'v2':
                            metaread_filename  = os.path.join(run['dir'],run['fastq_path'].format(split_affix=affix, read='R2', library_prefix=lib['library_prefix']))
                            bioread_filename = os.path.join(run['dir'],run['fastq_path'].format(split_affix=affix, read='R1', library_prefix=lib['library_prefix']))

                        filtered_part_filename = filtered_filename.format(run_name=run['name'], split_affix=affix, library_name=lib_name)
                        filtered_part_path = os.path.join(self.project_dir, lib_name, 'filtered_parts', filtered_part_filename)
                        part = V1V2Filtering(filtered_fastq_filename=filtered_part_path,
                            project=self, 
                            bioread_filename=bioread_filename,
                            metaread_filename=metaread_filename,
                            run_name=run['name'],
                            library_name=lib_name,
                            part_name=affix)

                        if run['name'] not in self.runs:
                            self.runs[run['name']] = []
                        self.runs[run['name']].append(part)
                        self.libraries[lib_name].parts.append(part)

            elif run['version'] == 'v3' or run['version'] == 'v3-miseq':
                for affix in split_affixes:
                    filtered_part_filename = filtered_filename.format(run_name=run['name'], split_affix=affix,
                        library_name='{library_name}', library_index='{library_index}')
                    part_filename = os.path.join(self.project_dir, '{library_name}', 'filtered_parts', filtered_part_filename)

                    input_filename = os.path.join(run['dir'], run['fastq_path'].format(split_affix=affix, read='{read}'))
                    part = V3Demultiplexer(run['libraries'], project=self, part_filename=part_filename, input_filename=input_filename, run_name=run['name'], part_name=affix,
                        run_version_details=run['version'])

                    if run['name'] not in self.runs:
                        self.runs[run['name']] = []
                    self.runs[run['name']].append(part)

                    for lib in run_libraries:
                        lib_name = lib['library_name']
                        lib_index = lib['library_index']
                        if lib_name not in self.libraries:
                            self.libraries[lib_name] = IndropsLibrary(name=lib_name, project=self, version=run['version'])
                        self.libraries[lib_name].parts.append(part.libraries[lib_index])


    @property
    def paths(self):
        if not hasattr(self, '_paths'):
            script_dir = os.path.dirname(os.path.realpath(__file__))
            #Read defaults
            with open(os.path.join(script_dir, 'default_parameters.yaml'), 'r') as f:
                paths = yaml.load(f)['paths']
            # Update with user provided values
            paths.update(self.yaml['paths'])

            paths['python'] = os.path.join(paths['python_dir'], 'python')
            paths['java'] = os.path.join(paths['java_dir'], 'java')
            paths['bowtie'] = os.path.join(paths['bowtie_dir'], 'bowtie')
            paths['samtools'] = os.path.join(paths['samtools_dir'], 'samtools')
            paths['trimmomatic_jar'] = os.path.join(script_dir, 'bins', 'trimmomatic-0.33.jar')
            paths['rsem_tbam2gbam'] = os.path.join(paths['rsem_dir'], 'rsem-tbam2gbam')
            paths['rsem_prepare_reference'] = os.path.join(paths['rsem_dir'], 'rsem-prepare-reference')

            self._paths = type('Paths_anonymous_object',(object,),paths)()
            self._paths.trim_polyA_and_filter_low_complexity_reads_py = os.path.join(script_dir, 'trim_polyA_and_filter_low_complexity_reads.py')
            self._paths.quantify_umifm_from_alignments_py = os.path.join(script_dir, 'quantify_umifm_from_alignments.py')
            self._paths.count_barcode_distribution_py = os.path.join(script_dir, 'count_barcode_distribution.py')
            self._paths.gel_barcode1_list = os.path.join(script_dir, 'ref/barcode_lists/gel_barcode1_list.txt')
            self._paths.gel_barcode2_list = os.path.join(script_dir, 'ref/barcode_lists/gel_barcode2_list.txt')
        return self._paths

    @property
    def parameters(self):
        if not hasattr(self, '_parameters'):
            #Read defaults
            with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'default_parameters.yaml'), 'r') as f:
                self._parameters = yaml.load(f)['parameters']
            # Update with user provided values
            if 'parameters' in self.yaml:
                for k, d in self.yaml['parameters'].items():
                    self._parameters[k].update(d)

        return self._parameters

    @property
    def gel_barcode1_revcomp_list_neighborhood(self):
        if not hasattr(self, '_gel_barcode1_list_neighborhood'):
            self._gel_barcode1_revcomp_list_neighborhood = build_barcode_neighborhoods(self.paths.gel_barcode1_list, True)
        return self._gel_barcode1_revcomp_list_neighborhood
    
    @property
    def gel_barcode2_revcomp_list_neighborhood(self):
        if not hasattr(self, '_gel_barcode2_revcomp_list_neighborhood'):
            self._gel_barcode2_revcomp_list_neighborhood = build_barcode_neighborhoods(self.paths.gel_barcode2_list, True)
        return self._gel_barcode2_revcomp_list_neighborhood

    @property
    def gel_barcode2_list_neighborhood(self):
        if not hasattr(self, '_gel_barcode2_list_neighborhood'):
            self._gel_barcode2_list_neighborhood = build_barcode_neighborhoods(self.paths.gel_barcode2_list, False)
        return self._gel_barcode2_list_neighborhood

    @property
    def stable_barcode_names(self):
        if not hasattr(self, '_stable_barcode_names'):
            with open(self.paths.gel_barcode1_list) as f:
                rev_bc1s = [rev_comp(line.rstrip()) for line in f]
            with open(self.paths.gel_barcode2_list) as f:
                bc2s = [line.rstrip() for line in f]
                rev_bc2s = [rev_comp(bc2) for bc2 in bc2s]

            # V1, V2 names:
            v1v2_names = {}
            barcode_iter = product(rev_bc1s, rev_bc2s)
            name_iter = product(string.ascii_uppercase, repeat=4)
            for barcode, name in zip(barcode_iter, name_iter):
                v1v2_names['-'.join(barcode)] = 'bc' + ''.join(name)

            # V3 names:
            v3_names = {}
            barcode_iter = product(bc2s, rev_bc2s)
            name_iter = product(string.ascii_uppercase, repeat=4)
            for barcode, name in zip(barcode_iter, name_iter):
                v3_names['-'.join(barcode)] = 'bc' + ''.join(name)


            self._stable_barcode_names = {
                'v1' : v1v2_names,
                'v2' : v1v2_names,
                'v3': v3_names,
                'v3-miseq':v3_names,
            }
        return self._stable_barcode_names

    def project_check_dir(self, path):
        if not self.read_only:
            check_dir(path)

    def filter_gtf(self, gzipped_transcriptome_gtf, gtf_with_genenames_in_transcript_id):

        # A small number of gene are flagged as having two different biotypes.
        gene_biotype_dict = defaultdict(set)

        # Read through GTF file once to get all gene names
        for line in subprocess.Popen(["gzip", "--stdout", "-d", gzipped_transcriptome_gtf], stdout=subprocess.PIPE).stdout:
            # Skip non-gene feature lines.
            if '\tgene\t' not in line:
                continue
                
            gene_biotype_match = re.search(r'gene_biotype \"(.*?)\";', line)
            gene_name_match = re.search(r'gene_name \"(.*?)\";', line)
            if gene_name_match and gene_biotype_match:
                gene_name = gene_name_match.group(1)
                gene_biotype = gene_biotype_match.group(1)
                
                # Record biotype.
                gene_biotype_dict[gene_name].add(gene_biotype)

        # Detect read-through genes by name. Name must be a fusion of two other gene names 'G1-G2'.
        readthrough_genes = set()
        for gene in gene_biotype_dict.keys():
            if '-' in gene and len(gene.split('-')) == 2:
                g1, g2 = gene.split('-')
                if g1 in gene_biotype_dict and g2 in gene_biotype_dict:
                    readthrough_genes.add(gene)


        # Detect pseudogenes: genes where all associated biotypes have 'pseudogene' in name
        pseudogenes = set()
        for gene, biotypes in gene_biotype_dict.items():
            if all('pseudogene' in b for b in biotypes):
                pseudogenes.add(gene)

        all_genes = set(gene_biotype_dict.keys())
        valid_genes = all_genes.difference(pseudogenes).difference(readthrough_genes)

        transcripts_counter = defaultdict(int)


        # Go through GTF file again, annotating each transcript_id with the gene and outputting to a new GTF file.
        output_gtf = open(gtf_with_genenames_in_transcript_id, 'w')
        for line in subprocess.Popen(["gzip", "--stdout", "-d", gzipped_transcriptome_gtf], stdout=subprocess.PIPE).stdout:
            # Skip non-transcript feature lines.
            if 'transcript_id' not in line:
                continue
                
            gene_name_match = re.search(r'gene_name \"(.*?)\";', line)
            if gene_name_match:
                gene_name = gene_name_match.group(1)
                if gene_name in valid_genes:
                    
                    # An unusual edgecase in the GTF for Danio Rerio rel89
                    if ' ' in gene_name:
                        gene_name = gene_name.replace(' ', '_')

                    out_line = re.sub(r'(?<=transcript_id ")(.*?)(?=";)', r'\1|'+gene_name, line)
                    output_gtf.write(out_line)
                    if '\ttranscript\t' in line:
                        transcripts_counter['valid'] += 1
                elif gene_name in pseudogenes and '\ttranscript\t' in line:
                    transcripts_counter['pseudogenes'] += 1
                elif gene_name in readthrough_genes and '\ttranscript\t' in line:
                    transcripts_counter['readthrough_genes'] += 1
        output_gtf.close()

        print_to_stderr('Filtered GTF contains %d transcripts (%d genes)' % (transcripts_counter['valid'], len(valid_genes)))
        print_to_stderr('   - ignored %d transcripts from %d pseudogenes)' % (transcripts_counter['pseudogenes'], len(pseudogenes)))
        print_to_stderr('   - ignored %d read-through transcripts (%d genes)' % (transcripts_counter['readthrough_genes'], len(readthrough_genes)))

    def build_transcriptome(self, gzipped_genome_softmasked_fasta_filename, gzipped_transcriptome_gtf,
            mode='strict'):
        import pyfasta
        
        index_dir = os.path.dirname(self.paths.bowtie_index)
        self.project_check_dir(index_dir)

        genome_filename = os.path.join(index_dir, '.'.join(gzipped_genome_softmasked_fasta_filename.split('.')[:-1]))

        gtf_filename = os.path.join(index_dir, gzipped_transcriptome_gtf.split('/')[-1])
        gtf_prefix = '.'.join(gtf_filename.split('.')[:-2])
        # gtf_with_genenames_in_transcript_id = gtf_prefix + '.annotated.gtf'
        gtf_with_genenames_in_transcript_id = self.paths.bowtie_index + '.gtf'

        print_to_stderr('Filtering GTF')
        self.filter_gtf(gzipped_transcriptome_gtf, gtf_with_genenames_in_transcript_id)
        # accepted_gene_biotypes_for_NA_transcripts = set(["protein_coding","IG_V_gene","IG_J_gene","TR_J_gene","TR_D_gene","TR_V_gene","IG_C_gene","IG_D_gene","TR_C_gene"])
        # tsl1_or_tsl2_strings = ['transcript_support_level "1"', 'transcript_support_level "1 ', 'transcript_support_level "2"', 'transcript_support_level "2 ']
        # tsl_NA =  'transcript_support_level "NA'

        # def filter_ensembl_transcript(transcript_line):

        #     line_valid_for_output = False
        #     if mode == 'strict':
        #         for string in tsl1_or_tsl2_strings:
        #             if string in line:
        #                 line_valid_for_output = True
        #                 break
        #         if tsl_NA in line:
        #             gene_biotype = re.search(r'gene_biotype \"(.*?)\";', line)
        #             if gene_biotype and gene_biotype.group(1) in accepted_gene_biotypes_for_NA_transcripts:
        #                 line_valid_for_output = True
        #         return line_valid_for_output

        #     elif mode == 'all_ensembl':
        #         line_valid_for_output = True
        #         return line_valid_for_output



        # print_to_stderr('Filtering GTF')
        # output_gtf = open(gtf_with_genenames_in_transcript_id, 'w')
        # for line in subprocess.Popen(["gzip", "--stdout", "-d", gzipped_transcriptome_gtf], stdout=subprocess.PIPE).stdout:
        #     if 'transcript_id' not in line:
        #         continue

        #     if filter_ensembl_transcript(line):
        #         gene_name = re.search(r'gene_name \"(.*?)\";', line)
        #         if gene_name:
        #             gene_name = gene_name.group(1)
        #             out_line = re.sub(r'(?<=transcript_id ")(.*?)(?=";)', r'\1|'+gene_name, line)
        #             output_gtf.write(out_line)
        # output_gtf.close()

        print_to_stderr('Gunzipping Genome')
        p_gzip = subprocess.Popen(["gzip", "-dfc", gzipped_genome_softmasked_fasta_filename], stdout=open(genome_filename, 'wb'))
        if p_gzip.wait() != 0:
            raise Exception(" Error in rsem-prepare reference ")

        p_rsem = subprocess.Popen([self.paths.rsem_prepare_reference, '--bowtie', '--bowtie-path', self.paths.bowtie_dir,
                            '--gtf', gtf_with_genenames_in_transcript_id, 
                            '--polyA', '--polyA-length', '5', genome_filename, self.paths.bowtie_index])

        if p_rsem.wait() != 0:
            raise Exception(" Error in rsem-prepare reference ")

        print_to_stderr('Finding soft masked regions in transcriptome')
        
        transcripts_fasta = pyfasta.Fasta(self.paths.bowtie_index + '.transcripts.fa')
        soft_mask = {}
        for tx, seq in transcripts_fasta.items():
            seq = str(seq)
            soft_mask[tx] = set((m.start(), m.end()) for m in re.finditer(r'[atcgn]+', seq))
        with open(self.paths.bowtie_index + '.soft_masked_regions.pickle', 'w') as out:
            pickle.dump(soft_mask, out)

class IndropsLibrary():

    def __init__(self, name='', project=None, version=''):
        self.project = project
        self.name = name
        self.parts = []
        self.version = version

        self.paths = {}
        for lib_dir in ['filtered_parts', 'quant_dir']:
            dir_path = os.path.join(self.project.project_dir, self.name, lib_dir)
            self.project.project_check_dir(dir_path)
            self.paths[lib_dir] = dir_path
        self.paths = type('Paths_anonymous_object',(object,),self.paths)()

        self.paths.abundant_barcodes_names_filename = os.path.join(self.project.project_dir, self.name, 'abundant_barcodes.pickle')
        self.paths.filtering_statistics_filename = os.path.join(self.project.project_dir, self.name, self.name+'.filtering_stats.csv')
        self.paths.barcode_abundance_histogram_filename = os.path.join(self.project.project_dir, self.name, self.name+'.barcode_abundance.png')
        self.paths.barcode_abundance_by_barcode_filename = os.path.join(self.project.project_dir, self.name, self.name+'.barcode_abundance_by_barcode.png')
        self.paths.missing_quants_filename = os.path.join(self.project.project_dir, self.name, self.name+'.missing_barcodes.pickle')

    @property
    def barcode_counts(self):
        if not hasattr(self, '_barcode_counts'):
            self._barcode_counts = defaultdict(int)
            for part in self.parts:
                for k, v in part.part_barcode_counts.items():
                    self._barcode_counts[k] += v

        return self._barcode_counts

    @property
    def abundant_barcodes(self):
        if not hasattr(self, '_abundant_barcodes'):
            with open(self.paths.abundant_barcodes_names_filename) as f:
                self._abundant_barcodes = pickle.load(f)
        return self._abundant_barcodes

    def sorted_barcode_names(self, min_reads=0, max_reads=10**10):
        return [name for bc,(name,abun) in sorted(self.abundant_barcodes.items(), key=lambda i:-i[1][1]) if (abun>min_reads) & (abun<max_reads)]

    def identify_abundant_barcodes(self, make_histogram=True, absolute_min_reads=250):
        """
        Identify which barcodes are above the absolute minimal abundance, 
        and make a histogram summarizing the barcode distribution
        """
        keep_barcodes = []
        for k, v in self.barcode_counts.items():
            if v > absolute_min_reads:
                keep_barcodes.append(k)

        abundant_barcodes = {}
        print_to_stderr(" %d barcodes above absolute minimum threshold" % len(keep_barcodes))
        for bc in keep_barcodes:
            abundant_barcodes[bc] = (self.project.stable_barcode_names[self.version][bc], self.barcode_counts[bc])

        self._abundant_barcodes = abundant_barcodes
        with open(self.paths.abundant_barcodes_names_filename, 'w') as f:
            pickle.dump(abundant_barcodes, f)

        # Create table about the filtering process
        with open(self.paths.filtering_statistics_filename, 'w') as filtering_stats:

            header = ['Run', 'Part', 'Input Reads', 'Valid Structure', 'Surviving Trimmomatic', 'Surviving polyA trim and complexity filter']

            if self.version == 'v1' or self.version == 'v2':
                structure_parts = ['W1_in_R2', 'empty_read',  'No_W1', 'No_polyT', 'BC1', 'BC2', 'Umi_error']
                header += ['W1 in R2', 'empty read',  'No W1 in R1', 'No polyT', 'BC1', 'BC2', 'UMI_contains_N']
            elif self.version == 'v3' or self.version == 'v3-miseq':
                structure_parts = ['Invalid_BC1', 'Invalid_BC2', 'UMI_contains_N']
                header += ['Invalid BC1', 'Invalid BC2', 'UMI_contains_N']

            trimmomatic_parts = ['dropped']
            header += ['Dropped by Trimmomatic']

            complexity_filter_parts = ['rejected_because_too_short', 'rejected_because_complexity_too_low']
            header += ['Too short after polyA trim', 'Read complexity too low']

            filtering_stats.write(','.join(header)+'\n')

            for part in self.parts:
                with open(part.filtering_metrics_filename) as f:
                    part_stats = yaml.load(f)
                    line = [part.run_name, part.part_name, part_stats['read_structure']['Total'], part_stats['read_structure']['Valid'], part_stats['trimmomatic']['output'], part_stats['complexity_filter']['output']]
                    line += [part_stats['read_structure'][k] if k in part_stats['read_structure'] else 0 for k in structure_parts]
                    line += [part_stats['trimmomatic'][k] if k in part_stats['trimmomatic'] else 0 for k in trimmomatic_parts]
                    line += [part_stats['complexity_filter'][k] if k in part_stats['complexity_filter'] else 0 for k in complexity_filter_parts]
                    line = [str(l) for l in line]
                    filtering_stats.write(','.join(line)+'\n')

        print_to_stderr("Created Library filtering summary:")
        print_to_stderr("  " + self.paths.filtering_statistics_filename)
 
        # Make the histogram figure
        if not make_histogram:
            return

        count_freq = defaultdict(int)
        for bc, count in self.barcode_counts.items():
            count_freq[count] += 1

        x = np.array(count_freq.keys())
        y = np.array(count_freq.values())
        w = x*y

        # need to use non-intenactive Agg backend
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(x, bins=np.logspace(0, 6, 50), weights=w, color='green')
        ax.set_xscale('log')
        ax.set_xlabel('Reads per barcode')
        ax.set_ylabel('#reads coming from bin')
        fig.savefig(self.paths.barcode_abundance_histogram_filename)

        print_to_stderr("Created Barcode Abundance Histogram at:")
        print_to_stderr("  " + self.paths.barcode_abundance_histogram_filename)


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(list(self.barcode_counts.values()), bins=np.logspace(2, 6, 50), color='green')
        ax.set_xlim((1, 10**6))
        ax.set_xscale('log')
        ax.set_xlabel('Reads per barcode')
        ax.set_ylabel('# of barcodes')
        fig.savefig(self.paths.barcode_abundance_by_barcode_filename)
        print_to_stderr("Created Barcode Abundance Histogram by barcodes at:")
        print_to_stderr("  " + self.paths.barcode_abundance_by_barcode_filename)

    def sort_reads_by_barcode(self, index=0):
        self.parts[index].sort_reads_by_barcode(self.abundant_barcodes)

    def get_reads_for_barcode(self, barcode, run_filter=[]):
        for part in self.parts:
            if (not run_filter) or (part.run_name in run_filter):
                for line in part.get_reads_for_barcode(barcode):
                    yield line

    def output_barcode_fastq(self, analysis_prefix='', min_reads=750, max_reads=10**10, total_workers=1, worker_index=0, run_filter=[]):
        if analysis_prefix:
            analysis_prefix = analysis_prefix + '.'

        output_dir_path = os.path.join(self.project.project_dir, self.name, 'barcode_fastq')
        self.project.project_check_dir(output_dir_path)

        sorted_barcode_names = self.sorted_barcode_names(min_reads=min_reads, max_reads=max_reads)

        # Identify which barcodes belong to this worker
        barcodes_for_this_worker = []
        i = worker_index
        while i < len(sorted_barcode_names):
            barcodes_for_this_worker.append(sorted_barcode_names[i])
            i += total_workers

        print_to_stderr("""[%s] This worker assigned %d out of %d total barcodes.""" % (self.name, len(barcodes_for_this_worker), len(sorted_barcode_names)))        

        for barcode in barcodes_for_this_worker:
            barcode_fastq_filename = analysis_prefix+'%s.%s.fastq' % (self.name, barcode)
            print_to_stderr("  "+barcode_fastq_filename)
            with open(os.path.join(output_dir_path, barcode_fastq_filename), 'w') as f:
                for line in self.get_reads_for_barcode(barcode, run_filter):
                    f.write(line)

    def quantify_expression(self, analysis_prefix='', max_reads=10**10, min_reads=750, min_counts=0, total_workers=1, worker_index=0, no_bam=False, run_filter=[]):
        if analysis_prefix:
            analysis_prefix = analysis_prefix + '.'

        sorted_barcode_names = self.sorted_barcode_names(min_reads=min_reads, max_reads=max_reads)
        #print_to_stderr("   min_reads: %d sorted_barcode_names counts: %d" % (min_reads, len(sorted_barcode_names)))

        # Identify which barcodes belong to this worker
        barcodes_for_this_worker = []
        i = worker_index
        while i < len(sorted_barcode_names):
            barcodes_for_this_worker.append(sorted_barcode_names[i])
            i += total_workers

        counts_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.counts.tsv' % (analysis_prefix, worker_index, total_workers))
        ambig_counts_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.ambig.counts.tsv' % (analysis_prefix, worker_index, total_workers))
        ambig_partners_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.ambig.partners' % (analysis_prefix, worker_index, total_workers))
        metrics_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.metrics.tsv' % (analysis_prefix, worker_index, total_workers))
        ignored_for_output_filename = counts_output_filename+'.ignored'

        merged_bam_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.bam'% (analysis_prefix, worker_index, total_workers))
        merged_bam_index_filename = merged_bam_filename + '.bai'

        get_barcode_genomic_bam_filename = lambda bc: os.path.join(self.paths.quant_dir, '%s%s.genomic.sorted.bam' % (analysis_prefix, bc))

        # If we wanted BAM output, and the merge BAM and merged BAM index are present, then we are done
        if (not no_bam) and (os.path.isfile(merged_bam_filename) and os.path.isfile(merged_bam_index_filename)):
            print_to_stderr('Indexed, merged BAM file detected for this worker. Done.')
            return 

        # Otherwise, we have to check what we need to quantify

        
        """
        Function to determine which barcodes this quantification worker might have already quantified.
        This tries to handle interruption during any step of the process.

        The worker is assigned some list of barcodes L. For every barcode:
            - It could have been quantified
                - but have less than min_counts ---> so it got written to `ignored` file.
                - and quantification succeeded, meaning
                    1. there is a line (ending in \n) in the `metrics` file. 
                    2. there is a line (ending in \n) in the `quantification` file.
                    3. there (could) be a line (ending in \n) in the `ambiguous quantification` file.
                    4. there (could) be a line (ending in \n) in the `ambiguous quantification partners` file.
                        [If any line doesn't end in \n, then likely the output of that line was interrupted!]
                    5. (If BAM output is desired) There should be a sorted genomic BAM
                    6. (If BAM output is desired) There should be a sorted genomic BAM index
        """
        succesfully_previously_quantified = set()
        previously_ignored = set()
        header_written = False

        if os.path.isfile(counts_output_filename) and os.path.isfile(metrics_output_filename):
            # Load in list of ignored barcodes
            if os.path.isfile(ignored_for_output_filename):
                with open(ignored_for_output_filename, 'r') as f:
                    previously_ignored = set([line.rstrip().split('\t')[0] for line in f])

            # Load the metrics data into memory
            # (It should be fairly small, this is fast and safe)
            existing_metrics_data = {}
            with open(metrics_output_filename, 'r') as f:
                existing_metrics_data = dict((line.partition('\t')[0], line) for line in f if line[-1]=='\n')


            # Quantification data could be large, read it line by line and output it back for barcodes that have a matching metrics line.
            with open(counts_output_filename, 'r') as in_counts, \
                     open(counts_output_filename+'.tmp', 'w') as tmp_counts, \
                     open(metrics_output_filename+'.tmp', 'w') as tmp_metrics:

                for line in in_counts:
                    # The first worker is reponsible for written the header.
                    # Make sure we carry that over
                    if (not header_written) and (worker_index==0):
                        tmp_counts.write(line)
                        tmp_metrics.write(existing_metrics_data['Barcode'])
                        header_written = True
                        continue

                    # This line has incomplete output, skip it.
                    # (This can only happen with the last line)
                    if line[-1] != '\n':
                        continue

                    barcode = line.partition('\t')[0]

                    # Skip barcode if we don't have existing metrics data
                    if barcode not in existing_metrics_data:
                        continue

                    # Check if we BAM required BAM files exist
                    barcode_genomic_bam_filename = get_barcode_genomic_bam_filename(barcode)
                    bam_files_required_and_present = no_bam or (os.path.isfile(barcode_genomic_bam_filename) and os.path.isfile(barcode_genomic_bam_filename+'.bai'))
                    if not bam_files_required_and_present:
                        continue

                    # This passed all the required checks, write the line to the temporary output files
                    tmp_counts.write(line)
                    tmp_metrics.write(existing_metrics_data[barcode])
                    succesfully_previously_quantified.add(barcode)

            shutil.move(counts_output_filename+'.tmp', counts_output_filename)
            shutil.move(metrics_output_filename+'.tmp', metrics_output_filename)

            # For any 'already quantified' barcode, make sure we also copy over the ambiguity data
            with open(ambig_counts_output_filename, 'r') as in_f, \
                 open(ambig_counts_output_filename+'.tmp', 'w') as tmp_f:
                 f_first_line = (worker_index == 0)
                 for line in in_f:
                    if f_first_line:
                        tmp_f.write(line)
                        f_first_line = False
                        continue
                    if (line.partition('\t')[0] in succesfully_previously_quantified) and (line[-1]=='\n'):
                        tmp_f.write(line)
            shutil.move(ambig_counts_output_filename+'.tmp', ambig_counts_output_filename)

            with open(ambig_partners_output_filename, 'r') as in_f, \
                 open(ambig_partners_output_filename+'.tmp', 'w') as tmp_f:
                 for line in in_f:
                    if (line.partition('\t')[0] in succesfully_previously_quantified) and (line[-1]=='\n'):
                        tmp_f.write(line)
            shutil.move(ambig_partners_output_filename+'.tmp', ambig_partners_output_filename)

        barcodes_to_quantify = [bc for bc in barcodes_for_this_worker if (bc not in succesfully_previously_quantified and bc not in previously_ignored)]


        print_to_stderr("""[%s] This worker assigned %d out of %d total barcodes.""" % (self.name, len(barcodes_for_this_worker), len(sorted_barcode_names)))
        if len(barcodes_for_this_worker)-len(barcodes_to_quantify) > 0:
            print_to_stderr("""    %d previously quantified, %d previously ignored, %d left for this run.""" % (len(succesfully_previously_quantified), len(previously_ignored), len(barcodes_to_quantify)))
        


        print_to_stderr(('{0:<14.12}'.format('Prefix') if analysis_prefix else '') + '{0:<14.12}{1:<9}'.format("Library", "Barcode"), False)
        print_to_stderr("{0:<8s}{1:<8s}{2:<10s}".format("Reads", "Counts", "Ambigs"))
        for barcode in barcodes_to_quantify:
            self.quantify_expression_for_barcode(barcode,
                counts_output_filename, metrics_output_filename,
                ambig_counts_output_filename, ambig_partners_output_filename,
                no_bam=no_bam, write_header=(not header_written) and (worker_index==0), analysis_prefix=analysis_prefix,
                min_counts = min_counts, run_filter=run_filter)
            header_written = True
        print_to_stderr("Per barcode quantification completed.")

        if no_bam:
            return

        #Gather list of barcodes with output from the metrics file
        genomic_bams = []
        with open(metrics_output_filename, 'r') as f:
            for line in f:
                bc = line.partition('\t')[0]
                if bc == 'Barcode': #This is the line in the header
                    continue
                genomic_bams.append(get_barcode_genomic_bam_filename(bc))

        print_to_stderr("Merging BAM output.")
        try:
            subprocess.check_output([self.project.paths.samtools, 'merge', '-f', merged_bam_filename]+genomic_bams, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, err:
            print_to_stderr("   CMD: %s" % str(err.cmd)[:400])
            print_to_stderr("   stdout/stderr:")
            print_to_stderr(err.output)
            raise Exception(" === Error in samtools merge === ")

        print_to_stderr("Indexing merged BAM output.")
        try:
            subprocess.check_output([self.project.paths.samtools, 'index', merged_bam_filename], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, err:
            print_to_stderr("   CMD: %s" % str(err.cmd)[:400])
            print_to_stderr("   stdout/stderr:")
            print_to_stderr(err.output)
            raise Exception(" === Error in samtools index === ")

        print(genomic_bams)
        for filename in genomic_bams:
            os.remove(filename)
            os.remove(filename + '.bai')

    def quantify_expression_for_barcode(self, barcode, counts_output_filename, metrics_output_filename,
            ambig_counts_output_filename, ambig_partners_output_filename,
            min_counts=0, analysis_prefix='', no_bam=False, write_header=False, run_filter=[]):
        print_to_stderr(('{0:<14.12}'.format(analysis_prefix) if analysis_prefix else '') + '{0:<14.12}{1:<9}'.format(self.name, barcode), False)

        unaligned_reads_output = os.path.join(self.paths.quant_dir, '%s%s.unaligned.fastq' % (analysis_prefix,barcode))
        aligned_bam = os.path.join(self.paths.quant_dir, '%s%s.aligned.bam' % (analysis_prefix,barcode))

        # Bowtie command
        bowtie_cmd = [self.project.paths.bowtie, self.project.paths.bowtie_index, '-q', '-',
            '-p', '1', '-a', '--best', '--strata', '--chunkmbs', '1000', '--norc', '--sam',
            '-shmem', #should sometimes reduce memory usage...?
            '-m', str(self.project.parameters['bowtie_arguments']['m']),
            '-n', str(self.project.parameters['bowtie_arguments']['n']),
            '-l', str(self.project.parameters['bowtie_arguments']['l']),
            '-e', str(self.project.parameters['bowtie_arguments']['e']),
            ]
        if self.project.parameters['output_arguments']['output_unaligned_reads_to_other_fastq']:
            bowtie_cmd += ['--un', unaligned_reads_output]

        # Quantification command
        script_dir = os.path.dirname(os.path.realpath(__file__))
        quant_cmd = [self.project.paths.python, self.project.paths.quantify_umifm_from_alignments_py,
            '-m', str(self.project.parameters['umi_quantification_arguments']['m']),
            '-u', str(self.project.parameters['umi_quantification_arguments']['u']),
            '-d', str(self.project.parameters['umi_quantification_arguments']['d']),
            '--min_non_polyA', str(self.project.parameters['umi_quantification_arguments']['min_non_polyA']),
            '--library', str(self.name),
            '--barcode', str(barcode),
            '--counts', counts_output_filename,
            '--metrics', metrics_output_filename,
            '--ambigs', ambig_counts_output_filename,
            '--ambig-partners', ambig_partners_output_filename,
            '--min-counts', str(min_counts),
        ]
        if not no_bam:
            quant_cmd += ['--bam', aligned_bam]
        if write_header:
            quant_cmd += ['--write-header']

        if self.project.parameters['umi_quantification_arguments']['split-ambigs']:
            quant_cmd.append('--split-ambig')
        if self.project.parameters['output_arguments']['filter_alignments_to_softmasked_regions']:
            quant_cmd += ['--soft-masked-regions', self.project.paths.bowtie_index + '.soft_masked_regions.pickle']

        # Spawn processes

        p1 = subprocess.Popen(bowtie_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(quant_cmd, stdin=p1.stdout, stderr=subprocess.PIPE)
        
                
        for line in self.get_reads_for_barcode(barcode, run_filter=run_filter):
            try:
                p1.stdin.write(line)
            except IOError as e:
                print_to_stderr('\n')
                print_to_stderr(p1.stderr.read())
                raise Exception('\n === Error on piping data to bowtie ===')


        p1.stdin.close()

        if p1.wait() != 0:
            print_to_stderr('\n')
            print_to_stderr(p1.stderr.read())
            raise Exception('\n === Error on bowtie ===')

        if p2.wait() != 0:
            print_to_stderr(p2.stderr.read())
            raise Exception('\n === Error on Quantification Script ===')
        print_to_stderr(p2.stderr.read(), False)

        if no_bam:
            # We are done here
            return False

        if not os.path.isfile(aligned_bam):
            raise Exception("\n === No aligned bam was output for barcode %s ===" % barcode)

        genomic_bam = os.path.join(self.paths.quant_dir, '%s%s.genomic.bam' % (analysis_prefix,barcode))
        sorted_bam = os.path.join(self.paths.quant_dir, '%s%s.genomic.sorted.bam' % (analysis_prefix,barcode))
        try:
            subprocess.check_output([self.project.paths.rsem_tbam2gbam, self.project.paths.bowtie_index, aligned_bam, genomic_bam], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, err:
            print_to_stderr("   CMD: %s" % str(err.cmd)[:100])
            print_to_stderr("   stdout/stderr:")
            print_to_stderr(err.output)
            raise Exception(" === Error in rsem-tbam2gbam === ")

        try:
            subprocess.check_output([self.project.paths.samtools, 'sort', '-o', sorted_bam, genomic_bam], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, err:
            print_to_stderr("   CMD: %s" % str(err.cmd)[:100])
            print_to_stderr("   stdout/stderr:")
            print_to_stderr(err.output)
            raise Exception(" === Error in samtools sort === ")

        try:
            subprocess.check_output([self.project.paths.samtools, 'index', sorted_bam], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, err:
            print_to_stderr("   CMD: %s" % str(err.cmd)[:100])
            print_to_stderr("   stdout/stderr:")
            print_to_stderr(err.output)
            raise Exception(" === Error in samtools index === ")

        os.remove(aligned_bam)
        os.remove(genomic_bam)


        return True

    def aggregate_counts(self, analysis_prefix='', process_ambiguity_data=False):
        if analysis_prefix:
            analysis_prefix = analysis_prefix + '.'
            quant_output_files = [fn[len(analysis_prefix):].split('.')[0] for fn in os.listdir(self.paths.quant_dir) if ('worker' in fn and fn[:len(analysis_prefix)]==analysis_prefix)]
        else:
            quant_output_files = [fn.split('.')[0] for fn in os.listdir(self.paths.quant_dir) if (fn[:6]=='worker')]
        
        worker_names = [w[6:] for w in quant_output_files]
        worker_indices = set(int(w.split('_')[0]) for w in worker_names)

        total_workers = set(int(w.split('_')[1]) for w in worker_names)
        if len(total_workers) > 1:
            raise Exception("""Quantification for library %s, prefix '%s' was run with different numbers of total_workers.""" % (self.name, analysis_prefix))
        total_workers = list(total_workers)[0]

        missing_workers = []
        for i in range(total_workers):
            if i not in worker_indices:
                missing_workers.append(i)
        if missing_workers:
            missing_workers = ','.join([str(i) for i in sorted(missing_workers)])
            raise Exception("""Output from workers %s (total %d) is missing. """ % (missing_workers, total_workers))

        aggregated_counts_filename = os.path.join(self.project.project_dir, self.name, self.name+'.'+analysis_prefix+'counts.tsv')
        aggregated_quant_metrics_filename = os.path.join(self.project.project_dir, self.name, self.name+'.'+analysis_prefix+'quant_metrics.tsv')
        aggregated_ignored_filename = os.path.join(self.project.project_dir, self.name, self.name+'.'+analysis_prefix+'ignored_barcodes.txt')
        aggregated_bam_output = os.path.join(self.project.project_dir, self.name, self.name+'.'+analysis_prefix+'bam')

        aggregated_ambig_counts_filename = os.path.join(self.project.project_dir, self.name, self.name+'.'+analysis_prefix+'ambig_counts.tsv')
        aggregated_ambig_partners_filename = os.path.join(self.project.project_dir, self.name, self.name+'.'+analysis_prefix+'ambig_partners.tsv')

        agg_counts = open(aggregated_counts_filename, mode='w')
        agg_metrics = open(aggregated_quant_metrics_filename, mode='w')
        agg_ignored = open(aggregated_ignored_filename, mode='w')
        if process_ambiguity_data:
            agg_ambigs = open(aggregated_ambig_counts_filename, mode='w')
            agg_ambig_partners = open(aggregated_ambig_partners_filename, mode='w')

        end_of_counts_header = 0
        end_of_metrics_header = 0
        end_of_ambigs_header = 0
        print_to_stderr('  Concatenating output from all workers.')
        for worker_index in range(total_workers):
            counts_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.counts.tsv' % (analysis_prefix, worker_index, total_workers))
            ambig_counts_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.ambig.counts.tsv' % (analysis_prefix, worker_index, total_workers))
            ambig_partners_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.ambig.partners' % (analysis_prefix, worker_index, total_workers))
            metrics_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.metrics.tsv' % (analysis_prefix, worker_index, total_workers))
            ignored_for_output_filename = counts_output_filename+'.ignored'

            # Counts
            with open(counts_output_filename, 'r') as f:
                shutil.copyfileobj(f, agg_counts)

            # Metrics
            with open(metrics_output_filename, 'r') as f:
                shutil.copyfileobj(f, agg_metrics)

            # Ignored
            if os.path.isfile(counts_output_filename+'.ignored'):
                with open(counts_output_filename+'.ignored', 'r') as f:
                    shutil.copyfileobj(f, agg_ignored)

            if process_ambiguity_data:
                with open(ambig_counts_output_filename, 'r') as f:
                    shutil.copyfileobj(f, agg_ambigs)

                with open(ambig_partners_output_filename, 'r') as f:
                    shutil.copyfileobj(f, agg_ambig_partners)

        print_to_stderr('  GZIPping concatenated output.')
        agg_counts.close()
        subprocess.Popen(['gzip', '-f', aggregated_counts_filename]).wait()
        agg_metrics.close()
        subprocess.Popen(['gzip', '-f', aggregated_quant_metrics_filename]).wait()
        print_to_stderr('Aggregation completed in %s.gz' % aggregated_counts_filename)

        if process_ambiguity_data:
            agg_ambigs.close()
            subprocess.Popen(['gzip', '-f', aggregated_ambig_counts_filename]).wait()
            agg_ambig_partners.close()
            subprocess.Popen(['gzip', '-f', aggregated_ambig_partners_filename]).wait()

        target_bams = [os.path.join(self.paths.quant_dir, '%sworker%d_%d.bam'% (analysis_prefix, worker_index, total_workers)) for worker_index in range(total_workers)]
        target_bams = [t for t in target_bams if os.path.isfile(t)]
        if target_bams:
            print_to_stderr('  Merging BAM files.')
            p1 = subprocess.Popen([self.project.paths.samtools, 'merge', '-f', aggregated_bam_output]+target_bams, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            if p1.wait() == 0:
                print_to_stderr('  Indexing merged BAM file.')
                p2 = subprocess.Popen([self.project.paths.samtools, 'index', aggregated_bam_output], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                if p2.wait() == 0:
                    for filename in target_bams:
                        os.remove(filename)
                        os.remove(filename + '.bai')
                else:
                    print_to_stderr(" === Error in samtools index ===")
                    print_to_stderr(p2.stderr.read())
            else:
                print_to_stderr(" === Error in samtools merge ===")
            print_to_stderr(p1.stderr.read())     

        # print_to_stderr('Deleting per-worker counts files.')
        # for worker_index in range(total_workers):
        #     counts_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.counts.tsv' % (analysis_prefix, worker_index, total_workers))
        #     os.remove(counts_output_filename)

        #     ambig_counts_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.ambig.counts.tsv' % (analysis_prefix, worker_index, total_workers))
        #     os.remove(ambig_counts_output_filename)

        #     ambig_partners_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.ambig.partners' % (analysis_prefix, worker_index, total_workers))
        #     os.remove(ambig_partners_output_filename)

        #     metrics_output_filename = os.path.join(self.paths.quant_dir, '%sworker%d_%d.metrics.tsv' % (analysis_prefix, worker_index, total_workers))
        #     os.remove(metrics_output_filename)

        #     ignored_for_output_filename = counts_output_filename+'.ignored'
        #     os.remove(ignored_for_output_filename)


class LibrarySequencingPart():
    def __init__(self, filtered_fastq_filename=None, project=None, run_name='', library_name='', part_name=''):
        self.project = project
        self.run_name = run_name
        self.part_name = part_name
        self.library_name = library_name
        self.filtered_fastq_filename = filtered_fastq_filename
        self.barcode_counts_pickle_filename = filtered_fastq_filename + '.counts.pickle'
        self.filtering_metrics_filename = '.'.join(filtered_fastq_filename.split('.')[:-1]) + 'metrics.yaml'

        self.sorted_gzipped_fastq_filename = filtered_fastq_filename + '.sorted.fastq.gz'
        self.sorted_gzipped_fastq_index_filename = filtered_fastq_filename + '.sorted.fastq.gz.index.pickle'

    @property
    def is_filtered(self):
        if not hasattr(self, '_is_filtered'):
            self._is_filtered = os.path.exists(self.filtered_fastq_filename) and os.path.exists(self.barcode_counts_pickle_filename)
        return self._is_filtered
    
    @property
    def is_sorted(self):
        if not hasattr(self, '_is_sorted'):
            self._is_sorted = os.path.exists(self.sorted_gzipped_fastq_filename) and os.path.exists(self.sorted_gzipped_fastq_index_filename)
        return self._is_sorted

    @property
    def part_barcode_counts(self):
        if not hasattr(self, '_part_barcode_counts'):
            with open(self.barcode_counts_pickle_filename, 'r') as f:
                self._part_barcode_counts = pickle.load(f)
        return self._part_barcode_counts

    @property
    def sorted_index(self):
        if not hasattr(self, '_sorted_index'):
            with open(self.sorted_gzipped_fastq_index_filename, 'r') as f:
                self._sorted_index = pickle.load(f)
        return self._sorted_index

    def contains_library_in_query(self, query_libraries):
        return self.library_name in query_libraries

    def sort_reads_by_barcode(self, abundant_barcodes={}):
        sorted_barcodes = [j for j,v in sorted(abundant_barcodes.items(), key=lambda i:-i[1][1])]
        sorted_barcodes = [j for j in sorted_barcodes if j in self.part_barcode_counts]

        barcode_buffers = {}
        barcode_gzippers = {}
        for bc in sorted_barcodes + ['ignored']:
            barcode_buffers[bc] = BytesIO()
            barcode_gzippers[bc] = gzip.GzipFile(fileobj=barcode_buffers[bc], mode='wb')

        total_processed_reads = 0
        total_ignored_reads = 0
        bcs_with_data = set()
        bcs_with_tmp_data = set()
        barcode_tmp_filename = lambda bc: '%s.%s.tmp.gz' % (self.sorted_gzipped_fastq_filename, bc)


        total_reads = sum(self.part_barcode_counts.values())
        print_to_stderr('Sorting %d reads from %d barcodes above absolute minimum threshold.' % (total_reads, len(abundant_barcodes)))
        with open(self.filtered_fastq_filename, 'r') as input_fastq:
            for name, seq, qual in from_fastq(input_fastq):
                total_processed_reads += 1
                bc = name.split(':')[0]

                if total_processed_reads%1000000 == 0:
                    print_to_stderr('Read in %.02f percent of all reads (%d)' % (100.*total_processed_reads/total_reads, total_processed_reads))
                
                if bc in abundant_barcodes:
                    barcode_gzippers[bc].write(to_fastq(name, seq, qual))
                    bcs_with_data.add(bc)
                else:
                    total_ignored_reads += 1
                    barcode_gzippers['ignored'].write(to_fastq(name, seq, qual))
                    bcs_with_data.add('ignored')


        sorted_output_index = {}
        with open(self.sorted_gzipped_fastq_filename, 'wb') as sorted_output:
            for original_bc in sorted_barcodes + ['ignored']:
                if original_bc != 'ignored':
                    new_bc_name = abundant_barcodes[original_bc][0]
                    barcode_reads_count = self.part_barcode_counts[original_bc]
                else:
                    new_bc_name = 'ignored'
                    barcode_reads_count = total_ignored_reads

                start_pos = sorted_output.tell()
                barcode_gzippers[original_bc].close()
                if original_bc in bcs_with_data:
                    barcode_buffers[original_bc].seek(0)
                    shutil.copyfileobj(barcode_buffers[original_bc], sorted_output)
                barcode_buffers[original_bc].close()
                end_pos = sorted_output.tell()

                if end_pos > start_pos:
                    sorted_output_index[new_bc_name] = (original_bc, start_pos, end_pos, end_pos-start_pos, barcode_reads_count)

        with open(self.sorted_gzipped_fastq_index_filename, 'w') as f:
            pickle.dump(sorted_output_index, f)      

    def get_reads_for_barcode(self, barcode):
        if barcode not in self.sorted_index:
            raise StopIteration

        original_barcode, start_byte_offset, end_byte_offset, byte_length, barcode_reads = self.sorted_index[barcode]

        with open(self.sorted_gzipped_fastq_filename, 'rb') as sorted_output:
            sorted_output.seek(start_byte_offset)
            byte_buffer = BytesIO(sorted_output.read(byte_length))
            ungzipper = gzip.GzipFile(fileobj=byte_buffer, mode='rb')
            while True:
                yield next(ungzipper)

    @contextmanager
    def trimmomatic_and_low_complexity_filter_process(self):
        """
        We start 3 processes that are connected with Unix pipes.

        Process 1 - Trimmomatic. Doesn't support stdin/stdout, so we instead use named pipes (FIFOs). It reads from FIFO1, and writes to FIFO2. 
        Process 2 - In line complexity filter, a python script. It reads from FIFO2 (Trimmomatic output) and writes to the ouput file. 
        Process 3 - Indexer that counts the number of reads for every barcode. This reads from stdin, writes the reads to stdout and writes the index as a pickle to stderr.

        When these are done, we start another process to count the results on the FastQ file.
        """
        filtered_dir = os.path.dirname(self.filtered_fastq_filename) #We will use the same directory for creating temporary FIFOs, assuming we have write access.
        
        self.filtering_statistics_counter = defaultdict(int)
        with FIFO(dir=filtered_dir) as fifo2, open(self.filtered_fastq_filename, 'w') as filtered_fastq_file, open(self.filtered_fastq_filename+'.counts.pickle', 'w') as filtered_index_file:
            
            low_complexity_filter_cmd = [self.project.paths.python, self.project.paths.trim_polyA_and_filter_low_complexity_reads_py,
                '-input', fifo2.filename, 
                '--min-post-trim-length', self.project.parameters['trimmomatic_arguments']['MINLEN'],
                '--max-low-complexity-fraction', str(self.project.parameters['low_complexity_filter_arguments']['max_low_complexity_fraction']),
                ]
            counter_cmd = [self.project.paths.python,  self.project.paths.count_barcode_distribution_py]

            p2 = subprocess.Popen(low_complexity_filter_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p3 = subprocess.Popen(counter_cmd, stdin=p2.stdout, stdout=filtered_fastq_file, stderr=filtered_index_file)

            with FIFO(dir=filtered_dir) as fifo1:

                trimmomatic_cmd = [self.project.paths.java, '-Xmx500m', '-jar', self.project.paths.trimmomatic_jar,
                        'SE', '-threads', "1", '-phred33', fifo1.filename, fifo2.filename]
                for arg in self.project.parameters['trimmomatic_arguments']['argument_order']:
                    val = self.project.parameters['trimmomatic_arguments'][arg]
                    trimmomatic_cmd.append('%s:%s' % (arg, val))

                p1 = subprocess.Popen(trimmomatic_cmd, stderr=subprocess.PIPE)

                fifo1_filehandle = open(fifo1.filename, 'w')
                yield fifo1_filehandle
                fifo1_filehandle.close()
                trimmomatic_stderr = p1.stderr.read().splitlines()
                if trimmomatic_stderr[2] != 'TrimmomaticSE: Completed successfully':
                    raise Exception('Trimmomatic did not complete succesfully on %s' % filtered_filename)
                trimmomatic_metrics = trimmomatic_stderr[1].split() 
                # ['Input', 'Reads:', #READS, 'Surviving:', #SURVIVING, (%SURVIVING), 'Dropped:', #DROPPED, (%DROPPED)]
                trimmomatic_metrics = {'input' : trimmomatic_metrics[2], 'output': trimmomatic_metrics[4], 'dropped': trimmomatic_metrics[7]}
                p1.wait()

            complexity_filter_metrics = pickle.load(p2.stderr)
            p2.wait()
            p3.wait()


        filtering_metrics = {
            'read_structure' : dict(self.filtering_statistics_counter),
            'trimmomatic' : trimmomatic_metrics,
            'complexity_filter': complexity_filter_metrics,
        }
        with open(self.filtering_metrics_filename, 'w') as f:
            yaml.dump(dict(filtering_metrics), f, default_flow_style=False)


class V1V2Filtering(LibrarySequencingPart):

    def __init__(self, bioread_filename=None, metaread_filename=None, *args, **kwargs):

        self.bioread_filename = bioread_filename
        self.metaread_filename = metaread_filename
        LibrarySequencingPart.__init__(self, *args, **kwargs)


    def filter_and_count_reads(self):
        """
        Input the two raw FastQ files
        Output: 
            - A single fastQ file that uses the read name to store the barcoding information
            - A pickle of the number of reads originating from each barcode 
        """
        # Relevant paths
        r1_filename, r2_filename = self.metaread_filename, self.bioread_filename

        #Get barcode neighborhoods
        bc1s = self.project.gel_barcode1_revcomp_list_neighborhood
        bc2s = self.project.gel_barcode2_revcomp_list_neighborhood 


        # This starts a Trimmomatic process, a low complexity filter process, and will 
        # upon closing, start the barcode distribution counting process.
        last_ping = time.time()
        ping_every_n_reads = 1000000
        ping_header = "{0:>12}{1:>16}{2:>12}{3:>10}{4:>10}{5:>10}{6:>10}{7:>10}{8:>10}{9:>10}"
        ping_header = ping_header.format("Total Reads", "", "Valid Reads", "W1 in R2", "Empty", "No W1", "No polyT", "No BC1", "No BC2", "No UMI")
        ping_template = "{total:12d}    {rate:5.1f} sec/M {Valid:12.1%}{W1_in_R2:10.1%}{empty_read:10.1%}{No_W1:10.1%}{No_polyT:10.1%}{BC1:10.1%}{BC2:10.1%}{Umi_error:10.1%}"
        def print_ping_to_log(last_ping):
            sec_per_mil = (time.time()-last_ping)/(ping_every_n_reads/10**6) if last_ping else 0.0
            total = self.filtering_statistics_counter['Total']
            if total > 0:
                ping_format_data = {k: float(self.filtering_statistics_counter[k])/total for k in ['Valid', 'W1_in_R2', 'empty_read',  'No_W1', 'No_polyT', 'BC1', 'BC2', 'Umi_error']}
                print_to_stderr(ping_template.format(total=total, rate=sec_per_mil, **ping_format_data))


        with self.trimmomatic_and_low_complexity_filter_process() as trim_process:
            #Iterate over the weaved reads
            for r_name, r1_seq, r1_qual, r2_seq, r2_qual in self._weave_fastqs(r1_filename, r2_filename):
                    
                # Check if they should be kept
                keep, result = self._process_reads(r1_seq, r2_seq, valid_bc1s=bc1s, valid_bc2s=bc2s)

                # Write the the reads worth keeping
                if keep:
                    bc, umi = result
                    trim_process.write(to_fastq_lines(bc, umi, r2_seq, r2_qual, r_name))
                    self.filtering_statistics_counter['Valid'] += 1
                else:
                    self.filtering_statistics_counter[result] += 1

                # Track speed per M reads
                self.filtering_statistics_counter['Total'] += 1
                if self.filtering_statistics_counter['Total']%(10*ping_every_n_reads) == 1:
                    print_to_stderr(ping_header)

                if self.filtering_statistics_counter['Total']%ping_every_n_reads == 0:
                    print_ping_to_log(last_ping)
                    last_ping = time.time()

            print_ping_to_log(False)

        print_to_stderr(self.filtering_statistics_counter)

    def _weave_fastqs(self, r1_fastq, r2_fastq):
        """
        Merge 2 FastQ files by returning paired reads for each.
        Returns only R1_seq, R2_seq and R2_qual.
        """

        is_gz_compressed = False
        is_bz_compressed = False
        if r1_fastq.split('.')[-1] == 'gz' and r2_fastq.split('.')[-1] == 'gz':
            is_gz_compressed = True
            
        #Added bz2 support VS
        if r1_fastq.split('.')[-1] == 'bz2' and r2_fastq.split('.')[-1] == 'bz2':
            is_bz_compressed = True

        # Decompress Gzips using subprocesses because python gzip is incredibly slow.
        if is_gz_compressed:    
            r1_gunzip = subprocess.Popen("gzip --stdout -d %s" % (r1_fastq), shell=True, stdout=subprocess.PIPE)
            r1_stream = r1_gunzip.stdout
            r2_gunzip = subprocess.Popen("gzip --stdout -d %s" % (r2_fastq), shell=True, stdout=subprocess.PIPE)
            r2_stream = r2_gunzip.stdout
        elif is_bz_compressed:
            r1_bunzip = subprocess.Popen("bzcat %s" % (r1_fastq), shell=True, stdout=subprocess.PIPE)
            r1_stream = r1_bunzip.stdout
            r2_bunzip = subprocess.Popen("bzcat %s" % (r2_fastq), shell=True, stdout=subprocess.PIPE)
            r2_stream = r2_bunzip.stdout
        else:
            r1_stream = open(r1_fastq, 'r')
            r2_stream = open(r2_fastq, 'r')

        while True:
            #Read 4 lines from each FastQ
            name = next(r1_stream).rstrip()[1:].split()[0] #Read name
            r1_seq = next(r1_stream).rstrip() #Read seq
            next(r1_stream) #+ line
            r1_qual = next(r1_stream).rstrip() #Read qual
            
            next(r2_stream) #Read name
            r2_seq = next(r2_stream).rstrip() #Read seq
            next(r2_stream) #+ line
            r2_qual = next(r2_stream).rstrip() #Read qual
            
            # changed to allow for empty reads (caused by adapter trimming)
            if name:
                yield name, r1_seq, r1_qual, r2_seq, r2_qual
            else:
            # if not r1_seq or not r2_seq:
                break

        r1_stream.close()
        r2_stream.close()

    def _process_reads(self, name, read, valid_bc1s={}, valid_bc2s={}):
        """
        Returns either:
            True, (barcode, umi)
                (if read passes filter)
            False, name of filter that failed
                (for stats collection)
        
        R1 anatomy: BBBBBBBB[BBB]WWWWWWWWWWWWWWWWWWWWWWCCCCCCCCUUUUUUTTTTTTTTTT______________
            B = Barcode1, can be 8, 9, 10 or 11 bases long.
            W = 'W1' sequence, specified below
            C = Barcode2, always 8 bases
            U = UMI, always 6 bases
            T = Beginning of polyT tail.
            _ = Either sequencing survives across the polyT tail, or signal starts dropping off
                (and start being anything, likely with poor quality)
        """

        minimal_polyT_len_on_R1 = 7
        hamming_threshold_for_W1_matching = 3

        w1 = "GAGTGATTGCTTGTGACGCCTT"
        rev_w1 = "AAGGCGTCACAAGCAATCACTC" #Hard-code so we don't recompute on every one of millions of calls
        # If R2 contains rev_W1, this is almost certainly empty library
        if rev_w1 in read:
            return False, 'W1_in_R2'

        # # With reads sufficiently long, we will often see a PolyA sequence in R2. 
        # if polyA in read:
        #     return False, 'PolyA_in_R2'

        # Check for polyT signal at 3' end.
        # 44 is the length of BC1+W1+BC2+UMI, given the longest PolyT
        #BC1: 8-11 bases
        #W1 : 22 bases
        #BC2: 8 bases
        #UMI: 6 bases

        # check for empty reads (due to adapter trimming)
        if not read:
            return False, 'empty_read'
        
        #Check for W1 adapter
        #Allow for up to hamming_threshold errors
        if w1 in name:
            w1_pos = name.find(w1)
            if not 7 < w1_pos < 12:
                return False, 'No_W1'
        else:
            #Try to find W1 adapter at start positions 8-11
            #by checking hamming distance to W1.
            for w1_pos in range(8, 12):
                if string_hamming_distance(w1, name[w1_pos:w1_pos+22]) <= hamming_threshold_for_W1_matching:
                    break
            else:
                return False, 'No_W1'
                
        bc2_pos=w1_pos+22
        umi_pos=bc2_pos+8
        polyTpos=umi_pos+6
        expected_poly_t = name[polyTpos:polyTpos+minimal_polyT_len_on_R1]
        if string_hamming_distance(expected_poly_t, 'T'*minimal_polyT_len_on_R1) > 3:
                 return False, 'No_polyT'
            
        bc1 = str(name[:w1_pos])
        bc2 = str(name[bc2_pos:umi_pos])
        umi = str(name[umi_pos:umi_pos+6])
        
        #Validate barcode (and try to correct when there is no ambiguity)
        if valid_bc1s and valid_bc2s:
            # Check if BC1 and BC2 can be mapped to expected barcodes
            if bc1 in valid_bc1s:
                # BC1 might be a neighboring BC, rather than a valid BC itself. 
                bc1 = valid_bc1s[bc1]
            else:
                return False, 'BC1'
            if bc2 in valid_bc2s:
                bc2 = valid_bc2s[bc2]
            else:
                return False, 'BC2'
            if 'N' in umi:
                return False, 'UMI_error'
        bc = '%s-%s'%(bc1, bc2)
        return True, (bc, umi)

class V3Demultiplexer():

    def __init__(self, library_indices, project=None, part_filename="", input_filename="", run_name="", part_name="", run_version_details="v3"):

        self.run_version_details = run_version_details
        self.input_filename = input_filename
        self.project = project
        self.run_name = run_name
        self.part_name = part_name
        self.libraries = {}
        for lib in library_indices:
            lib_index = lib['library_index']
            lib_name = lib['library_name']
            library_part_filename = part_filename.format(library_name=lib_name, library_index=lib_index)
            self.libraries[lib_index] = LibrarySequencingPart(filtered_fastq_filename=library_part_filename, project=project, run_name=run_name, library_name=lib_name, part_name=part_name)

    def _weave_fastqs(self, fastqs):
        last_extension = [fn.split('.')[-1] for fn in fastqs]
        if all(ext == 'gz' for ext in last_extension):
            processes = [subprocess.Popen("gzip --stdout -d %s" % (fn), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for fn in fastqs]
            streams = [r.stdout for r in processes]
        elif all(ext == 'bz2' for ext in last_extension):
            processes = [subprocess.Popen("bzcat %s" % (fn), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for fn in fastqs]
            streams = [r.stdout for r in processes]
        elif all(ext == 'fastq' for ext in last_extension):
            streams = [open(fn, 'r') for fn in fastqs]
        else:
            raise("ERROR: Different files are compressed differently. Check input.")

        while True:
            names = [next(s)[:-1].split()[0] for s in streams]
            seqs = [next(s)[:-1] for s in streams]
            blanks = [next(s)[:-1]  for s in streams]
            quals = [next(s)[:-1]  for s in streams]

            assert all(name==names[0] for name in names)
            yield names[0], seqs, quals

        for s in streams:
            s.close()


    def _process_reads(self, name, seqs, quals, valid_bc1s={}, valid_bc2s={}, valid_libs={}):
        """
        Returns either:
            True, (barcode, umi)
                (if read passes filter)
            False, name of filter that failed
                (for stats collection)
        """

        r1, r2, r3, r4 = seqs
        if self.run_version_details=='v3-miseq':
            r2 = rev_comp(r2)
            r4 = rev_comp(r4)

        if r3 in valid_libs:
            lib_index = valid_libs[r3]
        else:
            return False, r3, 'Invalid_library_index'

        if r2 in valid_bc1s:
            bc1 = valid_bc1s[r2]
        else:
            return False, lib_index, 'Invalid_BC1'

        orig_bc2 = r4[:8]
        umi = r4[8:8+6]
        polyA = r4[8+6:]

        if orig_bc2 in valid_bc2s:
            bc2 = valid_bc2s[orig_bc2]
        else:
            return False, lib_index, 'Invalid_BC2'

        if 'N' in umi:
            return False, lib_index, 'UMI_contains_N'

        final_bc = '%s-%s' % (bc1, bc2)
        return True, lib_index, (final_bc, umi)


    def filter_and_count_reads(self):
        # Prepare error corrected index sets
        self.sequence_to_index_mapping = {}
        libs = self.libraries.keys()
        self.sequence_to_index_mapping = dict(zip(libs, libs))
        index_neighborhoods = [set(seq_neighborhood(lib, 1)) for lib in libs]
        for lib, clibs in zip(libs, index_neighborhoods):
            # Quick check that error-correction maps to a single index
            for clib in clibs:
                if sum(clib in hood for hood in index_neighborhoods)==1:
                    self.sequence_to_index_mapping[clib] = lib

        # Prepare error corrected barcode sets
        error_corrected_barcodes = self.project.gel_barcode2_list_neighborhood
        error_corrected_rev_compl_barcodes = self.project.gel_barcode2_revcomp_list_neighborhood

        # Open up our context managers
        manager_order = [] #It's imperative to exit managers the opposite order than we open them!
        trim_processes = {}
        trim_processes_managers = {}

        for lib in self.libraries.keys():
            manager_order.append(lib)
            trim_processes_managers[lib] = self.libraries[lib].trimmomatic_and_low_complexity_filter_process()
            trim_processes[lib] = trim_processes_managers[lib].__enter__()

        overall_filtering_statistics = defaultdict(int)

        # Paths for the 4 expected FastQs
        input_fastqs = []
        for r in ['R1', 'R2', 'R3', 'R4']:
            input_fastqs.append(self.input_filename.format(read=r))

        last_ping = time.time()
        ping_every_n_reads = 1000000
        ping_header = "{0:>12}{1:>16}{2:>12}{3:>10}{4:>10}{5:>10}{6:>10}   |" + ''.join("{%d:>12.10}"%i for i in range(7,7+len(manager_order)))
        ping_header = ping_header.format("Total Reads", "", "Valid Reads", "No index", "No BC1", "No BC2", "No UMI", *[self.libraries[k].library_name for k in manager_order])
        ping_template = "{total:12d}    {rate:5.1f} sec/M {Valid:12.1%}{Invalid_library_index:10.1%}{Invalid_BC1:10.1%}{Invalid_BC2:10.1%}{UMI_contains_N:10.1%}   |{"+":>12.1%}{".join(manager_order)+":>12.1%}"
        
        def print_ping_to_log(last_ping):
            sec_per_mil = (time.time() - last_ping)/(float(ping_every_n_reads)/10**6) if last_ping else 0
            total = overall_filtering_statistics['Total']
            ping_format_data = {k: float(overall_filtering_statistics[k])/total for k in ['Valid', 'Invalid_library_index', 'Invalid_BC1',  'Invalid_BC2', 'UMI_contains_N']}
            if overall_filtering_statistics['Valid'] > 0:
                ping_format_data.update({k: float(self.libraries[k].filtering_statistics_counter['Valid'])/overall_filtering_statistics['Valid'] for k in manager_order})
            print_to_stderr(ping_template.format(total=total, rate=sec_per_mil, **ping_format_data))

        common__ = defaultdict(int)
        print_to_stderr('Filtering %s, file %s' % (self.run_name, self.input_filename))
        for r_name, seqs, quals in self._weave_fastqs(input_fastqs):

            # Python 3 compatibility in mind!
            seqs = [s.decode('utf-8') for s in seqs]

            keep, lib_index, result = self._process_reads(r_name, seqs, quals,
                                                    error_corrected_barcodes, error_corrected_rev_compl_barcodes, 
                                                    self.sequence_to_index_mapping)
            common__[seqs[1]] += 1
            if keep:
                bc, umi = result
                bio_read = seqs[0]
                bio_qual = quals[0]
                trim_processes[lib_index].write(to_fastq_lines(bc, umi, bio_read, bio_qual, r_name[1:]))
                self.libraries[lib_index].filtering_statistics_counter['Valid'] += 1
                self.libraries[lib_index].filtering_statistics_counter['Total'] += 1
                overall_filtering_statistics['Valid'] += 1

            else:
                if result != 'Invalid_library_index':
                    self.libraries[lib_index].filtering_statistics_counter[result] += 1
                    self.libraries[lib_index].filtering_statistics_counter['Total'] += 1
                overall_filtering_statistics[result] += 1

            # Track speed per M reads
            overall_filtering_statistics['Total'] += 1

            if overall_filtering_statistics['Total']%(ping_every_n_reads*10)==1:
                print_to_stderr(ping_header)
            
            if overall_filtering_statistics['Total']%ping_every_n_reads == 0:
                print_ping_to_log(last_ping)
                last_ping = time.time()
                
        print_ping_to_log(False)
        # Close up the context managers
        for lib in manager_order[::-1]:
            trim_processes_managers[lib].__exit__(None, None, None)

    def contains_library_in_query(self, query_libraries):
        for lib in self.libraries.values():
            if lib.contains_library_in_query(query_libraries):
                return True
        return False





if __name__=="__main__":

    import sys, argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('project', type=argparse.FileType('r'), help='Project YAML File.')
    parser.add_argument('-l', '--libraries', type=str, help='[all] Library name(s) to work on. If blank, will iterate over all libraries in project.', nargs='?', default='')
    parser.add_argument('-r', '--runs', type=str, help='[all] Run name(s) to work on. If blank, will iterate over all runs in project.', nargs='?', default='')
    parser.add_argument('command', type=str, choices=['info', 'filter', 'identify_abundant_barcodes', 'sort', 'quantify', 'aggregate', 'build_index', 'get_reads', 'output_barcode_fastq'])
    parser.add_argument('--total-workers', type=int, help='[all] Total workers that are working together. This takes precedence over barcodes-per-worker.', default=1)
    parser.add_argument('--worker-index', type=int, help='[all] Index of current worker (the first worker should have index 0).', default=0)
    parser.add_argument('--min-reads', type=int, help='[quantify] Minimum number of reads for barcode to be processed', nargs='?', default=750)
    parser.add_argument('--max-reads', type=int, help='[quantify] Maximum number of reads for barcode to be processed', nargs='?', default=100000000)
    parser.add_argument('--min-counts', type=int, help='[aggregate] Minimun number of UMIFM counts for barcode to be aggregated', nargs='?', default=0)
    parser.add_argument('--analysis-prefix', type=str, help='[quantify/aggregate/convert_bam/merge_bam] Prefix for analysis files.', nargs='?', default='')
    parser.add_argument('--no-bam', help='[quantify] Do not output alignments to bam file.', action='store_true')
    parser.add_argument('--genome-fasta-gz', help='[build_index] Path to gzipped soft-masked genomic FASTA file.')
    parser.add_argument('--ensembl-gtf-gz', help='[build_index] Path to gzipped ENSEMBL GTF file. ')
    parser.add_argument('--mode', help='[build_index] Stringency mode for transcriptome build. [strict|all_ensembl]', default='strict')
    parser.add_argument('--override-yaml', help="[all] Dictionnary to update project YAML with.. [You don't need this.]", nargs='?', default='')

    args = parser.parse_args()
    project = IndropsProject(args.project)
    if args.override_yaml:
        override = eval(args.override_yaml)
        if 'paths' in override:
            project.yaml['paths'].update(override['paths'])
        if 'parameters' in override:
            for k,v in override['parameters'].items():
                project.yaml['parameters'][k].update(v)
        if hasattr(project, '_paths'):
            del project._paths
        if hasattr(project, '_parameters'):
            del project._parameters

    target_libraries = []
    if args.libraries:
        for lib in args.libraries.split(','):
            assert lib in project.libraries
            if lib not in target_libraries:
                target_libraries.append(lib)
    else:
        target_libraries = project.libraries.keys()
    lib_query = set(target_libraries)

    target_runs = []
    if args.runs:
        for run in args.runs.split(','):
            assert run in project.runs
            target_runs.append(run)
    else:
        target_runs = project.runs.keys()

    target_library_parts = []
    for lib in target_libraries:
        for pi, part in enumerate(project.libraries[lib].parts):
            if part.run_name in target_runs:
                target_library_parts.append((lib, pi))

    if args.command == 'info':
        print_to_stderr('Project Name: ' + project.name)
        target_run_parts = []
        for run in target_runs:
            target_run_parts += [part for part in project.runs[run] if part.contains_library_in_query(lib_query)]
        print_to_stderr('Total library parts in search query: ' + str(len(target_run_parts)))

    elif args.command == 'filter':
        target_run_parts = []
        for run in target_runs:
            target_run_parts += [part for part in project.runs[run] if part.contains_library_in_query(lib_query)]

        for part in worker_filter(target_run_parts, args.worker_index, args.total_workers):
            print_to_stderr('Filtering run "%s", library "%s", part "%s"' % (part.run_name, part.library_name if hasattr(part, 'library_name') else 'N/A', part.part_name))
            part.filter_and_count_reads()

    elif args.command == 'identify_abundant_barcodes':
        for library in worker_filter(target_libraries, args.worker_index, args.total_workers):
            project.libraries[library].identify_abundant_barcodes()

    elif args.command == 'sort':
        for library, part_index in worker_filter(target_library_parts, args.worker_index, args.total_workers):
            print_to_stderr('Sorting %s, part "%s"' % (library, project.libraries[library].parts[part_index].filtered_fastq_filename))
            project.libraries[library].sort_reads_by_barcode(index=part_index)

    elif args.command == 'quantify':
        for library in target_libraries:
            project.libraries[library].quantify_expression(worker_index=args.worker_index, total_workers=args.total_workers,
                    min_reads=args.min_reads, max_reads=args.max_reads, min_counts=args.min_counts,
                    analysis_prefix=args.analysis_prefix,
                    no_bam=args.no_bam, run_filter=target_runs)

            for part in project.libraries[library].parts:
                if hasattr(part, '_sorted_index'):
                    del part._sorted_index

    elif args.command == 'aggregate':
        for library in target_libraries:
            project.libraries[library].aggregate_counts(analysis_prefix=args.analysis_prefix)

    elif args.command == 'build_index':
        project.build_transcriptome(args.genome_fasta_gz, args.ensembl_gtf_gz, mode=args.mode)

    elif args.command == 'get_reads':
        for library in target_libraries:
            sorted_barcode_names = project.libraries[library].sorted_barcode_names(min_reads=args.min_reads, max_reads=args.max_reads)
            for bc in sorted_barcode_names:
                for line in project.libraries[library].get_reads_for_barcode(bc, run_filter=target_runs):
                    sys.stdout.write(line)

            for part in project.libraries[library].parts:
                if hasattr(part, '_sorted_index'):
                    del part._sorted_index

    elif args.command == 'output_barcode_fastq':
        for library in target_libraries:
            project.libraries[library].output_barcode_fastq(worker_index=args.worker_index, total_workers=args.total_workers,
                    min_reads=args.min_reads, max_reads=args.max_reads, analysis_prefix=args.analysis_prefix, run_filter=target_runs)