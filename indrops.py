import os, subprocess
import itertools
import operator
from collections import defaultdict
import errno

# cPickle is a faster version of pickle that isn't installed in python3
# inserted try statement just in case
try:
   import cPickle as pickle
except:
   import pickle

import numpy as np
import re

# product: product(A, B) returns the same as ((x,y) for x in A for y in B).
# combination: Return r length subsequences of elements from the input iterable.
from itertools import product, combinations
import time
import matplotlib
import yaml

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

def rev_comp(seq):
    tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(tbl[s] for s in seq[::-1])

def to_fastq_lines(bc, umi, seq, qual):
    """
    Return string that can be written to fastQ file
    """
    return '@'+bc+':'+umi+'\n'+seq+'\n+\n'+qual+'\n'

def from_fastq(handle):
    while True:
        name = next(handle).rstrip()[1:] #Read name
        seq = next(handle).rstrip() #Read seq
        next(handle) #+ line
        qual = next(handle).rstrip() #Read qual
        if not name or not seq or not qual:
            break
        yield name, seq, qual

def build_barcode_neighborhoods(barcode_file):
    """
    Given a set of barcodes, produce sequences which can unambiguously be
    mapped to these barcodes, within 2 substitutions. If a sequence maps to 
    multiple barcodes, get rid of it. However, if a sequences maps to a bc1 with 
    1change and another with 2changes, keep the 1change mapping.
    """

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
    
    # contains all mutants that map uniquely to a barcode
    clean_mapping = dict()

    # contain single or double mutants 
    mapping1 = defaultdict(set)
    mapping2 = defaultdict(set)
    
    #Build the full neighborhood and iterate through barcodes
    with open(barcode_file, 'rU') as f:
        # iterate through each barcode (rstrip cleans string of whitespace)
        for line in f:
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

def print_to_log(msg):
    """
    Wrapper to eventually log in smart way, instead of using 'print()'
    """
    sys.stderr.write(str(msg)+'\n')

def parallelized_using_workers(original_func):

    def func_wrapper(self, barcodes_per_worker=0, worker_index=0, **kwargs):
        with open(self.output_paths['good_barcodes_with_names'], 'r') as f:
                sorted_barcode_names = sorted(pickle.load(f).values())

        if barcodes_per_worker == 0:
            barcodes_per_worker = len(sorted_barcode_names)+1

        for i in range(worker_index*barcodes_per_worker, (worker_index+1)*barcodes_per_worker):
            if i >= len(sorted_barcode_names):
                break
            chosen_barcode = sorted_barcode_names[i]

            # Actually execute the method being wrapped!
            original_func(self, chosen_barcode)

    return func_wrapper


# -----------------------
#
# Core objects
#
# -----------------------


class IndropsAnalysis():

    def __init__(self, parameters):
        self.parameters = parameters
        self.user_paths = parameters['project_paths']
        self.user_paths.update(parameters['general_paths'])

        self.output_paths = {
            'read_fail_counts': 'stats/filtering_metrics.yaml',
            'barcode_histogram': 'stats/barcode_abundance_histogram.png',
            'filtered_fastq': 'pre_split/filtered.fastq', #Fastq file after removal of bad reads, but before split
            'barcode_read_counts': 'pre_split/barcode_read_counts.pickle',
            'good_barcodes_with_names': 'pre_split/good_barcodes_with_names.pickle',
            'split_dir': 'post_split/',
            'split_fastq_dir': 'post_split/filtered_fastq/', #Directory where individual barcode fastqs will be placed
            'split_trimmed_fastq_dir': 'post_split/trimmed_fastq/', #Directory where individual barcode fastqs will be placed after trimming
            'split_quant_dir': 'post_split/quant_output/', #Directory where quantification output will be placed
            'aggregated_counts': 'aggregated_counts/',
        }

        # Update all these paths relative to the user-specified output dir.
        for k, rel_path in self.output_paths.items():
            abs_path = os.path.join(self.user_paths['output_dir'], rel_path)
            self.output_paths[k] = abs_path

            #Check, or create relevant output dirs
            check_dir(os.path.dirname(abs_path))

    def filter_and_count_reads(self):
        """
        Input the two raw FastQ files
        Output: 
            - A single fastQ file that uses the read name to store the barcoding information
            - A pickle of number of reads originating from each barcode 
        """
        #Prepare data collection
        barcode_read_counter = defaultdict(int)
        filter_fail_counter = defaultdict(int)
        kept_reads = 0

        #Get barcode neighborhoods
        bc1s = build_barcode_neighborhoods(self.user_paths['gel_barcode1_list'])
        bc2s = build_barcode_neighborhoods(self.user_paths['gel_barcode2_list'])

        i = 0
        start_time = time.time()
        
        with open(self.output_paths['filtered_fastq'], 'w') as output_fastq:
            for r1_seq, r1_qual, r2_seq, r2_qual in self._weave_fastqs(self.user_paths['raw_R1_fastq'], self.user_paths['raw_R2_fastq']):
                    
                # We currently ignore R1 qualities
                keep, result = self._process_reads(r1_seq, r2_seq, valid_bc1s=bc1s, valid_bc2s=bc2s)

                # why not look at quality of r1_seq as well to toss out extra reads here?
                i += 1
                if i%1000000 == 0:
                    sec_per_mil = (time.time()-start_time)/(float(i)/10**6)
                    print_to_log('%d reads parsed, kept %d reads, %.02f seconds per M reads.' % (i, kept_reads, sec_per_mil))

                if keep:
                    bc, umi = result
                    kept_reads += 1
                    barcode_read_counter[bc] += 1
                    output_fastq.write(to_fastq_lines(bc, umi, r2_seq, r2_qual))
                else:
                    filter_fail_counter[result] += 1

        #Save results
        with open(self.output_paths['barcode_read_counts'], 'w') as f:
            pickle.dump(dict(barcode_read_counter), f)


        filtering_statistics = {
            'Total Reads' : i,
            'Valid Reads' : kept_reads,
            'Rejected Reads' : i - kept_reads,
            'Valid Fraction' : float(kept_reads)/i,
            'Rejection Flags' : dict(filter_fail_counter)
        }
        with open(self.output_paths['read_fail_counts'], 'w') as f:
            yaml.dump(dict(filtering_statistics), f, default_flow_style=False)

        print_to_log('%d reads parsed, kept %d reads.' % (i, kept_reads))

    def _weave_fastqs(self, r1_fastq, r2_fastq):
        """
        Merge 2 FastQ files by returning paired reads for each.
        Returns only R1_seq, R2_seq and R2_qual.
        """

        is_gz_compressed = False
        if r1_fastq.split('.')[-1] == 'gz' and r2_fastq.split('.')[-1] == 'gz':
            is_gz_compressed = True

        # Decompress Gzips using subprocesses because python gzip is incredibly slow.
        if is_gz_compressed:    
            r1_gunzip = subprocess.Popen("gzip --stdout -d %s" % (r1_fastq), shell=True, stdout=subprocess.PIPE)
            r1_stream = r1_gunzip.stdout
            r2_gunzip = subprocess.Popen("gzip --stdout -d %s" % (r2_fastq), shell=True, stdout=subprocess.PIPE)
            r2_stream = r2_gunzip.stdout
        else:
            r1_stream = open(r1_fastq, 'r')
            r2_stream = open(r2_fastq, 'r')


        while True:
            #Read 4 lines from each FastQ
            name = next(r1_stream) #Read name
            r1_seq = next(r1_stream).rstrip() #Read seq
            next(r1_stream) #+ line
            r1_qual = next(r1_stream).rstrip() #Read qual
            
            next(r2_stream) #Read name
            r2_seq = next(r2_stream).rstrip() #Read seq
            next(r2_stream) #+ line
            r2_qual = next(r2_stream).rstrip() #Read qual
            
            # changed to allow for empty reads (caused by adapter trimming)
            if name:
                yield r1_seq, r1_qual, r2_seq, r2_qual
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
        # 47 is the length of BC1+W1+BC2+UMI, given the longest PolyT
        expected_poly_t = name[47:47+minimal_polyT_len_on_R1:]
        if string_hamming_distance(expected_poly_t, 'T'*minimal_polyT_len_on_R1) > 3:
            return False, 'No_polyT'

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
            
        bc1 = name[:w1_pos]
        bc2 = str(name[w1_pos+22:w1_pos+22+8])
        umi = str(name[w1_pos+22+8: w1_pos+22+8+6])
        
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
        bc = '%s-%s'%(bc1, bc2)
        return True, (bc, umi)


    def make_barcode_abundance_histogram(self):
        """
        Takes the read-count-by-barcode pickle and outputs a histogram used 
        to determine a treshold on the minimal number of reads coming from good barcodes
        """
        with open(self.output_paths['barcode_read_counts'], 'r') as f:
            barcode_read_counter = pickle.load(f)

        count_freq = defaultdict(int)
        for bc, count in barcode_read_counter.items():
            count_freq[count] += 1

        x = np.array(count_freq.keys())
        y = np.array(count_freq.values())
        w = x*y

        # need to use Agg backend (which is noninteractive)
        matplotlib.use('Agg')
            
        from matplotlib import pyplot as plt
        ax = plt.subplot(111)
        ax.hist(x, bins=np.logspace(0, 6, 50), weights=w)
        ax.set_xscale('log')
        ax.set_xlabel('Reads per barcode')
        ax.set_ylabel('#reads coming from bin')
        plt.savefig(self.output_paths['barcode_histogram'])

        print_to_log("Created Barcode Abundance Histogram at:")
        print_to_log("  " + self.output_paths['barcode_histogram'])

    def choose_good_barcodes(self, threshold=15000):
        """
        Takes the read-count-by-barcode pickle and a minimal threshold value, 
        Outputs a list of barcodes to the retained and assigns a 'bcN' name to each barcode
        """
        with open(self.output_paths['barcode_read_counts'], 'r') as f:
            barcode_read_counter = pickle.load(f)

        good_barcodes = []
        for bc, count in barcode_read_counter.items():
            if count >= threshold:
                good_barcodes.append(bc)

        
        str_length_of_largest_bc = len(str(len(good_barcodes)+1))
        name_fmt = 'bc%0' + str(str_length_of_largest_bc) + 'd'

        barcode_names = {}
        for i, bc in enumerate(sorted(good_barcodes, key=lambda b: barcode_read_counter[b], reverse=True)):
            barcode_names[bc] = name_fmt % (i+1)

        with open(self.output_paths['good_barcodes_with_names'], 'w') as f:
            pickle.dump(barcode_names, f)

        print_to_log('Keeping %d barcodes.' % len(barcode_names))

    def split_reads_by_barcode(self, output_unassigned_reads=False):
        """
        # TODO : Connect 'output_unassigned_reads' to outside parameters!
        Starts with the list of good barcodes and the filtered FastQ
        Splits the filtered FastQ into one FastQ file per read. 
        """
        with open(self.output_paths['good_barcodes_with_names'], 'r') as f:
            barcode_names = pickle.load(f)

        #Because we use 'a' to write our split fastqs, we need to make sure everythign is empty.
        for fn in os.listdir(self.output_paths['split_fastq_dir']):
            os.remove(os.path.join(self.output_paths['split_fastq_dir'], fn))

        total_processed_reads = 0
        start_time = time.time()

        pre_write_buffer_size = 0
        pre_write = defaultdict(list)
        unassigned_filename = os.path.join(self.output_paths['split_fastq_dir'], 'unassigned.fastq')

        with open(self.output_paths['filtered_fastq'], 'r') as input_fastq:
            for name, seq, qual in from_fastq(input_fastq):
                pre_write_buffer_size += 1
                total_processed_reads += 1
                bc, umi = name.split(':')
                
                if bc in barcode_names:
                    bc_name = barcode_names[bc]
                    filename = os.path.join(self.output_paths['split_fastq_dir'], '%s.fastq' % bc_name)
                    pre_write[filename].append(to_fastq_lines(bc, umi, seq, qual))

                elif output_unassigned_reads:
                    pre_write[unassigned_filename].append(to_fastq_lines(bc, umi, seq, qual))


                if pre_write_buffer_size % 1000000 == 0:
                    for fn, chunks in pre_write.items():
                        with open(fn, 'a') as out:
                            for chunk in chunks:
                                out.write(chunk)
                    pre_write_buffer_size = 0
                    pre_write = defaultdict(list)

                if total_processed_reads % 1000000 == 0:
                    sec_per_mil = (time.time()-start_time)/(float(total_processed_reads)/10**6)
                    print_to_log('%d reads processed, %.02f seconds per M reads.' % (total_processed_reads, sec_per_mil))

        #Make sure we write anything possibly left in 'pre_write'
        for fn, chunks in pre_write.items():
            with open(fn, 'a') as out:
                for chunk in chunks:
                    out.write(chunk)

    @parallelized_using_workers
    def quantify_expression_for_barcode(self, barcode):

        fastq_input = os.path.join(self.output_paths['split_fastq_dir'], '%s.fastq' % barcode)
        if not os.path.isfile(fastq_input):
            print_to_log(fastq_input)
            print_to_log("Does not exists.")
            return 

        trimmed_fastq = os.path.join(self.output_paths['split_trimmed_fastq_dir'], '%s.fastq' % barcode)

        counts_output = os.path.join(self.output_paths['split_quant_dir'], '%s.counts' % barcode)
        quant_metrics_output = os.path.join(self.output_paths['split_quant_dir'], '%s.quant_metrics' % barcode)
        oversequencing_metrics_output = os.path.join(self.output_paths['split_quant_dir'], '%s.oversequencing' % barcode)
        unaligned_reads_output = os.path.join(self.output_paths['split_quant_dir'], '%s.unaligned.fastq' % barcode)
        aligned_bam_output = os.path.join(self.output_paths['split_quant_dir'], '%s.aligned.bam' % barcode)

        # Build Trimmomatic Trim command
        trimmomatic_cmd = ['java', '-jar', self.user_paths['trimmomatic'], 'SE', '-threads', "1", '-phred33', fastq_input, trimmed_fastq]
        for arg, val in self.parameters['trimmomatic_arguments'].items():
            trimmomatic_cmd.append('%s:%s' % (arg, val))

        subprocess.call(trimmomatic_cmd)

        # Build Alignment and Quantification command
        bowtie_exec = os.path.join(self.user_paths['bowtie'], 'bowtie')

        bowtie_cmd = [bowtie_exec, self.user_paths['bowtie_index_prefix'], '-q', trimmed_fastq,
            '-p', '1', '-a', '--best', '--strata', '--chunkmbs', '1000', '--norc', '--sam',
            '-m', str(self.parameters['bowtie_arguments']['m']),
            '-n', str(self.parameters['bowtie_arguments']['n']),
            '-l', str(self.parameters['bowtie_arguments']['l']),
            '-e', str(self.parameters['bowtie_arguments']['e']),
            ]

        if self.parameters['output_arguments']['output_unaligned_reads_to_other_fastq']:
            bowtie_cmd += ['--un', unaligned_reads_output]


        script_dir = os.path.dirname(os.path.realpath(__file__))

        quant_cmd = ['python', os.path.join(script_dir, 'filter_alignments.py'),
            '-m', str(self.parameters['umi_quantification_arguments']['m']),
            '-u', str(self.parameters['umi_quantification_arguments']['u']),
            '-d', str(self.parameters['umi_quantification_arguments']['d']),
            '--min_non_polyA', str(self.parameters['umi_quantification_arguments']['min_non_polyA']),
            '--counts', counts_output,
        ]
        if self.parameters['umi_quantification_arguments']['split-ambigs']:
            quant_cmd.append('--split-ambig')
        if self.parameters['output_arguments']['output_oversequencing_metrics']:
            quant_cmd += ['--umifm_oversequencing', oversequencing_metrics_output]
        if self.parameters['output_arguments']['output_umifm_calculation_metrics']:
            quant_cmd += ['--metrics', output_umifm_calculation_metrics]


        final_pipe = aligned_bam_output if self.parameters['output_arguments']['output_alignment_to_bam'] else '/dev/null'
        final_cmd = ' '.join(bowtie_cmd) + ' | ' + ' '.join(quant_cmd) + ' > ' + final_pipe

        subprocess.call(final_cmd, shell=True)


    def aggregate_counts(self):
        with open(self.output_paths['good_barcodes_with_names'], 'r') as f:
            sorted_barcode_names = sorted(pickle.load(f).values())

        counts_data = defaultdict(lambda : defaultdict(float))
        ambig_counts_data = defaultdict(lambda : defaultdict(float))
        ambiguity_partners = defaultdict(set)

        full_gene_list = set()

        print_to_log('Started aggregating barcodes.')
        i = 0
        missing_barcodes = []
        for barcode in sorted_barcode_names:
            i += 1
            if (i % 100)== 0:
                print_to_log('Read %d barcodes.' % i)
            counts_filename = os.path.join(self.output_paths['split_quant_dir'], '%s.counts' % barcode)

            # Check that the file we expect is actually there
            if not os.path.isfile(counts_filename):
                missing_barcodes.append(barcode)
                continue

            with open(counts_filename, 'r') as f:
                
                # If we have an empty file, reading the header will raise 'StopIteration'
                # so we can catch it here
                try:
                    header = next(f).rstrip('\n').split('\t')
                except StopIteration:
                    print_to_log(str(int(barcode[2:])-1))
                    continue

                for line in f:
                    row = line.rstrip('\n').split('\t')
                    gene = row[0]
                    counts = float(row[1])
                    ambig_counts = float(row[2])
                    partners = set(row[3].split())

                    full_gene_list.add(gene)
                    if counts > 0:
                        counts_data[gene][barcode] = counts
                    if ambig_counts > 0:
                        ambig_counts_data[gene][barcode] = ambig_counts
                    ambiguity_partners[gene] = ambiguity_partners[gene].union(partners)

        print_to_log('Missed the following barcodes: '+ ','.join(missing_barcodes))
        print_to_log('Corresponding indices: ')
        print_to_log(' '.join([str(int(bc[2:])-1) for bc in missing_barcodes]))


        print_to_log('Finished processing')
        output_filename = os.path.join(self.output_paths['aggregated_counts'], 'full_counts.txt')
        output_file = open(output_filename, 'w')
        ambig_filename = os.path.join(self.output_paths['aggregated_counts'], 'ambig_counts.txt')
        ambig_file = open(ambig_filename, 'w')

        output_header = ['gene'] + ['Sum_counts', 'Sum_ambig', 'Ambigs'] + sorted_barcode_names
        to_output_line = lambda row: '%s\n' % '\t'.join([str(r) for r in row])

        output_file.write(to_output_line(output_header))
        ambig_file.write(to_output_line(output_header))

        print_to_log('Starting output')
        for gene in sorted(full_gene_list):
            per_sample_counts = [counts_data[gene][s] for s in sorted_barcode_names]
            per_sample_ambig_counts = [ambig_counts_data[gene][s] for s in sorted_barcode_names]

            counts_row = [gene] + [sum(per_sample_counts), sum(per_sample_ambig_counts), ' '.join(ambiguity_partners[gene])] + per_sample_counts
            ambig_counts_row = [gene] + [sum(per_sample_counts), sum(per_sample_ambig_counts), ' '.join(ambiguity_partners[gene])] + per_sample_ambig_counts

            output_file.write(to_output_line(counts_row))
            ambig_file.write(to_output_line(ambig_counts_row))

        print_to_log('Aggregation completed in %s' % output_filename)

    @parallelized_using_workers
    def sort_and_index_bam(self, chosen_barcode):

        aligned_bam = os.path.join(self.output_paths['split_quant_dir'], '%s.aligned.bam' % chosen_barcode)
        genomic_bam = os.path.join(self.output_paths['split_quant_dir'], '%s.genomic.bam' % chosen_barcode)
        sorted_bam = os.path.join(self.output_paths['split_quant_dir'], '%s.genomic.sorted.bam' % chosen_barcode)

        if os.path.isfile(aligned_bam):
            samtools_exec = os.path.join(self.user_paths['samtools'], 'samtools')
            rsem = os.path.join(self.user_paths['rsem'], 'rsem-tbam2gbam')

            rsem_transcriptome_to_bam_cmd = [rsem, self.user_paths['bowtie_index_prefix'], aligned_bam, genomic_bam]
            sort_cmd = [samtools_exec, 'sort', '-o', sorted_bam, '-O', 'bam', '-T', chosen_barcode+'.temp', genomic_bam]
            index_cmd = [samtools_exec, 'index', sorted_bam]

            subprocess.call(rsem_transcriptome_to_bam_cmd)
            subprocess.call(sort_cmd)
            subprocess.call(index_cmd)
        else:
            print_to_log('File not found: ' + aligned_bam)
def build_transcriptome_from_ENSEMBL_files(input_fasta_filename, bowtie_index_prefix,
    gtf_filename="", 
    polyA_tail_length=50,
    fasta_line_width=60,
    bowtie_path=""):

    output_fasta_filename = os.path.join(os.path.dirname(input_fasta_filename), '.'.join(input_fasta_filename.split('.')[:-1])+'.annotated.fa')

    print_to_log('Annotating transcriptome in %s\nSaving to %s' % (input_fasta_filename, output_fasta_filename))

    def parse_gtf_attributes(attr):
        attr_list = attr.split('; ')[:-1]
        attrs = dict([a[:-1].split(' "') for a in attr_list])
        return attrs

    def fasta_read_iter(filename):
        """
        Simple function to iterate over records in a FASTA file
        """
        curr_name, curr_seq = None, ""
        with open(filename, 'r') as f:
            for line in f:
                if line[0] == '>':
                    if curr_name is not None:
                        yield curr_name, curr_seq

                    curr_name = line.strip()[1:]
                    curr_seq = ""
                else:
                    curr_seq += line.rstrip()

    accepted_transcript_biotypes = set(["lincRNA",
        "translated_processed_pseudogene",
        "protein_coding",
        "antisense",
        "misc_RNA",
        "transcribed_processed_pseudogene",
        "processed_transcript",
        "processed_pseudogene",])

    # Use ENSEMBL GTF File to create mapping from 'ENSG***' IDs to gene names.

    gene_id_to_gene_name = {}
    if gtf_filename:

        print_to_log('Using GTF file %s to infer gene names.' % (gtf_filename))

        gtf_line_counter = 0
        for line in subprocess.Popen(["gzip", "--stdout", "-d", gtf_filename], stdout=subprocess.PIPE).stdout:
            if line[0] == '#':
                continue

            split_line = line.rstrip().split('\t')
            line_attrs = parse_gtf_attributes(split_line[-1])
           
            gtf_line_counter += 1
            if gtf_line_counter % 100000 == 0:
                print_to_log('Processed %d lines from GTF.' % gtf_line_counter)

            if 'gene_id' in line_attrs:
                gene_id_to_gene_name[line_attrs['gene_id']] = line_attrs['gene_name']

    output_fasta = open(output_fasta_filename, 'w')

    for name, seq in fasta_read_iter(input_fasta_filename):
        split_name = name.split()
        transcript_id = split_name[0]
        transcript_attrs = dict(s.split(':') for s in split_name[1:] if len(s.split(':')) == 2)

        
        if transcript_attrs['transcript_biotype'] in accepted_transcript_biotypes:
            # Always try to convert gene IDs to gene names. 
            # If no GTF was specified, then gene_id_to_gene_name is empty.
            if transcript_attrs['gene'] in gene_id_to_gene_name:
                gene_name = gene_id_to_gene_name[transcript_attrs['gene']]
            else:
                gene_name = transcript_attrs['gene']


            # Name is now 'Transcript_ID|Gene'
            new_name = '%s|%s' % (transcript_id, gene_name)
            # Append polyA tail of desired length.
            seq = seq + 'A'*polyA_tail_length

            # Output name
            output_fasta.write('>' + new_name + '\n')

            # Output sequence, with fixed line width
            for i in range(0, len(seq), fasta_line_width):
                output_fasta.write(seq[i:i+fasta_line_width] + '\n')

    output_fasta.close()


    print_to_log('Transcriptome Annotation completed.')
    print_to_log('Building index.')

    bowtie_build_exec = os.path.join(bowtie_path, 'bowtie-build')
    subprocess.call([bowtie_build_exec, output_fasta_filename, bowtie_index_prefix])



if __name__=="__main__":

    import sys, argparse
    parser = argparse.ArgumentParser()
    
    subparsers = parser.add_subparsers(dest='command', help='Choose among the possible commands.')

    parser_index = subparsers.add_parser('index')
    parser_index.add_argument('--polyA', type=int,
        help='Length of polyA tail to append to transcript sequences.', default=50)
    parser_index.add_argument('transcriptome', type=str,
        help='ENSEMBL Transcriptome Fasta. Example for human: ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz',)
    parser_index.add_argument('index_prefix', type=str,
        help='Index Prefix. Add this path to your parameters.yaml file.',)
    parser_index.add_argument('--gtf-gz', type=str,
        help='Matching ENSEMBL gzipped GTF Fasta. Example for human: ftp://ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz',
        default="")
    parser_index.add_argument('--bowtie-path', type=str,
        help='Directory containing Bowtie (1) executables',
        default="")

    parser_preprocess = subparsers.add_parser('preprocess')
    parser_preprocess.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    
    parser_split = subparsers.add_parser('split_barcodes')
    parser_split.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_split.add_argument('--read-count-threshold', type=int,
        help="Minimal read count to be considered 'abundant' barcode")

    parser_quantify = subparsers.add_parser('quantify')
    parser_quantify.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_quantify.add_argument('--barcodes-per-worker', type=int, help='Barcodes to be processed by each worker.', default=0)
    parser_quantify.add_argument('--worker-index', type=int, help='Index of current worker. (Starting at 0). Make sure max(worker-index)*(barcodes-per-worker) > total barcodes.', default=0)

    parser_aggregate = subparsers.add_parser('aggregate')
    parser_aggregate.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')

    parser_sort_bam = subparsers.add_parser('sort_bam')
    parser_sort_bam.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_sort_bam.add_argument('--barcodes-per-worker', type=int, help='Barcodes to be processed by each worker.', default=0)
    parser_sort_bam.add_argument('--worker-index', type=int, help='Index of current worker. (Starting at 0). Make sure max(worker-index)*(barcodes-per-worker) > total barcodes.', default=0)


    args = parser.parse_args()

    if args.command == 'index':
        build_transcriptome_from_ENSEMBL_files(args.transcriptome,
            args.index_prefix,
            gtf_filename=args.gtf_gz, 
            polyA_tail_length=args.polyA,
            bowtie_path=args.bowtie_path)

    else:
        parameters = yaml.load(args.parameters)
        analysis = IndropsAnalysis(parameters)
        if args.command == 'preprocess':
            analysis.filter_and_count_reads()
            analysis.make_barcode_abundance_histogram()
        elif args.command == 'split_barcodes':
            analysis.choose_good_barcodes(args.read_count_threshold)
            analysis.split_reads_by_barcode()
        elif args.command == 'quantify':
            analysis.quantify_expression_for_barcode(args.barcodes_per_worker, args.worker_index)
        elif args.command == 'aggregate':
            analysis.aggregate_counts()
        elif args.command == 'sort_bam':
            analysis.sort_and_index_bam(args.barcodes_per_worker, args.worker_index)
        