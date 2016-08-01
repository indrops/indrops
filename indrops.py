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
import shutil

# product: product(A, B) returns the same as ((x,y) for x in A for y in B).
# combination: Return r length subsequences of elements from the input iterable.
from itertools import product, combinations
import time
import matplotlib
import yaml

import tempfile

from contextlib import contextmanager

# -----------------------
#
# Default parameters
#
# -----------------------

import subprocess, threading

class ThreadedSubprocess(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            print 'Thread started'
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()
            print 'Thread finished'

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print 'Terminating process'
            self.process.terminate()
            thread.join()
        print self.process.returncode

# command = Command("echo 'Process started'; sleep 2; echo 'Process finished'")
# command.run(timeout=3)
# command.run(timeout=1)




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

def print_to_log(msg):
    """
    Wrapper to eventually log in smart way, instead of using 'print()'
    """
    sys.stderr.write(str(msg)+'\n')

def parallelized_using_workers(original_func):
    """
    Wrapper to help parallelize functions that independently operate over barcodes.
    For parallelization, specify a number of barcodes_per_worker, and a specific worker index. 
    For testing, specify a single barcode to get processed. 
    """

    def func_wrapper(self, barcodes_per_worker=0, worker_index=0, target_barcode=None, total_workers=0, missing='', **kwargs):

        if target_barcode is not None:
            original_func(self, target_barcode)

        else: 
            with open(self.output_paths['good_barcodes_with_names'], 'r') as f:
                sorted_barcode_names = sorted(pickle.load(f).values())

            if missing:
                # Chose barcodes to run based on a missing pattern. 
                barcodes_for_this_worker = []

                dirname = os.path.dirname(self.user_paths['output_dir'] + missing)
                basename = os.path.basename(self.user_paths['output_dir'] + missing)
                existing_files = set(os.listdir(dirname))

                for bc in sorted_barcode_names:
                    if basename % bc not in existing_files:
                        barcodes_for_this_worker.append(bc)

                print_to_log('Round-up of missing barcodes:')
                print_to_log(barcodes_for_this_worker)

            elif total_workers > 0: #Total workers over-rides barcodes per worker.
                barcodes_for_this_worker = []
                i = worker_index
                while i < len(sorted_barcode_names):
                    barcodes_for_this_worker.append(sorted_barcode_names[i])
                    i += total_workers
            elif barcodes_per_worker == 0 and total_workers == 0:
                barcodes_for_this_worker = sorted_barcode_names
            else:
                barcodes_for_this_worker = [sorted_barcode_names[i] for i in range(worker_index*barcodes_per_worker, (worker_index+1)*barcodes_per_worker)]

            for chosen_barcode in barcodes_for_this_worker:
                # Actually execute the method being wrapped!
                original_func(self, chosen_barcode)

    return func_wrapper

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
            os.rmdir(self.tmpdir)


# -----------------------
#
# Core objects
#
# -----------------------

class v2Demultiplexer():

    def __init__(self, yaml_parameters_file, split_file_index=0):
        #Pass the yaml parameters path so we can use it to try to find the library-specific yamls
        self.parameters = yaml.load(yaml_parameters_file)
        self.split_file_index = split_file_index
        self.libraries = {}
        self.sequence_to_index_mapping = {}

        for lib_index, lib_yaml in self.parameters['library_indices'].items():
            lib_yaml_path = os.path.abspath(os.path.join(os.path.dirname(yaml_parameters_file.name), lib_yaml))
            with open(lib_yaml_path) as f:
                lib_params = yaml.load(f)
            self.libraries[lib_index] = IndropsAnalysis(lib_params, split_file_index)

            new_affixes = ['%s_%s_%s' % (self.parameters['run_name'], lib_index, affix) for affix in self.parameters['split_affixes']]
            self.libraries[lib_index].build_filtered_reads_paths(new_affixes)

        libs = self.libraries.keys()
        self.sequence_to_index_mapping = dict(zip(libs, libs))
        if self.parameters['error_correct_library_indices']:
            index_neighborhoods = [set(seq_neighborhood(lib, 1)) for lib in libs]
            for lib, clibs in zip(libs, index_neighborhoods):
                # Quick check that error-correction maps to a single index
                for clib in clibs:
                    if sum(clib in hood for hood in index_neighborhoods)==1:
                        self.sequence_to_index_mapping[clib] = lib

    def _weave_fastqs(self, fastqs):
        last_extension = [fn.split('.')[-1] for fn in fastqs]
        if all(ext == 'gz' for ext in last_extension):
            processes = [subprocess.Popen("gzip --stdout -d %s" % (fn), shell=True, stdout=subprocess.PIPE) for fn in fastqs]
            streams = [r.stdout for r in processes]
        elif all(ext == 'bz2' for ext in last_extension):
            processes = [subprocess.Popen("bzcat %s" % (fn), shell=True, stdout=subprocess.PIPE) for fn in fastqs]
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


    def _v2_process_reads(self, name, seqs, quals, valid_bc1s={}, valid_bc2s={}, valid_libs={}):
        """
        Returns either:
            True, (barcode, umi)
                (if read passes filter)
            False, name of filter that failed
                (for stats collection)
        """

        r1, r2, r3, r4 = seqs

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

    def v2_filter_and_count_reads(self):

        # Prepare error corrected barcode sets
        error_corrected_barcodes = build_barcode_neighborhoods(self.parameters['barcode_list'], False)
        error_corrected_rev_compl_barcodes = build_barcode_neighborhoods(self.parameters['barcode_list'], True)

        # Open up our context managers
        trim_processes = {}
        trim_processes_managers = {}
        filtering_statistics = {}
        filtering_statistics_managers = {}

        for lib in self.libraries.keys():
            trim_processes_managers[lib] = self.libraries[lib].trimmomatic_and_low_complexity_filter_process()
            trim_processes[lib] = trim_processes_managers[lib].__enter__()

            filtering_statistics_managers[lib] = self.libraries[lib].filtering_statistics()
            filtering_statistics[lib] = filtering_statistics_managers[lib].__enter__()

        overall_filtering_statistics = defaultdict(int)

        # Paths for the 4 expected FastQs
        input_fastqs = []
        for r in [1, 2, 3, 4]:
            input_fastqs.append(self.parameters['fastq_path_base'] % (self.parameters['split_affixes'][self.split_file_index], r))

        last_ping = time.time()
        ping_every_n_reads = 1000000

        for r_name, seqs, quals in self._weave_fastqs(input_fastqs):

            # Python 3 compatibility in mind!
            seqs = [s.decode('utf-8') for s in seqs]

            keep, lib_index, result = self._v2_process_reads(r_name, seqs, quals,
                                                    error_corrected_barcodes, error_corrected_rev_compl_barcodes, 
                                                    self.sequence_to_index_mapping)
            if keep:
                bc, umi = result
                bio_read = seqs[0]
                bio_qual = quals[0]
                trim_processes[lib_index].write(to_fastq_lines(bc, umi, bio_read, bio_qual, r_name[1:]))
                # print_to_log(to_fastq_lines(bc, umi, bio_read, bio_qual, r_name))
                filtering_statistics[lib_index]['Valid'] += 1
                overall_filtering_statistics['Valid'] += 1

            else:
                if result != 'Invalid_library_index':
                    filtering_statistics[lib_index][result] += 1
                overall_filtering_statistics[result] += 1

            # Track speed per M reads
            overall_filtering_statistics['Total'] += 1
            
            
            if overall_filtering_statistics['Total']%ping_every_n_reads == 0:
                sec_per_mil = (time.time()-last_ping)/(float(ping_every_n_reads)/10**6)
                last_ping = time.time()
                print_to_log('%d reads parsed, kept %d reads, currently %.02f seconds/M reads.' % (overall_filtering_statistics['Total'], overall_filtering_statistics['Valid'], sec_per_mil))
                
                for k in ['Valid', 'Invalid_library_index', 'Invalid_BC1',  'Invalid_BC2', 'UMI_contains_N']:
                    print_to_log("     %.02fpct %s" % (100.*overall_filtering_statistics[k]/overall_filtering_statistics['Total'], k))
                print_to_log("     ----")
                for lib in self.libraries.keys():
                    print_to_log("     %.02fpct %s" % (100.*filtering_statistics[lib]['Valid']/overall_filtering_statistics['Valid'], lib))

        print_to_log('%d reads parsed, kept %d reads.' % (overall_filtering_statistics['Total'], overall_filtering_statistics['Valid']))

        # Close up the context managers
        for lib in self.libraries.keys():
            print_to_log(lib)
            print_to_log(filtering_statistics[lib])

            trim_processes_managers[lib].__exit__(None, None, None)
            filtering_statistics_managers[lib].__exit__(None, None, None)

        print_to_log('Closed up shop')


class IndropsAnalysis():

    def __init__(self, parameters, split_file_index=0):
        self.parameters = parameters
        self.user_paths = parameters['project_paths']
        self.user_paths.update(parameters['general_paths'])

        self.split_file_index = split_file_index

        # Add defaults for any newly added parameters.
        if 'min_non_polyA' not in self.parameters['umi_quantification_arguments']:
            self.parameters['umi_quantification_arguments']['min_non_polyA'] = 0

        if 'java' not in self.user_paths:
            self.user_paths['java'] = 'java'

        if 'python' not in self.user_paths:
            self.user_paths['python'] = 'python'

        self.output_paths = {
            'stats_dir' : 'stats/',
            'pre_split_dir' : 'pre_split/',
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


        # For the raw input, check if we have a single or multiple files
        # Update the output requirements accordingly
        if self.user_paths['raw_R1_fastq'] and self.user_paths['raw_R2_fastq']:

            self.output_paths['raw_file_pairs'] = [(self.user_paths['raw_R1_fastq'], self.user_paths['raw_R2_fastq'])]
            self.output_paths['filtered_fastq'] = [os.path.join(self.output_paths['pre_split_dir'], 'filtered.fastq')]
            self.output_paths['read_fail_counts'] = [os.path.join(self.output_paths['stats_dir'], 'filtering_metrics.png')]
            self.output_paths['barcode_read_counts'] = [p + '.counts.pickle' for p in self.output_paths['filtered_fastq']]
            self.output_paths['barcode_read_indices'] = [p + '.index.pickle' for p in self.output_paths['filtered_fastq']]
            

        if self.user_paths['split_raw_R1_prefix'] and self.user_paths['split_raw_R2_prefix'] and self.user_paths['split_suffixes']:
            suffixes = self.user_paths['split_suffixes']
            r1s = [self.user_paths['split_raw_R1_prefix']%suf for suf in suffixes]
            r2s = [self.user_paths['split_raw_R2_prefix']%suf for suf in suffixes]
            self.output_paths['raw_file_pairs'] = zip(r1s, r2s)

            self.build_filtered_reads_paths(suffixes)
        else:
            self.output_paths['raw_file_pairs'] = []

        if self.user_paths['split_suffixes']:
            self.build_filtered_reads_paths()

        self.output_paths['barcode_histogram'] = os.path.join(self.output_paths['stats_dir'], 'barcode_abundance_histogram.png')
        self.output_paths['barcode_sorted_reads'] = os.path.join(self.output_paths['split_dir'], 'barcoded_sorted.fastq')
        self.output_paths['barcode_sorted_reads_index'] = self.output_paths['barcode_sorted_reads'] + '.index.pickle'

    def build_filtered_reads_paths(self, suffixes):
        self.output_paths['filtered_fastq'] = [os.path.join(self.output_paths['pre_split_dir'], 'filtered_%s.fastq' % suf) for suf in suffixes]
        self.output_paths['read_fail_counts'] = [os.path.join(self.output_paths['stats_dir'], 'filtering_metrics_%s.png' % suf) for suf in suffixes]
        self.output_paths['barcode_read_counts'] = [p + '.counts.pickle' for p in self.output_paths['filtered_fastq']]
        self.output_paths['barcode_read_indices'] = [p + '.index.pickle' for p in self.output_paths['filtered_fastq']]

    @contextmanager
    def trimmomatic_and_low_complexity_filter_process(self):
        """
        We start 3 processes that are connected with Unix pipes.

        Process 1 - Trimmomatic. Doesn't support stdin/stdout, so we instead use named pipes (FIFOs). It reads from FIFO1, and writes to FIFO2. 
        Process 2 - In line complexity filter, a python script. It reads from FIFO2 (Trimmomatic output) and writes to the ouput file. 

        When these are done, we start another process to count the results on the FastQ file.
        """
        filtered_filename = self.output_paths['filtered_fastq'][self.split_file_index]
        with FIFO(dir=self.output_paths['pre_split_dir']) as fifo2:
            script_dir = os.path.dirname(os.path.realpath(__file__))
            low_complexity_filter_cmd = [self.user_paths['python'], os.path.join(script_dir, 'filter_low_complexity_reads.py'),
                '-input', fifo2.filename,
                '-output', filtered_filename,
                ]
            p2 = subprocess.Popen(low_complexity_filter_cmd)

            with FIFO(dir=self.output_paths['pre_split_dir']) as fifo1:

                trimmomatic_cmd = [self.user_paths['java'], '-Xmx500m', '-jar', self.user_paths['trimmomatic'], 'SE', '-threads', "1", '-phred33', fifo1.filename, fifo2.filename]
                for arg, val in self.parameters['trimmomatic_arguments'].items():
                    trimmomatic_cmd.append('%s:%s' % (arg, val))
                p1 = subprocess.Popen(trimmomatic_cmd)

            
                fifo1_filehandle = open(fifo1.filename, 'w')
                yield fifo1_filehandle
                fifo1_filehandle.write('\n') # THIS IS AWFUL, but for some reason, otherwise Trimmomatic won't detect that we are at EOF...!!!! This makes it crash.
                fifo1_filehandle.close()
                p1.wait()

            p2.wait()

        counter_cmd = [self.user_paths['python'], os.path.join(script_dir, 'count_barcode_distribution.py'), filtered_filename]
        p3 = subprocess.Popen(counter_cmd)
        p3.wait()

    @contextmanager
    def filtering_statistics(self):
        filtering_statistics_counter = defaultdict(int)
        yield filtering_statistics_counter

        filtering_statistics = {
            'Total Reads' : filtering_statistics_counter['Total'],
            'Valid Reads' : filtering_statistics_counter['Valid'],
            'Rejected Reads' : filtering_statistics_counter['Total'] - filtering_statistics_counter['Valid'],
            'Rejection Flags' : dict(filtering_statistics_counter)
        }

        failure_counts_filename = self.output_paths['read_fail_counts'][self.split_file_index]
        with open(failure_counts_filename, 'w') as f:
            yaml.dump(dict(filtering_statistics), f, default_flow_style=False)


    ########
    #
    # V1 Beads design
    #
    #######
    def v1_filter_and_count_reads(self):
        """
        Input the two raw FastQ files
        Output: 
            - A single fastQ file that uses the read name to store the barcoding information
            - A pickle of number of reads originating from each barcode 
        """
        # Relevant paths
        r1_filename, r2_filename = self.output_paths['raw_file_pairs'][self.split_file_index]


        #Get barcode neighborhoods
        bc1s = build_barcode_neighborhoods(self.user_paths['gel_barcode1_list'])
        bc2s = build_barcode_neighborhoods(self.user_paths['gel_barcode2_list'])


        # This starts a Trimmoatic process, a low complexity filter process, and will 
        # upon closing, start the barcode distribution counting process.
        start_time = time.time()
        with self.trimmomatic_and_low_complexity_filter_process() as trim_process, self.filtering_statistics() as filtering_statistics_counter :

            #Iterate over the weaved reads
            for r_name, r1_seq, r1_qual, r2_seq, r2_qual in self._v1_weave_fastqs(r1_filename, r2_filename):
                    
                # Check if they should be kept
                keep, result = self._v1_process_reads(r1_seq, r2_seq, valid_bc1s=bc1s, valid_bc2s=bc2s)

                # Write the the reads worth keeping
                if keep:
                    bc, umi = result
                    trim_process.write(to_fastq_lines(bc, umi, r2_seq, r2_qual, r_name))

                    filtering_statistics_counter['Valid'] += 1
                else:
                    filtering_statistics_counter[result] += 1

                # Track speed per M reads
                filtering_statistics_counter['Total'] += 1
                if filtering_statistics_counter['Total']%1000000 == 0:
                    sec_per_mil = (time.time()-start_time)/(float(filtering_statistics_counter['Total'])/10**6)
                    print_to_log('%d reads parsed, kept %d reads, %.02f seconds per M reads.' % (filtering_statistics_counter['Total'], filtering_statistics_counter['Valid'], sec_per_mil))

        print_to_log('%d reads parsed, kept %d reads.' % (filtering_statistics_counter['Total'], filtering_statistics_counter['Valid']))

    def _v1_weave_fastqs(self, r1_fastq, r2_fastq):
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

    def _v1_process_reads(self, name, read, valid_bc1s={}, valid_bc2s={}):
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
#             Nstart=0
#             if name.startswith('N')
#               Nstart=1
            if not 7 < w1_pos < 12:
                # print_to_log(name)
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
        #was:expected_poly_t = name[44:44+minimal_polyT_len_on_R1:]
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

    ########
    #
    # Design independent analysis
    #
    #######

    def get_barcode_read_counts_from_pickle(self):
        barcode_read_counter = defaultdict(int)
        for filename in self.output_paths['barcode_read_counts']:
            with open(filename, 'r') as f:
                partial_counts = pickle.load(f)
                for bc, c in partial_counts.items():
                    barcode_read_counter[bc] += c
        return barcode_read_counter


    def make_barcode_abundance_histogram(self):
        """
        Takes the read-count-by-barcode pickle and outputs a histogram used 
        to determine a treshold on the minimal number of reads coming from good barcodes
        """
        barcode_read_counter = self.get_barcode_read_counts_from_pickle()

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
        barcode_read_counter = self.get_barcode_read_counts_from_pickle()

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

        for filtered_fn in self.output_paths['filtered_fastq']:
            print_to_log('Now on file %s' % filtered_fn)
            with open(filtered_fn, 'r') as input_fastq:
                for name, seq, qual in from_fastq(input_fastq):
                    pre_write_buffer_size += 1
                    total_processed_reads += 1
                    bc = name.split(':')[0]
                    
                    if bc in barcode_names:
                        bc_name = barcode_names[bc]
                        filename = os.path.join(self.output_paths['split_fastq_dir'], '%s.fastq' % bc_name)
                        pre_write[filename].append(to_fastq(name, seq, qual))

                    elif output_unassigned_reads:
                        pre_write[unassigned_filename].append(to_fastq(name, seq, qual))


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

        # Merge individual files, keeping an index.
        barcode_read_counter = self.get_barcode_read_counts_from_pickle()

        sorted_output_index = {}
        sorted_output = open(self.output_paths['barcode_sorted_reads'], 'w')
        for bc, bc_name in sorted(barcode_names.items(), key=lambda i: i[1]):
            filename = os.path.join(self.output_paths['split_fastq_dir'], '%s.fastq' % bc_name)

            sorted_output_index[bc_name] = (bc, sorted_output.tell(), 4 * barcode_read_counter[bc])
            # print_to_log(str((bc, sorted_output.tell(), 4 * barcode_read_counter[bc])))
            with open(filename, 'r') as single_barcode:
                shutil.copyfileobj(single_barcode, sorted_output)
            os.remove(filename)

        if output_unassigned_reads:
            with open(unassigned_filename) as f:
                shutil.copyfileobj(f, sorted_output)
            os.remove(unassigned_filename)

        sorted_output.close()
        with open(self.output_paths['barcode_sorted_reads_index'], 'w') as f:
            pickle.dump(sorted_output_index, f)

    def get_reads_for_barcode(self, barcode):
        if not hasattr(self, '_barcode_sorted_reads_index'):
            with open(self.output_paths['barcode_sorted_reads_index'], 'r') as f:
                self._barcode_sorted_reads_index = pickle.load(f)

        original_barcode, byte_offset, total_lines = self._barcode_sorted_reads_index[barcode]
        with open(self.output_paths['barcode_sorted_reads'], 'r') as barcode_sorted_reads:
            barcode_sorted_reads.seek(byte_offset)
            for i in range(total_lines):
                yield next(barcode_sorted_reads)


    @parallelized_using_workers
    def quantify_expression_for_barcode(self, barcode):

        counts_output = os.path.join(self.output_paths['split_quant_dir'], '%s.counts' % barcode)
        quant_metrics_output = os.path.join(self.output_paths['split_quant_dir'], '%s.quant_metrics' % barcode)
        oversequencing_metrics_output = os.path.join(self.output_paths['split_quant_dir'], '%s.oversequencing' % barcode)
        output_umifm_calculation_metrics = os.path.join(self.output_paths['split_quant_dir'], '%s.umifm_stats' % barcode)
        unaligned_reads_output = os.path.join(self.output_paths['split_quant_dir'], '%s.unaligned.fastq' % barcode)
        aligned_bam_output = os.path.join(self.output_paths['split_quant_dir'], '%s.aligned.bam' % barcode)

        # Bowtie command

        bowtie_exec = os.path.join(self.user_paths['bowtie'], 'bowtie')
        bowtie_cmd = [bowtie_exec, self.user_paths['bowtie_index_prefix'], '-q', '-',
            '-p', '1', '-a', '--best', '--strata', '--chunkmbs', '1000', '--norc', '--sam',
            '-shmem', #should sometimes reduce memory usage...?
            '-m', str(self.parameters['bowtie_arguments']['m']),
            '-n', str(self.parameters['bowtie_arguments']['n']),
            '-l', str(self.parameters['bowtie_arguments']['l']),
            '-e', str(self.parameters['bowtie_arguments']['e']),
            ]
        if self.parameters['output_arguments']['output_unaligned_reads_to_other_fastq']:
            bowtie_cmd += ['--un', unaligned_reads_output]

        
        # Quantification command

        script_dir = os.path.dirname(os.path.realpath(__file__))
        quant_cmd = [self.user_paths['python'], os.path.join(script_dir, 'filter_alignments.py'),
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
        if self.parameters['output_arguments']['low_complexity_mask']:
            low_complexity_mask_pickle_path = self.user_paths['bowtie_index_prefix'] + '.low_complexity_regions.pickle'
            quant_cmd += ['--low_complexity_mask', low_complexity_mask_pickle_path]

        if self.parameters['output_arguments']['output_alignment_to_bam']:
            p2_output = open(aligned_bam_output, 'w')
        else:
            p2_output = open(os.devnull, 'w')


        # Spawn processes
        p1 = subprocess.Popen(bowtie_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(quant_cmd, stdin=p1.stdout, stdout=p2_output)

        for line in self.get_reads_for_barcode(barcode):
            p1.stdin.write(line)
        p1.stdin.close()
        p2.wait() #Wait to finish before moving on


    def aggregate_counts(self, process_ambiguity_data=False, minimal_counts=0):
        with open(self.output_paths['good_barcodes_with_names'], 'r') as f:
            sorted_barcode_names = sorted(pickle.load(f).values())


        existing_barcodes = os.listdir(self.output_paths['split_quant_dir'])
        missing_barcodes = []
        for barcode in sorted_barcode_names:
            if not ('%s.counts' % barcode) in existing_barcodes:
                missing_barcodes.append(barcode)
                continue

        print_to_log('Missing the following barcodes: '+ ','.join(missing_barcodes))
        print_to_log('Corresponding indices: ')
        print_to_log(','.join([str(int(bc[2:])-1) for bc in missing_barcodes]))

        missing_barcodes = set(missing_barcodes)
        i = 0

        output_filename = os.path.join(self.output_paths['aggregated_counts'], 'full_counts.txt')
        output_file = open(output_filename, 'w')

        i = 0
        gene_names = []

        all_counts = None
        for i, barcode in enumerate(sorted_barcode_names):
            if (i % 100)== 0:
                print_to_log('Read %d barcodes (of %d).' % (i, len(sorted_barcode_names)))
            if barcode in missing_barcodes:
                continue

            counts_filename = os.path.join(self.output_paths['split_quant_dir'], '%s.counts' % barcode)
            with open(counts_filename, 'r') as f:
                try:
                    header = next(f)
                except StopIteration:
                    print_to_log(str(int(barcode[2:])-1))
                    continue

                if i == 0:
                    first_file_counts = []
                    for line in f:
                        row = line.rstrip('\n').split('\t')
                        first_file_counts.append(float(row[1]))
                        gene_names.append(row[0])

                    all_counts = np.zeros((len(gene_names), len(sorted_barcode_names)))
                    all_counts[:, 0] = np.array(first_file_counts)
                    is_first_file = False

                else:
                    positive_counts_mask = []
                    file_counts = []
                    for j, line in enumerate(f):
                        row = line.rstrip('\n').split('\t')
                        count = float(row[1])
                        if count > 0:
                            file_counts.append(count)
                            positive_counts_mask.append(j)
                    all_counts[positive_counts_mask, i] = np.array(file_counts)



        barcode_name_array = np.array(sorted_barcode_names)
        barcodes_above_threshold = all_counts.sum(axis=0) > minimal_counts
        gene_total_counts = all_counts.sum(axis=1)
        print_to_log('Output only %d barcodes above threshold.' % sum(barcodes_above_threshold))


        to_output_line = lambda row: ('\t'.join([str(r) for r in row]))
        output_header = ['gene'] + ['Sum_counts', 'Sum_ambig', 'Ambigs'] + list(barcode_name_array[barcodes_above_threshold])
        output_file.write(to_output_line(output_header))
        for i in range(all_counts.shape[0]):
            line = [gene_names[i], gene_total_counts[i], '0.0', ''] + list(all_counts[i, barcodes_above_threshold])
            output_file.write('\n' + to_output_line(line))
            if (i % 1000)== 0:
                print_to_log('Wrote %d genes.' % i)

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


    parser_v2preprocess = subparsers.add_parser('v2_preprocess')
    parser_v2preprocess.add_argument('parameters', type=argparse.FileType('r'), help='Sequencing run YAML File.')
    parser_v2preprocess.add_argument('--split-file-index', type=int, help='Index for split files.', default=0)

    parser_preprocess = subparsers.add_parser('preprocess')
    parser_preprocess.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_preprocess.add_argument('--split-file-index', type=int, help='Split file to preprocess, leave as 0 if there is a single file.', default=0)

    parser_histogram = subparsers.add_parser('histogram')
    parser_histogram.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')

    parser_split = subparsers.add_parser('split_barcodes')
    parser_split.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_split.add_argument('--read-count-threshold', type=int,
        help="Minimal read count to be considered 'abundant' barcode")

    parser_quantify = subparsers.add_parser('quantify')
    parser_quantify.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_quantify.add_argument('--barcodes-per-worker', type=int, help='Barcodes to be processed by each worker.', default=0)
    parser_quantify.add_argument('--worker-index', type=int, help='Index of current worker. (Starting at 0). Make sure max(worker-index)*(barcodes-per-worker) > total barcodes.', default=0)
    parser_quantify.add_argument('--total-workers', type=int, help='Total workers that are working together. This takes precedence over barcodes-per-worker.', default=0)
    parser_quantify.add_argument('--missing', type=str, help='Pattern for file name to search', default='')

    parser_aggregate = subparsers.add_parser('aggregate')
    parser_aggregate.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_aggregate.add_argument('--minimal-counts', type=int, help='Minimal counts for barcodes that get written to output.', default=0)

    parser_sort_bam = subparsers.add_parser('sort_bam')
    parser_sort_bam.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_sort_bam.add_argument('--barcodes-per-worker', type=int, help='Barcodes to be processed by each worker.', default=0)
    parser_sort_bam.add_argument('--worker-index', type=int, help='Index of current worker. (Starting at 0). Make sure max(worker-index)*(barcodes-per-worker) > total barcodes.', default=0)

    parser_get_reads = subparsers.add_parser('get_reads')
    parser_get_reads.add_argument('parameters', type=argparse.FileType('r'), help='Project Parameters YAML File.')
    parser_get_reads.add_argument('barcode', type=str, help='Barcode ID. (i.e. bc0020)')
    parser_get_reads.add_argument('output', type=argparse.FileType('w'), nargs='?', default=sys.stdout)

    args = parser.parse_args()

    if args.command == 'index':
        build_transcriptome_from_ENSEMBL_files(args.transcriptome,
            args.index_prefix,
            gtf_filename=args.gtf_gz, 
            polyA_tail_length=args.polyA,
            bowtie_path=args.bowtie_path)

    else:
        
        if args.command == 'v2_preprocess':
            demultiplexer = v2Demultiplexer(args.parameters, args.split_file_index)
            demultiplexer.v2_filter_and_count_reads()

        elif args.command in ['preprocess', 'histogram', 'split_barcodes', 'quantify', 'aggregate', 'sort_bam', 'get_reads']:
            parameters = yaml.load(args.parameters)
            analysis = IndropsAnalysis(parameters, args.split_file_index if hasattr(args, 'split_file_index') else 0)
            if args.command == 'preprocess':
                analysis.v1_filter_and_count_reads()
                if len(analysis.output_paths['raw_file_pairs']) == 1:
                    analysis.make_barcode_abundance_histogram()

            elif args.command == 'histogram':
                analysis.make_barcode_abundance_histogram()

            elif args.command == 'split_barcodes':
                analysis.choose_good_barcodes(args.read_count_threshold)
                analysis.split_reads_by_barcode()

            elif args.command == 'quantify':
                analysis.quantify_expression_for_barcode(args.barcodes_per_worker, args.worker_index, total_workers=args.total_workers, missing=args.missing)
            
            elif args.command == 'aggregate':
                analysis.aggregate_counts(minimal_counts=args.minimal_counts)
            
            elif args.command == 'sort_bam':
                analysis.sort_and_index_bam(args.barcodes_per_worker, args.worker_index)
            
            elif args.command == 'get_reads':
                for line in analysis.get_reads_for_barcode(args.barcode):
                    args.output.write(line)

