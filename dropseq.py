import os, subprocess
import itertools
import operator
from collections import defaultdict
import cPickle as pickle

from itertools import product, combinations
from future.builtins import dict #So dict.keys/values/items are memory efficient
import time

def string_hamming_distance(str1, str2):
    """
    Fast hamming distance over 2 strings known to be of same length.
    """
    return sum(itertools.imap(operator.ne, str1, str2))

def rev_comp(seq):
    tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(tbl[s] for s in seq[::-1])

def seq_neighborhood(seq, n_subs=1):
    """
    Given a sequence, yield all sequences within n_subs substitutions of that sequence.
    """
    for positions in combinations(range(len(seq)), n_subs):
        for subs in product(*("ATGCN",)*n_subs):
            seq_copy = list(seq)
            for p, s in zip(positions, subs):
                seq_copy[p] = s
            yield ''.join(seq_copy)
            
def build_neighborhoods(barcode_file):
    """
    Given a set of barcodes, produce sequences which can unambiguously be mapped to these barcodes,
    within 2 substitutions. If a sequence maps to multiple barcodes, get rid of it.
    However, if a sequences maps to a bc1 with 1change and another with 2changes, keep the 1change mapping.
    """
    
    clean_mapping = dict()
    mapping1 = defaultdict(set)
    mapping2 = defaultdict(set)
    
    #Build the full neighborhood
    with open(barcode_file, 'rU') as f:
        for line in f:
            barcode = rev_comp(line.rstrip())
            clean_mapping[barcode] = barcode
            
            for n in seq_neighborhood(barcode, 1):
                mapping1[n].add(barcode)
            for n in seq_neighborhood(barcode, 2):
                mapping2[n].add(barcode)   
            
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

def process_reads(name, read, polyT_len=7, hamming_threshold=3,
                      valid_bc1s={}, valid_bc2s={}):
    """
    Returns either:
        True, (barcode, umi)
            (if read passes filter)
        False, name of filter that failed
            (for stats collection)
    """
    w1 = "GAGTGATTGCTTGTGACGCCTT"
    rev_w1 = "AAGGCGTCACAAGCAATCACTC"
    polyA = "A"*8
    
    if rev_w1 in read:
        return False, 'W1_in_R2'

    if polyA in read:
        return False, 'PolyA_in_R2'

    #Check for polyT signal at 3' end.
    poly_t = name[-polyT_len:]
    if string_hamming_distance(poly_t, 'T'*polyT_len) > 3:        
        return False, 'No_polyT'
    
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
            if string_hamming_distance(w1, name[w1_pos:w1_pos+22]) <= hamming_threshold:
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
            bc1 = valid_bc1s[bc1]
        else:
            return False, 'BC1'
        if bc2 in valid_bc2s:
            bc2 = valid_bc2s[bc2]
        else:
            return False, 'BC2'
    bc = '%s-%s'%(bc1, bc2)
    return True, (bc, umi)

def weave_fastqs(r1_fastq, r2_fastq, filetype='gz'):
    """
    Merge 2 FastQ files by returning paired reads for each.
    Returns only R1_seq, R2_seq and R2_qual.
    """
    # Decompress Gzips using subprocesses because python gzip is incredibly slow.
    if filetype == 'gz':    
        r1_gunzip = subprocess.Popen("gzip --stdout -d %s" % (r1_fastq), shell=True, stdout=subprocess.PIPE)
        r1_stream = r1_gunzip.stdout
        r2_gunzip = subprocess.Popen("gzip --stdout -d %s" % (r2_fastq), shell=True, stdout=subprocess.PIPE)
        r2_stream = r2_gunzip.stdout
    elif filetype == 'fq' or filetype == 'fastq':
        r1_stream = open(r1_fastq, 'r')
        r2_stream = open(r2_fastq, 'r')


    while True:
        #Read 4 lines from each FastQ
        next(r1_stream) #Read name
        r1_seq = next(r1_stream).rstrip() #Read seq
        next(r1_stream) #+ line
        next(r1_stream) #Read qual
        
        next(r2_stream) #Read name
        r2_seq = next(r2_stream).rstrip() #Read seq
        next(r2_stream) #+ line
        r2_qual = next(r2_stream).rstrip() #Read qual
        
        if not r1_seq or not r2_seq:
            break
        yield r1_seq, r2_seq, r2_qual
    r1_stream.stdout.close()
    r2_stream.stdout.close()

def to_fastq_lines(bc, umi, seq, qual):
    """
    Return string that can be written to fastQ file
    """
    return '@'+bc+':'+umi+'\n'+seq+'\n+\n'+qual+'\n'

def run(paths, keep_barcodes=[]):
    #Prepare data collection
    barcode_read_counter = defaultdict(int)
    filter_fail_counter = defaultdict(int)
    kept_reads = 0

    #Get barcode neighborhoods
    bc1s = build_neighborhoods(paths['bc1s'])
    bc2s = build_neighborhoods(paths['bc2s'])

    i = 0
    start_time = time.time()
    
    with open(paths['filtered_fastq'], 'w') as output_fastq:
        for r1_seq, r2_seq, r2_qual in weave_fastqs(paths['r1_input'], paths['r2_input'], paths['input_filetype']):
            keep, result = process_reads(r1_seq, r2_seq, valid_bc1s=bc1s, valid_bc2s=bc2s)

            i += 1
            if i%1000000 == 0:
                sec_per_mil = (time.time()-start_time)/(float(i)/10**6)
                print('%d reads parsed, kept %d reads, %.02f seconds per M reads.' % (i, kept_reads, sec_per_mil))

            if keep:
                bc, umi = result
                kept_reads += 1
                barcode_read_counter[bc] += 1
                output_fastq.write(to_fastq_lines(bc, umi, r2_seq, r2_qual))
            else:
                filter_fail_counter[result]

    #Do something with results
    with open(paths['barcode_read_counts'], 'w') as f:
        pickle.dump(dict(barcode_read_counter), f)


def run_multithreaded(paths):
    import multiprocessing, multiprocessing.queues
    from multiprocessing import Process, current_process
    import string
    import cPickle as pickle

    class SharedCounter(object):
        """ A synchronized shared counter.

        The locking done by multiprocessing.Value ensures that only a single
        process or thread may read or write the in-memory ctypes object. However,
        in order to do n += 1, Python performs a read followed by a write, so a
        second process may read the old value before the new one is written by the
        first process. The solution is to use a multiprocessing.Lock to guarantee
        the atomicity of the modifications to Value.

        This class comes almost entirely from Eli Bendersky's blog:
        http://eli.thegreenplace.net/2012/01/04/shared-counter-with-pythons-multiprocessing/

        """

        def __init__(self, n = 0):
            self.count = multiprocessing.Value('i', n)

        def increment(self, n = 1):
            """ Increment the counter by n (default = 1) """
            with self.count.get_lock():
                self.count.value += n

        @property
        def value(self):
            """ Return the value of the counter """
            return self.count.value


    class Queue(multiprocessing.queues.Queue):
        """ A portable implementation of multiprocessing.Queue.


        Because of multithreading / multiprocessing semantics, Queue.qsize() may
        raise the NotImplementedError exception on Unix platforms like Mac OS X
        where sem_getvalue() is not implemented. This subclass addresses this
        problem by using a synchronized shared counter (initialized to zero) and
        increasing / decreasing its value every time the put() and get() methods
        are called, respectively. This not only prevents NotImplementedError from
        being raised, but also allows us to implement a reliable version of both
        qsize() and empty().

        Taken wholesale from https://github.com/vterron/lemon/
        """

        def __init__(self, *args, **kwargs):
            super(Queue, self).__init__(*args, **kwargs)
            self.size = SharedCounter(0)

        def put(self, *args, **kwargs):
            self.size.increment(1)
            super(Queue, self).put(*args, **kwargs)

        def get(self, *args, **kwargs):
            self.size.increment(-1)
            return super(Queue, self).get(*args, **kwargs)

        def qsize(self):
            """ Reliable implementation of multiprocessing.Queue.qsize() """
            return self.size.value

        def empty(self):
            """ Reliable implementation of multiprocessing.Queue.empty() """
            return not self.qsize()

    #Get barcode neighborhoods
    bc1s = build_neighborhoods(paths['bc1s'])
    bc2s = build_neighborhoods(paths['bc2s'])

    #Setup queues
    read_queue = Queue(2**15-1)
    stats_queue = Queue(100)


    def worker_task(read_queue, stats_queue, paths, bc1s, bc2s):

        #Prepare stats collection
        barcode_read_counter = defaultdict(int)
        filter_fail_counter = defaultdict(int)
        total_reads = 0
        kept_reads = 0

        worker_output_filename = paths['filtered_fastq'] + '.part_' + current_process().name
        worker_stats_filename = paths['filtered_fastq'] + '.part_' + current_process().name + '.stats'

        with open(worker_output_filename, 'w') as output_fastq:
            job = read_queue.get()
            while job is not None:
                (r1_seq, r2_seq, r2_qual) = job
                total_reads += 1
                keep, result = process_reads(r1_seq, r2_seq, valid_bc1s=bc1s, valid_bc2s=bc2s)
                if keep:
                    bc, umi = result
                    output_fastq.write(to_fastq_lines(bc, umi, r2_seq, r2_qual))

                    kept_reads += 1
                    barcode_read_counter[bc] += 1
                else:
                    filter_fail_counter[result] += 1
                job = read_queue.get()
            #Job was none, add another None to the queue so other workers can pick it up
            read_queue.put(None)

        with open(worker_stats_filename, 'w') as output_stats:
            stats = {
                'total_reads': total_reads,
                'kept_reads': kept_reads,
                'barcode_read_count': dict(barcode_read_counter),
                'filter_fail_count': dict(filter_fail_counter),
            }
            pickle.dump(stats, output_stats)

    def generator_task(read_queue, paths):
        #Keep pumping reads in the read queue...
        i = 0
        start_time = time.time()
        for r1_seq, r2_seq, r2_qual in weave_fastqs(paths['r1_input'], paths['r2_input'], paths['input_filetype']):
            i += 1
            read_queue.put((r1_seq, r2_seq, r2_qual))
            if i%10000 == 0:
                sec_per_mil = (time.time()-start_time)/(float(i)/10**6)
                print('%d reads parsed, %d reads in queue, %.02f' % (i, read_queue.qsize(), sec_per_mil))
        #Send the queue sentinel message
        read_queue.put(None)

    process_list = []
    process_list.append(Process(target=generator_task, name='manager', 
                        args=(read_queue, paths)))

    worker_count = 4
    worker_names = ["%02d"%i for i in range(1, worker_count+1)]
    for name in worker_names:
        process_list.append(Process(target=worker_task, name=name,
                            args=(read_queue, stats_queue, paths, bc1s, bc2s)))

    for proc in process_list:
        proc.start()

    for proc in process_list:
        proc.join()
        print "worker", proc.name, "finished!"
    print('well...')

    #Prepare for stats aggregation
    barcode_read_counter = defaultdict(int)
    filter_fail_counter = defaultdict(int)
    total_reads = 0
    kept_reads = 0

    for name in worker_names:
        with open(paths['filtered_fastq'] + '.part_' + name + '.stats', 'r') as f:
            stats = pickle.load(f)
        total_reads += stats['total_reads']
        kept_reads += stats['kept_reads']

        for k, v in stats['barcode_read_count'].items():
            barcode_read_counter[k] += v

        for k, v in stats['filter_fail_count'].items():
            filter_fail_counter[k] += v

    print(total_reads, kept_reads, len(barcode_read_counter))




if __name__=="__main__":
    import sys

    if len(sys.argv)==1 or sys.argv[1] == 'local':

        base_dir = '/Users/averes/Projects/Melton/temp_dropseq/'
        barcode_dir = '/Users/averes/Projects/Melton/temp_dropseq/'
        paths = {
            'r1_input': os.path.join(base_dir, 'test.R1.fastq.aa.gz'),
            'r2_input': os.path.join(base_dir, 'test.R2.fastq.aa.gz'),
            'input_filetype': 'gz',
            'barcode_read_counts': os.path.join(base_dir, 'stats', 'barcode_read_counts.pickle'),
            'filtered_fastq': os.path.join(base_dir, 'test.filtered.fastq'),
            'bc1s': os.path.join(barcode_dir, 'gel_barcode1_list.txt'),
            'bc2s': os.path.join(barcode_dir, 'gel_barcode2_list.txt'),
        }

    elif sys.argv[1] == 's6d13':
        base_dir = '/n/regal/melton_lab/adrianveres/datasets/S6D13_cells/data/'
        barcode_dir = '/n/beta_cell/Users/adrianveres/dropseq_data/'
        paths = {
            'r1_input': os.path.join(base_dir, 'S6D13-100_S0.R1.fastq.gz'),
            'r2_input': os.path.join(base_dir, 'S6D13-100_S0.R2.fastq.gz'),
            'input_filetype': 'gz',
            'barcode_read_counts': os.path.join(base_dir, 'stats', 'barcode_read_counts.pickle'),
            'filtered_fastq': os.path.join(base_dir, 'S6D13-100.filtered.fastq'),
            'bc1s': os.path.join(barcode_dir, 'gel_barcode1_list.txt'),
            'bc2s': os.path.join(barcode_dir, 'gel_barcode2_list.txt'),
        }

    # paths = {
    #     'r1_input': os.path.join(base_dir, 'S6D13-100.R1.fastq.aa.gz'),
    #     'r2_input': os.path.join(base_dir, 'S6D13-100.R2.fastq.aa.gz'),
    #     'input_filetype': 'gz',
    #     'filtered_fastq': os.path.join(base_dir, 'test.filtered.fastq'),
    #     'bc1s': os.path.join(barcode_dir, 'gel_barcode1_list.txt'),
    #     'bc2s': os.path.join(barcode_dir, 'gel_barcode2_list.txt'),
    # }

    run(paths)
    # run_multiprocess(paths)

