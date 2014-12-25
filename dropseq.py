import os, subprocess
import itertools
import operator
from collections import defaultdict
import cPickle as pickle

import numpy as np
import re

from itertools import product, combinations
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

def from_fastq(handle):
    while True:
        name = next(handle).rstrip()[1:] #Read name
        seq = next(handle).rstrip() #Read seq
        next(handle) #+ line
        qual = next(handle).rstrip() #Read qual
        if not name or not seq or not qual:
            break
        yield name, seq, qual


def filter_and_count_reads(paths):
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

def barcode_histogram(paths):
    """
    Takes the read-count-by-barcode pickle and outputs a histogram used 
    to determine a treshold on the minimal number of reads coming from good barcodes
    """
    with open(paths['barcode_read_counts'], 'r') as f:
        barcode_read_counter = pickle.load(f)

    count_freq = defaultdict(int)
    for bc, count in barcode_read_counter.items():
        count_freq[count] += 1

    x = np.array(count_freq.keys())
    y = np.array(count_freq.values())
    w = x*y

    import matplotlib.pyplot as plt
    ax = plt.subplot(111)
    ax.hist(x, bins=np.logspace(0, 6, 50), weights=w)
    ax.set_xscale('log')
    ax.set_xlabel('Reads per barcode')
    ax.set_ylabel('#reads coming from bin')
    plt.savefig(paths['barcode_histogram'])

def choose_good_barcodes(paths, threshold=15000):
    """
    Takes the read-count-by-barcode pickle and a minimal threshold value, 
    Outputs a list of barcodes to the retained and assigns a 'bcN' name to each barcode
    """
    with open(paths['barcode_read_counts'], 'r') as f:
        barcode_read_counter = pickle.load(f)

    good_barcodes = []
    for bc, count in barcode_read_counter.items():
        if count >= threshold:
            good_barcodes.append(bc)

    print(len(good_barcodes))
    barcode_names = {}
    for i, bc in enumerate(sorted(good_barcodes, key=lambda b: barcode_read_counter[bc], reverse=True)):
        barcode_names[bc] = 'bc%d' % (i+1)

    with open(paths['good_barcodes_with_names'], 'w') as f:
        pickle.dump(barcode_names, f)

def split_reads_by_barcode(paths):
    """
    Starts with the list of good barcodes and the filtered FastQ
    Splits the filtered FastQ into one FastQ file per read. 
    """
    with open(paths['good_barcodes_with_names'], 'r') as f:
        barcode_names = pickle.load(f)

    i = 0
    j = 0
    pre_write = defaultdict(list)

    with open(paths['filtered_fastq'], 'r') as input_fastq:
        for name, seq, qual in from_fastq(input_fastq):
            i += 1
            bc, umi = name.split(':')
            if bc in barcode_names:
                j += 1
                bc_name = barcode_names[bc]
                filename = os.path.join(paths['split_barcodes_dir'], '%s.fastq' % bc_name)
                pre_write[filename].append(to_fastq_lines(bc, umi, seq, qual))
                if j % 1000000 == 0:
                    for fn, chunks in pre_write.items():
                        with open(fn, 'a') as out:
                            for chunk in chunks:
                                out.write(chunk)
                    j = 0
                    pre_write = defaultdict(list)

    for fn, chunks in pre_write.items():
        with open(fn, 'a') as out:
            for chunk in chunks:
                out.write(chunk)

def prepare_transcriptome_index():
    in_genes = '/Users/averes/Projects/Melton/temp_dropseq/genes.gtf'
    out_genes = '/Users/averes/Projects/Melton/temp_dropseq/genes.annotated.gtf'
    with open(in_genes, 'r') as in_f, open(out_genes, 'w') as out_f:
        for line in in_f:
            chr_name = line.rstrip().split('\t')[0]
            if '_hap' in chr_name or 'Un_gl' in chr_name or '_' in chr_name:
                continue
            gene_name = re.search(r'gene_id \"(.*?)\";', line).group(1)
            out_line = re.sub(r'(?<=transcript_id ")(.*?)(?=";)', r'\1|'+gene_name, line)
            out_f.write(out_line)

if __name__=="__main__":
    
    #Change to relevant directory
    base_dir = '/n/regal/melton_lab/adrianveres/datasets/S6D13_cells/data/'

    #Where you have the two barcode lists
    barcode_dir = '/n/beta_cell/Users/adrianveres/dropseq_data/' 
    
    paths = {
        'r1_input': os.path.join(base_dir, 'S6D13-100_S0.R1.fastq.gz'),
        'r2_input': os.path.join(base_dir, 'S6D13-100_S0.R2.fastq.gz'),
        'input_filetype': 'gz', #Either 'gz' (so the script deals with compression) or 'fq'/'fastq' so it doesn't
        'barcode_read_counts': os.path.join(base_dir, 'stats', 'barcode_read_counts.pickle'), #Temp file
        'good_barcodes_with_names': os.path.join(base_dir, 'stats', 'good_barcodes_with_names.pickle'), #Temp file
        'filtered_fastq': os.path.join(base_dir, 'S6D13-100.filtered.fastq'), #Fastq file after removal of bad reads, but before split
        'barcode_histogram': os.path.join(base_dir, 'reads_from_barcodes.png')
        'split_barcodes_dir': os.path.join(base_dir, 'barcodes'), #Directory where individual barcode fastqs will be placed
        'bc1s': os.path.join(barcode_dir, 'gel_barcode1_list.txt'),
        'bc2s': os.path.join(barcode_dir, 'gel_barcode2_list.txt'),
    

    #See the functiond description for these steps (I suggest running them one at a time)
    import sys
    if len(sys.argv)>1:
        if sys.argv[1] == 'filter':
            filter_and_count_reads(paths)
        elif sys.argv[1] == 'histo':
            barcode_histogram(paths) #Inspect this histogram to set the threshold below
        elif sys.argv[1] == 'choose_barcodes':
            choose_good_barcodes(paths, int(sys.argv[2]))
        elif sys.argv[1] == 'split':
            split_reads_by_barcode(paths)

