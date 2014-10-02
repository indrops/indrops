import subprocess
import os

hg19 = '/Users/averes/Projects/Melton/hg19/hg19.transcripts.annotated'
mm10 = '/Users/averes/Projects/Melton/mm10/mm10_refseq_annotated_rsem'
mixed = '/Users/averes/Projects/Melton/hg19_mm10/mm10_hg19_mixed'

# with open('/Users/averes/Projects/Melton/temp_seq_data/good_barcodes.txt') as f:
#     for i,line in enumerate(f):
#         print('Running %d' % i)
#         bc = line.rstrip()
#         input_filename = '/Users/averes/Projects/Melton/temp_seq_data/by_barcode_unique/%s.fastq' % bc
#         output_filename = '/Users/averes/Projects/Melton/temp_seq_data/by_barcode_counts/%s.txt' % bc

#         p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -k 200 --sam' % (ref, input_filename), stdout=subprocess.PIPE, shell=True)
#         p2 = subprocess.Popen('/Users/averes/miniconda3/envs/py27/bin/python alignments_to_counts.py --gtf genes.gtf --out %s' % output_filename,  stdin=p1.stdout, shell=True)
#         p1.stdout.close()
#         p2.wait()


# RUN FOR TEST FILES


# input_filename = '/Users/averes/Projects/Melton/temp_dropseq/killme.fq'

input_filename = '/Users/averes/Projects/Melton/temp_dropseq/reads_369_299.fq'
counts_filename = '/Users/averes/Projects/Melton/temp_dropseq/reads_369_299.mixed.counts'
raw_sam_filename = '/Users/averes/Projects/Melton/temp_dropseq/reads_369_299.raw.sam'
ref = mm10

input_filename = '/Users/averes/Projects/Melton/temp_dropseq/bugtest.fq'
counts_filename = '/Users/averes/Projects/Melton/temp_dropseq/bugtest.mixed.counts'
raw_sam_filename = '/Users/averes/Projects/Melton/temp_dropseq/reads_369_299.raw.sam'
ref = mm10

# input_filename = '/Users/averes/Projects/Melton/temp_dropseq/bc1.trimmed.fastq'
# counts_filename = '/Users/averes/Projects/Melton/temp_dropseq/bc1.trimmed.counts'
# raw_sam_filename = '/Users/averes/Projects/Melton/temp_dropseq/bc1.trimmed.sam'
# ref = hg19
# input_filename = '/Users/averes/Projects/Melton/temp_dropseq/head_to_head/reads_100_157.fq'
# counts_filename = '/Users/averes/Projects/Melton/temp_dropseq/hth.mixed.counts'


# input_filename = '/Users/averes/Projects/Melton/temp_dropseq/av1.fastq'
# raw_sam_filename = '/Users/averes/Projects/Melton/temp_dropseq/killme2.raw.sam'
# counts_filename = '/Users/averes/Projects/Melton/temp_dropseq/av1.mixed.counts'
# raw_bam_filename = '/Users/averes/Projects/Melton/temp_dropseq/av1.mixed.raw.bam'


# input_filename = '/Users/averes/Projects/Melton/temp_dropseq/av1.fastq'
# p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -a --best --strata --sam --norc' % (mixed, input_filename), stdout=subprocess.PIPE, shell=True)
p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -S --threads 1 -q --phred33-quals --norc -n 1 -l 15 -e 200 -a -m 200 --best --strata' % (ref, input_filename), stdout=subprocess.PIPE, shell=True)
# p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -k 200 --sam' % (mm10, input_filename), stdout=subprocess.PIPE, shell=True)
p2 = subprocess.Popen('/Users/averes/miniconda3/envs/py27/bin/python filter_alignments.py -m 4 --counts %s --split_ambi > /dev/null' % (counts_filename),  stdin=p1.stdout, shell=True)

# MIXED REF run
# p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -a --best --strata --sam --norc -n 2 --seedlen 15 --chunkmbs 300' % (mixed, input_filename), stdout=subprocess.PIPE, shell=True)
# p2 = subprocess.Popen('/Users/averes/miniconda3/envs/py27/bin/python filter_alignments.py -m 4 --counts %s --split_ambi --mixed_ref > /dev/null' % (counts_filename),  stdin=p1.stdout, shell=True)
# p1.stdout.close()
# p2.wait()


# p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s --phred33 -m 200 -a --best --strata --sam --norc -n 2 --seedlen 15 --chunkmbs 300 > %s' % (ref, input_filename, raw_sam_filename), shell=True)
# p1.wait()





# with file(counts_filename, 'r') as f:
#     h = next(f).rstrip().split('\t')
#     from collections import defaultdict
#     cc = defaultdict(float)
#     for line in f:
#         d = line.rstrip().split('\t')
#         ref = d[0].split(':')[1]
#         cc[ref] += float(d[1])

# print(cc)


# SILLY 2 READ EXAMPLE

# p1 = subprocess.Popen('/usr/local/bin/bowtie %s - -m 200 -k 200 --sam' % ref, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
# p2 = subprocess.Popen('/Users/averes/miniconda3/envs/py27/bin/python alignments_to_counts_v_allon.py --gtf genes.gtf',  stdin=p1.stdout, shell=True)
# for i,r in enumerate(['GCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC', 'GGGCTCGTGAAGCATGTGGGGGTGAGCCCAGGGGCC']):
#     print(r)
#     fastq = '@AAA-TTT:ABT%d\n%s\n+\n%s\n' % (i, r, 'F'*len(r))
#     p1.stdin.write(fastq)
# p1.stdin.close()
# p2.wait()
