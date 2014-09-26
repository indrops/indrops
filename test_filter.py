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
input_filename = '/Users/averes/Projects/Melton/Dropseq/av1.fastq'
# p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -a --best --strata --sam --norc' % (mixed, input_filename), stdout=subprocess.PIPE, shell=True)
p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -a --best --strata --sam --norc -n 2 --seedlen 15 --chunkmbs 300' % (mixed, input_filename), stdout=subprocess.PIPE, shell=True)
# p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -k 200 --sam' % (ref, input_filename), stdout=subprocess.PIPE, shell=True)
p2 = subprocess.Popen('/Users/averes/miniconda3/envs/py27/bin/python filter_alignments.py -m 1 --counts killme3.counts.txt --split_ambi > killme.bam',  stdin=p1.stdout, shell=True)
p1.stdout.close()
p2.wait()




# SILLY 2 READ EXAMPLE

# p1 = subprocess.Popen('/usr/local/bin/bowtie %s - -m 200 -k 200 --sam' % ref, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
# p2 = subprocess.Popen('/Users/averes/miniconda3/envs/py27/bin/python alignments_to_counts_v_allon.py --gtf genes.gtf',  stdin=p1.stdout, shell=True)
# for i,r in enumerate(['GCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC', 'GGGCTCGTGAAGCATGTGGGGGTGAGCCCAGGGGCC']):
#     print(r)
#     fastq = '@AAA-TTT:ABT%d\n%s\n+\n%s\n' % (i, r, 'F'*len(r))
#     p1.stdin.write(fastq)
# p1.stdin.close()
# p2.wait()
