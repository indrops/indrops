import subprocess
import os

ref = '/Users/averes/Dropbox/RSEM_index/rsem'
ref = '/Users/averes/Projects/Melton/mm10_transcriptome_reindex/mm10_refseq_annotated_rsem'


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
input_filename = 'killme.fq'
p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -a --best --strata --sam --norc' % (ref, input_filename), stdout=subprocess.PIPE, shell=True)
# p1 = subprocess.Popen('/usr/local/bin/bowtie %s %s -m 200 -k 200 --sam' % (ref, input_filename), stdout=subprocess.PIPE, shell=True)
p2 = subprocess.Popen('/Users/averes/miniconda3/envs/py27/bin/python filter_alignments.py --counts killme.counts > killme.bam',  stdin=p1.stdout, shell=True)
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
