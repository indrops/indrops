import subprocess
import os

if __name__=="__main__":
    import sys, argparse
    if 'SLURM_ARRAY_TASK_ID' in os.environ:
        array_position = int(os.environ['SLURM_ARRAY_TASK_ID'])

    else:
        array_position = 1
        pass
        # raise Exception
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help='Number of barcodes to align per job.', type=int, default=20)
    args = parser.parse_args()

    base_dir = '/n/regal/melton_lab/adrianveres/datasets/S6D13_cells/data/'
    fastq_dir = os.path.join(base_dir, 'barcodes', 'fastq')
    trimmed_fastq_dir = os.path.join(base_dir, 'barcodes', 'trimmed')
    bam_dir = os.path.join(base_dir, 'barcodes', 'bam')
    counts_dir = os.path.join(base_dir, 'barcodes', 'counts')
    unal_dir = os.path.join(base_dir, 'barcodes', 'unal')

    for i in range(1+(args.n*array_position), 1+(args.n*(array_position+1))):
        print('')
        print('-----')
        print('  BC%d'%i)
        print('-----')
        cmd_params = {

            'index': '/n/regal/melton_lab/adrianveres/datasets/references/rsem/hg19.transcripts.annotated',
            'fastq_file': os.path.join(fastq_dir, 'bc%d.fastq' % i),
            'trimmed_fastq_file': os.path.join(trimmed_fastq_dir, 'bc%d.trimmed.fastq' % i),
            'bam_file': os.path.join(bam_dir, 'bc%d.bam' % i),
            'counts_file': os.path.join(counts_dir, 'bc%d.counts' % i),
            'unal_file': os.path.join(unal_dir, 'bc%d.unal.fastq' % i),
            'python': '/n/home15/adrianveres/envs/main/bin/python',
            'filter': '/n/beta_cell/Users/adrianveres/Dropseq/filter_alignments.py',
            'bowtie': '/n/sw/fasrcsw/apps/Core/bib/2014.05.19-fasrc01/active/bin/bowtie',
            'trimmomatic' : '/n/sw/centos6/Trimmomatic-0.32/trimmomatic-0.32.jar',

        }

        trim_cmd = "java -jar %(trimmomatic)s SE -threads 1 -phred33  %(fastq_file)s  %(trimmed_fastq_file)s LEADING:28 SLIDINGWINDOW:4:20 MINLEN:16" % cmd_params
        count_cmd = "%(bowtie)s %(index)s %(fastq_file)s --un %(unal_file)s -p 1 -m 20 -a --best --strata --sam --norc -n 2 --seedlen 15 --chunkmbs 300 | %(python)s %(filter)s -d 400 -m 10 -u 2 --counts %(counts_file)s > %(bam_file)s" % cmd_params
        if os.path.isfile(cmd_params['fastq_file']):
            subprocess.call(trim_cmd, shell=True)
            subprocess.call(count_cmd, shell=True)
        else:
            print(' FASTQ NOT FOUND ')
