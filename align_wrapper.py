import subprocess
import os
import errno
import sys

def check_dir(path):
    """
    Checks if directory already exists or not and creates it if it doesn't
    """
    try:
        os.mkdir(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

if __name__=="__main__":
    import sys, argparse

    # for reporting
    print_err = lambda x: sys.stderr.write(str(x)+'\n')

    # get data from command line call of align_wrapper.py
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help='Number of barcodes to align per job.', 
                        type=int, default=20)
    parser.add_argument('filename', metavar='F', type=str, 
                        help='File prefix to be analyzed.')
    parser.add_argument('-c', help='Counter for which chunk of barcodes to submit.',
                        type=int, default=0)
    args = parser.parse_args()

    # define paths
    base_dir = os.path.join('/n/regal/melton_lab/adrianveres/sequencing_libraries/islet_indrops/analyses/', 
                            args.filename)
    fastq_dir = os.path.join(base_dir, 'barcodes', 'fastq')
    trimmed_fastq_dir = os.path.join(base_dir, 'barcodes', 'trimmed')
    bam_dir = os.path.join(base_dir, 'barcodes', 'bam')
    counts_dir = os.path.join(base_dir, 'barcodes', 'counts')
    unal_dir = os.path.join(base_dir, 'barcodes', 'unal')
    frag_dir = os.path.join(base_dir, 'barcodes', 'fragments')

    for path in [base_dir, fastq_dir, trimmed_fastq_dir, 
                 bam_dir, counts_dir, unal_dir, frag_dir]:
        check_dir(path)

    # Maybe we have an array ID, hacky and only work on SLURM/Odyssey
    array_size = args.n
    if 'SLURM_ARRAY_TASK_ID' in os.environ:
        array_position = int(os.environ['SLURM_ARRAY_TASK_ID'])
    else:
        array_position = args.c

    # work on user-defined chunk of barcodes
    for i in range(1+(array_size*array_position), 1+(array_size*(array_position+1))):
        # logging what barcodes have been processed so far
        print_err('')
        print_err('-----')
        print_err('  BC%d' % i)
        print_err('-----')
        cmd_params = {

            'index': '/n/beta_cell/Users/adrianveres/dropseq_data/hg19/hg19.transcripts.annotated',
            'fastq_file': os.path.join(fastq_dir, 'bc%d.fastq' % i),
            'trimmed_fastq_file': os.path.join(trimmed_fastq_dir, 
                                               'bc%d.trimmed.fastq' % i),
            'bam_file': os.path.join(bam_dir, 'bc%d.bam' % i),
            'counts_file': os.path.join(counts_dir, 'bc%d.counts' % i),
            'unal_file': os.path.join(unal_dir, 'bc%d.unal.fastq' % i),
            'frags_file': os.path.join(frag_dir, 'bc%d.frag' % i),
            'python': '/n/home15/adrianveres/envs/py27/bin/python',
            'filter': '/n/beta_cell/Users/adrianveres/Dropseq/' + \
                      'filter_alignments.py',
            'bowtie': '/n/sw/fasrcsw/apps/Core/bowtie/1.1.1-fasrc01/bowtie',
            'trimmomatic' : '/n/sw/centos6/Trimmomatic-0.32/trimmomatic-0.32.jar',

        }

        trim_cmd = "java -jar %(trimmomatic)s SE -threads 1 -phred33  %(fastq_file)s  %(trimmed_fastq_file)s LEADING:28 SLIDINGWINDOW:4:20 MINLEN:16" % cmd_params
        """
        SE: single end mode
        LEADING: minimum quality to keep a base_dir
        SLIDINGWINDOW: number of bases to average across 
            average quality score required
        MINLEN: minimum length of reads to be kept
        """

        count_cmd = "%(bowtie)s %(index)s %(trimmed_fastq_file)s --un %(unal_file)s -p 1 -m 20 -a --best --strata --sam --norc -n 2 --seedlen 15 --chunkmbs 300 | %(python)s %(filter)s -d 400 -m 10 -u 2 --counts %(counts_file)s --fragments %(frags_file)s > %(bam_file)s" % cmd_params
        """
        bowtie
            --un: unaligned fastqs
            -p: threads
            -m: max number of unique alignments allowed for a read
            -a: report all valid alignments
            --strata: report only alignments for a read that fall into best stratum
            --best: required for --strata option
            --sam: print alignments in SAM format
            --norc: no reverse complement alignments
            -n: maximum number of permitted mismatches
        filter_alignments.py
            -d: distance from end of transcript
            -m: multiple alignment threshold
            -u: ambig_count_threshold
            --counts:
        """

        if os.path.isfile(cmd_params['fastq_file']):
            subprocess.call(trim_cmd, shell=True)
            subprocess.call(count_cmd, shell=True)
        else:
            print('FASTQ FILE %s NOT FOUND ' % i)
            print(cmd_params['fastq_file'])
