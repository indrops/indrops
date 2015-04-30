import sh
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', help='File prefix to analyze.', type=str)
parser.add_argument('-n', help='Number of fastqs to analyze per job.', type=int, default=50)
args = parser.parse_args()

file_name = args.f
n = args.n
scripts_dir = "/home/man36/dawnchorus/scripts/dbseq"


# fix this so that in future this directory gets passed to align_wrapper.py so I dont have to change it in two places
base_dir = os.path.join("/groups/neuroduo/Aurel/dawnchorus/dawnchorus_data/dropseq/150421_b2/analyses", file_name, "barcodes")

sh.cd(base_dir)

# number of fastq files to process
fq_count = len([name for name in os.listdir('fastq') if os.path.isfile('fastq/' + name)])
# chunk size (isn't this redundant with the args.n?)
n = 50

for i in range(fq_count / n + 1):
    log_id = "align_wrapper_" + str(i*n + 1) + "_" + str(i*n + n) + "__LOG.txt"
    bsub = sh.bsub.bake(W='12:00', q="short", 
        R='select[mem>=16000] && rusage[mem=16000]', 
        o=os.path.join(base_dir, log_id)
    )

    bsub("python %s %s -n %d -c %d" % (os.path.join(scripts_dir, "align_wrapper.py"), file_name, n, i))
       