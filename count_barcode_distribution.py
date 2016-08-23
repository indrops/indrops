import re
from collections import defaultdict
try:
   import cPickle as pickle
except:
   import pickle

from indrops import from_fastq, to_fastq

def count():

    barcode_read_counter = defaultdict(int)

    for name, seq, qual in from_fastq(sys.stdin):
        split_name = name.split(':')
        cell_name = split_name[0]

        barcode_read_counter[cell_name] += 1
        sys.stdout.write(to_fastq(name, seq, qual))
    pickle.dump(dict(barcode_read_counter), sys.stderr)

if __name__=="__main__":
    import sys, argparse
    count()
