import re
from collections import defaultdict
try:
   import cPickle as pickle
except:
   import pickle

from indrops import from_fastq, to_fastq

def count(args):

    barcode_read_counter = defaultdict(int)

    args.counts = open(args.filename.name + '.counts.pickle', 'w')

    for name, seq, qual in from_fastq(args.filename):
        split_name = name.split(':')
        cell_name = split_name[0]

        barcode_read_counter[cell_name] += 1

    pickle.dump(dict(barcode_read_counter), args.counts)

    args.counts.close()

if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=argparse.FileType('r'), nargs='?', default=sys.stdin)

    args = parser.parse_args()
    count(args)
