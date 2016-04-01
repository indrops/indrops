import re
from collections import defaultdict
try:
   import cPickle as pickle
except:
   import pickle

from indrops import from_fastq, to_fastq

def index_and_count(args):

    barcode_read_counter = defaultdict(int)
    barcode_read_index = defaultdict(list)

    args.index = open(args.output.name + '.index.pickle', 'w')
    args.counts = open(args.output.name + '.counts.pickle', 'w')

    for name, seq, qual in from_fastq(args.input):
        split_name = name.split(':')
        cell_name = split_name[0]

        barcode_read_index[cell_name].append(args.output.tell())
        barcode_read_counter[cell_name] += 1

        args.output.write(to_fastq(name, seq, qual))

    pickle.dump(dict(barcode_read_index), args.index)
    pickle.dump(dict(barcode_read_counter), args.counts)

    args.index.close()
    args.counts.close()
    args.output.close()

if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=argparse.FileType('r'), nargs='?', default=sys.stdin)
    parser.add_argument('-output', type=argparse.FileType('w'), nargs='?')

    args = parser.parse_args()
    index_and_count(args)
