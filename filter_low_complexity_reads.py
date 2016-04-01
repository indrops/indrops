import re
try:
   import cPickle as pickle
except:
   import pickle

from indrops import from_fastq, to_fastq

def low_complexity_filter(args):

    low_complexity_fraction_threshold = 0.45
    min_trimmed_length = 20

    total_reads = 0
    kept_reads = 0
    for name, seq, qual in from_fastq(args.input):

        total_reads += 1
        keep_read = True

        low_complexity_bases = sum([m.end()-m.start() for m in re.finditer('A{5,}|T{5,}|C{5,}|G{5,}', seq)])
        low_complexity_fraction = float(low_complexity_bases)/len(seq)
        if low_complexity_fraction > low_complexity_fraction_threshold:
            keep_read = False
        else:
            #Identify length of polyA tail.
            polyA_length = 0
            for s in seq[::-1]:
                if s!='A':
                    break
                polyA_length += 1

            read_length = len(seq)
            trim_at_position = read_length - polyA_length + 4

            if trim_at_position < min_trimmed_length:
                keep_read = False
            else:
                new_seq = seq[:trim_at_position]
                new_qual = qual[:trim_at_position]

        if keep_read:
            args.output.write(to_fastq(name, new_seq, new_qual))
            kept_reads += 1

        elif args.rejected:
            # Output to rejected read pile
            args.rejected.write(to_fastq(name, seq, qual))

    sys.stderr.write('Kept %d out of %d.\n' % (kept_reads, total_reads))

if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=argparse.FileType('r'), nargs='?', default=sys.stdin)
    parser.add_argument('-output', type=argparse.FileType('w'), nargs='?', default=sys.stdout)
    parser.add_argument('-rejected', type=argparse.FileType('w'), nargs='?', default=False)
    args = parser.parse_args()
    low_complexity_filter(args)
