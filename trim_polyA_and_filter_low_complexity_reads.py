import re
try:
   import cPickle as pickle
except:
   import pickle

from indrops import from_fastq, to_fastq

def low_complexity_filter(args):
    total_reads = 0
    kept_reads = 0
    rejected_because_complexity_too_low = 0
    rejected_because_too_short = 0

    keep_polyA_length = 4
    single_base_runs_regex = '|'.join(['%s{%d,}'%(b, keep_polyA_length+1) for b in 'ATCG'])

    for name, seq, qual in from_fastq(args.input):

        total_reads += 1
        keep_read = True

        #Identify length of polyA tail.
        polyA_length = 0
        for s in seq[::-1]:
            if s!='A':
                break
            polyA_length += 1

        read_length = len(seq)
        trim_at_position = read_length - max(polyA_length - keep_polyA_length, 0)  

        if trim_at_position < args.min_post_trim_length:
            keep_read = False
            rejected_because_too_short += 1
        else:
            new_seq = seq[:trim_at_position]
            new_qual = qual[:trim_at_position]


        low_complexity_bases = sum([m.end()-m.start() for m in re.finditer(single_base_runs_regex, new_seq)])
        low_complexity_fraction = float(low_complexity_bases)/len(new_seq)

        if low_complexity_fraction > args.max_low_complexity_fraction:
            keep_read = False
            rejected_because_complexity_too_low += 1

        if keep_read:
            output_lines = to_fastq(name, new_seq, new_qual)
            args.output.write(output_lines)
            kept_reads += 1

        elif args.rejected:
            args.rejected.write(to_fastq(name, seq, qual))

    if args.metrics:
        pickle.dump({'input': total_reads, 'output': kept_reads, 'rejected_because_complexity_too_low': rejected_because_complexity_too_low, 'rejected_because_too_short': rejected_because_too_short}, args.metrics)

    sys.stderr.write('Kept %d out of %d.\n' % (kept_reads, total_reads))

if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=argparse.FileType('r'), nargs='?', default=sys.stdin)
    parser.add_argument('-output', type=argparse.FileType('w'), nargs='?', default=sys.stdout)
    parser.add_argument('-rejected', type=argparse.FileType('w'), nargs='?', default=False)
    parser.add_argument('-metrics', type=argparse.FileType('w'), nargs='?', default=sys.stderr)
    parser.add_argument('--max-low-complexity-fraction', type=float, nargs='?', default=1.0)
    parser.add_argument('--min-post-trim-length', type=int, nargs='?', default=20)
    args = parser.parse_args()
    low_complexity_filter(args)
