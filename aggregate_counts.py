import os.path
from collections import defaultdict


def aggregate(args):
    filenames = args.filenames
    output_handle = args.out
    if args.ambigs:
        ambigs_output_handle = args.ambigs
    else:
        ambigs_output_handle = None
    samples = []
    counts_data = defaultdict(dict)
    ambig_counts_data = defaultdict(dict)
    ambiguity_partners = defaultdict(set)

    print('Aggregating %d files' % len(filenames))
    for filename in filenames:
        print('Processing %s' % filename)
        sample_name = os.path.basename(filename).split('.')[0]
        samples.append(sample_name)
        with open(filename, 'r') as f:
            header = next(f).rstrip('\n').split('\t')
            for line in f:
                row = line.rstrip('\n').split('\t')
                gene = row[0]
                counts = float(row[1])
                ambig_counts = float(row[2])
                partners = set(row[3].split())

                counts_data[gene][sample_name] = counts
                ambig_counts_data[gene][sample_name] = ambig_counts
                ambiguity_partners[gene] = ambiguity_partners[gene].union(partners)


    # print(set(ambiguity_counts.values()))
    out_header = ['gene'] + ['Sum_counts', 'Sum_ambig', 'Ambigs'] + samples
    output_handle.write('%s\n' % '\t'.join(out_header))

    to_output_line = lambda row: '%s\n' % '\t'.join([str(r) for r in row])

    for gene in sorted(counts_data.keys()):
        per_sample_counts = [counts_data[gene][s] for s in samples]
        per_sample_ambig_counts = [ambig_counts_data[gene][s] for s in samples]

        counts_row = [gene] + [sum(per_sample_counts), sum(per_sample_ambig_counts), ' '.join(ambiguity_partners[gene])] + per_sample_counts
        ambig_counts_row = [gene] + [sum(per_sample_counts), sum(per_sample_ambig_counts), ' '.join(ambiguity_partners[gene])] + per_sample_ambig_counts

        output_handle.write(to_output_line(counts_row))
        if ambigs_output_handle:
            ambigs_output_handle.write(to_output_line(ambig_counts_row))

if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=argparse.FileType('w'), help='Output filename')
    parser.add_argument('--ambigs', type=argparse.FileType('w'), default=None)
    parser.add_argument('filenames', type=str, nargs='+')
    args = parser.parse_args()
    aggregate(args)