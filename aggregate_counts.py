import os.path
from collections import defaultdict


def aggregate(filenames, output_handle):
    data = defaultdict(dict)
    samples = []
    ambiguity_counts = defaultdict(int)
    ambiguity_partners = defaultdict(set)

    for filename in filenames:
        sample_name = os.path.basename(filename).split('.')[0]
        samples.append(sample_name)
        with open(filename, 'r') as f:
            header = next(f).rstrip('\n').split('\t')
            for line in f:
                row = line.rstrip('\n').split('\t')
                gene = row[0]
                counts = int(row[1])
                ambig_counts = int(row[2])
                partners = set(row[3].split())

                data[gene][sample_name] = counts
                ambiguity_counts[gene] += ambig_counts
                ambiguity_partners[gene] = ambiguity_partners[gene].union(partners)


    print(set(ambiguity_counts.values()))
    out_header = ['gene'] + ['Sum_counts', 'Sum_ambig', 'Ambigs'] + samples
    output_handle.write('%s\n' % '\t'.join(out_header))
    for gene in sorted(data.keys()):
        per_sample_counts = [data[gene][s] for s in samples]

        row = [gene] + [sum(per_sample_counts), ambiguity_counts[gene], ' '.join(ambiguity_partners[gene])] + per_sample_counts 
        row = [str(r) for r in row]
        output_handle.write('%s\n' % '\t'.join(row))


if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=argparse.FileType('w'), help='Output filename')
    parser.add_argument('filenames', type=str, nargs='+')
    args = parser.parse_args()
    aggregate(args.filenames, args.out)