import csv
f = file('/Users/averes/Projects/Melton/temp_dropseq/killme.mixed.counts', 'r')
h = next(f).rstrip().split('\t')
h
from collections import defaultdict
cc = defaultdict(float)
for line in f:
    d = line.rstrip().split('\t')
    ref = d[0].split(':')[1]
    print(float(d[1]))
    cc[ref] += float(d[1])

print(cc)