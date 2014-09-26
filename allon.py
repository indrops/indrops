import csv
f = file('killme3.counts.txt', 'r')
h = next(f).rstrip().split('\t')
h
from collections import defaultdict
cc = defaultdict(float)
for line in f:
    d = line.rstrip().split('\t')
    if 'mm10' in d[0]:
        ref = 'mm10'
    elif 'hg19' in d[0]:
        ref = 'hg19'
    cc[ref] += float(d[1])

print(cc)