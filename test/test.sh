#!/bin/bash
mkdir -p ref
cd ref
wget ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa
wget ftp://ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz
cd ..

# This step took ~30 minutes on my laptop
python indrops.py index ref/Homo_sapiens.GRCh38.cdna.all.fa ref/GRCh38.81.index --gtf-gz ref/Homo_sapiens.GRCh38.81.gtf.gz

python indrops.py preprocess test/test.yaml 
python indrops.py split_barcodes test/test.yaml --read-count-threshold 50 
python indrops.py quantify
python indrops.py aggregate