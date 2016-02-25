import re
in_genes="Homo_sapiens.GRCh38.82.with_tid.gtf"
# out_genes="Homo_sapiens.GRCh38.82.annotated.gtf"
out_genes="Homo_sapiens.GRCh38.82.annotated.gtf"
accepted_gene_biotypes_for_NA_transcripts = set(["IG_V_gene","IG_J_gene","protein_coding","TR_J_gene","TR_D_gene","TR_V_gene","IG_C_gene","IG_D_gene","TR_C_gene"])
with open(in_genes, 'r') as in_f, open(out_genes, 'w') as out_f:
    for line in in_f:
        chr_name = line.rstrip().split('\t')[0]
        # Check the transcript_support level
        # This should be faster than a regex
        # We need to support the case where see these two types of annotations:
        #   transcript_support_level "1"
        #   transcript_support_level "1 (assigned to previous version X)"
        #   transcript_support_level "2" <- Clear example of a gene like this is NKX6.1
        #   transcript_support_level "2 (assigned to previous version X)"
        line_valid_for_output = False
        if 'transcript_support_level "1"' in line or 'transcript_support_level "1 ' in line or 'transcript_support_level "2"' in line or 'transcript_support_level "2 ' in line:
            line_valid_for_output = True
        elif 'transcript_support_level "NA' in line:
            # Transcript Support Level Not Analysed. Pseudogenes, single exon transcripts, HLA, T-cell receptor and Ig transcripts are not analysed and therefore not given any of the TSL categories.
            # Keep only a few ones annotated as "IG_V_gene","IG_J_gene","protein_coding","TR_J_gene","TR_D_gene","TR_V_gene","IG_C_gene","IG_D_gene","TR_C_gene"
            gene_biotype = re.search(r'gene_biotype \"(.*?)\";', line)
            if gene_biotype and gene_biotype.group(1) in accepted_gene_biotypes_for_NA_transcripts:
                line_valid_for_output = True
        if line_valid_for_output:
            gene_name = re.search(r'gene_name \"(.*?)\";', line)
            if gene_name:
                gene_name = gene_name.group(1)
                out_line = re.sub(r'(?<=transcript_id ")(.*?)(?=";)', r'\1|'+gene_name, line)
                out_f.write(out_line)