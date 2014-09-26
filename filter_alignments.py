import pysam
from collections import defaultdict
import cPickle
from copy import copy
from itertools import combinations

print_err = lambda x: sys.stderr.write(str(x)+'\n')

def run(args):
    multiple_alignment_threshold = args.m
    distance_from_tx_end = args.d
    counts_output_handle = args.counts
    split_ambiguities = args.split_ambi
    print_err('Going for it with %d' % distance_from_tx_end)
    #Assume that references are name 'transcript_name|gene_name'
    tx_to_gid = lambda tx: tx.split('|')[1] 

    umis_for_geneset = defaultdict(set)
    sam_input = pysam.Samfile("-", "r" )

    ref_lengths = copy(sam_input.lengths)
    sam_output = pysam.Samfile('-', "wb", template=sam_input)

    def process_read_alignments(alignments):
        #Transcript IDs
        tx_ids = [sam_input.getrname(a.tid) for a in alignments]
        tx_id_lens = [len(tx_id) for tx_id in tx_ids]

        #Map to Gene IDs
        g_ids = [tx_to_gid(tx_id) for tx_id in tx_ids]
        genes = set(g_ids)

        unique = True
        rescued_non_unique = False
        failed_m_threshold = False

        if 1 < len(genes):
            unique = False

            close_alignments = [a for a in alignments if (ref_lengths[a.tid] - a.aend)<distance_from_tx_end]
            close_tx_ids = [sam_input.getrname(a.tid) for a in close_alignments]
            close_g_ids = [tx_to_gid(tx_id) for tx_id in close_tx_ids]
            close_genes = set(close_g_ids)

            if 0 < len(close_genes) < len(genes):
                alignments = close_alignments
                genes = close_genes
                if len(close_genes) == 1:
                    rescued_non_unique = True

        #Choose 1 alignment per gene to output.
        chosen_alignments = {}
        if 0 < len(genes) <= multiple_alignment_threshold:
            for gene in genes:
                gene_alignments = [a for a in alignments if tx_to_gid(sam_input.getrname(a.tid)) == gene]
                chosen_alignment = sorted(gene_alignments, key=lambda a: ref_lengths[a.tid], reverse=True)[0]
                chosen_alignments[gene] = chosen_alignment
        else:
            failed_m_threshold = True

        read_filter_status = (unique, rescued_non_unique, failed_m_threshold)
        return chosen_alignments, read_filter_status

    # --------------------------
    # Process SAM input
    # (loading everything into memory)
    # --------------------------

    uniq_count = 0
    rescued_count = 0
    non_uniq_count = 0
    failed_m_count = 0

    current_read = None
    read_alignments = []

    reads_by_umi = defaultdict(dict)


    rev = 0
    non_rev = 0
    for alignment in sam_input:
        #Skip alignments that failed to align...
        if alignment.tid == -1:
            continue

        if not current_read == alignment.qname: #Bowtie is giving info about a different read, so let's process the last one before we proceed.
            #We process the data from the previous read
            if read_alignments: #Check that our read has any alignments

                chosen_alignments, processing_stats = process_read_alignments(read_alignments)
                if chosen_alignments:
                    split_name = current_read.split(':')
                    if len(split_name) == 2:
                        umi = split_name[1] #Adrian format
                    else:
                        umi = split_name[4] #Allon format
                    seq = read_alignments[0].seq
                    reads_by_umi[umi][seq] = chosen_alignments

                uniq_count += processing_stats[0]
                non_uniq_count += not(processing_stats[0] or processing_stats[1] or processing_stats[2])
                rescued_count += processing_stats[1]
                failed_m_count += processing_stats[2]

            #We reset the current read info
            current_read = alignment.qname
            read_alignments = []

        read_alignments.append(alignment)

    else:
        #We might have left-over alignments when we are done 
        chosen_alignments, processing_stats = process_read_alignments(read_alignments)
        if chosen_alignments:
            split_name = current_read.split(':')
            if len(split_name) == 2:
                umi = split_name[1] #Adrian format
            else:
                umi = split_name[4] #Allon format
            seq = read_alignments[0].seq
            reads_by_umi[umi][seq] = chosen_alignments

        uniq_count += processing_stats[0]
        non_uniq_count += not(processing_stats[0] or processing_stats[1] or processing_stats[2])
        rescued_count += processing_stats[1]
        failed_m_count += processing_stats[2]
    # -----------------------------
    # Time to filter based on UMIs
    # (and output)
    # --------------------------
    
    umi_counts = defaultdict(float)
    ambig_umi_counts = defaultdict(float)
    ambig_gene_partners = defaultdict(set)
    ambig_clique_count = defaultdict(list)
    for umi, umi_reads in reads_by_umi.items():
        
        #Invert the (read, gene) mapping
        aligns_by_gene = defaultdict(lambda: defaultdict(set))
        for read, read_genes in umi_reads.items():
            for gene, alignment in read_genes.items():
                aligns_by_gene[gene][len(read_genes)].add(alignment)

        #Pick the best alignment for each gene:
        # - least other alignments
        # - highest alignment quality 
        # - longest read
        best_alignment_for_gene = {}
        for gene, alignments in aligns_by_gene.items():
            min_ambiguity_alignments = alignments[min(alignments.keys())]
            max_qual = max(a.mapq for a in min_ambiguity_alignments)
            max_qual_alignments = filter(lambda a: a.mapq==max_qual, min_ambiguity_alignments)
            best_alignment_for_gene[gene] = max(max_qual_alignments, key=lambda a: a.qlen)

        #Compute hitting set
        g0 = set.union(*(set(gs) for gs in umi_reads.values())) #Union of the gene sets of all reads from that UMI
        r0 = set(umi_reads.keys())
        gene_read_mapping = dict()
        for g in g0:
            for r in r0:
                gene_read_mapping[(g, r)] = float(g in umi_reads[r])/(len(umi_reads[r])**2)
                # gene_read_mapping[(g, r)] = float(g in umi_reads[r])

        target_genes = dict()
        #Keys are genes, values are the number of ambiguous partner of each gene

        while len(r0) > 0:
            #For each gene in g0, compute how many reads point ot it
            gene_contrib = dict((gi, sum(gene_read_mapping[(gi, r)] for r in r0)) for gi in g0)

            #Maximum value
            max_contrib = max(gene_contrib.values())

            #Gene with max contrib
            max_contrib_genes = filter(lambda g: gene_contrib[g]==max_contrib, gene_contrib.keys())

            #Pick a gene. Which doesn't matter until the last step
            g = max_contrib_genes[0]
            # target_genes.add(g)
            
            for r in copy(r0): #Take a copy of r0 doesn't change as we iterate through it
                if gene_read_mapping[(g, r)]: #Remove any reads from r0 that contributed to the picked gene.
                    r0.remove(r)

            # If we had equivalent picks, 
            # and their gene contrib value is now 0
            # they were ambiguity partners
            if len(max_contrib_genes) > 1:

                # Update the gene contribs based on the new r0, but on the 'old' g0.
                # That is why we remove g from g0 after this step only
                gene_contrib = dict((gi, sum(gene_read_mapping[(gi, r)] for r in r0)) for gi in g0)
                ambig_partners = filter(lambda g: gene_contrib[g]==0, max_contrib_genes)

                ambig_clique_count[len(ambig_partners)].append(umi)
                
                if len(ambig_partners) > 1:
                    for g_alt in ambig_partners:
                        ambig_gene_partners[g_alt].add(frozenset(ambig_partners))
                        target_genes[g_alt] = float(len(ambig_partners))    
            else:
                target_genes[g] = 1.
                ambig_clique_count[1].append(umi)

            #Remove g here, so that g is part of the updated gene_contrib, when necessary
            g0.remove(g)

        else: 
            #Debug code
            # if umi == 'GCCTTT':
            #     print_err(umi)
            #     for r, gs in umi_reads.items():
            #         print_err(' '+r[:10])
            #         for g in gs:
            #             print_err('  '+g)
            pass

        #For each target gene, output the best alignment
        #and record umi count
        for gene, ambigs in target_genes.items():
            sam_output.write(best_alignment_for_gene[gene])

            split_between = ambigs if split_ambiguities else 1.
            umi_counts[gene] += 1./split_between
            ambig_umi_counts[gene] += (1./split_between if ambigs>1 else 0)


    #Output the counts per gene
    all_genes = set()
    for ref in sam_input.references:
        gene = ref.split('|')[1]
        all_genes.add(gene)
    
    counts_output_handle.write('%s\n' % '\t'.join(['Gene Symbol', 'Counts', 'Ambig counts', 'Ambig partners']))
    for gene in sorted(all_genes):
        if ambig_gene_partners[gene]:
            ambig_partners = frozenset.union(*ambig_gene_partners[gene])-frozenset((gene,))
        else:
            ambig_partners = []
        data_row = [gene, str(umi_counts[gene]), str(ambig_umi_counts[gene]), ' '.join(ambig_partners)]
        counts_output_handle.write('%s\n' % '\t'.join(data_row))
    counts_output_handle.close()

    # #Also debug code for now
    # with open('sets.py', 'w') as f:
    #     for k, v in ambig_clique_count.items():
    #         f.write('s%d = %s\n' % (k, str(set(v))))


    #Output the fixing metrics
    sys.stderr.write("\nReads with unique mapping (%d), with unique after rescue (%d), with non unique (%d), failed m threshold (%d)" % (uniq_count, rescued_count, non_uniq_count, failed_m_count))
    sys.stderr.write("\n")
    sys.stderr.write("Number of counts with degree of ambiguity:\n")
    for k, v in ambig_clique_count.items():
        sys.stderr.write(' %d : %d\n' % (k-1, len(v)))


if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', help='Ignore reads with more than M alignments, after filtering on distance from transcript end.', type=int, default=4)
    parser.add_argument('-d', help='Maximal distance from transcript end.', type=int, default=525)
    parser.add_argument('--split_ambi', help="If umi is assigned to m genes, add 1/m to each gene's count (instead of 1)", action='store_true', default=False)
    parser.add_argument('--counts', type=argparse.FileType('w'))
    args = parser.parse_args()
    run(args)
