import pysam
from collections import defaultdict
import cPickle
from copy import copy
from itertools import combinations

def print_to_log(msg):
    """
    Wrapper to eventually log in smart way, instead of using 'print()'
    """
    sys.stderr.write(str(msg)+'\n')

def quant(args):
    #Convert arg to more explicit names
    multiple_alignment_threshold = args.m
    distance_from_tx_end = args.d
    split_ambiguities = args.split_ambi
    ambig_count_threshold = args.u
    using_mixed_ref = args.mixed_ref

    #Assume that references are named 'transcript_name|gene_name'
    tx_to_gid = lambda tx: tx.split('|')[1] 

    umis_for_geneset = defaultdict(set)
    sam_input = pysam.Samfile("-", "r" )

    # Tuple containing lengths of reference sequences
    ref_lengths = copy(sam_input.lengths)

    # Bam file to be generated
    sam_output = pysam.Samfile('-', "wb", template=sam_input)

    def process_read_alignments(alignments):
        """input: one-element list of a single alignment from a bam file 
        corresponding to a given barcode"""

        # Remove any alignments that aren't supported by a certain number of non-poly A bases!
        dependent_on_polyA_tail = False
        if args.min_non_polyA > 0:
            polyA_independent_alignments = []
            for a in alignments:
                start_of_polyA = ref_lengths[a.tid] - args.polyA
                if a.aend < start_of_polyA:
                    # The alignment doesn't overlap the polyA tail. 
                    polyA_independent_alignments.append(a)
                else:
                    non_polyA_part = start_of_polyA - a.pos
                    if non_polyA_part > args.min_non_polyA:
                        polyA_independent_alignments.append(a)

            dependent_on_polyA_tail = len(polyA_independent_alignments) == 0
            alignments = polyA_independent_alignments



        # We need to obtain Transcript IDs in terms of reference names (Transcrupt_ID|Gene_ID)
        # as opposed to the arbitrary 'a.tid' number
        tx_ids = [sam_input.getrname(a.tid) for a in alignments]


        #Map to Gene IDs
        g_ids = [tx_to_gid(tx_id) for tx_id in tx_ids]
        # finally remove all copies to get a comprehensive unique list of genes
        # found for this barcode
        genes = set(g_ids)

        # Does the alignment map to multiple genes or just one?
        unique = True
        # Was the alignment non-unique, but then rescued to being unique?
        rescued_non_unique = False
        # Even after rescue, was the alignment mapping to more than M genes?
        failed_m_threshold = False

        # The same read could align to transcripts from different genes. 
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

        #Choose 1 alignment per gene, that we will write to the output BAM.
        chosen_alignments = {}
        keep_read = 0 < len(genes) <= multiple_alignment_threshold

        # We need different logic if we are using a mixed organism reference
        if using_mixed_ref:
            refs = set(g.split(':')[1] for g in genes)
            keep_read = (len(refs) == 1) and (0 < len(genes) <= multiple_alignment_threshold)
            

        if keep_read:
            for gene in genes:
                gene_alignments = [a for a in alignments if tx_to_gid(sam_input.getrname(a.tid)) == gene]
                chosen_alignment = sorted(gene_alignments, key=lambda a: ref_lengths[a.tid], reverse=True)[0]
                chosen_alignments[gene] = chosen_alignment
            
        else:
            failed_m_threshold = True

        read_filter_status = (unique, rescued_non_unique, failed_m_threshold, dependent_on_polyA_tail)
        return chosen_alignments, read_filter_status

    # --------------------------
    # Process SAM input
    # (we load everything into memory, so if a single barcode has truly very deep sequencing, we could get into trouble
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

        # The If statements detects that Bowtie is giving info about a different read,
        # so let's process the last one before proceeding
        if not current_read == alignment.qname: 
            #Check that our read has any alignments
            if read_alignments: 
                chosen_alignments, processing_stats = process_read_alignments(read_alignments)
                if chosen_alignments:
                    split_name = current_read.split(':')
                    if len(split_name) == 2:
                        umi = split_name[1] #Adrian format
                    else:
                        umi = split_name[4] #Old Allon format
                    seq = read_alignments[0].seq
                    reads_by_umi[umi][seq] = chosen_alignments

                uniq_count += processing_stats[0]
                non_uniq_count += not(processing_stats[0] or processing_stats[1] or processing_stats[2])
                rescued_count += processing_stats[1]
                failed_m_count += processing_stats[2]

            # We reset the current read info
            current_read = alignment.qname
            read_alignments = []

        read_alignments.append(alignment)

    # Only runs if preceding for loop terminated without break
    # This is not very DRY...
    else:
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

    oversequencing = []
    distance_from_transcript_end = []

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

        # Compute hitting set
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

            #Maximum value of how many reads poitn to any gene
            max_contrib = max(gene_contrib.values())

            #Gene with max contrib
            max_contrib_genes = filter(lambda g: gene_contrib[g]==max_contrib, gene_contrib.keys())

            #Pick a gene among those with the highest value. Which doesn't matter until the last step
            g = max_contrib_genes[0]
            
            read_count_for_umifm = 0
            # umifm_reads = []
            umifm_assigned_unambiguously = False


            for r in copy(r0): #Take a copy of r0 doesn't change as we iterate through it
                if gene_read_mapping[(g, r)]: #Remove any reads from r0 that contributed to the picked gene.
                    r0.remove(r)

                    #Count how many reads we are removing (this is the degree of over-sequencing)
                    read_count_for_umifm += 1
                    # umifm_reads.append(r)


            # If we had equivalent picks, 
            # and their gene contrib value is now 0
            # they were ambiguity partners
            if len(max_contrib_genes) > 1:

                # Update the gene contribs based on the new r0, but on the 'old' g0.
                # That is why we remove g from g0 after this step only
                gene_contrib = dict((gi, sum(gene_read_mapping[(gi, r)] for r in r0)) for gi in g0)
                ambig_partners = filter(lambda g: gene_contrib[g]==0, max_contrib_genes)

                ambig_clique_count[len(ambig_partners)].append(umi)
        
                #Ambig partners will often be a 1-element set. That's ok.
                #Then it will be equivalent to "target_genes[g] = 1."
                if len(ambig_partners) <= ambig_count_threshold:

                    if len(ambig_partners) == 1:
                        umifm_assigned_unambiguously = True
                    
                    for g_alt in ambig_partners:
                        ambig_gene_partners[g_alt].add(frozenset(ambig_partners))
                        target_genes[g_alt] = float(len(ambig_partners))

            else:
                umifm_assigned_unambiguously = True
                target_genes[g] = 1.
                ambig_clique_count[1].append(umi)

            if bool(args.umifm_oversequencing) & umifm_assigned_unambiguously:

                alignment = best_alignment_for_gene[g]

                ref_length = ref_lengths[alignment.tid]
                alignment_start = alignment.pos
                alignment_end = alignment.aend
                dist_from_end = ref_lengths[alignment.tid] - alignment.aend
                oversequencing.append((read_count_for_umifm, g, umi, ref_length, alignment_start, alignment_end))
                # oversequencing.append((read_count_for_umifm, g, umi, dist_from_end, '-'.join(umifm_reads)))


            #Remove g here, so that g is part of the updated gene_contrib, when necessary
            g0.remove(g)

        #For each target gene, output the best alignment
        #and record umi count
        for gene, ambigs in target_genes.items():
            sam_output.write(best_alignment_for_gene[gene])
            
            split_between = ambigs if split_ambiguities else 1.
            umi_counts[gene] += 1./split_between
            ambig_umi_counts[gene] += (1./split_between if ambigs>1 else 0)

    if args.umifm_oversequencing:
        args.umifm_oversequencing.write('Reads in UMIFM,Gene,Reference Length,Alignment Start,Alignment End\n') 
        for metadata in oversequencing:
            args.umifm_oversequencing.write(','.join([str(x) for x in metadata]) + '\n')

    #Output the counts per gene
    all_genes = set()
    for ref in sam_input.references:
        gene = ref.split('|')[1]
        all_genes.add(gene)
    
    args.counts.write('%s\n' % '\t'.join(['Gene Symbol', 'Counts', 'Ambig counts', 'Ambig partners']))
    for gene in sorted(all_genes):
        if ambig_gene_partners[gene]:
            ambig_partners = frozenset.union(*ambig_gene_partners[gene])-frozenset((gene,))
        else:
            ambig_partners = []
        data_row = [gene, str(umi_counts[gene]), str(ambig_umi_counts[gene]), ' '.join(ambig_partners)]
        args.counts.write('%s\n' % '\t'.join(data_row))
    args.counts.close()

    #Output the fixing metrics
    if args.metrics:
        args.metrics.write('Reads with a single alignment\t%d\n' %  uniq_count)
        args.metrics.write('Reads with a single alignment after  filtering on distance from transcript end \t%d\n' %  rescued_count)
        args.metrics.write('Reads with M=>, >1 alignments\t%d\n' %  non_uniq_count)
        args.metrics.write('(Rejected) Reads with >M alignments\t%d\n' %  failed_m_count)
        args.metrics.write('#Numbers of UMIFM with degree of ambiguity')
        for k, v in ambig_clique_count.items():
            args.metrics.write('%d\t%d\n' % (k-1, len(v)))


    
    print_to_log("\nReads with unique mapping (%d), with unique after rescue (%d), with non unique (%d), failed m threshold (%d)" % (uniq_count, rescued_count, non_uniq_count, failed_m_count))
    print_to_log("\n")
    print_to_log("Number of counts with degree of ambiguity:\n")
    for k, v in ambig_clique_count.items():
        sys.stderr.write(' %d : %d\n' % (k-1, len(v)))


if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', help='Ignore reads with more than M alignments, after filtering on distance from transcript end.', type=int, default=4)
    parser.add_argument('-u', help='Ignore counts from UMI that should be split among more than U genes.', type=int, default=4)
    parser.add_argument('-d', help='Maximal distance from transcript end.', type=int, default=525)
    parser.add_argument('--polyA', help='Length of polyA tail in reference transcriptome.', type=int, default=125)
    parser.add_argument('--split_ambi', help="If umi is assigned to m genes, add 1/m to each gene's count (instead of 1)", action='store_true', default=False)
    parser.add_argument('--mixed_ref', help="Reference is mixed, with records named 'gene:ref', should only keep reads that align to one ref.", action='store_true', default=False)
    parser.add_argument('--counts', type=argparse.FileType('w'))
    parser.add_argument('--min_non_polyA', type=int, default=0)
    parser.add_argument('--umifm_oversequencing', type=argparse.FileType('w'), nargs='?')
    parser.add_argument('--metrics', type=argparse.FileType('w'), nargs='?')
    args = parser.parse_args()
    quant(args)

