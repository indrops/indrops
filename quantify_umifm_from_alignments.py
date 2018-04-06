import pysam
from collections import defaultdict
try:
   import cPickle as pickle
except:
   import pickle

from copy import copy
from itertools import combinations
from numpy import memmap
# from indrops import load_indexed_memmapped_array

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
    sam_input = pysam.AlignmentFile("-", "r" )

    # Tuple containing lengths of reference sequences
    ref_lengths = copy(sam_input.lengths)

    # Bam file to be generated
    if args.bam:
        sam_output = pysam.AlignmentFile(args.bam, "wb", template=sam_input)


    # Load cache of low complexity regions
    soft_masked_regions = None
    if args.soft_masked_regions:
        low_complexity_regions = pickle.load(args.soft_masked_regions)
        soft_masked_regions = defaultdict(set)
        for tx, regions in low_complexity_regions.items():
            if regions:
                soft_masked_regions[tx] = set.union(*[set(range(a,b)) for a,b in regions])
    soft_masked_fraction_threshold = 0.5

    def process_read_alignments(alignments):
        """input: one-element list of a single alignment from a bam file 
        corresponding to a given barcode"""

        # Remove any alignments that aren't supported by a certain number of non-poly A bases.
        dependent_on_polyA_tail = False
        if args.min_non_polyA > 0:
            polyA_independent_alignments = []
            for a in alignments:
                start_of_polyA = ref_lengths[a.reference_id] - args.polyA
                if a.reference_end < start_of_polyA:
                    # The alignment doesn't overlap the polyA tail. 
                    polyA_independent_alignments.append(a)
                else:
                    non_polyA_part = start_of_polyA - a.reference_start
                    if non_polyA_part > args.min_non_polyA:
                        polyA_independent_alignments.append(a)

            dependent_on_polyA_tail = len(polyA_independent_alignments) == 0
            alignments = polyA_independent_alignments

        # Remove any alignments that are mostly to low complexity regions
        if soft_masked_regions:
            for a in alignments:
                tx_id = sam_input.getrname(a.reference_id)
                soft_masked_bases = soft_masked_regions[tx_id].intersection(set(range(a.reference_start, a.reference_end)))
                soft_masked_fraction = float(len(soft_masked_bases))/(a.reference_end - a.reference_start)
                a.setTag('XC', '%.2f' % soft_masked_fraction)

            alignments = [a for a in alignments if float(a.opt('XC')) < soft_masked_fraction_threshold]

        # We need to obtain Transcript IDs in terms of reference names (Transcrupt_ID|Gene_ID)
        # as opposed to the arbitrary 'a.reference_id' number
        tx_ids = [sam_input.getrname(a.reference_id) for a in alignments]

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

            close_alignments = [a for a in alignments if (ref_lengths[a.reference_id] - a.reference_end)<distance_from_tx_end]
            close_tx_ids = [sam_input.getrname(a.reference_id) for a in close_alignments]
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
                gene_alignments = [a for a in alignments if tx_to_gid(sam_input.getrname(a.reference_id)) == gene]
                chosen_alignment = sorted(gene_alignments, key=lambda a: ref_lengths[a.reference_id], reverse=True)[0]
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
    not_aligned_count = 0

    current_read = None
    read_alignments = []

    reads_by_umi = defaultdict(dict)

    rev = 0
    non_rev = 0
    for alignment in sam_input:
        #Skip alignments that failed to align...
        if alignment.reference_id == -1:
            not_aligned_count += 1
            # if args.bam:
            #     sam_output.write(alignment)
            continue

        # The If statements detects that Bowtie is giving info about a different read,
        # so let's process the last one before proceeding
        if not current_read == alignment.query_name: 
            #Check that our read has any alignments
            if read_alignments: 
                chosen_alignments, processing_stats = process_read_alignments(read_alignments)
                if chosen_alignments:
                    split_name = current_read.split(':')
                    if len(split_name) == 2:
                        umi = split_name[1] #Old Adrian Format
                    elif len(split_name) == 3:
                        umi = split_name[1] #Adrian format
                    else:
                        umi = split_name[4] #Old Allon format
                    seq = read_alignments[0].seq
                    reads_by_umi[umi][alignment.query_name] = chosen_alignments

                uniq_count += processing_stats[0]
                non_uniq_count += not(processing_stats[0] or processing_stats[1] or processing_stats[2])
                rescued_count += processing_stats[1]
                failed_m_count += processing_stats[2]

            # We reset the current read info
            current_read = alignment.query_name
            read_alignments = []

        read_alignments.append(alignment)

    # Only runs if preceding for loop terminated without break
    # This is not very DRY...
    else:
        if read_alignments:
            chosen_alignments, processing_stats = process_read_alignments(read_alignments)
            if chosen_alignments:
                split_name = current_read.split(':')
                if len(split_name) == 2:
                    umi = split_name[1] #Old Adrian Format
                elif len(split_name) == 3:
                    umi = split_name[1] #Adrian format
                else:
                    umi = split_name[4] #Allon format
                seq = read_alignments[0].seq
                reads_by_umi[umi][alignment.query_name] = chosen_alignments

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

    temp_sam_output = []

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
            # min_ambiguity_alignments = alignments[min(alignments.keys())]
            # max_qual = max(a.mapq for a in min_ambiguity_alignments)
            # max_qual_alignments = filter(lambda a: a.mapq==max_qual, min_ambiguity_alignments)
            # best_alignment_for_gene[gene] = max(max_qual_alignments, key=lambda a: a.qlen)
            best_alignment_for_gene[gene] = alignments[min(alignments.keys())]

        # Compute hitting set
        g0 = set.union(*(set(gs) for gs in umi_reads.values())) #Union of the gene sets of all reads from that UMI
        r0 = set(umi_reads.keys())
        gene_read_mapping = dict()
        for g in g0:
            for r in r0:
                gene_read_mapping[(g, r)] = float(g in umi_reads[r])/(len(umi_reads[r])**2)

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
               
        
                #Ambig partners will often be a 1-element set. That's ok.
                #Then it will be equivalent to "target_genes[g] = 1."
                if len(ambig_partners) <= ambig_count_threshold:
                    if len(ambig_partners) == 1:
                        umifm_assigned_unambiguously = True
                        ambig_clique_count[0].append(umi)
                    
                    for g_alt in ambig_partners:
                        ambig_gene_partners[g_alt].add(frozenset(ambig_partners))
                        target_genes[g_alt] = float(len(ambig_partners))
                        if len(ambig_partners) != 1:
                            ambig_clique_count[len(ambig_partners)].append(umi)

            else:
                umifm_assigned_unambiguously = True
                target_genes[g] = 1.
                ambig_clique_count[1].append(umi)

            #Remove g here, so that g is part of the updated gene_contrib, when necessary
            g0.remove(g)

        #For each target gene, output the best alignment
        #and record umi count
        for gene, ambigs in target_genes.items():
            supporting_alignments = best_alignment_for_gene[gene]
            if args.bam:
                for alignment_for_output in best_alignment_for_gene[gene]:
                    # Add the following tags to aligned reads:
                    # XB - Library Name
                    # XB - Barcode Name
                    # XU - UMI sequence
                    # XO - Oversequencing number (how many reads with the same UMI are assigned to this gene)
                    # YG - Gene identity
                    # YK - Start of the alignment, relative to the transcriptome
                    # YL - End of the alignment, relative to the transcriptome
                    # YT - Length of alignment transcript
                    alignment_for_output.setTag('XL', args.library)
                    alignment_for_output.setTag('XB', args.barcode)
                    alignment_for_output.setTag('XU', umi)
                    alignment_for_output.setTag('XO', len(supporting_alignments))
                    alignment_for_output.setTag('YG', gene)
                    alignment_for_output.setTag('YK', int(alignment_for_output.pos))
                    alignment_for_output.setTag('YL', int(alignment_for_output.reference_end))
                    alignment_for_output.setTag('YT', int(ref_lengths[alignment.reference_id]))
                    temp_sam_output.append(alignment_for_output)
            
            split_between = ambigs if split_ambiguities else 1.
            umi_counts[gene] += 1./split_between
            ambig_umi_counts[gene] += (1./split_between if ambigs>1 else 0)

    #Output the counts per gene
    all_genes = set()
    for ref in sam_input.references:
        gene = ref.split('|')[1]
        all_genes.add(gene)


    sorted_all_genes = sorted(all_genes)
    sorted_metric_columns = ['total_input_reads','single_alignment','rescued_single_alignment','non_unique_less_than_m','non_unique_more_than_m','not_aligned','unambiguous_umifm','umifm_degrees_of_ambiguity_2','umifm_degrees_of_ambiguity_3','umifm_degrees_of_ambiguity_>3']
    output_umi_counts = [umi_counts[gene] for gene in sorted_all_genes]

    if args.write_header:
        args.counts.write('\t'.join(['barcode'] + sorted_all_genes) + '\n')
        args.ambigs.write('\t'.join(['barcode'] + sorted_all_genes) + '\n')
        args.metrics.write('\t'.join(["Barcode","Reads","Reads with unique alignment","Reads with unique alignment within shorter distance of 3'-end","Reads with less than `m` multiple alignments","Reads with more than than `m` multiple alignments","Reads with no alignments", "UMIFM","Ambig UMIFM (between 2 genes)","Ambig UMIFM (between 3 genes)","Ambig UMIFM (between more than 3 genes)",]) + '\n')


    if sum(output_umi_counts) >= args.min_counts:
        ignored = False
        args.counts.write('\t'.join([args.barcode] + [str(int(u)) for u in output_umi_counts]) + '\n')

        # Output sam data
        if args.bam:
            for alignment in temp_sam_output:
                sam_output.write(alignment)
            sam_output.close()

        # Output ambig data
        output_ambig_counts = [ambig_umi_counts[gene] for gene in sorted_all_genes]
        if sum(output_ambig_counts) > 0:

            args.ambigs.write('\t'.join([args.barcode] + [str(int(u)) for u in output_ambig_counts]) + '\n') 
            output_ambig_partners = {}
            for gene in sorted_all_genes:
                if ambig_gene_partners[gene]:
                    gene_partners = frozenset.union(*ambig_gene_partners[gene])-frozenset((gene,))
                    if gene_partners:
                        output_ambig_partners[gene] = gene_partners
            args.ambig_partners.write(args.barcode + '\t'+ str(output_ambig_partners) + '\n')
    else:
        ignored = True
        with open(args.counts.name + '.ignored', 'a') as f:
            f.write(args.barcode + '\n')

    args.counts.close()
    args.ambigs.close()
    args.ambig_partners.close()
    

    #Output the fixing metrics
    total_input_reads = uniq_count + rescued_count + non_uniq_count + failed_m_count + not_aligned_count
    metrics_data = {
        'total_input_reads': total_input_reads,
        'single_alignment': uniq_count,
        'rescued_single_alignment': rescued_count,
        'non_unique_less_than_m': non_uniq_count,
        'non_unique_more_than_m': failed_m_count,
        'not_aligned': not_aligned_count,
        'unambiguous_umifm' : 0,
        'umifm_degrees_of_ambiguity_2' : 0,
        'umifm_degrees_of_ambiguity_3' : 0,
        'umifm_degrees_of_ambiguity_>3' : 0,
    }

    for k, v in ambig_clique_count.items():
        if k == 0:
            metrics_data['unambiguous_umifm'] += len(v)
        elif k == 1:
            metrics_data['unambiguous_umifm'] += len(v)
        elif k == 2:
            metrics_data['umifm_degrees_of_ambiguity_2'] += len(v)
        elif k == 3:
            metrics_data['umifm_degrees_of_ambiguity_3'] += len(v)
        elif k > 3:
            metrics_data['umifm_degrees_of_ambiguity_>3'] += len(v)


    args.metrics.write('\t'.join([args.barcode] + [str(metrics_data[c]) for c in sorted_metric_columns]) + '\n')
    log_output_line = "{0:<8d}{1:<8d}{2:<10d}".format(total_input_reads, metrics_data['unambiguous_umifm'],
        metrics_data['umifm_degrees_of_ambiguity_2']+metrics_data['umifm_degrees_of_ambiguity_3']+metrics_data['umifm_degrees_of_ambiguity_>3'])
    if ignored:
        log_output_line += '  [Ignored from output]'
    print_to_log(log_output_line)

if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', help='Ignore reads with more than M alignments, after filtering on distance from transcript end.', type=int, default=4)
    parser.add_argument('-u', help='Ignore counts from UMI that should be split among more than U genes.', type=int, default=4)
    parser.add_argument('-d', help='Maximal distance from transcript end.', type=int, default=525)
    parser.add_argument('--polyA', help='Length of polyA tail in reference transcriptome.', type=int, default=5)
    parser.add_argument('--split_ambi', help="If umi is assigned to m genes, add 1/m to each gene's count (instead of 1)", action='store_true', default=False)
    parser.add_argument('--mixed_ref', help="Reference is mixed, with records named 'gene:ref', should only keep reads that align to one ref.", action='store_true', default=False)
    parser.add_argument('--min_non_polyA', type=int, default=0)

    # parser.add_argument('--counts', type=argparse.FileType('w'))
    # parser.add_argument('--metrics', type=argparse.FileType('w'))

    parser.add_argument('--counts', type=argparse.FileType('a'))
    parser.add_argument('--metrics', type=argparse.FileType('a'))
    parser.add_argument('--ambigs', type=argparse.FileType('a'))
    parser.add_argument('--ambig-partners', type=argparse.FileType('a'))

    parser.add_argument('--barcode', type=str)
    parser.add_argument('--library', type=str, default='')
    parser.add_argument('--min-counts', type=int, default=0)
    parser.add_argument('--write-header', action='store_true')
    
    
    parser.add_argument('--bam', type=str, nargs='?', default='')
    parser.add_argument('--soft-masked-regions', type=argparse.FileType('r'), nargs='?')
    args = parser.parse_args()
    quant(args)

