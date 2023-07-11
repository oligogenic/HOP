

def order_genes_rvis(gene_A, gene_B):
    rvisA = gene_A['rvis']
    rvisB = gene_B['rvis']

    if (rvisA is None or rvisA == 'None') and rvisB is not None and rvisB != 'None':
        final_geneA, final_geneB = gene_B, gene_A

    elif (rvisA is not None and rvisA != 'None') and (rvisB is None or rvisB == 'None'):
        final_geneA, final_geneB = gene_A, gene_B

    elif (rvisA is None or rvisA == 'None') and (rvisB is None or rvisB == 'None'):
        final_geneA, final_geneB = gene_A, gene_B

    else:
        if min(float(rvisA), float(rvisB)) == float(rvisA):
            final_geneA, final_geneB = gene_A, gene_B
        else:
            final_geneA, final_geneB = gene_B, gene_A
    return final_geneA['gene_name'] + ';' + final_geneB['gene_name']


def get_unique_genes(list_genes):
    unique_genes = []
    for g in list_genes:
        if g not in unique_genes:
            unique_genes.append(g)
    return unique_genes


def create_olida_patient_pairs(olida_comb_variants, patient_variants, annotated_pairs):
    """Create pairs of genes of the patient variants with the OLIDA variants"""
    pairs_with_olida = []
    genes_patients_variants = get_unique_genes([v['gene'] for v in patient_variants.values()])
    genes_positive_instances = get_unique_genes([v['gene'] for v in olida_comb_variants.values()])
    for g1 in genes_positive_instances:
        for g2 in genes_patients_variants:
            if g1['gene_name'] != g2['gene_name']:
                pairs_with_olida.append(order_genes_rvis(g1, g2))

    pairs_with_olida.append(
        order_genes_rvis(genes_positive_instances[0], genes_positive_instances[1]))  # Add the actual OLIDA comb
    pairs_with_olida_annotated = {k: annotated_pairs[k] for k in pairs_with_olida}
    return pairs_with_olida_annotated