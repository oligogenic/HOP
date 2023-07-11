from __future__ import division
import itertools
import json



def create_pair_dict(pairdata):
    """ Takes a list of GenePairData objects and creates a dictionary:
    [(geneA,geneB)]=GenePairData object"""
    pairdict = {}
    for pair in pairdata:
        pairdict[(pair.geneA.gene_name, pair.geneB.gene_name)] = pair

    return pairdict


def create_variant_gene_dict(variants):
    """
    Function that takes a list with variant dicts and returns
    a variant-gene dictionary with the format
    {gene:{hom:[variant_object],het:[variant_object],hemi:[variant_object]}}
    """

    # initiate dictionary
    genevar_dict = {}

    # parse the variants
    for var in variants:

        gene = var['gene']['gene_name']
        zyg = var['gene_var']['zygosity']

        # append information to both dictionaries
        if gene not in genevar_dict:

            genevar_dict.update({gene: {'hom': [], 'het': [], 'hemi': []}})

            # append in the heterozygous/homozygous/hemizygous dictionary of the gene

            if zyg in ['hom', 'Hom', 'homozygous', 'Homozygous', 'Homo']:
                genevar_dict[gene]['hom'] += [var]

            elif zyg in ['het', 'Het', 'Heterozygous', 'heterozygous']:
                genevar_dict[gene]['het'] += [var]
            elif zyg in ['hemi', 'Hemi', 'Hemizygous', 'hemizygous', 'hem', 'Hem']:
                genevar_dict[gene]['hemi'] += [var]

        else:

            # append in the heterozygous/homozygous/hemizygous dictionary of the gene
            if zyg in ['hom', 'Hom', 'homozygous', 'Homozygous', 'Homo']:
                genevar_dict[gene]['hom'] += [var]

            elif zyg in ['het', 'Het', 'Heterozygous', 'heterozygous']:
                genevar_dict[gene]['het'] += [var]

            elif zyg in ['hemi', 'Hemi', 'Hemizygous', 'hemizygous', 'hem', 'Hem']:
                genevar_dict[gene]['hemi'] += [var]

    return genevar_dict


def get_compvars(var_list):
    '''
    Function that takes a list of heterozygous variants and returns a list
    of lists with all possible combinations of those variants.

    It orders the variants inside the combination based on CADD score, so that
    the variant with the highest score is appended first.'''

    combs_list = []

    # iterate over the variants
    for varA, varB in itertools.combinations(var_list, 2):

        # find CADD score
        varA_pred = varA['gene_var']['cadd_raw']
        varB_pred = varB['gene_var']['cadd_raw']

        try:
            comb = find_dangerous_allele(varA, varB, varA_pred, varB_pred)
        except ValueError:
            print('##### Check that all CADD scores have a value.')
            print('##### Script is terminated')
            quit()

        combs_list += [comb]

    return combs_list


def find_dangerous_allele(var1, var2, var1_pred, var2_pred):
    ''' Function that takes two variants and re-orders their appearance
    based on their pathogenicity prediction.
    The first variant is the one the highest score.
    '''

    # first check values availability
    if var1_pred in ['NULL', 'NaN', 'NA', None] and \
            var2_pred not in ['NULL', 'NaN', 'NA', None]:
        comp_var = [var2, var1]

    elif var2_pred in ['NULL', 'NaN', 'NA', None] and \
            var1_pred not in ['NULL', 'NaN', 'NA', None]:
        comp_var = [var1, var2]

    elif var2_pred in ['NULL', 'NaN', 'NA', None] and \
            var1_pred in ['NULL', 'NaN', 'NA', None]:
        comp_var = [var1, var2]

    else:
        # find most dangerous
        var1_pred = float(var1_pred)
        var2_pred = float(var2_pred)

        most_dangerous = max(var1_pred, var2_pred)

        # create heterozygous compound combination
        if most_dangerous == var1_pred:
            comp_var = [var1, var2]
        else:
            comp_var = [var2, var1]

    return comp_var


def alleles_to_vector(comb, varfeats, genefeats, pairfeats):
    '''
    Function that takes a CombinationData object and returns the vectorised version of it.
    '''

    # load default values from config JSON
    with open("data/varcopp_data/imputation_values.json", "r") as predictors_config_file:
        predictors_config = json.load(predictors_config_file)
        defaults_imputation = predictors_config["imputation_values"]

    # initiate the vector
    vector = []
    h = 0  # vector index

    # initiate dictionary with the vector start position of each feature
    feat_pos_dict = {}

    # collect variant information
    # iterate over the features to vectorize them
    for i in range(len(varfeats)):

        # vectorize for the aminoacid property change (Flex)
        if varfeats[i] in ['CADD1']:

            try:
                vector += [comb['varA'][0]['gene_var']['cadd_raw']]
            except ValueError:
                print('##### CADD1 raw score value is invalid.' + \
                      ' Please provide a score for CADD for all variants.')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_name = varfeats[i]
            feat_pos_dict[feat_name] = [h, len(vector) - 1]
            h = len(vector)

        # vectorize for CADD3
        elif varfeats[i] in ['CADD3']:

            try:
                vector += [comb['varB'][0]['gene_var']['cadd_raw']]
            except ValueError:
                print('##### CADD3 raw score value is invalid.' + \
                      ' Please provide a score for CADD for all variants.')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_name = varfeats[i]
            feat_pos_dict[feat_name] = [h, len(vector) - 1]
            h = len(vector)

        # vectorize for CADD2
        elif varfeats[i] in ['CADD2']:

            if len(comb['varA']) == 1:
                if comb['varA'][0]['gene_var']['zygosity'] == 'Heterozygous':
                    vector += [-3.0]
                if comb['varA'][0]['gene_var']['zygosity'] in ['Homozygous', 'Hemizygous']:
                    try:
                        vector += [comb['varA'][0]['gene_var']['cadd_raw']]
                    except ValueError:
                        print('##### CADD2 raw score value is invalid.' + \
                              ' Please provide a score for CADD for all variants.')
                        print('##### Script is terminated')
                        quit()

            if len(comb['varA']) == 2:
                try:
                    vector += [comb['varA'][1]['gene_var']['cadd_raw']]
                except ValueError:
                    print('##### CADD2 raw score value is invalid.' + \
                          ' Please provide a score for CADD for all variants.')
                    print('##### Script is terminated')
                    quit()

            # save vector position of the feature
            feat_name = varfeats[i]
            feat_pos_dict[feat_name] = [h, len(vector) - 1]
            h = len(vector)

        # vectorize for CADD4
        elif varfeats[i] in ['CADD4']:

            if len(comb['varB']) == 1:
                if comb['varB'][0]['gene_var']['zygosity'] == 'Heterozygous':
                    vector += [-3.0]
                if comb['varB'][0]['gene_var']['zygosity'] in ['Homozygous', 'Hemizygous']:
                    try:
                        vector += [comb['varB'][0]['gene_var']['cadd_raw']]
                    except ValueError:
                        print('##### CADD4 raw score value is invalid.' + \
                              ' Please provide a score for CADD for all variants.')
                        print('##### Script is terminated')
                        quit()

            if len(comb['varB']) == 2:
                try:
                    vector += [comb['varB'][1]['gene_var']['cadd_raw']]
                except ValueError:
                    print('##### CADD4 raw score value is invalid.' + \
                          ' Please provide a score for CADD for all variants.')
                    print('##### Script is terminated')
                    quit()

            # save vector position of the feature
            feat_name = varfeats[i]
            feat_pos_dict[feat_name] = [h, len(vector) - 1]
            h = len(vector)

    # collect information on genes info

    for i in range(len(genefeats)):

        # Haploinsufficiency probability A
        if genefeats[i] == 'ISPP_AR_A':
            try:
                if comb['varA'][0]['gene']['ISPP_AR'] in ['N/A', 'NA', 'NaN', 'nan', 'None'] or comb['varA'][0]['gene']['ISPP_AR'] is None:
                    vector += [defaults_imputation["ISPP_AR"]]
                else:
                    vector += [comb['varA'][0]['gene']['ISPP_AR']]
            except ValueError:
                print('##### ISPP_AR_A is invalid.' + \
                      ' Please provide a haploinsufficiency value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_pos_dict[genefeats[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)

        # Haploinsufficiency probability B
        elif genefeats[i] == 'ISPP_AR_B':
            try:
                if comb['varB'][0]['gene']['ISPP_AR'] in ['N/A', 'NA', 'NaN', 'nan', 'None'] or comb['varB'][0]['gene']['ISPP_AR'] is None:
                    vector += [defaults_imputation["ISPP_AR"]]
                else:
                    vector += [comb['varB'][0]['gene']['ISPP_AR']]
            except ValueError:
                print('##### HI_B is invalid.' + \
                      ' Please provide a haploinsufficiency value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_pos_dict[genefeats[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)

        # Recessiveness probability A
        elif genefeats[i] == 'ISPP_AD_A':

            try:
                if comb['varA'][0]['gene']['ISPP_AD'] in ['N/A', 'NA', 'NaN', 'nan', 'None'] or comb['varA'][0]['gene']['ISPP_AD'] is None:
                    vector += [defaults_imputation["ISPP_AD"]]
                else:
                    vector += [comb['varA'][0]['gene']['ISPP_AD']]
            except ValueError:
                print('##### ISPP_AD_A is invalid.' + \
                      ' Please provide a recessiveness value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_pos_dict[genefeats[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)

        # Recessiveness probability B
        elif genefeats[i] == 'ISPP_AD_B':

            try:
                if comb['varB'][0]['gene']['ISPP_AD'] in ['N/A', 'NA', 'NaN', 'nan', 'None'] or comb['varB'][0]['gene']['ISPP_AD'] is None:
                    vector += [defaults_imputation["ISPP_AD"]]
                else:
                    vector += [comb['varB'][0]['gene']['ISPP_AD']]
            except ValueError:
                print('##### ISPP_AD_B is invalid.' + \
                      ' Please provide a recessiveness value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_pos_dict[genefeats[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)

        elif genefeats[i] == 'ISPP_XL_A':
            if comb['varA'][0]['gene_var']['chr'] == 'X':
                try:
                    if comb['varA'][0]['gene']['ISPP_XL'] in ['N/A', 'NA', 'NaN', 'nan', 'None'] or comb['varA'][0]['gene']['ISPP_XL'] is None:
                        vector += [defaults_imputation["ISPP_XL"]]
                    else:
                        vector += [comb['varA'][0]['gene']['ISPP_XL']]
                except ValueError:
                    print('##### ISPP_XL_B is invalid.' + \
                          ' Please provide a recessiveness value or write' + \
                          ' "NA","NaN", or "nan".')
                    print('##### Script is terminated')
                    quit()

            # save vector position of the feature
            feat_pos_dict[genefeats[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)

        elif genefeats[i] == 'HIPred_A':

            try:
                if comb['varA'][0]['gene']['HI_pred'] in ['N/A', 'NA', 'NaN', 'nan', 'None'] or comb['varA'][0]['gene']['HI_pred'] is None:
                    vector += [defaults_imputation["HIPRED"]]
                else:
                    vector += [comb['varA'][0]['gene']['HI_pred']]
            except ValueError:
                print('##### HI_PredA is invalid.' + \
                      ' Please provide a recessiveness value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_pos_dict[genefeats[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)
        elif genefeats[i] == 'HIPred_B':

            try:
                if comb['varB'][0]['gene']['HI_pred'] in ['N/A', 'NA', 'NaN', 'nan', 'None'] or comb['varB'][0]['gene']['HI_pred'] is None:
                    vector += [defaults_imputation["HIPRED"]]
                else:
                    vector += [comb['varB'][0]['gene']['HI_pred']]
            except ValueError:
                print('##### HI_Pred_A is invalid.' + \
                      ' Please provide a recessiveness value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_pos_dict[genefeats[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)
        elif genefeats[i] == 'dn_ds_A':

            try:
                if comb['varA'][0]['gene']['dn_ds'] in ['N/A', 'NA', 'NaN', 'nan', None] or comb['varA'][0]['gene']['dn_ds'] is None:
                    vector += [defaults_imputation["DN_DS"]]
                else:
                    vector += [comb['varA'][0]['gene']['dn_ds']]
            except ValueError:
                print('##### ISPP_AD_B is invalid.' + \
                      ' Please provide a recessiveness value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

            # save vector position of the feature
            feat_pos_dict[genefeats[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)


    # collect information on pairs info
    for i in range(len(pairfeats)):

        # Biological distance
        if pairfeats[i] == 'BiolDist':

            try:
                if comb['BiolDist'] in ['N/A', 'NA', 'NaN', 'nan', 'None', None] or comb['BiolDist'] == [
                    None] or comb['BiolDist'] is None:
                    vector += [defaults_imputation["BIOL_DIST"]]
                else:
                    vector += [comb['BiolDist']]
            except ValueError:
                print('##### BiolDist is invalid.' + \
                      ' Please provide a value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

        elif pairfeats[i] == 'BPO_sim':

            try:
                if comb['BP_sim'] in ['N/A', 'NA', 'NaN', 'nan', 'None', None] or comb['BP_sim'] == [None] \
                        or comb['BP_sim'] is None:
                    vector += [defaults_imputation["BP_SIM"]]
                else:
                    vector += [comb['BP_sim']]
            except ValueError:
                print('##### BiolDist is invalid.' + \
                      ' Please provide a value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()

        elif pairfeats[i] == 'KG_dist':

            try:
                if comb['KG_dist'] in ['N/A', 'NA', 'NaN', 'nan', 'None', None] or comb['KG_dist'] == [
                    None] or comb['KG_dist'] is None:
                    vector += [defaults_imputation["KG_DIST"]]
                else:
                    vector += [comb['KG_dist']]
            except ValueError:
                print('##### BiolDist is invalid.' + \
                      ' Please provide a value or write' + \
                      ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()
        # save vector position of the feature
        feat_pos_dict[pairfeats[i]] = [h, len(vector) - 1]
        # initiate vector index
        h = len(vector)
    return vector, feat_pos_dict


def same_variant(varA, varB, comp_varsA, comp_varsB):
    """
    Check if two variants (or more precisely two sets of variants, to handle compound heterozygote) are the same.
    :param varA: a set of variant A
    :param varB: a set of variant B
    :param comp_varsA: list of coumpound het A
    :param comp_varsA: list of coumpound het B
    :return: True if a variant in varA is also found in varB, False otherwise.
    """
    vsA = varA if varA in comp_varsA else [varA]
    vsB = varB if varB in comp_varsB else [varB]

    for vA in vsA:
        for vB in vsB:
            if vA == vB:
                return True
    return False


def order_features(feats_order, feats_vector, feats_dict):
    new_vector = []
    for f in feats_order:
        index = feats_dict[f][0]
        new_vector.append(feats_vector[index])
    return new_vector


def chunked_transform_data(pairs, variants, feats_order, chunk_length,  olida_pair=None, file=None, index = None):
    varfeats = ['CADD1', 'CADD2', 'CADD3', 'CADD4']
    genefeats = ['ISPP_AR_A', 'ISPP_AR_B', 'ISPP_AD_A', 'ISPP_AD_B', 'ISPP_XL_A', 'HIPred_A', 'HIPred_B', 'dn_ds_A']
    pairfeats = ['BiolDist', 'BPO_sim', 'KG_dist']

    ''' Creates all possible digenic combinations between variants from the same individual'''

    combinations = {}
    combination_vectors = [] # this removes reiterating the combinations to get the vectors for predictions afterwards
    combinations_ids = []
    # create gene - variant dictionary
    d = create_variant_gene_dict([v for v in variants.values()])

    # create pair dictionary
    # pairdict = create_pair_dict(pairs)

    # track ID of digenic combination
    comb_id = 1
    counter = 0
    # start iterating over the pairs
    for pair in pairs:
        pair = pair.split(';')
        geneA, geneB = pair[0], pair[1]
        if d.get(geneA) is not None and d.get(geneB) is not None:
            # categorize and collect variants
            homvarsA = []
            homvarsB = []
            hetvarsA = []
            hetvarsB = []
            hemivarsA = []
            hemivarsB = []
            comp_varsA = []
            comp_varsB = []

            if len(d[geneA]['hom']) > 0:
                for el in d[geneA]['hom']:
                    homvarsA += [el]
            if len(d[geneA]['het']) > 0:
                for el in d[geneA]['het']:
                    hetvarsA += [el]
            if len(d[geneA]['hemi']) > 0:
                for el in d[geneA]['hemi']:
                    hemivarsA += [el]
            if len(d[geneB]['hom']) > 0:
                for el in d[geneB]['hom']:
                    homvarsB += [el]
            if len(d[geneB]['het']) > 0:
                for el in d[geneB]['het']:
                    hetvarsB += [el]
            if len(d[geneB]['hemi']) > 0:
                for el in d[geneB]['hemi']:
                    hemivarsB += [el]

            # create heterozygous compound variants and order based on CADD
            if len(d[geneA]['het']) > 1:
                comp_varsA = get_compvars(d[geneA]['het'])
            if len(d[geneB]['het']) > 1:
                comp_varsB = get_compvars(d[geneB]['het'])

            # combine all variants A
            if len(d[geneA]['het']) > 1 and set(pair) != olida_pair:
                varsA = homvarsA + comp_varsA + hemivarsA
            else:
                varsA = homvarsA + hetvarsA + hemivarsA + comp_varsA

            # iterate variants A
            for varA in varsA:
                if len(d[geneB]['het']) > 1 and set(pair) != olida_pair:
                    varsB = homvarsB + comp_varsB + hemivarsB
                else:
                    varsB = homvarsB + hetvarsB + hemivarsB + comp_varsB

                for varB in varsB:

                    # Handle case where a variant is duplicated because of association with 2 different canonical transcripts (and two different gene names)
                    # This happens for example when the variant is at a position where there is 1 gene on strand + and another gene on strand -.
                    # We don't want to create a pair with the same variant in this case.
                    if same_variant(varA, varB, comp_varsA, comp_varsB):
                        #print(f"Skipping {varA} - {varB}")
                        continue

                    comb = {'geneA':pair[0], 'geneB': pair[1], 'id':comb_id}

                    if varA in comp_varsA:
                        comb['varA'] = varA
                    else:
                        comb['varA'] = [varA]

                    if varB in comp_varsB:
                        comb['varB'] = varB
                    else:
                        comb['varB'] = [varB]


                    # append pair features
                    comb['BiolDist'] = pairs[';'.join([geneA, geneB])]['biol_dist']
                    comb['KG_dist'] = pairs[';'.join([geneA, geneB])]['kg_dist']
                    comb['BP_sim'] = pairs[';'.join([geneA, geneB])]['bpo_sim']

                    # transform combination to vector
                    vector_comb, feat_vp = alleles_to_vector(comb,
                                                             varfeats,
                                                             genefeats,
                                                             pairfeats)
                    vector_comb = order_features(feats_order, vector_comb, feat_vp)
                    varA_comb_key = [v['gene_var']['composite_key'] for v in comb['varA']]
                    varB_comb_key = [v['gene_var']['composite_key'] for v in comb['varB']]
                    if file is not None:
                        file.write('\t'.join([str(comb_id), '/'.join(varA_comb_key), '/'.join(varB_comb_key), geneA, geneB,
                                              ','.join([str(e) for e in vector_comb])]) + '\n')
                    else:

                        comb['vectorcomb'] = vector_comb
                        # add combination in the dictionaries for later use
                        if index is not None:
                            index_name = index + str(comb_id)
                        else:
                            index_name = comb_id
                        combinations[index_name] = ['/'.join(varA_comb_key), '/'.join(varB_comb_key), comb['geneA'], comb['geneB']]
                        combination_vectors.append(vector_comb)
                        combinations_ids.append(index_name)

                    # increase ID number
                    comb_id += 1
                    counter += 1
                    if counter == chunk_length:
                        counter = 0
                        yield combinations, combination_vectors, combinations_ids
                        combinations = {}
                        combination_vectors = []
                        combinations_ids = []

    yield combinations, combination_vectors, combinations_ids
