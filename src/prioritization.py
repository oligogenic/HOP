from .transform_data import chunked_transform_data
from .run_varcopp import run_varcopp
from .get_disease_relevance_score import get_rwr_gene_scores
from .load_data import load_combs
import os
import pickle5
import configparser
import gzip
import pickle

config = configparser.ConfigParser()
config.read('config.ini')
PATH_TO_DATA = config['BASE']['PATH_TO_DATA']


def full_varcopp_prediction(variants, pairs, fileName, varcopp_model, olida_pair=None):
    """Predict full exome with VarCoPP model"""
    feats_order = ['CADD1', 'CADD2', 'CADD3', 'CADD4', 'BiolDist', 'ISPP_AR_A', 'ISPP_AR_B', 'ISPP_AD_A', 'ISPP_AD_B',
                   'ISPP_XL_A', 'HIPred_A', 'HIPred_B', 'dn_ds_A', 'BPO_sim', 'KG_dist']
    outfile = open(PATH_TO_DATA + fileName + '_varcopp_combinations.txt', 'w+')

    # Initialize minmax for varcopp score
    best_varcopp_score = 0
    worst_varcopp_score = 1

    # Load varcopp model
    if varcopp_model is None:
        with open(PATH_TO_DATA + 'varcopp_data/varcopp_models/VarCoPP_400BRF.pickle', 'rb') as file:
            varcopp_model = pickle5.load(file)

    # Create combinations, get feature vectors, predict combinations and write results in chunks
    for combinations, combinations_vectors, combinations_ids in chunked_transform_data(pairs, variants, feats_order, 100000, olida_pair):
        final_results = run_varcopp(combinations_vectors, combinations_ids, varcopp_model)
        to_write = ''
        for id, res in final_results.items():
            to_write += '\t'.join([str(id)] + combinations[id] + [str(r) for r in res]) + '\n'
            if res[0] > best_varcopp_score:
                best_varcopp_score = res[0]
            if res[0] < worst_varcopp_score:
                worst_varcopp_score = res[0]
        outfile.write(to_write)
    outfile.close()
    # gzip the file with the varcopp predictions to make it smaller
    with open(PATH_TO_DATA + fileName + '_varcopp_combinations.txt', 'rb') as f_in, gzip.open(PATH_TO_DATA + fileName + '_varcopp_combinations.txt' + '.gz', 'wb') as f_out:
        f_out.writelines(f_in)
    os.remove(PATH_TO_DATA + fileName + '_varcopp_combinations.txt')
    with open(PATH_TO_DATA + fileName + 'minmax.p', 'wb') as file:
        pickle.dump((best_varcopp_score, worst_varcopp_score), file)
    return best_varcopp_score, worst_varcopp_score


def write_combinations_results(comb_results, output_file):
    """Write the results to an output file"""
    to_write = ''
    for c in comb_results:
        to_write += '\t'.join([str(round(it, 3)) if type(it) != str else str(it) for it in c]) + '\n'
    output_file.write(to_write)
    return None


def compute_pheno_relevance_and_normalize(combinations_results, rwr_results, max_pheno, min_pheno, max_varcopp, min_varcopp):
    """
    :param combinations_results:
    :param rwr_results:
    :param max_pheno:
    :param min_pheno:
    :return:
    """
    for id, comb in combinations_results.items():
        comb.pheno_relevance_score = (rwr_results.get(comb.geneA, 0) + rwr_results.get(comb.geneB, 0))/2
        comb.varcopp_score_normalized = (comb.varcopp_score - min_varcopp) / (max_varcopp - min_varcopp)
        if (max_pheno - min_pheno) != 0:
            comb.pheno_relevance_score_normalized = (comb.pheno_relevance_score - min_pheno) / (max_pheno - min_pheno)
        else:
            comb.pheno_relevance_score_normalized = 0
        comb.final_score = (comb.varcopp_score + comb.pheno_relevance_score) / 2
        comb.final_score_normalized = (comb.varcopp_score_normalized + comb.pheno_relevance_score_normalized) / 2
    return combinations_results


def compute_disease_relevance(row, rwr_results, max_rwr, min_rwr):
    if max_rwr-min_rwr != 0:
        return ((rwr_results.get(row[3],0) + rwr_results.get(row[4],0))/2 -min_rwr)/(max_rwr-min_rwr)
    else:
        return 0


def compute_disease_relevance_and_final_scores(comb_results, rwr_results_list, max_rwr_list, min_rwr_list, max_varcopp, min_varcopp, all_operators):
    varcopp_diff = max_varcopp - min_varcopp
    if type(rwr_results_list) == list:
        # If the rwr results are passed as a list it is the results for different types of seeds
        for row in comb_results:
            disease_and_final_scores = []
            row[5] = (float(row[5]) - min_varcopp) / varcopp_diff  # MinMax Normalization of VarCoPP
            for i in range(len(rwr_results_list)):  # For each type of seeds
                rwr_results = rwr_results_list[i]
                max_rwr = max_rwr_list[i]
                min_rwr = min_rwr_list[i]
                disease_relevance_score = compute_disease_relevance(row, rwr_results, max_rwr, min_rwr)
                disease_and_final_scores.append(disease_relevance_score)
                disease_and_final_scores.append((float(row[5]) + disease_relevance_score)/2)  # add the average
                if all_operators:
                    disease_and_final_scores.append(max([float(row[5]), disease_relevance_score]))
                    disease_and_final_scores.append(min([float(row[5]), disease_relevance_score]))
                    disease_and_final_scores.append(float(row[5]) * disease_relevance_score)
            row += disease_and_final_scores
    else:  # Else there is only one type of seed used
        for row in comb_results:
            disease_relevance_score = compute_disease_relevance(row, rwr_results_list, max_rwr_list, min_rwr_list)
            row[5] = (float(row[5]) - min_varcopp) / varcopp_diff
            row.append(disease_relevance_score)
            row.append((float(row[5]) + disease_relevance_score) / 2)
            if all_operators:
                row.append(max([float(row[5]), disease_relevance_score]))
                row.append(min([float(row[5]), disease_relevance_score]))
                row.append(float(row[5]) * disease_relevance_score)
    return [[round(it, 3) if type(it) != str else it for it in c] for c in comb_results]


def update_top_combs(top_combs, combinations_results, N, index):
    top_combs = sorted(top_combs + combinations_results, key=lambda i: i[index], reverse=True)[0:N]
    return top_combs


def combine_pheno_and_varcopp_scores(rwr_results, varcopp_file_name, output_file, max_rwr, min_rwr, max_varcopp, min_varcopp, chunksize=1000,
                                     all=False, nb_tops=100, index_col=-1):
    file = gzip.open(PATH_TO_DATA + varcopp_file_name, 'rt')
    top_combs = []
    for comb_results in load_combs(file, chunksize=chunksize):
        combinations_results = compute_disease_relevance_and_final_scores(comb_results, rwr_results, max_rwr, min_rwr, max_varcopp, min_varcopp, all)
        if output_file is not None:
            write_combinations_results(combinations_results, output_file)
        top_combs = update_top_combs(top_combs, combinations_results, nb_tops, index_col)
    return top_combs


def full_prioritization(variants, pairs, seeds, output_file_name, restart=0.3, write_raw=False, nb_tops=100, rank_col=-1):
    """
    :param patient_name:
    :param seeds: List of list of seeds with the first list being HPOs, second list being genes in the gene panel and third list being HPOs + gene_panel
    :param output_file_name:
    :return:
    """
    header = 'Comb_id\tVar_A\tVar_B\tGeneA\tGeneB\tPathogenicityScore\tVarCoPP_class\tDiseaseScore\tFinalScore\n'
    if write_raw:
        output_file = open(PATH_TO_DATA + output_file_name + '_raw.txt', 'w+')  # This file will contain the raw results not sorted
        output_file.write(header)
    else:
        output_file = None
    print('### Predicting with VarCoPP2.0')

    max_varcopp, min_varcopp = full_varcopp_prediction(variants, pairs, output_file_name, None, None)
    print('### Getting RWR results')
    exome_genes = list(set([v['gene']['gene_name'] for v in variants.values()]))  # Get all unique mutated genes
    rwr_results = get_rwr_gene_scores(seeds,restart)
    rwr_results_patient = {r: rwr_results.get(r, 0) for r in exome_genes}  # Keep only results for the patient's genes
    sorted_rwr_res = sorted(list(rwr_results_patient.values()), reverse=True)
    max_rwr = (sorted_rwr_res[0] + sorted_rwr_res[1]) / 2
    min_rwr = (sorted_rwr_res[-1] + sorted_rwr_res[-2]) / 2
    print('### Combining scores in final scores')
    top_combs = combine_pheno_and_varcopp_scores(rwr_results,
                                     output_file_name + '_varcopp_combinations.txt.gz',
                                     output_file,
                                     max_rwr, min_rwr, max_varcopp,min_varcopp, nb_tops=nb_tops,
                                                 index_col=rank_col)
    if write_raw:
        output_file.close()
        with open(PATH_TO_DATA + output_file_name + '_raw.txt', 'rb') as f_in, gzip.open(
                PATH_TO_DATA + output_file_name + '_raw.txt' + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
        os.remove(PATH_TO_DATA + output_file_name + '_raw.txt')

    if len(top_combs) > 0:
        with open(PATH_TO_DATA + output_file_name + '_top' + str(nb_tops) +'.txt', 'w+') as file:
            file.write(header)
            for c in top_combs:
                file.write('\t'.join([str(round(it, 3)) if type(it) != str else str(it) for it in c]) + '\n')
    os.remove(PATH_TO_DATA + output_file_name + '_varcopp_combinations.txt.gz')

    os.remove(PATH_TO_DATA + output_file_name + 'minmax.p')
    return None


