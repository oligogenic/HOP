"""Full pipeline to reproduce the results for the prioritization in the patients with the
OLIDA combinations from the independent set"""

from .prioritization import full_varcopp_prediction, combine_pheno_and_varcopp_scores
from.load_data import load_minmax_varcopp
from .create_synthetic_exomes import create_olida_patient_pairs
from .get_disease_relevance_score import get_rwr_gene_scores
from .preprocessing import input_preprocessing_all, input_preprocessing_crossval_indep
import os
import gzip
import configparser

config = configparser.ConfigParser()
config.read('config.ini')
PATH_TO_DATA = config['BASE']['PATH_TO_DATA']
PATH_TO_ANNOTATED_DATA = config['BASE']['PATH_TO_ANNOTATED_DATA']



def run_indep_priorization_synthetic_exome_from_template(patient_name, results_folder, restart=0.3, seed_type='All',
                                                         write_raw=False, nb_tops=100, final_score='HPO+Panel_AVGE'):
    """
    For one exome template run HOP for all the synthetic exomes generated with that exome template and the OLIDA combinations
    from the independent set.
    :param patient_name: Name of the patient to use as synthetic template
    :param results_folder: Folder where to save the results
    :param restart: value of the restart parameter, by default set to 0.3
    :return:
    """
    restart, nb_tops, patient_variants, patient_pairs, olida_instances = input_preprocessing_all(restart,patient_name, 'olida_testset_annotated.json', nb_tops)
    seeds_dictionaries, output_file_header, index_col = input_preprocessing_crossval_indep(seed_type, False, final_score, False)

    if patient_name + '_varcopp_combinations.txt.gz' not in os.listdir(PATH_TO_DATA + 'varcopp_data/varcopp_results'):
        print('### Predicting ' + patient_name + ' combinations with VarCoPP2.0')
        full_varcopp_prediction(patient_variants, patient_pairs, 'varcopp_data/varcopp_results/' + patient_name, None,None)

    if results_folder not in os.listdir(PATH_TO_DATA + 'results'):
        os.mkdir(PATH_TO_DATA + 'results/' + results_folder)

    for oli_id in olida_instances.keys():
        print('### Predicting ' + patient_name + '_' + oli_id + ' with HOP')
        # Get annotated positive instances, create pairs with the patient variants and predict these with varcopp
        comb_variants = olida_instances[oli_id]
        olida_pair = set([v['gene']['gene_name'] for v in comb_variants.values()])
        pairs_with_olida = create_olida_patient_pairs(comb_variants, patient_variants, patient_pairs)   # Create pairs based on the patient variants and the combination variants
        variants_with_olida = {v: vars for v, vars in patient_variants.items()}  # Create a copy of the patient variants
        variants_with_olida.update(comb_variants)  # Add the OLIDA variants to this copy

        # Predict the combinations generated with the patient variants and the OLIDA comb variants
        max_varcopp_olida, min_varcopp_olida = full_varcopp_prediction(variants_with_olida, pairs_with_olida, 'results/' + results_folder + '/' + patient_name + '_' + oli_id, None, olida_pair)
        synthetic_exome_genes = list(set([v['gene']['gene_name'] for v in variants_with_olida.values()]))

        # Based on the hpo terms associated with the combination, run rwr, score combinations and write the final results
        all_rwr_results = []
        all_max_rwr = []
        all_min_rwr = []
        for type in ['hpos', 'panel', 'hpopanel']:
            if seeds_dictionaries.get(type) is not None:
                seeds_dict = seeds_dictionaries[type]
                rwr_results = get_rwr_gene_scores(seeds_dict[oli_id], restart)
                rwr_results_patient = {r: rwr_results.get(r, 0) for r in synthetic_exome_genes}
                sorted_rwr_res = sorted(list(rwr_results_patient.values()), reverse=True)
                max_rwr = (sorted_rwr_res[0] + sorted_rwr_res[1]) / 2
                min_rwr = (sorted_rwr_res[-1] + sorted_rwr_res[-2]) / 2
                all_rwr_results.append(rwr_results)
                all_max_rwr.append(max_rwr)
                all_min_rwr.append(min_rwr)
        max_varcopp, min_varcopp = load_minmax_varcopp(patient_name, None)
        if write_raw:
            output_file = open(
                PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + oli_id + '_results.txt',
                'w+')
            output_file.write('\t'.join(output_file_header) + '\n')
        else:
            output_file = None
        # Adding the pheno score is done for both the combinations based on the patient variants alone
        # and the combinations of the patient variants with the olida variants
        top_combs_no_olida = combine_pheno_and_varcopp_scores(all_rwr_results, 'varcopp_data/varcopp_results/' + patient_name + '_varcopp_combinations.txt.gz',
                                         output_file, all_max_rwr, all_min_rwr, max([max_varcopp, max_varcopp_olida]),
                                         min([min_varcopp, min_varcopp_olida]),all=False, nb_tops=nb_tops, index_col=index_col)
        top_combs_olida = combine_pheno_and_varcopp_scores(all_rwr_results,
                                         'results/' + results_folder + '/' + patient_name + '_' + oli_id + '_varcopp_combinations.txt.gz',
                                         output_file,
                                         all_max_rwr, all_min_rwr, max([max_varcopp, max_varcopp_olida]),
                                         min([min_varcopp, min_varcopp_olida]), all=False, nb_tops=nb_tops, index_col=index_col)
        if write_raw:
            output_file.close()
            with open(
                    PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + oli_id + '_results.txt',
                    'rb') as f_in, gzip.open(
                PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + oli_id + '_results.txt' + '.gz',
                'wb') as f_out:
                f_out.writelines(f_in)
            os.remove(
                PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + oli_id + '_results.txt')

        top_combs = sorted(top_combs_no_olida + top_combs_olida, key=lambda i: i[index_col], reverse=True)[
                    0:nb_tops + 1]
        with open(PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + oli_id + '_top' + str(
                nb_tops) + '.txt', 'w+') as file:
            file.write('\t'.join(output_file_header) + '\n')
            for c in top_combs:
                file.write('\t'.join([str(round(it, 3)) if not isinstance(it, str) else str(it) for it in c]) + '\n')

        os.remove(PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + oli_id + '_varcopp_combinations.txt.gz')
        os.remove(PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + oli_id + 'minmax.p')
    return None
