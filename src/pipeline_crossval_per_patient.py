from .train_varcopp_models_crossval import train_crossval_models
from .prioritization import full_varcopp_prediction, combine_pheno_and_varcopp_scores
from.load_data import load_minmax_varcopp
from .create_synthetic_exomes import create_olida_patient_pairs
from .get_disease_relevance_score import get_rwr_gene_scores
from .preprocessing import input_preprocessing_all, input_preprocessing_crossval_indep
import os
import gzip
import configparser
import pickle
import pickle5

config = configparser.ConfigParser()
config.read('config.ini')
PATH_TO_DATA = config['BASE']['PATH_TO_DATA']
PATH_TO_ANNOTATED_DATA = config['BASE']['PATH_TO_DATA']


def load_combinations_fold(i):
    """Load OLIDA combinations ID from fold i"""
    with open(PATH_TO_DATA + 'crossval_folds/olida_combs_fold_' + str(i) + '.p', 'rb') as file:
        combinations_id = pickle.load(file)
    return combinations_id


def load_model(i):
    """Load VarCoPP model for fold i (i.e. varcopp model trained with all folds except i)"""
    with open(PATH_TO_DATA + 'varcopp_data/varcopp_models/varcopp_fold_' + str(i) + '.p', 'rb') as file:
        model = pickle5.load(file)
    return model


def all_files_present_crossval(patient_name, results_folder):
    """Check whether all crossval results files for one patient are present"""
    files = os.listdir(PATH_TO_DATA + results_folder)
    files_patient = [f for f in files if patient_name in f]
    return len(files_patient) == 301


def check_files_present(combs_fold, results_folder, patient_name):
    """Check whether all results files for one patient and one CV fold are present"""
    return all([f in os.listdir(PATH_TO_DATA + 'results' + results_folder) for f in [patient_name + '_' + pos_instance + '_results.txt' + '.gz' for pos_instance in combs_fold]])


def run_crossval_priorization_synthetic_exome_from_template(patient_name, results_folder, restart, compute_all_operators,
                                                            seed_type='All',
                                                            write_raw=False, nb_tops=100, final_score='HPO+Panel_AVGE'):
    """For a given patient, run the full crossvalidation by inserting sequentially one combination in the exome template
    and predicting the exome
    patientName: ID of one of the patient for which the annotated_data is in the 'patients_annotated_data' folder
    results_folder: Name of the folder where the results will be written
    restart: Value between 0 and 1 which represents the restart parameter of the RWR algorithm
    compute_all_operators: Boolean value of whether the combination of the final score should be done with all the operators tested or only the avge"""

    restart, nb_tops, patient_variants, patient_pairs, olida_instances = input_preprocessing_all(restart, patient_name,
                                                                                        'olida_trainset_annotated.json', nb_tops)
    seeds_dictionaries, output_file_header, index_col = input_preprocessing_crossval_indep(seed_type,compute_all_operators,
                                                                                           final_score, True)
    if 'varcopp_results' not in os.listdir(PATH_TO_DATA + 'varcopp_data'):
        os.mkdir(PATH_TO_DATA + 'varcopp_data/varcopp_results')
    # if the varcopp models for the different folds of the CV have not been trained, train them
    if not all(['varcopp_fold_' + str(i) + '.p' in os.listdir('data/varcopp_data/varcopp_models') for i in range(10)]):
        print('### Splitting training set and training VarCoPP2.0 models')
        train_crossval_models()

        if results_folder not in os.listdir(PATH_TO_DATA + 'results'):
            os.mkdir(PATH_TO_DATA + 'results/' + results_folder)

    for fold_nb in range(10):
        print('### Predicting fold ' + str(fold_nb))
        # Create folder for saving the varcopp predictions if the folder does not exist
        if 'fold_' + str(fold_nb) not in os.listdir(PATH_TO_DATA + 'varcopp_data/varcopp_results'):
            os.mkdir(PATH_TO_DATA + 'varcopp_data/varcopp_results/fold_' + str(fold_nb))


        varcopp_model = load_model(fold_nb)  # Load the model for the particular fold

        # Check if the patient was already predicted with this model and run prediction if this was not the case
        if patient_name + '_varcopp_combinations.txt.gz' not in os.listdir(PATH_TO_DATA + 'varcopp_data/varcopp_results/fold_' + str(fold_nb)):
            print('### Predicting ' + patient_name  + ' combinations with VarCoPP2.0')
            full_varcopp_prediction(patient_variants, patient_pairs, 'varcopp_data/varcopp_results/fold_' + str(fold_nb) + '/' + patient_name, varcopp_model)

        combs_fold = load_combinations_fold(fold_nb)
        for pos_instance in combs_fold:
            if patient_name + '_' + pos_instance + '_results.txt' + '.gz' not in os.listdir(PATH_TO_DATA + 'results/' + results_folder):
                # Get annotated OLIDA instances, create pairs with the patient variants and predict these with varcopp
                print('### Predicting ' + patient_name + '_' + pos_instance + ' with HOP')
                comb_variants = olida_instances[pos_instance]
                olida_pair = set([v['gene']['gene_name'] for v in comb_variants.values()])
                pairs_with_olida = create_olida_patient_pairs(comb_variants, patient_variants, patient_pairs)  # Create pairs based on the patient variants and the combination variants
                variants_with_olida = {v: vars for v, vars in patient_variants.items()}
                variants_with_olida.update(comb_variants)

                max_varcopp_with_olida, min_varcopp_with_olida = full_varcopp_prediction(variants_with_olida, pairs_with_olida, 'results/' + results_folder + '/' + patient_name + '_' + pos_instance, varcopp_model, olida_pair)
                synthetic_exome_genes = list(set([v['gene']['gene_name'] for v in variants_with_olida.values()]))

                all_rwr_results = []
                all_max_rwr = []
                all_min_rwr = []
                for type in ['hpos', 'panel', 'hpopanel']:
                    # Get the RWR results for the different seed types and record the max and min possible scores
                    if seeds_dictionaries.get(type) is not None:
                        seeds_dict = seeds_dictionaries[type]
                        rwr_results = get_rwr_gene_scores(seeds_dict[pos_instance], restart)
                        rwr_results_patient = {r: rwr_results.get(r, 0) for r in synthetic_exome_genes}
                        sorted_rwr_res = sorted(list(rwr_results_patient.values()), reverse=True)
                        max_rwr = (sorted_rwr_res[0] + sorted_rwr_res[1]) / 2
                        min_rwr = (sorted_rwr_res[-1] + sorted_rwr_res[-2]) / 2
                        all_rwr_results.append(rwr_results)
                        all_max_rwr.append(max_rwr)
                        all_min_rwr.append(min_rwr)



                max_varcopp, min_varcopp = load_minmax_varcopp(patient_name, fold_nb)
                if write_raw:
                    output_file = open(
                        PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + pos_instance + '_results.txt',
                        'w+')
                    output_file.write('\t'.join(output_file_header) + '\n')
                else:
                    output_file = None

                # Adding the pheno score is done for both the combinations based on the patient variants alone
                # and the combinations of the patient variants with the olida variants

                top_combs_no_olida = combine_pheno_and_varcopp_scores(all_rwr_results, 'varcopp_data/varcopp_results/fold_' + str(fold_nb) + '/' + patient_name + '_varcopp_combinations.txt.gz',
                                                 output_file, all_max_rwr, all_min_rwr, max([max_varcopp, max_varcopp_with_olida]),
                                                 min([min_varcopp, min_varcopp_with_olida]), all=compute_all_operators, nb_tops=nb_tops, index_col=index_col)
                top_combs_olida = combine_pheno_and_varcopp_scores(all_rwr_results, 'results/' + results_folder + '/' + patient_name + '_' + pos_instance + '_varcopp_combinations.txt.gz', output_file,
                                                 all_max_rwr, all_min_rwr, max([max_varcopp, max_varcopp_with_olida]),
                                                 min([min_varcopp, min_varcopp_with_olida]), all=compute_all_operators,nb_tops=nb_tops, index_col=index_col)
                if write_raw:
                    output_file.close()
                    with open(PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + pos_instance + '_results.txt',
                            'rb') as f_in, gzip.open(
                            PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + pos_instance + '_results.txt' + '.gz',
                            'wb') as f_out:
                        f_out.writelines(f_in)
                    os.remove(PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + pos_instance + '_results.txt')

                top_combs = sorted(top_combs_no_olida + top_combs_olida, key=lambda i: i[index_col], reverse=True)[0:nb_tops+1]
                with open(PATH_TO_DATA + 'results/' + results_folder + '/' +  patient_name + '_' + pos_instance + '_top' + str(nb_tops) + '.txt', 'w+') as file:
                    file.write('\t'.join(output_file_header) + '\n')
                    for c in top_combs:
                        file.write('\t'.join([str(round(it, 3)) if not isinstance(it, str) else str(it) for it in c]) + '\n')

                os.remove(PATH_TO_DATA + 'results/' + results_folder + '/' + patient_name + '_' + pos_instance + '_varcopp_combinations.txt.gz')
                os.remove(PATH_TO_DATA +'results/' + results_folder + '/' + patient_name + '_' + pos_instance + 'minmax.p')
        print('### Fold ' + str(fold_nb) + ' done!')
    return None

