from .prioritization import full_prioritization
from .load_data import load_olida_annotated_instances, load_patient_annotated_data
from .compute_RWR_matrix import *
from .preprocessing import input_preprocessing_all
import configparser
import os


config = configparser.ConfigParser()
config.read('config.ini')
PATH_TO_DATA = config['BASE']['PATH_TO_DATA']
PATH_TO_ANNOTATED_DATA = config['BASE']['PATH_TO_ANNOTATED_DATA']

def prioritize_single_exome(patient_name, results_folder, restart, seeds, oli_id, write_raw, nb_tops):
    # Parameter checks
    restart, nb_tops, patient_variants, patient_pairs, olida_instances = input_preprocessing_all(restart, patient_name,
                                                                                                 'olida_testset_annotated.json',
                                                                                                 nb_tops)
    seeds = seeds.split(',')
    if results_folder not in os.listdir(PATH_TO_DATA + 'results'):
        os.mkdir(PATH_TO_DATA + 'results/' + results_folder)

    output_filename = 'results/' + results_folder + '/' + patient_name

    if oli_id is not None:
        olida_instances = load_olida_annotated_instances(PATH_TO_ANNOTATED_DATA + 'olida_annotated_data/olida_testset_annotated.json')
        olida_instances.update(load_olida_annotated_instances(PATH_TO_ANNOTATED_DATA + 'olida_annotated_data/olida_trainset_annotated.json'))
        if oli_id in olida_instances.keys():
            oli_instance = olida_instances[oli_id]
            patient_variants.update(oli_instance)
            output_filename += '_' + oli_id
        else:
            print('The OLIDA ID you entered is not in the annotated set.')
            print('Script is terminated')
            quit()

    full_prioritization(patient_variants, patient_pairs, seeds, output_filename, restart, write_raw, nb_tops)
    print('### Done!')
    return None
