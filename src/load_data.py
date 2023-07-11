import configparser
import csv
import pickle
import json

config = configparser.ConfigParser()
config.read('config.ini')
PATH_TO_DATA = config['BASE']['PATH_TO_DATA']
PATH_TO_ANNOTATED_DATA = config['BASE']['PATH_TO_ANNOTATED_DATA']


def load_combs(file, chunksize):
    """Load combinations predicted with VarCoPP in chunks"""
    reader = csv.reader(file, delimiter='\t')
    comb_chunk = []
    for i, line in enumerate(reader):
        if i % chunksize == 0 and i > 0:
            yield comb_chunk
            comb_chunk = []
        comb_chunk.append(line)
    yield comb_chunk


def load_patient_annotated_data(patient_name, folder):
    """Load pre-annotated variants and pairs for a patient"""
    with open(PATH_TO_ANNOTATED_DATA + folder + '/' + patient_name + '_annotated_variants.json', 'r') as file:
        variants = json.load(file)
    with open(PATH_TO_ANNOTATED_DATA + folder + '/' + patient_name + '_annotated_pairs.json', 'r') as file:
        pairs = json.load(file)
    return variants, pairs


def load_olida_annotated_instances(fileName):
    """Load the dict of annotated positive instances """
    with open(fileName, 'r') as file:
        olida_annotated = json.load(file)
    return olida_annotated


def load_seeds(fileName):
    combinations_seeds = {}
    with open(PATH_TO_DATA + fileName, 'r') as file:
        file.readline()
        for line in file:
            line = line.strip('\n').split('\t')
            if len(line) > 1:
                combinations_seeds[line[0]] = line[1].split(',')
            else:
                combinations_seeds[line[0]] = []
    return combinations_seeds


def load_hpo_panel_seeds(file1, file2):
    combinations_seeds = {}
    with open(PATH_TO_DATA + file1, 'r') as file:
        file.readline()
        for line in file:
            line = line.strip('\n').split('\t')
            combinations_seeds[line[0]] = line[1].split(',')
    with open(PATH_TO_DATA + file2, 'r') as file:
        file.readline()
        for line in file:
            line = line.strip('\n').split('\t')
            if len(line) > 1:
                combinations_seeds[line[0]] += line[1].split(',')
            else:
                combinations_seeds[line[0]] += []
    combinations_seeds = {k:[i for i in v if i != ''] for k,v in combinations_seeds.items()}
    return combinations_seeds


def load_minmax_varcopp(patient_name, fold_nb):
    path = PATH_TO_DATA + 'varcopp_data/varcopp_results/'
    if fold_nb is None:
        with open(path + patient_name + 'minmax.p', 'rb') as file:
            return pickle.load(file)
    else:
        with open(path + 'fold_' + str(fold_nb) + '/' + patient_name + 'minmax.p', 'rb') as file:
            return pickle.load(file)
