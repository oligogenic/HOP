"""Computes disease relevance score for all pairs in a patient's exome.
Assumes a pre-computed matrix present in the 'data' folder"""
import configparser
import pickle
import gzip
import os
from .rwr import *
from .compute_RWR_matrix import compute_rw_matrix

config = configparser.ConfigParser()
config.read('config.ini')
PATH_TO_DATA = config['BASE']['PATH_TO_DATA']

def run_rwr(matrix, seed_vector, c, epsilon, max_iters):
    rwr = RWR(matrix)
    results = rwr.compute(seed_vector, c, epsilon, max_iters)
    return results

def load_matrix():
    if 'rwr_data' not in os.listdir(PATH_TO_DATA) or 'final_wholeKG_matrix.p' not in os.listdir(PATH_TO_DATA + 'rwr_data'):
        compute_rw_matrix()
    with open(PATH_TO_DATA +'rwr_data/final_wholeKG_matrix.p', 'rb') as file:
        matrix = pickle.load(file)
    return matrix

def load_instances_to_index():
    with open(PATH_TO_DATA + 'rwr_data/final_instances_to_index.p', 'rb') as file:
        instances_to_index = pickle.load(file)
    with open(PATH_TO_DATA + 'rwr_data/all_genes.p', 'rb') as file:
        all_genes = pickle.load(file)
    return all_genes, instances_to_index


def create_seed_vector(seeds_list, instances_to_index):
    seeds_vector = np.zeros(len(instances_to_index.keys()))
    seed_names = [s for s in seeds_list if s in instances_to_index]
    for s in seed_names:
        seeds_vector[instances_to_index[s]] = 1 / len(seed_names)
    return seeds_vector


def get_gene_combinations(patient_name):
    with gzip.open(PATH_TO_DATA + patient_name +'_predicted_combinations.p.gz', 'rb') as f:
        results = pickle.load(f)
    gene_combinations = [c.geneA + ',' + c.geneB for c in results.values()]
    return list(set(gene_combinations))


def get_rwr_results(matrix, instances_to_index,c, epsilon, max_iters, seeds_list, all_genes):
    seeds_vector = create_seed_vector(seeds_list, instances_to_index)
    rwr = RWR(matrix)
    results = rwr.compute(seeds_vector, c, epsilon, max_iters)
    list_instances_to_index_keys = list(instances_to_index.keys())
    results = {list_instances_to_index_keys[i]: results[i] for i in range(len(instances_to_index.keys())) if
               all_genes.get(list_instances_to_index_keys[i]) is not None}
    return results


def convert_to_geneName(results):
    with open(PATH_TO_DATA +'rwr_data/ensg_to_genename.p', 'rb') as file:
        ensg_genename = pickle.load(file)
    new_results = {ensg_genename.get(k, ''): r for k,r in results.items()}
    return new_results


def get_rwr_gene_scores(seeds, restart):
    matrix = load_matrix()
    all_genes, instances_to_index = load_instances_to_index()
    if seeds != ['']:
        results = get_rwr_results(matrix, instances_to_index, float(restart), 1e-9, 100, seeds, all_genes)
        max_rwr = max(results.values())
        min_rwr = min(results.values())
        if max_rwr != min_rwr:
            normalized_results = {g: (r-min_rwr)/(max_rwr-min_rwr) for g, r in results.items()}
        else:
            normalized_results = {g: 0 for g in all_genes.keys()}
    else:
        normalized_results = {g: 0 for g in all_genes.keys()}
    return convert_to_geneName(normalized_results)


