
import random
from imblearn.ensemble import BalancedRandomForestClassifier
from .prioritization import *


def get_annotated_instances(neutral_ids, positive_ids):
    neutral_ids = {n:1 for n in neutral_ids}
    neutral_vectors = []
    with open(PATH_TO_DATA + 'varcopp_data/training_data/1KGP_training_set_imputed_annotated.txt', 'r') as file:
        for line in file:
            line = line.strip('\n').split('\t')
            if neutral_ids.get(line[0]) is not None:
                neutral_vectors.append([float(el) for el in line[1::]])
    positive_vectors = []
    with open(PATH_TO_DATA + 'varcopp_data/training_data/OLIDA_training_set_annotated.txt', 'r') as file:
        for line in file:
            line = line.strip('\n').split('\t')
            if line[0] in positive_ids:
                positive_vectors.append([float(el) for el in line[1::]])

    return positive_vectors, neutral_vectors


def invert_dict(my_dict):
    new_dict = {}
    for k, v in my_dict.items():
        if isinstance(v, list):
            v = tuple(v)
        new_dict[v] = new_dict.get(v, []) + [k]
    return new_dict


def split_in_folds(nb_folds, comb_ids_to_gp):
    random.seed(505)  # For splitting the dataset consistently
    folds = [[] for _ in range(nb_folds)]
    gp_comb_ids = invert_dict(comb_ids_to_gp)
    gp_comb_ids = {gp: gp_comb_ids[gp] for _, gp in sorted(zip([len(comb_ids) for comb_ids in gp_comb_ids.values()], list(gp_comb_ids.keys())), reverse=True)}
    i = 0
    for gp in gp_comb_ids:
        folds[i] += gp_comb_ids[gp]
        i += 1
        if i % 10 == 0:
            folds = [f for _,f in sorted(zip([len(f2) for f2 in folds], [f3 for f3 in folds]))]
            i = 0
    return folds


def read_comb_ids_to_gp(fileName):
    res = {}
    with open(fileName, 'r') as file:
        file.readline()
        for line in file:
            line = line.strip('\n').split('\t')
            if len(line) == 3:
                res[line[0]] = (line[1], line[2])
            else:
                print(line)

    return res


def train_varcopp_model(positive_train, neutral_train):
    positive_vectors, neutral_vectors = get_annotated_instances(neutral_train, positive_train)
    brf = BalancedRandomForestClassifier(n_estimators=400,
                                         # sampling_strategy= "all",
                                         random_state=1)
    varcopp_model = brf.fit(positive_vectors + neutral_vectors, [1 for i in range(len(positive_train))] +
                            [0 for i in range(len(neutral_train))])

    return varcopp_model


def train_crossval_models():
    positives = read_comb_ids_to_gp(PATH_TO_DATA + 'varcopp_data/training_data/OLIDA_train_gene_pairs.txt')
    negatives = read_comb_ids_to_gp(PATH_TO_DATA + 'varcopp_data/training_data/1KGP_train_gene_pairs.txt')
    folds_positives = split_in_folds(10, positives)
    folds_negatives = split_in_folds(10, negatives)
    for i in range(10):
        pos_test = folds_positives[i]
        pos_train = sum([folds_positives[j] for j in range(10) if j != i], [])
        neut_train = sum([folds_negatives[j] for j in range(10) if j != i], [])
        varcopp_model = train_varcopp_model(pos_train, neut_train)
        with open(PATH_TO_DATA + 'varcopp_data/varcopp_models/varcopp_fold_' + str(i) + '.p', 'wb') as file:
            pickle5.dump(varcopp_model, file)
        if 'crossval_folds' not in os.listdir('data'):
            os.mkdir(PATH_TO_DATA + 'crossval_folds')
        with open(PATH_TO_DATA + 'crossval_folds/olida_combs_fold_' + str(i) + '.p', 'wb') as file:
            pickle.dump(pos_test, file)

    return None

