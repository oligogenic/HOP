from .load_data import load_patient_annotated_data, load_olida_annotated_instances, load_seeds, load_hpo_panel_seeds
import configparser

config = configparser.ConfigParser()
config.read('config.ini')
PATH_TO_DATA = config['BASE']['PATH_TO_DATA']
PATH_TO_ANNOTATED_DATA = config['BASE']['PATH_TO_ANNOTATED_DATA']



def input_preprocessing_all(restart, patient_name, olida_data, nb_tops):
    try:
        restart = float(restart)
    except ValueError:
        print('#### The restart value was not correctly specified. Please enter a value between 0 and 1')
        quit()
    if not (restart >= 0 and restart <= 1):
        print('#### The restart value was not correctly specified. Please enter a value between 0 and 1')
        quit()

    try:
        nb_tops = int(nb_tops)
    except ValueError:
        print('#### The nb of top combinations to return was not correctly specified. Please enter an integer value above 1.')
        quit()
    if nb_tops < 1:
        print(
            '#### The nb of top combinations to return was not correctly specified. Please enter an integer value above 1.')
        quit()
    # Load patient variants and pairs
    try:
        variants, pairs = load_patient_annotated_data(patient_name, folder='patients_annotated_data')
    except FileNotFoundError:
        print('##### Annotated patient file does not exist ' + patient_name)
        print('##### Check that the path in config.ini is correct and that you have downloaded the annotated data')
        quit()

    # Load annotated OLIDA combinations
    try:
        positive_instances = load_olida_annotated_instances(
            PATH_TO_ANNOTATED_DATA + 'olida_annotated_data/' + olida_data)
    except FileNotFoundError:
        print('### Annotated OLIDA data was not found.')
        print(
            '### Check that the path in config.ini is correct and that you have downloaded the annotated data')
        quit()
    return restart, nb_tops, variants, pairs, positive_instances


def input_preprocessing_crossval_indep(seed_type, compute_all_operators, final_score, crossval):
    # Load different types of seeds
    if crossval:
        hpopanel_annotations = load_hpo_panel_seeds('rwr_seeds/HPO_annotation_trainset.txt',
                                                    'rwr_seeds/Panel_annotation_trainset.txt')
        panel_annotations = load_seeds('rwr_seeds/Panel_annotation_trainset.txt')
        hpo_annotations = load_seeds('rwr_seeds/HPO_annotation_trainset.txt')
    else:
        hpopanel_annotations = load_hpo_panel_seeds('rwr_seeds/HPO_annotation_testset.txt',
                                                    'rwr_seeds/Panel_annotation_testset.txt')
        panel_annotations = load_seeds('rwr_seeds/Panel_annotation_testset.txt')
        hpo_annotations = load_seeds('rwr_seeds/HPO_annotation_testset.txt')
    if seed_type == 'All':
        seeds_dictionaries = {'hpos': hpo_annotations, 'panel': panel_annotations,
                              'hpopanel': hpopanel_annotations}
    elif seed_type == 'HPO':
        seeds_dictionaries = {'hpos': hpo_annotations}
    elif seed_type == 'Panel':
        seeds_dictionaries = {'panel': panel_annotations}
    elif seed_type == 'HPO+Panel':
        seeds_dictionaries = {'hpopanel': hpopanel_annotations}
    else:
        print('The seed type not correctly specified, options are "All", "HPO", "Panel" and "HPO+Panel" ')
        quit()
    ## Create file header
    output_file_header = ['Comb_id', 'Var_A', 'Var_B', 'GeneA', 'GeneB', 'PathogenicityScore', 'VarCoPP_class']
    for type in ['hpos', 'panel', 'hpopanel']:
        # Get the RWR results for the different seed types and record the max and min possible scores
        if seeds_dictionaries.get(type) is not None:
            output_file_header += ['DS_' + type, 'FS_avge_' + type]
            if compute_all_operators:
                output_file_header += ['FS_max_' + type, 'FS_min_' + type, 'FS_mult_' + type]
    input_to_colnum = {'Panel_MAX': 'FS_max_panel', 'Panel_MIN': 'FS_min_panel', 'Panel_AVGE': 'FS_avge_panel',
                       'Panel_MULT': 'FS_mult_panel', 'HPO+Panel_MAX': 'FS_max_hpopanel',
                       'HPO+Panel_MIN': 'FS_min_hpopanel', 'HPO+Panel_AVGE': 'FS_avge_hpopanel',
                       'HPO+Panel_MULT': 'FS_mult_hpopanel', 'HPO_MAX': 'FS_max_hpos', 'HPO_MIN': 'FS_min_hpos',
                       'HPO_AVGE': 'FS_avge_hpos', 'HPO_MULT': 'FS_mult_hpos'}
    try:
        index_col = output_file_header.index(input_to_colnum[final_score])
    except KeyError:
        print('The seeds and operator "' + final_score + '" is not correct.')
        print('Please write it in the form <seed_type>_<operator>.')
        quit()
    except ValueError:
        print('The seeds or operator "' + final_score + '" are not used for computation.')
        print('Please change input seed type, operators, or change the seeds/operators used for ranking.')
        quit()

    return seeds_dictionaries, output_file_header, index_col


