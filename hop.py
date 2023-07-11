import os

from src import pipeline_single_exome, pipeline_crossval_per_patient, pipeline_indepval_per_patient
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description='Type of usage of the HOP predictor')
    parser.add_argument('action', nargs='?', action="store")
    ### Arguments for all options
    parser.add_argument('-p', '--patient_name', help='Name of the template patient')
    parser.add_argument('-rf', '--results_folder', help='results folder')
    parser.add_argument('-r', '--restart', help='Value of the restart parameter',default=0.3)
    parser.add_argument('-n', '--nb_top', help='Number of top combinations you want to keep as top', default=100)
    parser.add_argument('-wr', '--write_raw', help='Whether or not to write the raw results', action='store_true', default=False)

    ### Arguments for crossvalidation and indep
    parser.add_argument('-o', '--all_operators',
                        help="Whether the different operators for combining the final score should be computed", action='store_true', default=False)
    parser.add_argument('-st', '--seed_type', choices=['All', 'HPO', 'Panel', 'HPO+Panel'],
                        help='Which type of seeds to use to compute the disease-relevance score, specify either HPO, Panel, HPO+Panel or All',
                        default='All')
    parser.add_argument('-fs', '--final_score',
                        choices=[a + '_' + b for a in ['HPO', 'Panel', 'HPO+Panel'] for b in ['MAX', 'MIN', 'AVGE', 'MULT']],
                        default='HPO+Panel_AVGE',
                        help='Which type of seeds to use to compute the disease-relevance score, specify either HPO, Panel, HPO+Panel or All')

    ### Arguments for prioritization of a single exome
    parser.add_argument('-c', '--combination', help='OLIDA combination to insert in the patient exome')
    parser.add_argument('-s', '--seeds', help='List of seeds separated by a comma')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if args.action == "crossval":
        pipeline_crossval_per_patient.run_crossval_priorization_synthetic_exome_from_template(args.patient_name, args.results_folder,
                                                                                              args.restart, args.all_operators,
                                                                                              args.seed_type, args.write_raw, args.nb_top,
                                                                                              args.final_score)
    elif args.action == "indep":
        pipeline_indepval_per_patient.run_indep_priorization_synthetic_exome_from_template(args.patient_name, args.results_folder,
                                                                                           args.restart,
                                                                                           args.seed_type,
                                                                                           args.write_raw, args.nb_top,
                                                                                           args.final_score)
    elif args.action == "predict":
        pipeline_single_exome.prioritize_single_exome(args.patient_name, args.results_folder, args.restart, args.seeds,
                                                      args.combination,args.write_raw, args.nb_top)
    else:
        raise Exception(
            f"Unknown action {args.action}. Command should take the form: predictor.py <crossval independent predict> [OPTIONS]")

    return None


if __name__ == "__main__":
    main()



