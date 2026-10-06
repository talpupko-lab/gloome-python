import argparse
import traceback

from sys import exit
from typing import Dict

from gloome.utils import *


class Config:
    def __init__(self, **attributes):
        self.MODE = MODE[3:4]
        self.COMMAND_LINE = COMMAND_LINE

        self.RESULTS_DIR = RESULTS_DIR
        self.LOGS_DIR = LOGS_DIR

        self.IN_DIR = IN_DIR
        self.OUT_DIR = OUT_DIR

        self.ACTIONS = ACTIONS
        self.VALIDATION_ACTIONS = VALIDATION_ACTIONS
        self.DEFAULT_ACTIONS = DEFAULT_ACTIONS
        self.MAIN_ACTIONS = MAIN_ACTIONS

        self.CURRENT_ARGS = DEFAULT_ARGUMENTS

        self.CALCULATED_ARGS = CALCULATED_ARGS

        self.WEBSERVER_NAME_CAPITAL = WEBSERVER_NAME_CAPITAL

        self.PROCESS_ID = None
        self.MSA_FILE = None
        self.TREE_FILE = None
        self.JOB_LOGGER = None

        self.USAGE = USAGE

        if attributes:
            for key, value in attributes.items():
                if key == 'PROCESS_ID':
                    self.change_process_id(value)
                else:
                    setattr(self, key, value)
        if len(self.COMMAND_LINE) > 4 and self.COMMAND_LINE[1].startswith('-') and self.COMMAND_LINE[3].startswith(
                '-'):
            self.parse_arguments()

        if not self.PROCESS_ID:
            self.change_process_id(get_new_process_id())

    def change_process_id(self, process_id: str):
        self.PROCESS_ID = process_id

        self.RESULTS_DIR = RESULTS_DIR
        self.LOGS_DIR = LOGS_DIR

        self.OUT_DIR = self.OUT_DIR.joinpath(self.PROCESS_ID)
        check_dir(self.OUT_DIR)
        self.IN_DIR = self.IN_DIR.joinpath(self.PROCESS_ID)

        self.CALCULATED_ARGS.file_path = self.OUT_DIR

        self.JOB_LOGGER = get_job_logger(f'{process_id}', self.LOGS_DIR)

    def check_and_set_input_and_output_variables(self):
        """get variables from input arguments and fill out the Variable Class properties"""
        if len(self.COMMAND_LINE) < 5:
            print('\tAt least two required parameters --msa_file --tree_file', self.USAGE, sep='\n')
            exit()

        if len(self.COMMAND_LINE) > 4 and self.COMMAND_LINE[1].startswith('-') and self.COMMAND_LINE[3].startswith('-'):
            if not self.check_arguments_for_errors():
                for error_type, error in self.CALCULATED_ARGS.err_list:
                    print(f'{error_type}: {error}')
                print(self.USAGE)
                exit()

    def get_selected_files(self) -> Dict[str, bool]:
        selected_files = {'file_interactive_tree_html': self.CURRENT_ARGS.file_interactive_tree_html,
                          'file_newick_tree_png': self.CURRENT_ARGS.file_newick_tree_png,
                          'file_table_of_coevolution_tsv': self.CURRENT_ARGS.file_table_of_coevolution_tsv,
                          'file_simulated_datasets_fastas': self.CURRENT_ARGS.file_simulated_datasets_fastas,
                          'file_barplot_of_correlation_svg': self.CURRENT_ARGS.file_barplot_of_correlation_svg,
                          'file_plot_distribution_of_correlation_svg':
                              self.CURRENT_ARGS.file_plot_distribution_of_correlation_svg,
                          'file_plot_correlation_by_rate_bin_svg':
                              self.CURRENT_ARGS.file_plot_correlation_by_rate_bin_svg,
                          'file_table_of_posterior_rates_tsv': self.CURRENT_ARGS.file_table_of_posterior_rates_tsv,
                          'file_table_of_pearson_correlation_tsv':
                              self.CURRENT_ARGS.file_table_of_pearson_correlation_tsv,
                          'file_table_of_nodes_tsv': self.CURRENT_ARGS.file_table_of_nodes_tsv,
                          'file_branch_position_probabilities_tsv':
                              self.CURRENT_ARGS.file_branch_position_probabilities_tsv,
                          'file_table_of_branches_tsv': self.CURRENT_ARGS.file_table_of_branches_tsv,
                          'file_log_likelihood_tsv': self.CURRENT_ARGS.file_log_likelihood_tsv,
                          'file_table_of_attributes_tsv': self.CURRENT_ARGS.file_table_of_attributes_tsv,
                          'file_table_of_parsimony_and_homoplasy_scores_tsv':
                              self.CURRENT_ARGS.file_table_of_parsimony_and_homoplasy_scores_tsv,
                          'file_phylogenetic_tree_nwk': self.CURRENT_ARGS.file_phylogenetic_tree_nwk}
        return selected_files

    def execute_action(self, func, *args, **kwargs):
        try:
            val = func(*args, **kwargs)
            if val:
                self.JOB_LOGGER.info(f'\n\tSuccessfully Command \'{func.__name__}\' executed successfully. -> {val}\n')
        except ValueError:
            format_exc = f'{traceback.format_exc()}'
            self.JOB_LOGGER.info(f'Failed to execute the command \'{func.__name__}\' -> {format_exc}')
            self.CALCULATED_ARGS.err_list.append((f'Failed to execute the command \'{func.__name__}\'', format_exc))

    def execute_calculation(self):
        if not self.CALCULATED_ARGS.err_list and self.DEFAULT_ACTIONS.get('calculate_tree', False):
            self.execute_action(self.ACTIONS.calculate_tree, self.CALCULATED_ARGS.newick_tree)
        if not self.CALCULATED_ARGS.err_list and self.DEFAULT_ACTIONS.get('calculate_ancestral_sequence', False):
            self.execute_action(self.ACTIONS.calculate_ancestral_sequence, self.CALCULATED_ARGS.newick_tree)
        if not self.CALCULATED_ARGS.err_list and self.DEFAULT_ACTIONS.get('calculate_correlation', False):
            self.execute_action(self.ACTIONS.calculate_correlation, self.CALCULATED_ARGS.newick_tree,
                                probability_lg=self.CURRENT_ARGS.probability_lg,
                                number_lg=self.CURRENT_ARGS.number_lg)
        if not self.CALCULATED_ARGS.err_list and self.DEFAULT_ACTIONS.get('execute_all_actions', False):
            self.execute_action(self.ACTIONS.execute_all_actions,
                                file_path=self.OUT_DIR,
                                newick_tree=self.CALCULATED_ARGS.newick_tree,
                                with_internal_nodes=self.CURRENT_ARGS.with_internal_nodes,
                                log_file=self.JOB_LOGGER.handlers[-1].baseFilename,
                                actions=self.MAIN_ACTIONS,
                                selected_files=self.get_selected_files(),
                                use_copap=not self.CURRENT_ARGS.is_do_not_use_copap,
                                probability_lg=self.CURRENT_ARGS.probability_lg,
                                number_lg=self.CURRENT_ARGS.number_lg,
                                number_datasets=self.CURRENT_ARGS.number_datasets)

    def check_arguments_for_errors(self) -> bool:
        if self.TREE_FILE.is_file():
            with open(self.TREE_FILE, 'r') as f:
                if self.VALIDATION_ACTIONS.get('del_bootstrap_values', False):
                    self.CALCULATED_ARGS.newick_text = self.ACTIONS.del_bootstrap_values(f.read().strip())
                else:
                    self.CALCULATED_ARGS.newick_text = f.read().strip()
        else:
            self.CALCULATED_ARGS.err_list.append((f'The File does not exist',
                                                  f'File "{self.TREE_FILE}" does not exist '))

        if self.MSA_FILE.is_file():
            with open(self.MSA_FILE, 'r') as f:
                self.CALCULATED_ARGS.msa = f.read().strip()
        else:
            self.CALCULATED_ARGS.err_list.append((f'The File does not exist',
                                                  f'File "{self.MSA_FILE}" does not exist '))

        if not self.CALCULATED_ARGS.err_list and self.VALIDATION_ACTIONS.get('check_data', False):
            self.CALCULATED_ARGS.err_list += (
                self.ACTIONS.check_data(self.CALCULATED_ARGS.newick_text,
                                        self.CALCULATED_ARGS.msa,
                                        self.CURRENT_ARGS.categories_quantity,
                                        self.CURRENT_ARGS.alpha,
                                        self.CURRENT_ARGS.pi_1,
                                        self.CURRENT_ARGS.coefficient_bl,
                                        self.CURRENT_ARGS.probability_lg,
                                        self.CURRENT_ARGS.number_lg,
                                        self.CURRENT_ARGS.number_datasets,
                                        self.CURRENT_ARGS.is_optimize_pi,
                                        self.CURRENT_ARGS.is_optimize_pi_average,
                                        self.CURRENT_ARGS.is_optimize_alpha,
                                        self.CURRENT_ARGS.is_optimize_bl,
                                        self.CURRENT_ARGS.is_do_not_use_copap,
                                        self.CURRENT_ARGS.file_interactive_tree_html,
                                        self.CURRENT_ARGS.file_newick_tree_png,
                                        self.CURRENT_ARGS.file_table_of_coevolution_tsv,
                                        self.CURRENT_ARGS.file_simulated_datasets_fastas,
                                        self.CURRENT_ARGS.file_barplot_of_correlation_svg,
                                        self.CURRENT_ARGS.file_plot_distribution_of_correlation_svg,
                                        self.CURRENT_ARGS.file_plot_correlation_by_rate_bin_svg,
                                        self.CURRENT_ARGS.file_table_of_posterior_rates_tsv,
                                        self.CURRENT_ARGS.file_table_of_pearson_correlation_tsv,
                                        self.CURRENT_ARGS.file_table_of_nodes_tsv,
                                        self.CURRENT_ARGS.file_branch_position_probabilities_tsv,
                                        self.CURRENT_ARGS.file_table_of_branches_tsv,
                                        self.CURRENT_ARGS.file_log_likelihood_tsv,
                                        self.CURRENT_ARGS.file_table_of_attributes_tsv,
                                        self.CURRENT_ARGS.file_table_of_parsimony_and_homoplasy_scores_tsv,
                                        self.CURRENT_ARGS.file_phylogenetic_tree_nwk,
                                        self.CURRENT_ARGS.rooting_method,
                                        self.CURRENT_ARGS.leaf))

        if not self.CALCULATED_ARGS.err_list and self.VALIDATION_ACTIONS.get('set_root', False):
            try:
                self.CALCULATED_ARGS.newick_text = self.ACTIONS.set_root(self.CALCULATED_ARGS.newick_text,
                                                                         self.CURRENT_ARGS.rooting_method,
                                                                         self.CURRENT_ARGS.leaf)
            except ValueError:
                self.CALCULATED_ARGS.err_list.append((f'TREE error', f'Failed to set root.'))

        if not self.CALCULATED_ARGS.err_list and self.VALIDATION_ACTIONS.get('check_tree', False):
            try:
                self.CALCULATED_ARGS.newick_tree = self.ACTIONS.check_tree(self.CALCULATED_ARGS.newick_text)
            except ValueError:
                self.CALCULATED_ARGS.err_list.append((f'TREE error', f'Wrong Phylogenetic tree format.'))

        if not self.CALCULATED_ARGS.err_list and self.DEFAULT_ACTIONS.get('set_tree_data', False):
            try:
                self.ACTIONS.set_tree_data(self.CALCULATED_ARGS.newick_tree, msa=self.CALCULATED_ARGS.msa,
                                           categories_quantity=self.CURRENT_ARGS.categories_quantity,
                                           alpha=self.CURRENT_ARGS.alpha, pi_1=self.CURRENT_ARGS.pi_1,
                                           coefficient_bl=self.CURRENT_ARGS.coefficient_bl,
                                           is_optimize_pi=self.CURRENT_ARGS.is_optimize_pi,
                                           is_optimize_pi_average=self.CURRENT_ARGS.is_optimize_pi_average,
                                           is_optimize_alpha=self.CURRENT_ARGS.is_optimize_alpha,
                                           is_optimize_bl=self.CURRENT_ARGS.is_optimize_bl)
            except ValueError:
                self.CALCULATED_ARGS.err_list.append((f'MSA error',
                                                      f'Wrong MSA format. Please provide MSA in FASTA format.'))

        if not self.CALCULATED_ARGS.err_list and (self.CURRENT_ARGS.is_optimize_pi_average or
                                                  self.CURRENT_ARGS.is_optimize_pi):
            try:
                self.CURRENT_ARGS.pi_1 = self.CALCULATED_ARGS.newick_tree.pi_1
            except ValueError:
                self.CALCULATED_ARGS.err_list.append((f'Strange error',
                                                      f'Strange error.'))

        if not self.CALCULATED_ARGS.err_list and self.CURRENT_ARGS.is_optimize_alpha:
            try:
                self.CURRENT_ARGS.alpha = self.CALCULATED_ARGS.newick_tree.alpha
            except ValueError:
                self.CALCULATED_ARGS.err_list.append((f'Strange error',
                                                      f'Strange error.'))

        if not self.CALCULATED_ARGS.err_list and self.CURRENT_ARGS.is_optimize_bl:
            try:
                self.CURRENT_ARGS.coefficient_bl = self.CALCULATED_ARGS.newick_tree.coefficient_bl
            except ValueError:
                self.CALCULATED_ARGS.err_list.append((f'Strange error',
                                                      f'Strange error.'))

        if self.CALCULATED_ARGS.err_list:
            self.JOB_LOGGER.info(f'Error list: \n{self.CALCULATED_ARGS.err_list}')
        else:
            self.JOB_LOGGER.info(f'Verification completed successfully\n')

        return not self.CALCULATED_ARGS.err_list

    def enable_default_actions(self):
        self.DEFAULT_ACTIONS.update({
            'calculate_tree': False,
            'calculate_ancestral_sequence': False,
            'calculate_correlation': False,
            'execute_all_actions': False
        })
        self.MAIN_ACTIONS.update({
            'compute_likelihood_of_tree': False,
            'draw_tree': False,
            'create_all_file_types': False
        })

        if 'compute_likelihood_of_tree' in self.MODE:
            self.DEFAULT_ACTIONS.update({'execute_all_actions': True})
            self.MAIN_ACTIONS.update({'compute_likelihood_of_tree': True})
        if 'draw_tree' in self.MODE:
            self.DEFAULT_ACTIONS.update({'calculate_tree': True,
                                         'calculate_ancestral_sequence': True,
                                         'execute_all_actions': True})
            self.MAIN_ACTIONS.update({'draw_tree': True})
        if 'create_all_file_types' in self.MODE or 'execute_all_actions' in self.MODE:
            self.DEFAULT_ACTIONS.update({'calculate_tree': True,
                                         'calculate_ancestral_sequence': True,
                                         'execute_all_actions': True})
            self.MAIN_ACTIONS.update({'create_all_file_types': True})
        if 'execute_all_actions' in self.MODE:
            self.DEFAULT_ACTIONS.update({'calculate_tree': True,
                                         'calculate_ancestral_sequence': True,
                                         'execute_all_actions': True})
            self.MAIN_ACTIONS.update({'compute_likelihood_of_tree': True,
                                      'draw_tree': True,
                                      'create_all_file_types': True})
        if not self.CURRENT_ARGS.is_do_not_use_copap:
            self.DEFAULT_ACTIONS.update({'calculate_correlation': True})

    def parse_arguments(self):
        """parse arguments and fill out the relevant Variable Class properties"""
        parser = argparse.ArgumentParser(prog=self.WEBSERVER_NAME_CAPITAL, description='GLOOME',
                                         usage='%(prog)s [options]')
        parser.add_argument('--msa_file', dest='msa_file', type=str, required=True,
                            help=f'Specify the msa filepath (required).')
        parser.add_argument('--tree_file', dest='tree_file', type=str, required=True,
                            help=f'Specify the newick filepath (required).')
        parser.add_argument('--out_dir', dest='out_dir', type=str, required=False, default=self.OUT_DIR,
                            help=f'Specify the outdir path (optional).')
        parser.add_argument('--process_id', dest='process_id', type=str, required=False, default=self.PROCESS_ID,
                            help=f'Specify a process ID or it will be generated automatically (optional).')
        parser.add_argument('--mode', dest='mode', required=False, action="extend", nargs="+", type=str,
                            help=f'Specify execution mode (optional). Possible options: ("draw_tree", '
                            f'"compute_likelihood_of_tree", "create_all_file_types", "execute_all_actions"). '
                            f'Default is {self.MODE[3:4]}.')
        parser.add_argument('--with_internal_nodes', dest='with_internal_nodes', type=int, required=False,
                            default=int(self.CURRENT_ARGS.with_internal_nodes), help=f'Specify the tree has internal '
                            f'nodes (optional). Default is {int(self.CURRENT_ARGS.with_internal_nodes)}.')
        parser.add_argument('--categories_quantity', dest='categories_quantity', type=int, required=False,
                            default=self.CURRENT_ARGS.categories_quantity, help=f'Specify categories quantity '
                            f'(optional). Default is {self.CURRENT_ARGS.categories_quantity}.')
        parser.add_argument('--alpha', dest='alpha', type=float, required=False, default=self.CURRENT_ARGS.alpha,
                            help=f'Specify alpha (optional). Default is {self.CURRENT_ARGS.alpha}.')
        parser.add_argument('--pi_1', dest='pi_1', type=float, required=False, default=self.CURRENT_ARGS.pi_1,
                            help=f'Specify pi_1 (optional). Default is {self.CURRENT_ARGS.pi_1}.')
        parser.add_argument('--coefficient_bl', dest='coefficient_bl', type=float, required=False,
                            help=f'Specify coefficient_bl (optional). Default is {self.CURRENT_ARGS.coefficient_bl}.',
                            default=self.CURRENT_ARGS.coefficient_bl)
        parser.add_argument('--probability_lg', dest='probability_lg', type=float, required=False,
                            help=f'Specify probability_lg (optional). Default is {self.CURRENT_ARGS.probability_lg}.',
                            default=self.CURRENT_ARGS.probability_lg)
        parser.add_argument('--number_lg', dest='number_lg', type=int, required=False,
                            help=f'Specify number_lg (optional). Default is {self.CURRENT_ARGS.number_lg}.',
                            default=self.CURRENT_ARGS.number_lg)
        parser.add_argument('--number_datasets', dest='number_datasets', type=int, required=False,
                            help=f'Specify number_datasets (optional). Default is {self.CURRENT_ARGS.number_datasets}.',
                            default=self.CURRENT_ARGS.number_datasets)
        parser.add_argument('--rooting_method', dest='rooting_method', type=str, required=False,
                            default=self.CURRENT_ARGS.rooting_method, help=f'Specify tree rooting method (optional). '
                            f'Possible options: ("mad", "mvr", "midpoint", "outgroup"). '
                            f'mad - Minimal Ancestor Deviation; '
                            f'mvr - Minimum Variance Rooting; '
                            f'midpoint - Midpoint Rooting; '
                            f'outgroup - Outgroup Rooting; '
                            f'Default is "{self.CURRENT_ARGS.rooting_method}"')
        parser.add_argument('--leaf', dest='leaf', type=str, required=False, default=self.CURRENT_ARGS.leaf,
                            help=f'Specify leaf for outgroup rooting (optional). Default is "{self.CURRENT_ARGS.leaf}"')
        parser.add_argument('--is_optimize_pi', dest='is_optimize_pi', type=int, required=False,
                            help=f'Specify is_optimize_pi (optional). Default is '
                            f'{int(self.CURRENT_ARGS.is_optimize_pi)}.', default=int(self.CURRENT_ARGS.is_optimize_pi))
        parser.add_argument('--is_optimize_pi_average', dest='is_optimize_pi_average', type=int, required=False,
                            help=f'Specify is_optimize_pi_average (optional). Default is '
                            f'{int(self.CURRENT_ARGS.is_optimize_pi_average)}.',
                            default=int(self.CURRENT_ARGS.is_optimize_pi_average))
        parser.add_argument('--is_optimize_alpha', dest='is_optimize_alpha', type=int, required=False,
                            help=f'Specify is_optimize_alpha (optional). Default is '
                            f'{int(self.CURRENT_ARGS.is_optimize_alpha)}.',
                            default=int(self.CURRENT_ARGS.is_optimize_alpha))
        parser.add_argument('--is_optimize_bl', dest='is_optimize_bl', type=int, required=False,
                            help=f'Specify is_optimize_bl (optional). Default is '
                            f'{int(self.CURRENT_ARGS.is_optimize_bl)}.', default=int(self.CURRENT_ARGS.is_optimize_bl))
        parser.add_argument('--is_do_not_use_copap', dest='is_do_not_use_copap', type=int, required=False,
                            help=f'Specify is_do_not_use_copap (optional). Default is '
                            f'{int(self.CURRENT_ARGS.is_do_not_use_copap)}.',
                            default=int(self.CURRENT_ARGS.is_do_not_use_copap))
        parser.add_argument('--file_interactive_tree_html', dest='file_interactive_tree_html', type=int, required=False,
                            help=f'Specify file_interactive_tree_html (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_interactive_tree_html)}.',
                            default=int(self.CURRENT_ARGS.file_interactive_tree_html))
        parser.add_argument('--file_newick_tree_png', dest='file_newick_tree_png', type=int, required=False,
                            help=f'Specify file_newick_tree_png (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_newick_tree_png)}.',
                            default=int(self.CURRENT_ARGS.file_newick_tree_png))
        parser.add_argument('--file_table_of_coevolution_tsv',
                            dest='file_table_of_coevolution_tsv', type=int, required=False,
                            help=f'Specify file_table_of_coevolution_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_table_of_coevolution_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_table_of_coevolution_tsv))
        parser.add_argument('--file_simulated_datasets_fastas',
                            dest='file_simulated_datasets_fastas', type=int, required=False,
                            help=f'Specify file_simulated_datasets_fastas (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_simulated_datasets_fastas)}.',
                            default=int(self.CURRENT_ARGS.file_simulated_datasets_fastas))
        parser.add_argument('--file_barplot_of_correlation_svg',
                            dest='file_barplot_of_correlation_svg', type=int, required=False,
                            help=f'Specify file_barplot_of_correlation_svg (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_barplot_of_correlation_svg)}.',
                            default=int(self.CURRENT_ARGS.file_barplot_of_correlation_svg))
        parser.add_argument('--file_plot_distribution_of_correlation_svg',
                            dest='file_plot_distribution_of_correlation_svg', type=int, required=False,
                            help=f'Specify file_plot_distribution_of_correlation_svg (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_plot_distribution_of_correlation_svg)}.',
                            default=int(self.CURRENT_ARGS.file_plot_distribution_of_correlation_svg))
        parser.add_argument('--file_plot_correlation_by_rate_bin_svg',
                            dest='file_plot_correlation_by_rate_bin_svg', type=int, required=False,
                            help=f'Specify file_plot_correlation_by_rate_bin_svg (optional). Default is'
                            f' {int(self.CURRENT_ARGS.file_plot_correlation_by_rate_bin_svg)}.',
                            default=int(self.CURRENT_ARGS.file_plot_correlation_by_rate_bin_svg))
        parser.add_argument('--file_table_of_posterior_rates_tsv',
                            dest='file_table_of_posterior_rates_tsv', type=int, required=False,
                            help=f'Specify file_table_of_posterior_rates_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_table_of_posterior_rates_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_table_of_posterior_rates_tsv))
        parser.add_argument('--file_table_of_pearson_correlation_tsv',
                            dest='file_table_of_pearson_correlation_tsv', type=int, required=False,
                            help=f'Specify file_table_of_pearson_correlation_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_table_of_pearson_correlation_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_table_of_pearson_correlation_tsv))
        parser.add_argument('--file_table_of_nodes_tsv', dest='file_table_of_nodes_tsv', type=int, required=False,
                            help=f'Specify file_table_of_nodes_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_table_of_nodes_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_table_of_nodes_tsv))
        parser.add_argument('--file_branch_position_probabilities_tsv', type=int, required=False,
                            dest='file_branch_position_probabilities_tsv', help=f'Specify '
                            f'file_branch_position_probabilities_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_branch_position_probabilities_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_branch_position_probabilities_tsv))
        parser.add_argument('--file_table_of_branches_tsv', dest='file_table_of_branches_tsv', type=int, required=False,
                            help=f'Specify file_table_of_branches_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_table_of_branches_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_table_of_branches_tsv))
        parser.add_argument('--file_log_likelihood_tsv', dest='file_log_likelihood_tsv', type=int, required=False,
                            help=f'Specify file_log_likelihood_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_log_likelihood_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_log_likelihood_tsv))
        parser.add_argument('--file_table_of_attributes_tsv', dest='file_table_of_attributes_tsv', type=int,
                            required=False, help=f'Specify file_table_of_attributes_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_table_of_attributes_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_table_of_attributes_tsv))
        parser.add_argument('--file_table_of_parsimony_and_homoplasy_scores_tsv', type=int, required=False,
                            dest='file_table_of_parsimony_and_homoplasy_scores_tsv',
                            help=f'Specify file_table_of_parsimony_and_homoplasy_scores_tsv (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_table_of_parsimony_and_homoplasy_scores_tsv)}.',
                            default=int(self.CURRENT_ARGS.file_table_of_parsimony_and_homoplasy_scores_tsv))
        parser.add_argument('--file_phylogenetic_tree_nwk', dest='file_phylogenetic_tree_nwk', type=int,
                            required=False, help=f'Specify file_phylogenetic_tree_nwk (optional). Default is '
                            f'{int(self.CURRENT_ARGS.file_phylogenetic_tree_nwk)}.',
                            default=int(self.CURRENT_ARGS.file_phylogenetic_tree_nwk))

        args = parser.parse_args()

        for arg_name, arg_value in vars(args).items():
            if arg_value is not None:
                if arg_name == 'process_id':
                    if arg_value != self.PROCESS_ID:
                        self.change_process_id(arg_value)
                elif arg_name == 'mode':
                    setattr(self, arg_name.upper(), arg_value)
                elif arg_name in ('msa_file', 'tree_file', 'out_dir'):
                    setattr(self, arg_name.upper(), Path(arg_value))
                elif arg_name in ('with_internal_nodes',
                                  'is_optimize_pi',
                                  'is_optimize_pi_average',
                                  'is_optimize_alpha',
                                  'is_optimize_bl',
                                  'is_do_not_use_copap',
                                  'file_interactive_tree_html',
                                  'file_newick_tree_png',
                                  'file_table_of_coevolution_tsv',
                                  'file_simulated_datasets_fastas',
                                  'file_barplot_of_correlation_svg',
                                  'file_plot_distribution_of_correlation_svg',
                                  'file_plot_correlation_by_rate_bin_svg',
                                  'file_table_of_posterior_rates_tsv',
                                  'file_table_of_pearson_correlation_tsv',
                                  'file_table_of_nodes_tsv',
                                  'file_branch_position_probabilities_tsv',
                                  'file_table_of_branches_tsv',
                                  'file_log_likelihood_tsv',
                                  'file_table_of_attributes_tsv',
                                  'file_table_of_parsimony_and_homoplasy_scores_tsv',
                                  'file_phylogenetic_tree_nwk'):
                    if hasattr(self.CURRENT_ARGS, arg_name):
                        setattr(self.CURRENT_ARGS, arg_name, bool(arg_value))
                else:
                    if hasattr(self, arg_name.upper()):
                        setattr(self, arg_name.upper(), arg_value)
                    if hasattr(self.CURRENT_ARGS, arg_name):
                        setattr(self.CURRENT_ARGS, arg_name, arg_value)

        self.enable_default_actions()

        return args
