from sys import argv
from typing import List, Tuple, Union
from types import FunctionType, MethodType
from pathlib import Path
from importlib.resources import files

from gloome.tree.tree import Tree
from gloome.services.service_functions import check_data, execute_all_actions

MODE = ['draw_tree', 'compute_likelihood_of_tree', 'create_all_file_types', 'execute_all_actions']
ROOTING_METHODS = [('mad', 'Minimal Ancestor Deviation'), ('mvr', 'Minimum Variance Rooting'),
                   ('midpoint', 'Midpoint Rooting'), ('outgroup', 'Outgroup Rooting')]
IS_PRODUCTION = True
MAX_CONTENT_LENGTH = 16 * 1000 * 1000 * 1000
PREFIX = '/'
APPLICATION_ROOT = PREFIX
DEBUG = not IS_PRODUCTION

WEBSERVER_NAME_CAPITAL = 'Gloome'

GLOOME = Path('/gloome')
BIN_DIR = GLOOME if GLOOME.exists() else Path.cwd()
RESULTS_DIR = BIN_DIR.joinpath('results')
IN_DIR = RESULTS_DIR.joinpath('in')
OUT_DIR = RESULTS_DIR.joinpath('out')
LOGS_DIR = BIN_DIR.joinpath('logs')
TMP_DIR = BIN_DIR.joinpath('tmp')

GLOOME_DIR = files('gloome')
DATA_DIR = GLOOME_DIR.joinpath('data')
INITIAL_DATA_DIR = DATA_DIR.joinpath('initial_data')

MSA_FILE_NAME = 'msa_file.msa'
TREE_FILE_NAME = 'tree_file.tree'


class Actions:

    def __init__(self, **attributes):
        if attributes:
            for key, value in attributes.items():
                if type(value) in (FunctionType, MethodType):
                    setattr(self, key, value)


class CalculatedArgs:
    err_list: List[Union[Tuple[str, ...], str]]

    def __init__(self, **attributes):
        self.err_list = []
        if attributes:
            for key, value in attributes.items():
                setattr(self, key, value)


class DefaultArgs:
    def __init__(self, **attributes):
        if attributes:
            for key, value in attributes.items():
                setattr(self, key, value)

    def get(self, attribute_name, default=None):
        if hasattr(self, attribute_name):
            return getattr(self, attribute_name)
        else:
            return default

    def update(self, *args, **kwargs) -> None:
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)
        if args:
            for arg in args:
                if isinstance(arg, dict):
                    for key, value in arg.items():
                        setattr(self, key, value)


COMMAND_LINE = argv

DEFAULT_FORM_ARGUMENTS = {
    'categories_quantity': 4,
    'alpha': 0.5,
    'pi_1': 0.5,
    'coefficient_bl': 1.0,
    'probability_lg': 0.5,
    'number_lg': 1,
    'number_datasets': 100,
    'rooting_method': ROOTING_METHODS[2][0],
    'rooting_methods': ROOTING_METHODS,
    'leaf': '',
    'leaves': [],
    'is_optimize_pi': True,
    'is_optimize_pi_average': False,
    'is_optimize_alpha': True,
    'is_optimize_bl': True,
    'is_do_not_use_copap': False,
    'file_interactive_tree_html': False,
    'file_newick_tree_png': False,
    'file_table_of_coevolution_tsv': True,
    'file_simulated_datasets_fastas': True,
    'file_barplot_of_correlation_svg': True,
    'file_plot_distribution_of_correlation_svg': True,
    'file_plot_correlation_by_rate_bin_svg': True,
    'file_table_of_posterior_rates_tsv': True,
    'file_table_of_pearson_correlation_tsv': True,
    'file_table_of_nodes_tsv': True,
    'file_branch_position_probabilities_tsv': True,
    'file_table_of_branches_tsv': True,
    'file_log_likelihood_tsv': True,
    'file_table_of_attributes_tsv': True,
    'file_table_of_parsimony_and_homoplasy_scores_tsv': True,
    'file_phylogenetic_tree_nwk': True
}

DEFAULT_ARGUMENTS = DefaultArgs(**{
    'with_internal_nodes': True,
    'sep': '\t'
    })

DEFAULT_ARGUMENTS.update(DEFAULT_FORM_ARGUMENTS)

ACTIONS = Actions(**{
                     'del_bootstrap_values': Tree.del_bootstrap_values,
                     'check_data': check_data,
                     'set_root': Tree.set_root,
                     'check_tree': Tree.rename_nodes,
                     'set_tree_data': Tree.set_tree_data,
                     'calculate_tree': Tree.calculate_tree,
                     'calculate_ancestral_sequence': Tree.calculate_ancestral_sequence,
                     'calculate_correlation': Tree.calculate_correlation,
                     'execute_all_actions': execute_all_actions
                     })

VALIDATION_ACTIONS = {
    'del_bootstrap_values': True,
    'check_data': True,
    'set_root': True,
    'check_tree': True
    }

DEFAULT_ACTIONS = {
    'set_tree_data': True,
    'calculate_tree': False,
    'calculate_ancestral_sequence': False,
    'calculate_correlation': False,
    'execute_all_actions': False
    }

MAIN_ACTIONS = {'compute_likelihood_of_tree': False,
                'draw_tree': False,
                'create_all_file_types': False}

CALCULATED_ARGS = CalculatedArgs(**{
                                    'file_path': None,
                                    'newick_text': None,
                                    'msa': None,
                                    'newick_tree': None
                                    })

USAGE = '''\tRequired parameters:
\t\t--msa_file <type=str>
\t\t\tSpecify the msa filepath.
\t\t--tree_file <type=str>
\t\t\tSpecify the newick filepath.
\tOptional parameters:
\t\t--out_dir <type=str>
\t\t\tSpecify the outdir path.
\t\t--process_id <type=str>
\t\t\tSpecify a process ID or it will be generated automatically.
\t\t--mode <type=str>
\t\t\tSpecify execution mode. Possible options: 
\t\t\t('draw_tree', 'compute_likelihood_of_tree', 'create_all_file_types', 'execute_all_actions'). 
\t\t\tDefault is {'execute_all_actions'}.
\t\t--with_internal_nodes <type=int> 
\t\t\tSpecify the tree has internal nodes. Default is 1.
\t\t--categories_quantity <type=int>
\t\t\tSpecify categories quantity. Default is 4.
\t\t--alpha <type=float>
\t\t\tSpecify alpha. Default is 0.5.
\t\t--pi_1 <type=float> 
\t\t\tSpecify pi_1. Default is 0.5.
\t\t--coefficient_bl <type=float> 
\t\t\tSpecify coefficient_bl. Default is 1.0.
\t\t--probability_lg <type=float> 
\t\t\tSpecify probability_lg. Default is 0.9.
\t\t--number_lg <type=float> 
\t\t\tSpecify number_lg. Default is 5.
\t\t--number_datasets <type=float> 
\t\t\tSpecify number_datasets. Default is 100.
\t\t--is_do_not_use_copap <type=int> 
\t\t\tSpecify is_do_not_use_copap. Default is 0.
\t\t--is_optimize_pi <type=int> 
\t\t\tSpecify is_optimize_pi. Default is 1.
\t\t--is_optimize_pi_average <type=int> 
\t\t\tSpecify is_optimize_pi_average. Default is 0.
\t\t--is_optimize_alpha <type=int> 
\t\t\tSpecify is_optimize_alpha. Default is 1.
\t\t--is_optimize_bl <type=int> 
\t\t\tSpecify is_optimize_bl. Default is 1.
\t\t--file_interactive_tree_html <type=int> 
\t\t\tSpecify file_interactive_tree_html. Default is 0.
\t\t--file_newick_tree_png <type=int> 
\t\t\tSpecify file_newick_tree_png. Default is 0.
\t\t--file_table_of_coevolution_tsv <type=int>
\t\t\tSpecify file_table_of_coevolution_tsv. Default is 1.
\t\t--file_simulated_datasets_fastas <type=int>
\t\t\tSpecify file_simulated_datasets_fastas. Default is 1.
\t\t--file_barplot_of_correlation_svg <type=int>
\t\t\tSpecify file_barplot_of_correlation_svg. Default is 1.
\t\t--file_plot_distribution_of_correlation_svg <type=int>
\t\t\tSpecify file_plot_distribution_of_correlation_svg. Default is 1.
\t\t--file_plot_correlation_by_rate_bin_svg <type=int>
\t\t\tSpecify file_plot_correlation_by_rate_bin_svg. Default is 1.
\t\t--file_table_of_posterior_rates_tsv <type=int>
\t\t\tSpecify file_table_of_posterior_rates_tsv. Default is 1.
\t\t--file_table_of_pearson_correlation_tsv <type=int>
\t\t\tSpecify file_table_of_pearson_correlation_tsv. Default is 1.
\t\t--file_table_of_nodes_tsv <type=int>
\t\t\tSpecify file_table_of_nodes_tsv. Default is 1.
\t\t--file_branch_position_probabilities_tsv 
\t\t\tSpecify file_branch_position_probabilities_tsv. Default is 1.
\t\t--file_table_of_branches_tsv <type=int> 
\t\t\tSpecify file_table_of_branches_tsv. Default is 1.
\t\t--file_log_likelihood_tsv <type=int> 
\t\t\tSpecify file_log_likelihood_tsv. Default is 1.
\t\t--file_table_of_attributes_tsv <type=int> 
\t\t\tSpecify file_table_of_attributes_tsv. Default is 1.
\t\t--file_table_of_parsimony_and_homoplasy_scores_tsv <type=int> 
\t\t\tSpecify file_table_of_parsimony_and_homoplasy_scores_tsv. Default is 1.
\t\t--file_phylogenetic_tree_nwk <type=int> 
\t\t\tSpecify file_phylogenetic_tree_nwk. Default is 1.
\t\t--rooting_method <type=str> 
\t\t\tSpecify tree rooting method. Possible options: ('mad', 'mvr', 'midpoint', 'outgroup').
\t\t\tmad - Minimal Ancestor Deviation
\t\t\tmvr - Minimum Variance Rooting
\t\t\tmidpoint - Midpoint Rooting
\t\t\toutgroup - Outgroup Rooting
\t\t\tDefault is 'midpoint'.
\t\t--leaf <type=str> 
\t\t\tSpecify leaf for outgroup rooting. Default is ''.'''
