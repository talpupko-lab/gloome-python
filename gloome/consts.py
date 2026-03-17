from os import getenv, chdir
from sys import argv
from typing import List, Tuple, Union
from dotenv import load_dotenv
from types import FunctionType, MethodType
from pathlib import Path
from urllib import parse
from importlib.resources import files

from gloome.tree.tree import Tree
from gloome.services.service_functions import check_data, execute_all_actions, recompile_json

load_dotenv()
MODE = ['draw_tree', 'compute_likelihood_of_tree', 'create_all_file_types', 'execute_all_actions']
IS_PRODUCTION = True
MAX_CONTENT_LENGTH = 16 * 1000 * 1000 * 1000
PREFIX = '/'
APPLICATION_ROOT = PREFIX
DEBUG = not IS_PRODUCTION

PREFERRED_URL_SCHEME = 'https'
WEBSERVER_NAME_CAPITAL = 'Gloome'
WEBSERVER_NAME = 'gloome.tau.ac.il'
WEBSERVER_URL = f'{PREFERRED_URL_SCHEME}://{WEBSERVER_NAME}'
RESULTS_URL = parse.urljoin(WEBSERVER_URL, 'results')
LOG_URL = WEBSERVER_URL

WEBSERVER_TITLE = '<b>GLOOME Server - Gain Loss Mapping Engine</b>'
MODULE_LOAD = 'module load mamba/mamba-1.5.8'

GLOOME = Path('/gloome')
BIN_DIR = GLOOME if GLOOME.exists() else Path.cwd()
RESULTS_DIR = BIN_DIR.joinpath('results')
IN_DIR = RESULTS_DIR.joinpath('in')
OUT_DIR = RESULTS_DIR.joinpath('out')
LOGS_DIR = BIN_DIR.joinpath('logs')
TMP_DIR = BIN_DIR.joinpath('tmp')
# APP_DIR = BIN_DIR.joinpath('app')
# TEMPLATES_DIR = APP_DIR.joinpath('templates')
# STATIC_DIR = APP_DIR.joinpath('static')
# ERROR_TEMPLATE = TEMPLATES_DIR.joinpath('404.html')
ENV = BIN_DIR.joinpath('.env')

ENVIRONMENT_DIR = BIN_DIR.joinpath('gloome_env2')
ENVIRONMENT_ACTIVATE = f'mamba activate {ENVIRONMENT_DIR}'

GLOOME_DIR = files('gloome')
DATA_DIR = GLOOME_DIR.joinpath('data')
INITIAL_DATA_DIR = DATA_DIR.joinpath('initial_data')

# HTTPDOCS_DIR = Path('/var/www/vhosts/gloome.tau.ac.il/httpdocs/')
# chdir(HTTPDOCS_DIR if HTTPDOCS_DIR.exists() else BIN_DIR)

MSA_FILE_NAME = 'msa_file.msa'
TREE_FILE_NAME = 'tree_file.tree'

REQUEST_WAITING_TIME = 20
REQUESTS_NUMBER = 24 * 60 * 60 * 3 / REQUEST_WAITING_TIME

if ENV.exists():
    UNDER_CONSTRUCTION = int(getenv('UNDER_CONSTRUCTION'))
    SECRET_KEY = getenv('SECRET_KEY')
    TOKEN = getenv('TOKEN')
    ACCOUNT = getenv('ACCOUNT')
    PARTITION = getenv('PARTITION')
    USE_OLD_SUBMITER = int(getenv('USE_OLD_SUBMITER'))

    LOGIN_NODE_URLS = getenv('LOGIN_NODE_URLS')
    USER_NAME = getenv('USER_NAME')
    USER_ID = getenv('USER_ID')
    USER_PASSWORD = getenv('USER_PASSWORD')
    ADMIN_EMAIL = getenv('ADMIN_EMAIL')
    SMTP_SERVER = getenv('SMTP_SERVER')
    SMTP_PORT = int(getenv('SMTP_PORT'))
    REPORT_RECEIVERS = getenv('REPORT_RECEIVERS').split()

    DEV_EMAIL = getenv('DEV_EMAIL')
    ADMIN_USER_NAME = getenv('ADMIN_USER_NAME')
    ADMIN_PASSWORD = getenv('ADMIN_PASSWORD')
    SEND_EMAIL_DIR_IBIS = getenv('SEND_EMAIL_DIR_IBIS')
    OWNER_EMAIL = getenv('OWNER_EMAIL')
else:
    UNDER_CONSTRUCTION = False
    SECRET_KEY = ''
    TOKEN = ''
    PARTITION = ''
    USE_OLD_SUBMITER = 0

    LOGIN_NODE_URLS = ''
    USER_NAME = ''
    USER_ID = ''
    USER_PASSWORD = ''
    ADMIN_EMAIL = ''
    SMTP_SERVER = ''
    SMTP_PORT = 0
    REPORT_RECEIVERS = []

    DEV_EMAIL = ''
    ADMIN_USER_NAME = ''
    ADMIN_PASSWORD = ''
    SEND_EMAIL_DIR_IBIS = ''
    OWNER_EMAIL = ''


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
    'e_mail': '',
    'is_optimize_pi': True,
    'is_optimize_pi_average': False,
    'is_optimize_alpha': True,
    'is_optimize_bl': True,
    'is_do_not_use_e_mail': True,
    'file_interactive_tree_html': False,
    'file_newick_tree_png': False,
    'file_table_of_nodes_tsv': True,
    'file_probability_per_pos_per_branches_tsv': True,
    'file_table_of_branches_tsv': True,
    'file_log_likelihood_tsv': True,
    'file_table_of_attributes_tsv': True,
    'file_phylogenetic_tree_nwk': True
}

DEFAULT_ARGUMENTS = DefaultArgs(**{
    'with_internal_nodes': True,
    'sep': '\t'
    })

DEFAULT_ARGUMENTS.update(DEFAULT_FORM_ARGUMENTS)

ACTIONS = Actions(**{
                     'check_data': check_data,
                     'check_tree': Tree.rename_nodes,
                     'set_tree_data': Tree.set_tree_data,
                     # 'compute_likelihood_of_tree': compute_likelihood_of_tree,
                     'calculate_tree': Tree.calculate_tree,
                     'calculate_ancestral_sequence': Tree.calculate_ancestral_sequence,
                     # 'draw_tree': draw_tree,
                     # 'create_all_file_types': create_all_file_types,
                     'execute_all_actions': execute_all_actions,
                     'recompile_json': recompile_json
                     })

VALIDATION_ACTIONS = {
    'check_data': True,
    'check_tree': True
    }

DEFAULT_ACTIONS = {
    'set_tree_data': True,
    'calculate_tree': False,
    # 'compute_likelihood_of_tree': False,
    'calculate_ancestral_sequence': False,
    # 'draw_tree': False,
    # 'create_all_file_types': False,
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
\t\t\tExecution mode style. Possible options: ('draw_tree', 'compute_likelihood_of_tree', 
\t\t\t'create_all_file_types', 'execute_all_actions'). Default is 'execute_all_actions'.
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
\t\t--file_table_of_nodes_tsv <type=int>
\t\t\tSpecify file_table_of_nodes_tsv. Default is 1.
\t\t--file_probability_per_pos_per_branches_tsv 
\t\t\tSpecify file_probability_per_pos_per_branches_tsv. Default is 1.
\t\t--file_table_of_branches_tsv <type=int> 
\t\t\tSpecify file_table_of_branches_tsv. Default is 1.
\t\t--file_log_likelihood_tsv <type=int> 
\t\t\tSpecify file_log_likelihood_tsv. Default is 1.
\t\t--file_table_of_attributes_tsv <type=int> 
\t\t\tSpecify file_table_of_attributes_tsv. Default is 1.
\t\t--file_phylogenetic_tree_nwk <type=int> 
\t\t\tSpecify file_phylogenetic_tree_nwk. Default is 1.
\t\t--e_mail <type=str> 
\t\t\tSpecify e_mail (technical parameter, do not change).
\t\t--is_do_not_use_e_mail <type=int> 
\t\t\tSpecify is_do_not_use_e_mail (technical parameter, do not change).'''

MENU = ({'name': 'Home', 'url': 'index',
         'submenu': ()
         },
        {'name': 'Overview', 'url': 'overview',
         'submenu': ()
         },
        {'name': 'Faq', 'url': 'faq',
         'submenu': ()
         },
        {'name': 'Gallery', 'url': 'gallery',
         'submenu': ()
         },
        {'name': 'Source code', 'url': 'source_code',
         'submenu': ()
         },
        {'name': 'Citing & credits', 'url': 'citing_and_credits',
         'submenu': ()
         }
        )
