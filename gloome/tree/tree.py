import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import re
import gc

from shutil import rmtree
from json import loads, dumps
from pathlib import Path
from d3blocks import D3Blocks
from typing import Optional, List, Union, Dict, Tuple, Set, Any, Callable
from Bio import Phylo
from Bio.Phylo.NewickIO import Writer
from numpy import ndarray, dtype, void
from scipy.stats import gamma, pearsonr, distributions, beta as sp_beta, kstest
from scipy.special import gammainc
from scipy.optimize import minimize_scalar
from io import StringIO

from .node import Node
from .npencoder import NpEncoder

eps = 5e-324
eps2 = 1e-10


class Tree:
    root: Optional[Node] = None
    alphabet: Optional[Tuple[str, ...]] = None
    msa: Optional[Dict[str, str]] = None
    rate_vector: Tuple[Union[float, np.float64, int], ...] = (1.0, )
    alpha: Optional[Union[float, np.float64, int]] = None
    categories_quantity: Optional[int] = None
    pi_1: Optional[Union[float, np.float64, int]] = None
    coefficient_bl: Optional[Union[float, np.float64, int]] = 1,
    log_likelihood_vector: Optional[np.ndarray] = None
    log_likelihood: Optional[Union[float, np.float64]] = None
    likelihood_vector: Optional[np.ndarray] = None
    likelihood: Optional[Union[float, np.float64]] = None
    posterior_rates: Optional[np.ndarray] = None
    correlation_vector: Optional[np.ndarray] = None
    calculated_ancestor_sequence: bool = False
    calculated_tree: bool = False
    calculated_likelihood: bool = False
    all_nodes: Dict[str, Node]
    all_nodes_objects: Optional[List[Node]] = None
    all_nodes_objects_post_order: Optional[List[Node]] = None
    nodes_objects: Optional[List[Node]] = None
    nodes_objects_post_order: Optional[List[Node]] = None
    leaves_objects: Optional[List[Node]] = None
    alphabet_length: int
    msa_length: int
    rate_vector_length: int

    def __init__(self, data: Optional[Union[str, Node]] = None, node_name: Optional[str] = None, **kwargs) -> None:
        """
         Args:
            data (str, Node, optional): `None` (default)
            node_name (str, optional): `None` (default)
            msa (Dict, str, optional): `None` (default)
            categories_quantity (float, optional): `None` (default)
            alpha (float, optional): `None` (default)
            beta (float, optional): `None` (default)
            pi_0 (float, np.float64, int, optional): `None` (default)
            pi_1 (float, np.float64, int, optional): `None` (default)
            coefficient_bl (float, np.float64, int, optional): `None` (default)
            is_optimize_pi (bool, optional): `None` (default)
            is_optimize_pi_average (bool, optional): `None` (default)
            is_optimize_alpha (bool, optional): `None` (default)
            is_optimize_bl (bool, optional): `None` (default)
        """
        available_parameters = {'data', 'node_name', 'msa', 'categories_quantity', 'alpha', 'beta', 'pi_0', 'pi_1',
                                'coefficient_bl', 'is_optimize_pi', 'is_optimize_pi_average', 'is_optimize_alpha',
                                'is_optimize_bl', 'seed'}
        invalid_parameters = set(kwargs.keys()) - available_parameters
        for key in invalid_parameters:
            del kwargs[key]

        if isinstance(data, str):
            data = self.del_bootstrap_values(data)
            self.newick_to_tree(data)
            if node_name and isinstance(node_name, str):
                self.rename_nodes(self, node_name)
        elif isinstance(data, Node):
            self.root = data
        else:
            self.root = Node('root')

        self.all_nodes_objects = self.get_all_nodes()
        self.all_nodes_objects_post_order = self.get_all_nodes(mode='post-order')
        self.all_nodes = {current_node.name: current_node for current_node in self.all_nodes_objects}
        self.nodes_objects = self.get_nodes()
        self.nodes_objects_post_order = self.get_nodes(mode='post-order')
        self.leaves_objects = self.get_leaves()

        self.msa = self.alphabet = self.categories_quantity = self.alpha = None
        self.rate_vector = (1.0, )
        self.rate_vector_length = 1
        self.pi_1, self.coefficient_bl = None, 1
        self.log_likelihood_vector = self.likelihood_vector = self.correlation_vector = self.posterior_rates = None
        self.log_likelihood = self.likelihood = 0.0

        self.calculated_ancestor_sequence = self.calculated_tree = self.calculated_likelihood = False

        if any(kwargs.values()):
            self.set_tree_data(**kwargs)
        elif invalid_parameters:
            print(f'There are invalid parameters: {", ".join(invalid_parameters)}')

    def __str__(self) -> str:

        return self.get_newick(with_internal_nodes=True)

    def __dir__(self) -> List[str]:

        return ['root', 'alphabet', 'msa', 'rate_vector', 'alpha', 'categories_quantity', 'pi_1', 'coefficient_bl',
                'log_likelihood_vector', 'log_likelihood', 'likelihood_vector', 'likelihood', 'posterior_rates',
                'correlation_vector', 'calculated_ancestor_sequence', 'calculated_tree', 'calculated_likelihood',
                'all_nodes', 'all_nodes_objects', 'nodes_objects', 'nodes_objects_post_order', 'leaves_objects',
                'alphabet_length', 'msa_length', 'rate_vector_length']

    def __len__(self) -> int:

        return self.get_node_count()

    def __eq__(self, other) -> bool:

        return str(self).lower() == str(other).lower()

    def __ne__(self, other) -> bool:

        return not self == other

    def __lt__(self, other) -> bool:

        return len(self) < len(other)

    def __le__(self, other) -> bool:

        return self < other or self == other or len(str(self)) < len(str(other))

    def __gt__(self, other) -> bool:

        return len(self) > len(other)

    def __ge__(self, other) -> bool:

        return self > other or self == other or len(str(self)) > len(str(other))

    def get_dict(self) -> Dict[str, Optional[Union[Node, float, np.float64, int, np.ndarray, bool, Tuple[str, ...],
                               Tuple[Union[float, np.float64, int], ...], Dict[str, str], List[Node]]]]:

        return {'root': self.root,
                'alphabet': self.alphabet,
                'msa': self.msa,
                'rate_vector': self.rate_vector,
                'alpha': self.alpha,
                'categories_quantity': self.categories_quantity,
                'pi_1': self.pi_1,
                'coefficient_bl': self.coefficient_bl,
                'log_likelihood_vector': self.log_likelihood_vector,
                'log_likelihood': self.log_likelihood,
                'likelihood_vector': self.likelihood_vector,
                'likelihood': self.likelihood,
                'posterior_rates': self.posterior_rates,
                'correlation_vector': self.correlation_vector,
                'calculated_ancestor_sequence': self.calculated_ancestor_sequence,
                'calculated_tree': self.calculated_tree,
                'calculated_likelihood': self.calculated_likelihood,
                'all_nodes': self.all_nodes,
                'all_nodes_objects': self.all_nodes_objects,
                'nodes_objects': self.nodes_objects,
                'nodes_objects_post_order': self.nodes_objects_post_order,
                'leaves_objects': self.leaves_objects,
                'alphabet_length': self.alphabet_length,
                'msa_length': self.msa_length,
                'rate_vector_length': self.rate_vector_length}

    def print_args(self, prefix_name: str = '', prefix: str = '', sort: bool = False) -> None:
        if all((prefix_name, prefix)):
            print(f'{prefix_name}\t\t>\t>\t>\t\t{prefix}')
        items = dict(sorted(self.get_dict().items())).items() if sort else self.get_dict().items()
        for key, value in items:
            print(f'{key}:\t{value}')

    def set_tree_data(self, msa: Optional[Union[Dict[str, str], str]] = None,
                      categories_quantity: Optional[int] = None,
                      alpha: Optional[float] = None,
                      beta: Optional[float] = None,
                      pi_0: Optional[Union[float, np.float64, int]] = None,
                      pi_1: Optional[Union[float, np.float64, int]] = None,
                      coefficient_bl: Optional[Union[float, np.float64, int]] = None,
                      is_optimize_pi: Optional[bool] = None,
                      is_optimize_pi_average: Optional[bool] = None,
                      is_optimize_alpha: Optional[bool] = None,
                      is_optimize_bl: Optional[bool] = None,
                      seed: Optional[int] = None) -> None:

        if seed is not None:
            np.random.seed(seed)

        if isinstance(msa, str):
            self.msa = self.get_msa_dict(msa)
        elif isinstance(msa, dict):
            self.msa = msa

        if isinstance(self.msa, dict) and self.msa:
            self.alphabet = self.get_alphabet_from_dict(self.msa)
        else:
            self.alphabet = self.get_alphabet()
            self.set_basic_msa()

        self.msa_length = len(next(iter(self.msa.values())))
        self.alphabet_length = len(self.alphabet)

        self.set_all(categories_quantity, alpha, beta, pi_0, pi_1, coefficient_bl)

        self.optimize_coefficient_bl(is_optimize_bl)
        self.optimize_pi(is_optimize_pi, is_optimize_pi_average)
        self.optimize_alpha(is_optimize_alpha)

        if (is_optimize_alpha or is_optimize_pi or is_optimize_pi_average) and is_optimize_bl:
            self.optimize_coefficient_bl(is_optimize_bl)

        self.set_all(categories_quantity=self.categories_quantity, alpha=self.alpha, pi_1=self.pi_1,
                     coefficient_bl=self.coefficient_bl)
        self.set_distance_taking_into_coefficient()

    def set_distance_taking_into_coefficient(self) -> None:
        for current_node in self.all_nodes_objects:
            current_node.distance_to_father_taking_into_coefficient = (current_node.distance_to_father *
                                                                       self.coefficient_bl)
            current_node.distance_to_nearest_taking_into_coefficient = (current_node.distance_to_nearest *
                                                                        self.coefficient_bl)
            current_node.distance_to_root_taking_into_coefficient = current_node.distance_to_root * self.coefficient_bl
            current_node.distance_to_root_vector_taking_into_coefficient = [i * self.coefficient_bl for i in
                                                                            current_node.distance_to_root_vector]

    def print_node_list(self, with_additional_details: bool = False, mode: Optional[str] = None,
                        filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None) -> None:
        """
        Print a list of nodes.

        This function prints a list of nodes.

        Args:
            with_additional_details (bool, optional): `False` (default)
            mode (str, optional): `None` (default), 'pre-order', 'in-order', 'post-order', 'level-order'
            filters (Dict, optional): `None` (default)

        Returns:
            None: This function does not return any value; it only prints the nodes to the standard output.
        """
        data_structure = self.get_list_nodes_info(with_additional_details, mode, filters)

        str_result = ''
        for i in data_structure:
            str_result = f'{str_result}\n{i}'
        print(str_result, '\n')

    def get_tree_info(self, filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None
                      ) -> pd.Series:
        nodes_info = self.get_list_nodes_info(True, filters=filters)

        return pd.Series([pd.Series(i) for i in nodes_info], index=[i.get('node') for i in nodes_info])

    def get_list_nodes_info(self,
                            with_additional_details: bool = False,
                            mode: Optional[str] = None,
                            filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None,
                            only_node_list: bool = False,
                            fields: Optional[Set[str]] = None
                            ) -> List[Union[Dict[str, Union[float, np.float64, bool, str, np.ndarray,
                                      List[float], List[np.float64]]], Node]]:
        """
        Retrieve a list of all nodes of the tree.

        This function collects all nodes of the tree. The function returns a list of nodes or a list of
        dictionaries with information about these nodes.

        Args:
            with_additional_details (bool, optional): `False` (default).
            mode (str, optional): `pre-order` (default), 'pre-order', 'in-order', 'post-order', 'level-order'.
            filters (Dict, optional):
            only_node_list (bool, optional): `False` (default).
            fields (Set[str], optional): `None` (default) builds every field. Restrict to only the keys actually
                needed -- see ``Node.get_node_info`` -- to skip building/serializing the big per-node tensors
                (up_vector/down_vector/marginal_vector/marginal_bl_vector/pmatrix) when they won't be used.

        Returns:
            list: A list of all nodes of the tree or a list of dictionaries with information about these nodes.
        """

        return self.root.get_list_nodes_info(with_additional_details, mode, filters, only_node_list, fields)

    def get_leaves(self, only_node_list: bool = True, mode: Optional[str] = None) -> List[Union[Node, str]]:

        return self.get_list_nodes_info(filters={'node_type': ['leaf']}, only_node_list=only_node_list, mode=mode)

    def get_nodes(self, only_node_list: bool = True, mode: Optional[str] = None) -> List[Union[Node, str]]:

        return self.get_list_nodes_info(filters={'node_type': ['node', 'root']}, only_node_list=only_node_list,
                                        mode=mode)

    def get_all_non_root_nodes(self, only_node_list: bool = True, mode: Optional[str] = None) -> List[Union[Node, str]]:

        return self.get_list_nodes_info(filters={'node_type': ['node', 'leaf']}, only_node_list=only_node_list,
                                        mode=mode)

    def get_all_nodes(self, only_node_list: bool = True, mode: Optional[str] = None) -> List[Union[Node, str]]:

        return self.get_list_nodes_info(only_node_list=only_node_list, mode=mode)

    def get_leaves_count(self) -> int:

        return self.get_node_count(filters={'node_type': ['leaf']})

    def get_nodes_count(self) -> int:

        return self.get_node_count(filters={'node_type': ['node', 'root']})

    def get_node_count(self, filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None) -> int:

        return len(self.get_list_nodes_info(filters=filters, only_node_list=True))

    def get_node_by_name(self, name: str) -> Optional[Node]:

        return self.root.get_node_by_name(name)

    def get_newick(self, with_internal_nodes: bool = False,
                   decimal_length: int = 0,
                   taking_into_coefficient: bool = False) -> str:

        """
        Convert the current tree structure to a Newick formatted string.

        This function serializes the tree into a Newick format, which is a standard format for representing
        tree structures.

        Args:
            with_internal_nodes (bool, optional):
            decimal_length (int, optional):
            taking_into_coefficient (bool, optional):

        Returns:
            str: A Newick formatted string representing the tree structure.
        """
        return f'{self.root.subtree_to_newick(with_internal_nodes, decimal_length, taking_into_coefficient)};'

    def find_node_by_name(self, name: str) -> bool:
        """
        Search for a node by its name in a tree structure.

        This function searches for a node with a specific name within a tree. If a root node is provided,
        the search starts from that node; otherwise, it searches from the default root of the tree.
        The function returns `True` if a node with the specified name is found, and `False` otherwise.

        Args:
            name (str): The name of the node to search for. This should be the exact name of the node
                        as a string.

        Returns:
            bool: `True` if a node with the specified name is found; `False` otherwise.
        """

        return name in self.get_list_nodes_info()

    def newick_to_tree(self, newick: str) -> Optional['Tree']:
        """
        Convert a Newick formatted string into a tree object.

        This function parses a Newick string, which represents a tree structure in a compact format,
        and constructs a tree object from it. The Newick format is often used in phylogenetics to
        describe evolutionary relationships among species.

        Args:
            newick (str): A string in Newick format representing the tree structure. The string
                              should be properly formatted according to Newick syntax.

        Returns:
            Tree: An object representing the tree structure parsed from the Newick string. The tree
                  object provides methods and properties to access and manipulate the tree structure.
        """
        newick = newick.replace(' ', '').strip()
        if newick.startswith('(') and newick.endswith(';'):

            len_newick = len(newick)
            list_end = [i for i in range(len_newick) if newick[i:i + 1] == ')']
            list_start = [i for i in range(len_newick) if newick[i:i + 1] == '(']
            list_children = []

            num = self.__counter()

            while list_start:
                int_start = list_start.pop(-1)
                int_end = min([i for i in list_end if i > int_start]) + 1
                list_end.pop(list_end.index(int_end - 1))
                node_name = newick[int_end: min([x for x in [newick.find(':', int_end), newick.find(',', int_end),
                                                 newick.find(';', int_end), newick.find(')', int_end)] if x >= 0])]
                distance_to_father = newick[int_end + len(node_name): min([x for x in [newick.find(',', int_end),
                                                                          newick.find(';', int_end), newick.find(')',
                                                                          int_end)] if x >= 0])]

                (visibility, node_name) = (True, node_name) if node_name else (False, 'nd' + str(num()).rjust(4, '0'))

                sub_str = newick[int_start:int_end]
                list_children.append({'children': sub_str, 'node': node_name, 'distance_to_father': distance_to_father,
                                      'visibility': visibility})

            list_children.sort(key=lambda x: len(x.get('children')), reverse=True)
            for i in range(len(list_children)):
                for j in range(i + 1, len(list_children)):
                    node_name = list_children[j].get('node') if list_children[j].get('visibility') else ''
                    list_children[i].update({'children': list_children[i].get('children').replace(
                        list_children[j].get('children') + node_name, list_children[j].get('node'))})
            for dict_children in list_children:
                if list_children.index(dict_children):
                    newick_node = self.get_node_by_name(dict_children.get('node'))
                else:
                    newick_node = self.__set_node(
                        f'{dict_children.get("node")}{dict_children.get("distance_to_father")}', num)
                    newick_node.distance_to_root_vector = [0.0]
                    newick_node.level = 1
                    self.root = newick_node
                self.__set_children_list_from_string(dict_children.get('children'), newick_node, num)
            for current_node in self.get_list_nodes_info(only_node_list=True):
                current_node.set_levels_and_distance_to_nearest()
                current_node.aliases = {'node': 'name', 'distance': 'distance_to_father'}
                if current_node.node_type in ('node', ) and self.is_bootstrap_value(current_node.name):
                    current_node.name = 'nd' + str(num()).rjust(4, '0')

            return self

    def get_html_tree(self, style: str = '', status: str = '') -> str:

        return self.structure_to_html_tree(self.tree_to_structure(), style, status)

    def tree_to_structure(self) -> Dict[str, str]:

        return self.subtree_to_structure(self.root)

    def add_distance_to_father(self, distance_to_father: float = 0) -> None:
        def add_distance(newick_node: Node) -> None:
            nonlocal distance_to_father
            newick_node.distance_to_father += distance_to_father
            newick_node.distance_to_father = round(newick_node.distance_to_father, 12)
            for child in newick_node.children:
                add_distance(child)

        add_distance(self.root)

    def get_edges_list(self) -> List[str]:
        list_result = []

        def get_list(newick_node: Node) -> None:
            nonlocal list_result
            if newick_node.father:
                list_result.append((newick_node.father.name, newick_node.name))
            for child in newick_node.children:
                get_list(child)

        get_list(self.root)

        return list_result

    def __set_children_list_from_string(self, str_children: str, father: Node, num) -> None:
        str_children = str_children[1:-1] if str_children.startswith('(') and str_children.endswith(
            ')') else str_children
        lst_nodes = str_children.split(',')
        father.node_type = 'node' if father.father else 'root'
        for str_node in lst_nodes:
            newick_node = self.__set_node(str_node.strip(), num)
            newick_node.node_type = 'leaf'
            newick_node.father = father
            newick_node.distance_to_root_vector = father.distance_to_root_vector.copy()
            newick_node.level = father.level + 1
            newick_node.distance_to_root_vector.append(newick_node.distance_to_father)
            newick_node.distance_to_root = round(sum(newick_node.distance_to_root_vector), 14)
            father.add_child(newick_node)

    def check_tree_for_binary(self) -> bool:
        nodes_list = self.get_list_nodes_info(True)
        for current_node in nodes_list:
            for key in current_node.keys():
                if key == 'children' and len(current_node.get(key)) > 2:
                    return False

        return True

    def tree_to_table(self, sort_values_by: Optional[Tuple[str, ...]] = None, decimal_length: int = 8, columns: Optional
                      [Dict[str, str]] = None, filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] =
                      None, distance_type: type = str, list_type: type = str, lists: Optional[Tuple[str, ...]] = None,
                      taking_into_coefficient: bool = True, decimals: int = 4) -> pd.DataFrame:
        suffix = '_taking_into_coefficient' if taking_into_coefficient else ''
        distance_name = f'distance{suffix}'
        full_distance_name = f'full_distance{suffix}'
        distance_to_nearest_name = f'distance_to_nearest{suffix}'

        columns = columns if columns else {'node': 'Name',
                                           'father_name': 'Parent',
                                           distance_name: 'Distance to parent',
                                           'children': 'Children',
                                           'level': 'Level',
                                           'node_type': 'Node type',
                                           distance_to_nearest_name: 'Distance to nearest leaf',
                                           'levels_to_nearest': 'Levels to nearest leaf',
                                           full_distance_name: 'Full distance',
                                           'up_vector': 'Up',
                                           'down_vector': 'Down',
                                           'likelihood': 'Likelihood',
                                           'likelihood_vector': 'Vector of likelihood',
                                           'sequence_likelihood': 'Likelihood of sequence',
                                           'log_likelihood': 'Log-likelihood',
                                           'log_likelihood_vector': 'Vector of log-likelihood',
                                           'marginal_vector': 'Marginal vector',
                                           'marginal_bl_vector': 'Marginal branch vector',
                                           'probability_vector': 'Probability vector',
                                           'sequence': 'Sequence',
                                           'ancestral_sequence': 'Ancestral Comparison',
                                           'probabilities_sequence_characters': 'character sequence probabilities',
                                           'probability_vector_gain': 'Gain probability',
                                           'probability_vector_loss': 'Loss probability'}
        lists = lists if lists else ('children', 'full_distance', 'full_distance_taking_into_coefficient', 'up_vector',
                                     'down_vector', 'marginal_vector', 'marginal_bl_vector', 'probability_vector',
                                     'probabilities_sequence_characters', 'log_likelihood_vector', 'likelihood_vector',
                                     'sequence', 'ancestral_sequence', 'probability_vector_gain',
                                     'probability_vector_loss')
        exceptions = ('sequence', 'ancestral_sequence')

        nodes_info = self.get_list_nodes_info(True, None, filters, fields=set(columns.keys()) | {'father_name'})

        tree_table = pd.DataFrame(nodes_info)
        if tree_table.empty:
            return pd.DataFrame(columns=np.asarray(columns.values()))

        existing_cols = [c for c in columns.keys() if c in tree_table.columns]
        tree_table = tree_table[existing_cols].copy()

        if 'father_name' in tree_table.columns:
            tree_table['father_name'] = tree_table['father_name'].fillna('root').replace('', 'root')

        if distance_name in tree_table.columns and distance_name in columns:
            if distance_type is str:
                blank_pad = ' ' * decimal_length
                tree_table[distance_name] = tree_table[distance_name].apply(
                    lambda x: f'{x:.10f}'.ljust(decimal_length, "0") if pd.notna(x) and x != '' else blank_pad
                )
            else:
                tree_table[distance_name] = tree_table[distance_name].astype(distance_type, errors='ignore')

        active_lists = [curr_col for curr_col in lists if curr_col in tree_table.columns and curr_col in columns]
        for col_name in active_lists:
            is_exception = col_name in exceptions
            tree_table[col_name] = [
                self.get_list_decimals(val, list_type, decimals, is_exception=is_exception, force_list=True)
                for val in tree_table[col_name]
            ]

        self.mark_unknown_characters(tree_table)

        tree_table = tree_table.rename(columns=columns)
        tree_table = tree_table.reindex(columns=columns.values())

        if sort_values_by:
            if isinstance(list_type, (list, tuple, set)):
                lists_names = {columns[k] for k in active_lists}
                sort_values_by = tuple(i for i in sort_values_by if i not in lists_names)

            if sort_values_by:
                tree_table = tree_table.sort_values(by=list(sort_values_by))

        return tree_table

    def mark_unknown_characters(self, tree_table: pd.DataFrame) -> None:
        if not self.msa or 'node' not in tree_table.columns:
            return
        unknown_characters = self.get_unknown_characters()
        for col_name in ('sequence', 'ancestral_sequence'):
            if col_name not in tree_table.columns:
                continue
            values = tree_table[col_name].tolist()
            for row, node_name in enumerate(tree_table['node']):
                observed = self.msa.get(node_name)
                if not observed or not any(char in observed for char in unknown_characters):
                    continue
                value = values[row]
                if isinstance(value, list) and len(value) == len(observed):
                    values[row] = [f'{char} ({j})' if char in unknown_characters else j
                                   for char, j in zip(observed, value)]
                elif isinstance(value, str) and len(value) == len(observed):
                    values[row] = ''.join(char if char in unknown_characters else j
                                          for char, j in zip(observed, value))
            tree_table[col_name] = values

    def calculate_ancestral_sequence(self, newick_node: Optional[Union[Node, str]] = None) -> str:
        if self.alphabet and not self.calculated_ancestor_sequence:
            node_list = []
            if not newick_node:
                node_list = self.all_nodes_objects[1:]
            else:
                node_list.append(newick_node)

            ancestral_alphabet = self.get_ancestral_alphabet()
            for current_node in node_list:
                current_node.ancestral_sequence = ''
                if current_node.father:
                    for i in range(len(current_node.sequence)):
                        if current_node.sequence[i] == current_node.father.sequence[i] == self.alphabet[0]:
                            current_node.ancestral_sequence += ancestral_alphabet[0]
                        elif ((current_node.sequence[i] != current_node.father.sequence[i])
                              and (current_node.sequence[i] == self.alphabet[0])):
                            current_node.ancestral_sequence += ancestral_alphabet[1]
                        elif ((current_node.sequence[i] != current_node.father.sequence[i])
                              and (current_node.sequence[i] == self.alphabet[1])):
                            current_node.ancestral_sequence += ancestral_alphabet[2]
                        elif current_node.sequence[i] == current_node.father.sequence[i] == self.alphabet[1]:
                            current_node.ancestral_sequence += ancestral_alphabet[3]
            self.calculated_ancestor_sequence = True

        return 'OK' if self.calculated_ancestor_sequence else ''

    def calculate_marginal(self) -> None:
        for current_node in self.all_nodes_objects:
            current_node.calculate_marginal(self.rate_vector_length, self.msa_length)

    def calculate_down(self) -> None:
        for current_node in self.all_nodes_objects:
            current_node.calculate_down(self.rate_vector_length, self.alphabet_length, self.msa_length)

    def calculate_up(self) -> None:
        self.alphabet_length = len(self.alphabet)
        self.msa_length = len(next(iter(self.msa.values())))
        self.rate_vector_length = len(self.rate_vector)

        self.initialize_leaf_up_vectors()
        self.initialize_node_up_vectors()

        self.likelihood_vector = self.root.likelihood_vector
        self.likelihood = self.root.likelihood
        self.log_likelihood_vector = self.root.log_likelihood_vector
        self.log_likelihood = self.root.log_likelihood

        self.calculated_likelihood = True

    def get_log_likelihood(self) -> Union[np.float64, float]:
        if self.msa and self.alphabet:
            self.calculate_up()

        return self.log_likelihood

    def initialize_node_up_vectors(self) -> None:
        for current_node in self.nodes_objects_post_order:
            current_node.calculate_up(self.rate_vector_length, self.alphabet_length, self.msa_length)

    def initialize_leaf_up_vectors(self) -> None:
        for leaf in self.leaves_objects:
            sequence = np.asarray(list(self.msa[leaf.name]))
            up_vector = np.zeros((self.rate_vector_length, self.alphabet_length, self.msa_length), dtype=np.float64)
            mask_0 = (sequence == '0')
            mask_1 = (sequence == '1')

            up_vector[:, 0, mask_0] = 1.0
            up_vector[:, 1, mask_1] = 1.0

            mask_gap = np.isin(sequence, ('-', ) + self.get_unknown_characters())
            up_vector[:, :, mask_gap] = 1.0

            leaf.up_vector = up_vector

    def get_msa_dict(self, msa: str, alphabet: Optional[Union[Tuple[str, ...], str]] = None, only_leaves: bool = True
                     ) -> Dict[str, Union[Tuple[int, ...], str]]:
        node_types = ['leaf'] if only_leaves else ['leaf', 'node', 'root']
        node_names = self.get_list_nodes_info(False, 'pre-order', {'node_type': node_types})
        msa_list = msa.strip().split()
        msa_list_size, msa_dict = len(msa_list), dict()
        if msa_list_size == 1:
            for i, node_name in enumerate(node_names):
                if alphabet:
                    value = [0] * len(alphabet)
                    value[alphabet.index(msa[i])] = 1
                    value = tuple(value)
                else:
                    value = msa[i]
                msa_dict.update({node_name: value})
        else:
            for j in range(msa_list_size // 2):
                node_name = msa_list[j + j][1::]
                if node_name in node_names:
                    value = msa_list[j + j + 1]
                    value = ''.join(value)
                    msa_dict.update({node_name: value})

        return msa_dict

    def calculate_correlation(self, prior: Optional[np.ndarray] = None, probability_lg: Union[float, np.float64] = 0.5,
                              number_lg: Union[float, np.float64, int] = 1) -> None:
        self.set_posterior_rates_vector(prior)
        self.set_pearson_correlation_vector(probability_lg, number_lg)

    def calculate_tree(self) -> Dict[str, Union[float, int, np.float64, np.ndarray]]:
        if self.msa and not self.calculated_tree:
            self.clean_all()

            self.calculate_up()
            self.calculate_down()
            self.calculate_marginal()

            self.calculated_tree = True

        return {'likelihood': self.likelihood,
                'likelihood_vector': self.likelihood_vector,
                'log_likelihood': self.log_likelihood,
                'log_likelihood_vector': self.log_likelihood_vector}

    def calculate_likelihood(self) -> None:
        if self.msa and self.alphabet and not self.calculated_likelihood:
            self.clean_all()
            self.calculate_up()

    def get_fasta_text(self, msa: Optional[Dict[str, str]] = None) -> str:

        return ''.join(f'>{k}\n{v}\n' for k, v in (self.msa if msa is None else msa).items()).strip()

    def get_json_structure(self, return_table: bool = False,
                           columns: Optional[Dict[str, str]] = None,
                           mode: str = 'node',
                           taking_into_coefficient: bool = True
                           ) -> Dict[str, Union[List[str], str]]:
        """
        Args:
            return_table (bool, optional): `False` (default).
            columns (dict, optional): `None` (default).
            mode (str, optional): 'node' (default), 'branch'.
            taking_into_coefficient (bool, optional): `True` (default).

        Returns:
            Dict: An dictionary representing the tree structure.
        """
        if return_table:
            columns, lists, decimals = self.get_columns(mode, columns, taking_into_coefficient)
            columns_names = {'node': 'Name', 'branch': 'Child node'}
            column_name = columns_names.get(mode, 'Name')

            table = self.tree_to_table(columns=columns, list_type=list, lists=lists, distance_type=float,
                                       taking_into_coefficient=taking_into_coefficient, decimals=decimals)
            dict_json = dict()
            for row in table.T:
                dict_row = dict()
                for key in columns.values():
                    dict_row.update({key: table[key][row]})
                dict_json.update({table[column_name][row]: dict_row})
        else:
            dict_json = self.root.node_to_json()

        return loads(dumps(dict_json, cls=NpEncoder).replace(f'\'', r'"'))

    def tree_to_fasta_file(self, file_name: str = 'file.fasta') -> str:
        fasta_text = self.get_fasta_text()

        return self.write_file(file_name, fasta_text)

    def probability_to_tsv(self, file_name: str = 'BranchPositionProbabilities.tsv', sep: str = '\t',
                           taking_into_coefficient: bool = True) -> str:
        self.make_dir(file_name)
        ancestral_comparison = ['absence', 'loss', 'gain', 'presence']
        probability_limit = 0.05
        rows = []

        suffix = '_taking_into_coefficient' if taking_into_coefficient else ''
        distance_to_father = f'distance_to_father{suffix}'
        distance_to_root = f'distance_to_root{suffix}'
        distance_to_nearest = f'distance_to_nearest{suffix}'

        for current_node in self.all_nodes_objects:
            branch_probability_vector = current_node.branch_probability_vector
            for pos, value in enumerate(branch_probability_vector):
                for i in range(1, 3):
                    row = {
                        'G/L': ancestral_comparison[i],
                        'POS': pos,
                        'branch': current_node.name,
                        'branchLength': getattr(current_node, distance_to_father),
                        'distance2root': getattr(current_node, distance_to_root),
                        'distance2NearestOTU': getattr(current_node, distance_to_nearest),
                        'numOfNodes2NearestOTU': current_node.levels_to_nearest,
                        'probability': value[i]
                    }
                    rows.append(row)

        df = pd.DataFrame(rows)
        df['G/L'] = df['G/L'].astype(pd.api.types.CategoricalDtype(categories=ancestral_comparison, ordered=True))
        df['POS'] = df['POS'].astype(int)
        df['branch'] = df['branch'].astype(pd.api.types.CategoricalDtype(categories=self.get_list_nodes_info(),
                                                                         ordered=True))
        df['branchLength'] = df['branchLength'].astype(float)
        df['distance2root'] = df['distance2root'].astype(float)
        df['distance2NearestOTU'] = df['distance2NearestOTU'].astype(float)
        df['numOfNodes2NearestOTU'] = df['numOfNodes2NearestOTU'].astype(int)
        df['probability'] = df['probability'].astype(float)
        df = df.sort_values(by=['POS', 'branch', 'G/L'])
        df = df.query(f'probability > {probability_limit} and `G/L` in {ancestral_comparison[1:3]}')
        df.to_csv(file_name, sep=sep, index=False, float_format='%.8f', chunksize=50000)

        return file_name

    def build_site_event_matrix(self) -> np.ndarray:
        """Return (msa_length, 2 * n_nodes) a matrix whose row for a site is the concatenation,
        over every non-root node, of [loss_probability, gain_probability] -- the exact per-site
        vector ``identify_event_candidates`` builds. Also returns the ordered node list (the root
        is excluded: it has no incoming branch, hence no gain/loss probability)."""
        nodes_list = self.all_nodes_objects[1:]
        loss = np.asarray([n.probability_vector_loss for n in nodes_list], dtype='float64')
        gain = np.asarray([n.probability_vector_gain for n in nodes_list], dtype='float64')

        site_matrix = np.empty((self.msa_length, 2 * len(nodes_list)), dtype='float64')
        site_matrix[:, 0::2], site_matrix[:, 1::2] = loss.T, gain.T

        return site_matrix

    def identify_event_candidates(self, event_threshold: Union[np.float64, float] = 0.5
                                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extracts gain/loss matrices, filters candidate sites, and categorizes evolutionary rates.

        Args:
            event_threshold: The minimum probability threshold to consider a site as an event.

        Returns:
            Tuple containing:
                - site_matrix (np.ndarray): Interleaved loss/gain probability matrix (MSA length x 2*Nodes).
                - candidates (np.ndarray): Indices of variable sites exceeding the threshold.
                - categories (np.ndarray): Assigned rate categories for each site.
        """
        site_matrix = self.build_site_event_matrix()

        candidates = np.where((site_matrix.max(axis=1) > event_threshold) & (site_matrix.var(axis=1) > 0))[0]
        categories = np.abs(self.posterior_rates[:, None] - np.array(self.rate_vector)[None, :]).argmin(axis=1)

        return site_matrix, candidates, categories

    def get_bins(self, site_matrix: np.ndarray, candidates: np.ndarray, categories: np.ndarray
                 ) -> Tuple[ndarray[Any, dtype[Any]], ndarray[Any, Any], List[Tuple[Union[Union[ndarray[Any, Any],
                            ndarray[Any, dtype[Any]], ndarray[Any, dtype[void]]], Any], ...]]]:
        correlations = self.get_row_correlations(site_matrix[candidates])
        i_idx, j_idx = np.triu_indices(len(candidates), k=1)
        pairs = np.column_stack((candidates[i_idx], candidates[j_idx]))
        r_values = correlations[i_idx, j_idx]
        bins = [tuple(sorted((categories[a], categories[b]))) for a, b in pairs]

        return pairs, r_values, bins

    def simulate_datasets(self, file_path: str = '',
                          sep: str = '\t',
                          number_datasets: int = 100,
                          event_threshold: Union[float, np.float64] = 0.5,
                          p_value_mode: str = 'empirical',
                          min_fit: int = 100,
                          use_simulated_datasets_file: bool = True,
                          use_coevolution_file: bool = False,
                          use_barplot_of_correlation_file: bool = False,
                          use_plot_distribution_of_correlation_file: bool = False,
                          use_plot_correlation_by_rate_bin_file: bool = False) -> Dict[str, str]:
        """
        Args (new/notable):
            event_threshold: the actual per-site candidate filter used for the coevolution table.
            p_value_mode: 'empirical' (default, rank-based against the raw per-bin null --
                resolution floored at 1/(n_in_bin + 1)) or 'beta' (fit a Beta(a, b) to each bin's
                null and read the tail off the fitted CDF, removing that floor; bins that don't
                fit/pass goodness-of-fit fall back to 'empirical' automatically -- see
                ``fit_beta_null``).
        """
        self.make_dir(file_path)
        if self.posterior_rates is None:
            self.set_posterior_rates_vector()

        result = {}
        file_simulated_datasets = f'{file_path}/SimulatedDatasets.fastas'
        file_coevolution = f'{file_path}/Coevolution.tsv'
        file_distribution_of_correlation = f'{file_path}/DistributionOfCorrelation.svg'
        file_barplot_of_correlation = f'{file_path}/BarplotOfCorrelation.svg'
        file_correlation_by_rate_bin = f'{file_path}/CorrelationByRateBin.svg'

        site_rate = np.asarray(self.posterior_rates, dtype=np.float64)
        branch_length = np.asarray([n.distance_to_father for n in self.all_nodes_objects[1:]], dtype=np.float64)

        a = 1.0 / (2 * (1 - self.pi_1))
        b = 1.0 / (2 * self.pi_1)
        mu = a + b
        t = branch_length[:, np.newaxis] * site_rate[np.newaxis, :] * self.coefficient_bl
        e = np.exp(-mu * t)
        p01 = a * (1 - e) / mu
        p11 = (a + b * e) / mu

        newick_text = self.get_newick(with_internal_nodes=True)

        null_pool, pairs, r_values, bins = {}, None, None, None

        if any((use_coevolution_file, use_barplot_of_correlation_file, use_plot_distribution_of_correlation_file,
                use_plot_correlation_by_rate_bin_file)):
            site_matrix, candidates, categories = self.identify_event_candidates(event_threshold)
            pairs, r_values, bins = self.get_bins(site_matrix, candidates, categories)

            null_pool = {current_bin: [] for current_bin in set(bins)}

        phylo_tree = Tree(newick_text, categories_quantity=self.categories_quantity, alpha=self.alpha,
                          pi_1=self.pi_1, coefficient_bl=self.coefficient_bl)

        for i in range(number_datasets):
            header = f'iterations = {i}'
            current_msa = self.generate_msa(msa_type=str,
                                            site_rate=self.posterior_rates,
                                            p01=p01, p11=p11,
                                            sites_quantity=self.msa_length,
                                            branch_length=branch_length,
                                            leaves=self.leaves_objects)

            phylo_tree.msa = phylo_tree.get_msa_dict(current_msa)
            phylo_tree.calculated_tree = phylo_tree.calculated_likelihood = False
            phylo_tree.calculate_tree()
            phylo_tree.set_posterior_rates_vector()

            if any((use_coevolution_file, use_barplot_of_correlation_file, use_plot_distribution_of_correlation_file,
                    use_plot_correlation_by_rate_bin_file)):
                current_site_matrix, current_candidates, current_categories = (
                    phylo_tree.identify_event_candidates(event_threshold))
                current_pairs, current_r_values, current_bins = phylo_tree.get_bins(current_site_matrix,
                                                                                    current_candidates,
                                                                                    current_categories)
                for current_bin, r_value in zip(current_bins, current_r_values):
                    null_pool.setdefault(current_bin, []).append(r_value)

            if use_simulated_datasets_file:
                current_content = f'{header}\n\n{current_msa}\n\n\n'
                with open(file_simulated_datasets, 'a', encoding='utf-8') as file:
                    file.write(current_content)

        if use_simulated_datasets_file:
            result.update({'Simulated datasets (fastas)': file_simulated_datasets})

        if any((use_coevolution_file, use_barplot_of_correlation_file, use_plot_distribution_of_correlation_file,
                use_plot_correlation_by_rate_bin_file)):
            pos1_list = []
            pos2_list = []
            rate_bin_list = []
            r_list = []
            p_value_list = []
            direction_list = []

            abs_null_by_bin = {k: np.abs(np.asarray(v, dtype=float)) for k, v in null_pool.items()}
            beta_by_bin: Dict[Any, Optional[Tuple[float, float]]] = {}
            pooled_beta = None
            if p_value_mode == 'beta' and abs_null_by_bin:
                pooled_signed = np.concatenate([np.asarray(v, dtype=float) for v in null_pool.values() if v])
                pooled_beta = self.fit_beta_null(pooled_signed, min_fit=min_fit) if pooled_signed.size else None

            for current_pair, current_r_value, current_bin in zip(pairs, r_values, bins):
                params = None
                if p_value_mode == 'beta':
                    if current_bin not in beta_by_bin:
                        beta_by_bin[current_bin] = self.fit_beta_null(
                            np.asarray(null_pool.get(current_bin, []), dtype=float), min_fit=min_fit)
                    params = beta_by_bin[current_bin] or pooled_beta
                if params is not None:
                    current_p_value = self.beta_p(current_r_value, params)
                else:
                    abs_null = abs_null_by_bin.get(current_bin, np.zeros(0))
                    n_bin = abs_null.size
                    current_p_value = ((1 + int(np.count_nonzero(abs_null >= abs(current_r_value)))) / (1 + n_bin)
                                       if n_bin else 1.0)
                pos1_list.append(current_pair[0])
                pos2_list.append(current_pair[1])
                r_list.append(current_r_value)
                rate_bin_list.append(current_bin)
                p_value_list.append(current_p_value)
                direction_list.append('co-occurrence' if current_r_value >= 0 else 'avoidance')

            p_values = np.asarray(p_value_list)
            q_values = self.benjamini_hochberg(p_values)
            df = pd.DataFrame({'POS1': np.asarray(pos1_list, np.int32),
                               'POS2': np.asarray(pos2_list, np.int32),
                               'r': np.round(np.asarray(r_list), decimals=14),
                               'rate-bin': rate_bin_list,
                               'p-value': np.round(p_values, decimals=14),
                               'q-value': np.round(q_values, decimals=14),
                               'direction': direction_list})
            df.sort_values(by=['q-value', 'p-value', 'r'], key=lambda x: x.abs() if x.name == 'r' else x, inplace=True,
                           ascending=[True, True, False])

            del pos1_list, pos2_list, rate_bin_list, r_list, p_value_list, direction_list
            gc.collect()

            r_data = df['r'].dropna().to_numpy()
            np.clip(r_data, -1 + eps2, 1 - eps2, out=r_data)
            alpha_est, beta_est, loc_est, scale_est = sp_beta.fit(r_data, floc=-1, fscale=2)
            info_text = f'Alpha: {alpha_est:.4f}\nBeta: {beta_est:.4f}'

            del r_data
            gc.collect()

            if use_coevolution_file:
                df.to_csv(file_coevolution, sep=sep, index=False, float_format='%.8f', chunksize=50000)
                result.update({'Table of coevolution (tsv)': file_coevolution})

            if use_plot_distribution_of_correlation_file:
                fig, ax = plt.subplots(figsize=(7, 5))
                ax.tick_params(axis='both', labelsize=8)
                sns.kdeplot(df['r'], fill=True, color='blue', ax=ax)
                ax.set_xlim(-1, 1)
                ax.set_xlabel('Correlation (r)')
                ax.set_ylabel('Density')
                ax.set_title('Distribution of original correlation coefficients')

                ax.text(0.05, 0.95, info_text,
                        transform=ax.transAxes,
                        fontsize=9,
                        verticalalignment='top',
                        bbox={'boxstyle': 'round,pad=0.5',
                              'facecolor': 'white',
                              'alpha': 0.6,
                              'edgecolor': 'gray'})

                fig.savefig(file_distribution_of_correlation, dpi=300, bbox_inches='tight')
                plt.close(fig)

                result.update({'Plot of distribution of correlation (svg)': file_distribution_of_correlation})

            if use_barplot_of_correlation_file:
                fig, ax = plt.subplots(figsize=(7, 5))
                ax.tick_params(axis='both', labelsize=8)
                sns.histplot(df['r'], bins=20, kde=True, color='green', ax=ax)
                ax.set_xlim(-1, 1)
                ax.set_xlabel('Correlation (r)')
                ax.set_ylabel('Count')
                ax.set_title('Barplot of original correlation coefficients')

                ax.text(0.05, 0.95, info_text,
                        transform=ax.transAxes,
                        fontsize=8,
                        verticalalignment='top',
                        bbox={'boxstyle': 'round,pad=0.5',
                              'facecolor': 'white',
                              'alpha': 0.6,
                              'edgecolor': 'gray'})

                fig.savefig(file_barplot_of_correlation, dpi=300, bbox_inches='tight')
                plt.close(fig)

                result.update({'Barplot of correlation (svg)': file_barplot_of_correlation})

            if use_plot_correlation_by_rate_bin_file:
                def make_clean_label(tup):
                    if isinstance(tup, tuple) and len(tup) == 2:
                        return f'({tup[0]:.0f}, {tup[1]:.0f})'
                    return str(tup)

                df['bin_clean'] = df['rate-bin'].apply(make_clean_label)
                unique_bins_sorted = sorted(df['rate-bin'].dropna().unique())
                categories_order = [make_clean_label(b) for b in unique_bins_sorted]
                bins_quantity = len(unique_bins_sorted)

                fig, ax = plt.subplots(figsize=(max(7, int(np.round(bins_quantity * 0.7))), 5))
                ax.tick_params(axis='both', labelsize=8)

                sns.stripplot(x='bin_clean',
                              y='r',
                              data=df,
                              order=categories_order,
                              palette='muted',
                              hue='bin_clean',
                              size=6,
                              jitter=0.15,
                              alpha=0.7,
                              ax=ax,
                              )

                ax.set_xlabel('Rate-bin categories')
                ax.set_ylabel('Correlation (r)')
                ax.set_title('Distribution of original correlation coefficients by rate-bin categories')

                fig.savefig(file_correlation_by_rate_bin, dpi=300, bbox_inches='tight')
                plt.close(fig)

                result.update({'Plot of correlation by rate-bin (svg)': file_correlation_by_rate_bin})

        return result

    def posterior_rates_to_tsv(self, file_name: str = 'PosteriorRates.tsv', sep: str = '\t') -> str:

        self.make_dir(file_name)
        if self.posterior_rates is None:
            self.set_posterior_rates_vector()

        df = pd.DataFrame({'POS': range(len(self.posterior_rates)),
                           'rate': self.posterior_rates})
        df.to_csv(file_name, sep=sep, index=False, float_format='%.8f', chunksize=50000)

        return file_name

    def pearson_correlation_to_tsv(self, file_name: str = 'PearsonCorrelation.tsv', sep: str = '\t',
                                   probability_lg: Union[float, np.float64] = 0.5,
                                   number_lg: Union[float, np.float64, int] = 1) -> str:

        self.make_dir(file_name)
        if self.correlation_vector is None:
            self.calculate_correlation(probability_lg=probability_lg, number_lg=number_lg)

        df = pd.DataFrame({'POS1': np.int32(self.correlation_vector[0]),
                           'POS2': np.int32(self.correlation_vector[1]),
                           'correlation': self.correlation_vector[2]})
        df.to_csv(file_name, sep=sep, index=False, float_format='%.8f', chunksize=50000)

        return file_name

    def attributes_to_tsv(self, file_name: str = 'TreeAttributes.tsv', sep: str = '\t') -> str:

        self.make_dir(file_name)
        data = loads(dumps({'π1 value': self.pi_1,
                            'Γ distribution α value': self.alpha,
                            'number of rate categories': self.categories_quantity,
                            'coefficient of branch lengths': self.coefficient_bl,
                            'rate vector': self.rate_vector,
                            'alphabet': self.alphabet,
                            'log likelihood': self.log_likelihood}, cls=NpEncoder))
        df = pd.DataFrame({k: ((v, ) if isinstance(v, (set, tuple, list)) else v) for k, v in data.items()
                           if v is not None}, index=[0])
        df.to_csv(file_name, sep=sep, index=False, float_format='%.8f', chunksize=50000)

        return file_name

    def parsimony_score_to_tsv(self, file_name: str = 'ParsimonyAndHomoplasyScores.tsv', sep: str = '\t') -> str:

        self.make_dir(file_name)
        df = self.get_parsimony_score()
        df.to_csv(file_name, sep=sep, index=False, float_format='%.8f', chunksize=50000)

        return file_name

    def likelihood_to_tsv(self, file_name: str = 'LogLikelihood.tsv', sep: str = '\t') -> str:

        self.make_dir(file_name)
        self.calculate_likelihood()
        df = pd.DataFrame({'POS': range(len(self.log_likelihood_vector)),
                           'log-likelihood': self.log_likelihood_vector})

        total_row = {'POS': 'TOTAL',
                     'log-likelihood': self.log_likelihood}

        df_with_total = pd.concat([df, pd.DataFrame([total_row])], ignore_index=True)
        df_with_total.to_csv(file_name, sep=sep, index=False, float_format='%.8f', chunksize=50000)

        return file_name

    def tree_to_tsv(self, file_name: str = 'Nodes.tsv', sep: str = '\t', mode: str = 'node',
                    taking_into_coefficient: bool = True, **kwargs) -> str:
        self.make_dir(file_name)
        columns, lists, decimals = kwargs.get('columns', None), kwargs.get('lists', None), kwargs.get('decimals', None)

        if columns is None or lists is None:
            columns, lists, decimals = self.get_columns(mode, columns, taking_into_coefficient)

        if kwargs.get('columns', None) is None:
            kwargs.update(columns=columns)
        if kwargs.get('lists', None) is None:
            kwargs.update(lists=lists)
        if kwargs.get('list_type', None) is None:
            kwargs.update(list_type=list)
        if kwargs.get('distance_type', None) is None:
            kwargs.update(distance_type=float)
        if kwargs.get('decimals', None) is None:
            kwargs.update(decimals=decimals)

        table = self.tree_to_table(taking_into_coefficient=taking_into_coefficient, **kwargs)

        table.to_csv(file_name, index=False, sep=sep, float_format='%.8f', chunksize=50000)

        return file_name

    def tree_to_newick_file(self, file_name: str = 'tree_file.tree', with_internal_nodes: bool = False,
                            decimal_length: int = 0, taking_into_coefficient: bool = True) -> str:
        newick_text = self.get_newick(with_internal_nodes, decimal_length, taking_into_coefficient)

        return self.write_file(file_name, newick_text)

    def tree_to_visual_format(self, file_name: str = 'VisualTree.svg', with_internal_nodes: bool = False,
                              file_extensions: Optional[Union[str, Tuple[str, ...]]] = None, show_axes: bool = False,
                              taking_into_coefficient: bool = True) -> Dict[str, str]:
        file_extensions = self.check_file_extensions_tuple(file_extensions, 'svg')

        self.make_dir(file_name)
        tmp_dir = Path(file_name).parent.joinpath('tmp')
        tmp_file = f'{tmp_dir.joinpath(f"{self.get_random_name()}.tree")}'
        self.make_dir(tmp_file)
        self.tree_to_newick_file(tmp_file, with_internal_nodes, taking_into_coefficient)
        phylogenetic_tree = Phylo.read(tmp_file, 'newick')

        j = file_name[::-1].find('.')
        file_names = dict()
        for file_extension in file_extensions:
            file_name = f'{file_name[:-(j + 1)]}.{file_extension}' if len(file_name) > j > -1 else (f'{file_name}.'
                                                                                                    f'{file_extension}')
            file_names.update({f'Newick tree ({file_extension})': file_name})
            if file_extension == 'txt':
                with open(file_name, 'w') as f:
                    Phylo.draw_ascii(phylogenetic_tree, f)
            else:
                Phylo.draw(phylogenetic_tree, do_show=False)
                plt.axis('on' if show_axes else 'off')
                kwargs = {'format': file_extension, 'bbox_inches': 'tight', 'dpi': 300} if (
                        file_extension == 'svg') else {'format': file_extension}
                plt.savefig(file_name, **kwargs)
                plt.close()
        rmtree(tmp_dir, ignore_errors=True)

        return file_names

    def tree_to_interactive_html(self, file_name: str = 'InteractiveTree.svg', taking_into_coefficient: bool = True
                                 ) -> str:
        self.make_dir(file_name)
        self.calculate_tree()
        self.calculate_ancestral_sequence()
        size_factor = min(1 + self.get_leaves_count() // 9, 6)
        columns, lists, decimals = self.get_columns(mode='tree_html', taking_into_coefficient=taking_into_coefficient)
        df = self.tree_to_table(columns=columns, distance_type=float, filters={'node_type': ['leaf', 'node', 'root']},
                                list_type=list, taking_into_coefficient=taking_into_coefficient, lists=lists,
                                decimals=decimals)
        df_copy = df.copy()
        del df['sequence'], df['node_type'], df['prob_characters']
        df = df.iloc[1:]

        d3 = D3Blocks(verbose='critical', chart='tree', frame=False)
        d3.set_node_properties(df)

        d3.font = {'size': 12}
        d3.hierarchy = [i for i in range(1, len(df_copy.T.count()) + 1)]
        d3.title = 'Phylogenetic tree'
        d3.filepath = file_name
        d3.figsize = (500, 500)
        d3.showfig, d3.overwrite, d3.reset_properties, d3.save_button = True, True, True, True
        d3.notebook = False
        d3.config = d3.chart.set_config(config=d3.config, filepath=d3.filepath, font=d3.font, title=d3.title,
                                        margin={'top': 20, 'right': 40, 'bottom': 20, 'left': 40},
                                        showfig=d3.showfig, overwrite=d3.overwrite, figsize=d3.figsize,
                                        reset_properties=d3.reset_properties, notebook=d3.notebook,
                                        hierarchy=d3.hierarchy, save_button=d3.save_button)

        colors = ['crimson', 'orangered', 'darkorange', 'gold', 'yellowgreen', 'forestgreen', 'mediumturquoise',
                  'dodgerblue', 'slateblue', 'darkviolet']
        colors_as = {'A': 'crimson', 'L': 'darkorange', 'G': 'forestgreen', 'P': 'slateblue'}
        for i in df_copy.T:
            ancestral_sequence = ''
            sequence = ''.join([self.draw_cell_html_table(colors[self.get_integer(self.get_shown_character(j))], j)
                                for j in df_copy['sequence'][i]])
            sequence = self.draw_row_html_table('Sequence', sequence)
            probability_coefficient = ''.join([self.draw_cell_html_table(colors[self.get_integer(j)], f'{j:.3f}')
                                              for j in df_copy['prob_characters'][i]])
            probability_coefficient = self.draw_row_html_table('Probability coefficient', probability_coefficient)
            if df_copy['node_type'][i] != 'root':
                ancestral_sequence = ''.join([self.draw_cell_html_table(colors_as[self.get_shown_character(j)], j)
                                              for j in df_copy['ancestral_sequence'][i]])
                ancestral_sequence = self.draw_row_html_table('Ancestral Comparison', ancestral_sequence)
            if df_copy['node_type'][i] != 'leaf':
                if df_copy['node_type'][i] == 'node':
                    d3.node_properties.get(df_copy['target'][i])['color'] = 'darkorange'
                    d3.node_properties.get(df_copy['target'][i])['size'] = 15 / size_factor
                if df_copy['node_type'][i] == 'root':
                    d3.node_properties.get(df_copy['target'][i])['color'] = 'firebrick'
                    d3.node_properties.get(df_copy['target'][i])['size'] = 20 / size_factor
            else:
                d3.node_properties.get(df_copy['target'][i])['color'] = 'forestgreen'
                d3.node_properties.get(df_copy['target'][i])['size'] = 10 / size_factor
            distance = f'<td class="h7 w-auto text-center">{df_copy["weight"][i]}</td>'
            info = (f'{self.draw_row_html_table("Distance", distance)}{sequence}{probability_coefficient}'
                    f'{ancestral_sequence}')
            d3.node_properties.get(df_copy['target'][i])['tooltip'] = self.draw_html_table(info)
            d3.font.update({'type': 'Anonymous Pro'})

        d3.set_edge_properties(df)
        d3.show()

        return file_name

    def tree_to_graph(self, file_name: str = 'graph.svg', file_extensions: Optional[Union[str, Tuple[str, ...]]] = None
                      ) -> Union[str, Dict[str, str]]:
        file_extensions = self.check_file_extensions_tuple(file_extensions, 'png')

        size_factor = min(1 + self.get_leaves_count() // 9, 6)
        self.make_dir(file_name)
        columns = {'node': 'Name', 'father_name': 'Parent', 'distance': 'Distance to parent'}
        table = self.tree_to_table(decimal_length=0, columns=columns)
        table = table.drop(0)
        j = file_name[::-1].find('.')
        file_name_sm = f'{file_name[:-j]}' if len(file_name) > j > -1 else f'{file_name}.'
        file_names = dict()
        for file_extension in file_extensions:
            file_name = f'{file_name_sm}{file_extension}'
            file_names.update({f'Graph ({file_extension})': file_name})
            graph = nx.Graph()
            for row in table.values:
                graph.add_edge(row[1], row[0], length=row[2] if row[2] else 0.0)
            if file_extension in ('svg', 'png'):
                nx.draw(graph, with_labels=True, font_color='Maroon', node_color='Gold', node_size=1000//size_factor,
                        font_size=12//size_factor, font_weight='bold')
                plt.savefig(file_name, **{'format': file_extension, 'bbox_inches': 'tight', 'dpi': 300})
                plt.close()
            if file_extension in ('dot', ):
                nx.drawing.nx_pydot.write_dot(graph, file_name)

        return file_names

    def optimize(self, func: Union[Callable, str], bracket: Tuple[Union[float, np.float64], ...] = (0.5,),
                 bounds: Tuple[Union[float, np.float64], ...] = (0.001, 0.999), args: Optional[Tuple[Any, ...]] = None,
                 result_fild: Optional[str] = None):
        func = self.__getattribute__(func) if isinstance(func, str) else func
        min_scalar = minimize_scalar(func, bracket=bracket, bounds=bounds) if args is None else (
            minimize_scalar(func, args=args, bracket=bracket, bounds=bounds))

        return min_scalar[result_fild] if result_fild else min_scalar

    def pi_optimization(self, pi: Union[float, np.float64], mode: int = 1) -> Union[float, np.float64]:
        current_pi = (pi, None)
        self.clean_all()
        self.set_pi(current_pi[mode], current_pi[::-1][mode])
        self.set_vars()

        return -self.get_log_likelihood()

    def alpha_optimization(self, alpha: Union[int, float, np.ndarray]) -> Union[float, np.float64]:
        self.clean_all()
        self.set_gamma_distribution_categories_vector(alpha)
        self.set_vars()

        return -self.get_log_likelihood()

    def coefficient_bl_optimization(self, coefficient_bl: Union[int, float, np.ndarray]) -> Union[float, np.float64]:
        self.clean_all()
        self.set_coefficient_bl(coefficient_bl)
        self.set_vars()

        return -self.get_log_likelihood()

    def optimize_coefficient_bl(self, is_optimize_bl: Optional[bool] = None) -> None:
        if is_optimize_bl:
            self.coefficient_bl = self.optimize(func=self.coefficient_bl_optimization, bracket=(1, ), bounds=(0.1, 10),
                                                result_fild='x')
            self.set_vars()

    def optimize_alpha(self, is_optimize_alpha: Optional[bool] = None) -> None:
        if is_optimize_alpha:
            self.alpha = self.optimize(func=self.alpha_optimization, bracket=(0.5, ), bounds=(0.1, 20), result_fild='x')
            self.set_vars()

    def optimize_pi(self, is_optimize_pi: Optional[bool] = None, is_optimize_pi_average: Optional[bool] = None,
                    mode: int = 1) -> None:
        if is_optimize_pi:
            self.pi_1 = self.optimize(func=self.pi_optimization, bracket=(0.5, ), bounds=(0.001, 0.999), args=(mode, ),
                                      result_fild='x')
            self.set_vars()

        elif is_optimize_pi_average:
            all_lines_list = list(self.msa.values())
            all_lines = ''.join(all_lines_list)
            known_characters_count = sum(all_lines.count(char) for char in self.alphabet)
            self.pi_1 = all_lines.count(self.alphabet[mode]) / known_characters_count
            self.set_vars()

    def clean_all(self):
        for current_node in self.all_nodes_objects:
            current_node.clean_all()
        self.log_likelihood_vector = self.likelihood_vector = self.correlation_vector = self.posterior_rates = None
        self.log_likelihood = self.likelihood = 0.0

    def set_all(self, categories_quantity: Optional[int] = None, alpha: Optional[float] = None,
                beta: Optional[float] = None, pi_0: Optional[Union[float, np.float64, int]] = None,
                pi_1: Optional[Union[float, np.float64, int]] = None,
                coefficient_bl: Optional[Union[float, np.float64, int]] = None) -> None:

        self.categories_quantity = categories_quantity
        self.set_alpha(alpha, beta)
        self.set_pi(pi_0, pi_1)
        self.set_coefficient_bl(coefficient_bl)
        self.set_gamma_distribution_categories_vector(self.alpha)
        self.set_vars()

    def set_gamma_distribution_categories_vector(self, alpha: Union[int, float, np.float64]) -> None:
        self.set_alpha(alpha)
        categories_vector = []
        gamma_percent_point = self.get_gamma_distribution_percent_point()
        for i in range(self.categories_quantity):
            lower_gamma_inc_1 = gammainc(self.alpha + 1, gamma_percent_point[i] * self.alpha)
            lower_gamma_inc_2 = gammainc(self.alpha + 1, gamma_percent_point[i + 1] * self.alpha)
            mean = (self.alpha / self.alpha) * (lower_gamma_inc_2 - lower_gamma_inc_1) / (1 / self.categories_quantity)
            categories_vector.append(mean)

        self.rate_vector = tuple(categories_vector)
        self.rate_vector_length = len(self.rate_vector)

    def set_coefficient_bl(self, coefficient_bl: Optional[Union[float, np.float64, int]] = None) -> None:
        self.coefficient_bl = 1.0 if coefficient_bl is None else coefficient_bl

    def set_alpha(self, alpha: Optional[float] = None, beta: Optional[float] = None) -> None:
        self.alpha = alpha if alpha else (beta if beta else 0.5)

    def set_pi(self, pi_0: Optional[Union[float, np.float64, int]] = None,
               pi_1: Optional[Union[float, np.float64, int]] = None) -> None:
        if pi_0:
            self.pi_1 = 1 - pi_0
        elif pi_1:
            self.pi_1 = pi_1
        else:
            self.pi_1 = 1 / self.alphabet_length

    def set_vars(self) -> None:
        if self.pi_1:
            frequency = (1 - self.pi_1, self.pi_1)
        else:
            frequency = (1 / self.alphabet_length, 1 / self.alphabet_length)
        for current_node in self.all_nodes_objects:
            current_node.alphabet = self.alphabet
            current_node.frequency = np.asarray(frequency, dtype=np.float64)
            current_node.pi_1 = self.pi_1
            current_node.coefficient_bl = self.coefficient_bl
            current_node.pmatrix = np.asarray([current_node.get_pmatrix(r) for r in self.rate_vector], dtype=np.float64)

    def get_gamma_distribution_percent_point(self) -> List[float]:
        probability_vector = np.linspace(0, 1, self.categories_quantity + 1)

        return gamma.ppf(probability_vector, a=self.alpha, scale=1/self.alpha)

    def set_basic_msa(self) -> None:
        self.msa = {leaf.name: self.alphabet[0] for leaf in self.leaves_objects}

    def generate_msa(self, msa_type: type = dict,
                     sites_quantity: int = 1,
                     site_rate: Optional[np.ndarray] = None,
                     p01: Optional[np.ndarray] = None,
                     p11: Optional[np.ndarray] = None,
                     branch_length: Optional[np.ndarray] = None,
                     branch_nodes: Optional[List[Node]] = None,
                     leaves: Optional[List[Node]] = None) -> Union[Dict[str, str], str]:

        if site_rate is None:
            site_rate = np.random.choice(self.rate_vector, 1 if sites_quantity is None else sites_quantity)
        if branch_nodes is None:
            branch_nodes = self.all_nodes_objects[1:]
        if branch_length is None:
            branch_length = np.asarray([n.distance_to_father for n in branch_nodes], dtype=np.float64)
        if p01 is None or p11 is None:
            a = 1.0 / (2 * (1 - self.pi_1))
            b = 1.0 / (2 * self.pi_1)
            mu = a + b
            t = branch_length[:, np.newaxis] * site_rate[np.newaxis, :] * self.coefficient_bl
            e = np.exp(-mu * t)
            p01 = a * (1 - e) / mu
            p11 = (a + b * e) / mu
        if leaves is None:
            leaves = self.leaves_objects

        states = {self.root.name: (np.random.random(sites_quantity) < self.pi_1).astype(np.int8)}
        for idx, node in enumerate(branch_nodes):
            parent_state = states[node.father.name]
            p = np.where(parent_state == 0, p01[idx], p11[idx])
            states[node.name] = (np.random.random(sites_quantity) < p).astype(np.int8)

        msa = {leaf.name: ''.join(self.alphabet[s] for s in states[leaf.name]) for leaf in leaves}

        return self.get_fasta_text(msa) if msa_type == str else msa

    def set_posterior_rates_vector(self, prior: Optional[np.ndarray] = None) -> None:
        if prior is None:
            prior = np.ones(self.rate_vector_length) / self.rate_vector_length
        prior = np.asarray(prior, dtype=np.float64)
        assert len(prior) == self.rate_vector_length, 'prior length must match number of rate categories'

        if not self.calculated_likelihood:
            self.calculate_up()

        likelihoods_per_rate = np.einsum('j,rji->ri', self.root.frequency, self.root.up_vector)
        invalid_mask = (likelihoods_per_rate <= 0.0) | np.isnan(likelihoods_per_rate)
        likelihoods_per_rate = np.where(invalid_mask, eps, likelihoods_per_rate)
        weighted = likelihoods_per_rate * prior[:, np.newaxis]
        weighted_sum = weighted.sum(axis=0)
        numerator = np.einsum('r,ri->i', self.rate_vector, weighted)
        posterior = np.divide(numerator, weighted_sum, where=(weighted_sum > 0), out=np.zeros_like(numerator))

        self.posterior_rates = posterior

    def set_pearson_correlation_vector(self, probability_lg: Union[float, np.float64] = 0.5,
                                       number_lg: Union[float, np.float64, int] = 1) -> None:
        if self.msa_length < 2:
            return None

        nodes_list = self.all_nodes_objects[1:]

        loss_vectors = np.array([node.probability_vector_loss for node in nodes_list])
        gain_vectors = np.array([node.probability_vector_gain for node in nodes_list])

        site_probs_matrix = np.stack((loss_vectors, gain_vectors), axis=1).reshape(-1, self.msa_length)

        condition_mask = site_probs_matrix >= probability_lg
        counts_per_site = np.sum(condition_mask, axis=0)
        unique_item = np.where(counts_per_site >= number_lg)[0]

        idx1, idx2 = np.triu_indices(len(unique_item), k=1)
        couples = np.column_stack((unique_item[idx1], unique_item[idx2]))

        filtered_probs = site_probs_matrix[:, unique_item]

        corr_matrix = np.corrcoef(filtered_probs.T)
        r_coefficients = corr_matrix[idx1, idx2]
        df = site_probs_matrix.shape[0] - 2
        valid_r_mask = np.abs(r_coefficients) < 1.0
        t_stat = np.zeros_like(r_coefficients)
        t_stat[valid_r_mask] = r_coefficients[valid_r_mask] * np.sqrt(df / (1.0 - r_coefficients[valid_r_mask] ** 2))

        p_values = np.zeros_like(r_coefficients)
        p_values[valid_r_mask] = distributions.t.sf(np.abs(t_stat[valid_r_mask]), df) * 2

        correlation_vector = np.zeros((4, len(couples)))
        correlation_vector[0] = unique_item[idx1]
        correlation_vector[1] = unique_item[idx2]
        correlation_vector[2] = r_coefficients
        correlation_vector[3] = p_values

        self.correlation_vector = correlation_vector

    def generate_site_rates(self, sites_quantity: int) -> np.ndarray:

        return np.random.choice(self.rate_vector, sites_quantity)

    def get_parsimony_score(self) -> pd.DataFrame:
        alphabet = self.alphabet
        leaves = self.leaves_objects
        n_leaves, msa_length = len(leaves), self.msa_length

        raw = np.frombuffer(''.join(self.msa[leaf.name] for leaf in leaves).encode('ascii'),
                            dtype=np.uint8).reshape(n_leaves, msa_length)

        has_state = [np.any(raw == ord(char), axis=0) for char in alphabet]
        m_vector = np.maximum(0, np.sum(has_state, axis=0) - 1)

        full_mask = (1 << len(alphabet)) - 1
        msa_matrix = np.full(raw.shape, full_mask, dtype=np.int32)
        for i, char in enumerate(alphabet):
            msa_matrix[raw == ord(char)] = 1 << i

        node_states: Dict[Node, np.ndarray] = {leaf: msa_matrix[i] for i, leaf in enumerate(leaves)}
        s_vector = np.zeros(self.msa_length, dtype=np.int32)

        for current_node in self.all_nodes_objects_post_order:
            children = current_node.children
            if not children:
                continue

            if len(children) == 2:
                first, second = node_states.pop(children[0]), node_states.pop(children[1])
                intersection, union = first & second, first | second
            else:
                children_matrix = np.array([node_states.pop(child) for child in children], dtype=np.int32)
                intersection = np.bitwise_and.reduce(children_matrix, axis=0)
                union = np.bitwise_or.reduce(children_matrix, axis=0)

            has_no_intersection = intersection == 0
            node_states[current_node] = np.where(has_no_intersection, union, intersection)
            s_vector += has_no_intersection.astype(np.int32)

        homoplasy_vector = s_vector - m_vector
        ci_vector = np.divide(m_vector, s_vector, out=np.ones_like(m_vector, dtype=np.float64),
                              where=(m_vector > 0) & (s_vector > 0))

        min_steps = np.sum(m_vector)
        parsimony_score = np.sum(s_vector)
        total_homoplasy = np.sum(homoplasy_vector)
        consistency_index = np.round(min_steps / parsimony_score, 4) if parsimony_score > 0 else 1.0

        df = pd.DataFrame({'POS': range(self.msa_length),
                           'Minimum Possible Number Of Steps (M)': m_vector,
                           'Parsimony Score (S)': s_vector,
                           'Homoplasy Score (H = S - M)': homoplasy_vector,
                           'Consistency Index (CI = M/S)': np.round(ci_vector, 4)})

        df['POS'] = df['POS'].astype(str)

        total_row = {'POS': 'TOTAL',
                     'Minimum Possible Number Of Steps (M)': min_steps,
                     'Parsimony Score (S)': parsimony_score,
                     'Homoplasy Score (H = S - M)': total_homoplasy,
                     'Consistency Index (CI = M/S)': consistency_index}

        df_with_total = pd.concat([df, pd.DataFrame([total_row])], ignore_index=True)

        return df_with_total

    @staticmethod
    def del_bootstrap_values(newick_text: str) -> str:
        pattern = r'(?<=\))([^:,()]+)(?=[:;])'

        raw_matches = re.findall(pattern, newick_text)
        if not raw_matches:
            return newick_text

        only_integers = [int(re.sub(r'^[\[\'](\d+)[\]\']$', r'\1', val)) for val in raw_matches if
                         val.isdigit() or re.sub(r'^[\[\'](\d+)[\]\']$', r'\1', val).isdigit()]

        is_random_node_ids = False
        if only_integers:
            list_length = len(only_integers)
            set_length = len(set(only_integers))
            if list_length == set_length:
                is_random_node_ids = True

        def replacer(match):
            val = match.group(1)

            if '/' in val or val in ['*', '#']:
                return ''

            if is_random_node_ids and val.isdigit():
                return val

            try:
                num = float(val)
                if num > 1000:
                    return val
                else:
                    return ''
            except ValueError:
                return val

        return re.sub(pattern, replacer, newick_text)

    @staticmethod
    def set_root_by_outgroup(newick_data: str, leaf_name: str) -> str:
        """
        Args:
            newick_data (str): A Newick formatted string representing the tree structure.
            leaf_name (str): leaf name.

        Returns:
            str: A Newick formatted string representing the tree structure.
        """
        phylo_tree = Phylo.read(StringIO(newick_data), 'newick')
        phylo_tree.root_with_outgroup(leaf_name)

        return ''.join(Writer((phylo_tree, )).to_strings(format_branch_length='%1.10f'))

    @staticmethod
    def calculate_mad_x(clade, dists: Dict[Any, Dict[Any, float]]) -> float:
        clade_leaves = clade.get_terminals()
        all_leaves = list(dists.keys())
        different_leaves = [Leaf for Leaf in all_leaves if Leaf not in clade_leaves]

        branch_length = clade.branch_length or 0

        numerator = 0
        denominator = 0

        for i in clade_leaves:
            d_in = clade.distance(i)

            for j in different_leaves:
                d_ij = dists[i][j]

                if d_ij != 0:
                    numerator += (d_ij - 2 * d_in) / (d_ij ** 2)
                    denominator += 2 / (d_ij ** 2)

        if denominator == 0:
            return branch_length / 2

        x = numerator / denominator

        return max(0, min(branch_length, x))

    @staticmethod
    def calculate_mad_score(clade, x: float, dists: Dict[Any, Dict[Any, float]]) -> float:
        clade_leaves = clade.get_terminals()
        all_leaves = list(dists.keys())
        different_leaves = [Leaf for Leaf in all_leaves if Leaf not in clade_leaves]

        deviations = []

        for i in clade_leaves:
            d_node_i = clade.distance(i)
            d_root_i = d_node_i + x

            for j in different_leaves:
                d_ij = dists[i][j]

                if d_ij != 0:
                    dev = (2 * d_root_i / d_ij) - 1
                    deviations.append(dev ** 2)

        if not deviations:
            return float('inf')

        return np.sqrt(np.mean(deviations))

    @staticmethod
    def calculate_mvr_x(phylo_tree, clade, dists: Dict[Any, Dict[Any, float]]) -> float:
        clade_leaves = clade.get_terminals()
        all_leaves = list(dists.keys())
        different_leaves = [Leaf for Leaf in all_leaves if Leaf not in clade_leaves]

        branch_length = clade.branch_length or 0

        n = len(all_leaves)
        n1 = len(clade_leaves)

        mu1 = np.mean([clade.distance(Leaf) for Leaf in clade_leaves])
        mu2 = np.mean([phylo_tree.distance(clade, Leaf) for Leaf in different_leaves])

        beta = n1 / n

        x = (mu2 - mu1 + branch_length * (1 - 2 * beta)) / (2 * (1 - beta))

        return max(0, min(branch_length, x))

    @staticmethod
    def calculate_mvr_score(phylo_tree, clade, x: float, dists: Dict[Any, Dict[Any, float]]) -> float:
        clade_leaves = clade.get_terminals()
        all_leaves = list(dists.keys())

        ref_leaf = clade_leaves[0]

        root_distances = []

        for leaf in all_leaves:
            if leaf in clade_leaves:
                d_root = clade.distance(leaf) + x
            else:
                d_root = phylo_tree.distance(ref_leaf, leaf) - x
            root_distances.append(d_root)

        return np.var(root_distances)

    @staticmethod
    def build_symmetric_newick(num_taxa: int, branch_length: Union[float, np.float64, int] = 0.5) -> str:
        """
        Recursively build a perfectly symmetric binary Newick tree.
        If num_taxa is less than 2 or not a power of 2,
        it will be automatically rounded down to the nearest power of 2 (minimum 2).
        Every branch (leaf and internal) has length branch_length
        """
        num_taxa = 2 if num_taxa < 2 else 1 << (num_taxa.bit_length() - 1)

        leaves = [f't{i + 1}' for i in range(num_taxa)]

        def build_subtree(taxa):
            taxa_length = len(taxa)
            if taxa_length == 1:
                return f'{taxa[0]}:{branch_length}'

            mid = taxa_length // 2

            return f'({build_subtree(taxa[:mid])},{build_subtree(taxa[mid:])}):{branch_length}'

        middle = num_taxa // 2

        return f'({build_subtree(leaves[:middle])},{build_subtree(leaves[middle:])});'

    @staticmethod
    def set_root_by_midpoint(newick_data: str) -> str:
        """
        Args:
            newick_data (str): A Newick formatted string representing the tree structure.

        Returns:
            str: A Newick formatted string representing the tree structure.
        """

        phylo_tree = Phylo.read(StringIO(newick_data), 'newick')
        if len(phylo_tree.root.clades) > 2:
            phylo_tree.root_at_midpoint()

        return ''.join(Writer((phylo_tree, )).to_strings(format_branch_length='%1.10f'))

    @staticmethod
    def get_round(obj: Union[int, float, np.float64, np.ndarray], decimals: int = 4) -> float:

        return float(np.round(obj, decimals))

    @staticmethod
    def get_list_decimals(obj: Any, list_type: type = str, decimals: int = 4,
                          is_exception: bool = False, force_list: bool = False) -> Any:
        if isinstance(obj, (int, float, np.float64)):
            return round(obj, decimals)

        if isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.number):
                obj = obj.round(decimals).tolist()
            elif obj.dtype == object:
                vectorized_round = np.vectorize(
                    lambda x: Tree.get_list_decimals(x, list_type, decimals, is_exception, force_list),
                    otypes=[object]
                )
                obj = vectorized_round(obj).tolist()
            else:
                obj = obj.tolist()

        if isinstance(obj, str):
            if len(obj) <= 1:
                return obj

            if force_list:
                obj = obj.split() if ' ' in obj.strip() else list(obj)

        if isinstance(obj, (list, tuple, set)):
            processed = [Tree.get_list_decimals(x, list_type, decimals, is_exception, force_list=False) for x in obj]

            target_type = list_type if list_type in (list, tuple, set) else list
            return target_type(processed)

        return obj

    @staticmethod
    def is_bootstrap_value(number_str: str, lower: Union[float, np.float64, int] = 0,
                           upper: Union[float, np.float64, int] = 100) -> bool:
        re_result = bool(re.fullmatch(r'^-?\d+(\.\d+)?$', number_str)) if len(number_str.strip()) else False

        return lower <= float(number_str) <= upper if re_result else False

    @staticmethod
    def get_columns(mode: str = 'node', columns: Optional[Dict[str, str]] = None,
                    taking_into_coefficient: bool = True) -> Tuple[Dict[str, str], Tuple[str, ...], int]:

        suffix = '_taking_into_coefficient' if taking_into_coefficient else ''
        distance_name = f'distance{suffix}'

        lists = ('children', 'full_distance', 'full_distance_taking_into_coefficient', 'up_vector', 'down_vector',
                 'marginal_vector', 'marginal_bl_vector', 'probability_vector', 'probabilities_sequence_characters',
                 'log_likelihood_vector', 'likelihood_vector', 'probability_vector_gain', 'probability_vector_loss',
                 'ancestral_sequence', 'sequence')
        decimals = 4
        if mode == 'node':
            columns = columns if columns else {'node': 'Name', 'node_type': 'Node type', distance_name:
                                               'Distance to parent', 'sequence': 'Sequence',
                                               'probabilities_sequence_characters': 'Probability coefficient',
                                               'ancestral_sequence': 'Ancestral Comparison'}
            lists = ('probabilities_sequence_characters', 'sequence', 'ancestral_sequence')
        elif mode == 'branch':
            columns = columns if columns else {'father_name': 'Parent node', 'node': 'Child node', distance_name:
                                               'Branch length', 'probability_vector_gain': 'Gain probability',
                                               'probability_vector_loss': 'Loss probability'}
            lists = ('probability_vector_gain', 'probability_vector_loss')
        elif mode == 'node_tsv':
            columns = columns if columns else {'node': 'Name', 'father_name': 'Parent', distance_name:
                                               'Distance to parent', 'children': 'Children', 'sequence': 'Sequence',
                                               'pmatrix': 'P matrix',
                                               'up_vector': 'Up vector',
                                               'down_vector': 'Down vector',
                                               'marginal_vector': 'Marginal vector',
                                               'probabilities_sequence_characters': 'Probability coefficient',
                                               'ancestral_sequence': 'Ancestral comparison',
                                               'sequence_likelihood':
                                               'Current conditional likelihood (subtree) of sequence',
                                               'log_likelihood': 'Current conditional log-likelihood (ln, subtree)',
                                               'log_likelihood_vector':
                                               'Vector of current conditional log-likelihood (ln, subtree)'}
            lists = ('children', 'pmatrix', 'up_vector', 'down_vector', 'marginal_vector',
                     'probabilities_sequence_characters')
            decimals = 8
        elif mode == 'branch_tsv':
            columns = columns if columns else {'father_name': 'Parent node', 'node': 'Child node', distance_name:
                                               'Branch length', 'probability_vector_gain': 'Gain probability',
                                               'probability_vector_loss': 'Loss probability',
                                               'branch_probability_vector': 'Branch probability vector'}
            lists = ('branch_probability_vector', 'probability_vector_gain', 'probability_vector_loss')
            decimals = 8
        elif mode == 'tree_html':
            columns = columns if columns else {'node': 'target', 'father_name': 'source', distance_name: 'weight',
                                               'sequence': 'sequence', 'probabilities_sequence_characters':
                                               'prob_characters', 'node_type': 'node_type', 'ancestral_sequence':
                                               'ancestral_sequence'}
            lists = ('probabilities_sequence_characters', 'sequence', 'ancestral_sequence')

        return columns, lists, decimals

    @staticmethod
    def write_file(file_name: str, file_text: str) -> str:
        try:
            Tree.make_dir(file_name)
            with open(file_name, 'w') as f:
                f.write(file_text)
        except Exception as e:
            print(f'An error occurred while saving the file: {e}')
            file_name = ''

        return file_name

    @staticmethod
    def get_unknown_characters() -> Tuple[str, ...]:

        return '?',

    @staticmethod
    def get_columns_list_for_sorting(mode: str = 'node') -> Dict[str, List[str]]:
        if mode == 'node':
            result = {'List for sorting': ['Name', 'Node type', 'Distance to parent', 'Sequence',
                                           'Probability coefficient', 'Ancestral Comparison']}
        else:
            result = {'List for sorting': ['Parent node', 'Child node', 'Branch length', 'Gain probability',
                                           'Loss probability']}

        return loads(dumps(result, cls=NpEncoder).replace(f'\'', r'"'))

    @staticmethod
    def get_random_name(lenght: int = 24) -> str:
        abc_list = [_ for _ in 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz1234567890']

        return ''.join(np.random.choice(abc_list, lenght))

    @staticmethod
    def get_ancestral_alphabet() -> Tuple[str, ...]:

        return 'A', 'L', 'G', 'P'

    @staticmethod
    def get_alphabet(search_argument: Optional[Union[Set[str], int, str]] = None) -> Tuple[str, ...]:
        alphabets = (('0', '1'), ('A', 'C', 'G', 'T'),
                     ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                      'V'))
        if not search_argument:
            return tuple(alphabets[0])
        if isinstance(search_argument, int):
            return tuple(alphabets[search_argument])
        if isinstance(search_argument, str):
            search_argument = set([i for i in search_argument])
        if isinstance(search_argument, set):
            for alphabet in alphabets:
                if not search_argument - set(alphabet):
                    return alphabet

    @staticmethod
    def find_dict_in_iterable(iterable: List[Union[Dict[str, Union[float, np.float64, bool, str, np.ndarray,
                                             List[float], List[np.float64]]], Node]], key: str,
                              value: Optional[Union[float, bool, str, List[float]]] = None
                              ) -> Dict[str, Union[float, np.float64, bool, str, np.ndarray, List[float],
                                        List[np.float64]]]:
        for index, dictionary in enumerate(iterable):
            if key in dictionary and (True if value is None else dictionary[key] == value):
                return dictionary

    @staticmethod
    def make_dir(file_path: str, **kwargs) -> None:
        dir_path = Path(file_path).parent
        if not dir_path.exists():
            dir_path.mkdir(mode=kwargs.get('mode', 0o777), parents=kwargs.get('parents', True),
                           exist_ok=kwargs.get('exist_ok', True))

    @staticmethod
    def check_file_extensions_tuple(file_extensions: Optional[Union[str, Tuple[str, ...]]] = None, default_value: str =
                                    'txt') -> Tuple[str, ...]:
        file_extensions = file_extensions if file_extensions else (default_value,)
        if isinstance(file_extensions, str):
            file_extensions = (file_extensions,)

        return file_extensions

    @staticmethod
    def check_newick(newick_text: str) -> bool:
        newick_text = newick_text.strip()

        return newick_text and newick_text.startswith('(') and newick_text.endswith(';')

    @staticmethod
    def __set_node(node_str: str, num: Callable) -> Node:
        if node_str.find(':') > -1:
            node_data: List[Union[str, int, float]] = node_str.split(':')
            node_data[0] = node_data[0] if node_data[0] else 'nd' + str(num()).rjust(4, '0')
            try:
                node_data[1] = float(node_data[1])
            except ValueError:
                node_data[1] = 0.0
        else:
            node_data = [node_str if node_str else 'nd' + str(num()).rjust(4, '0'), 0.0]

        newick_node = Node(node_data[0])
        newick_node.distance_to_father = float(node_data[1])

        return newick_node

    @staticmethod
    def __counter():
        count = 0

        def sub_function():
            nonlocal count
            count += 1
            return count

        return sub_function

    @staticmethod
    def get_row_correlations(matrix: np.ndarray) -> np.ndarray:
        centered = matrix - matrix.mean(axis=1, keepdims=True)
        norms = np.sqrt((centered ** 2).sum(axis=1))

        return (centered @ centered.T) / np.outer(norms, norms)

    @staticmethod
    def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
        order = np.argsort(p_values)
        ranked = p_values[order]
        scale = len(p_values) / np.arange(1, len(p_values) + 1)
        q_sorted = np.minimum.accumulate((ranked * scale)[::-1])[::-1]
        q = np.empty_like(q_sorted)
        q[order] = np.clip(q_sorted, 0, 1)

        return q

    @staticmethod
    def empirical_p(r: Any, bin_key: str, null_pool: Any) -> Union[np.float64, float]:
        if bin_key not in null_pool or len(null_pool[bin_key]) == 0:
            return 1.0

        null_distribution = np.asarray(null_pool[bin_key], dtype=float)
        n = null_distribution.size
        count_extreme = np.sum(np.abs(null_distribution) >= np.abs(r))
        p_value = (count_extreme + 1) / (n + 1)

        return p_value

    @staticmethod
    def fit_beta_null(null_r: np.ndarray, min_fit: int = 100, gof_alpha: float = 0.01) -> Optional[Tuple[float, float]]:
        """Returns the shape parameters, or `None` if the sample is too small to fit reliably, its
        variance is incompatible with any Beta on [-1, 1] (can happen for a degenerate/near-point
        sample), or the fit fails a goodness-of-fit check against its own bulk. A bad fit must
        never silently produce an overconfident tail p-value, so the caller falls back to
        `empirical_p` in that case."""
        null_r = np.asarray(null_r, dtype='float64')
        null_r = null_r[np.isfinite(null_r)]
        if null_r.size < min_fit:
            return None
        y = np.clip((null_r + 1.0) / 2.0, eps2, 1 - eps2)
        mean, var = y.mean(), y.var(ddof=1)
        max_var = mean * (1 - mean)
        if not (0.0 < var < max_var):
            return None
        common = max_var / var - 1.0
        a, b = mean * common, (1 - mean) * common
        if not (np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0):
            return None
        _, p_gof = kstest(null_r, 'beta', args=(a, b, -1, 2))

        return None if p_gof < gof_alpha else (a, b)

    @staticmethod
    def beta_p(r: Any, params: Tuple[float, float]) -> Union[np.float64, float]:
        """P(|R| >= |r|) under a null Beta(a, b, loc=-1, scale=2) fitted by ``fit_beta_null`` --
        the continuous analogue of ``empirical_p``'s exceedance count, same [-1, 1] support and
        the same |r| definition of 'at least as extreme'."""
        a, b = params
        r_abs = min(abs(r), 1 - eps2)

        return float(sp_beta.cdf(-r_abs, a, b, loc=-1, scale=2) + sp_beta.sf(r_abs, a, b, loc=-1, scale=2))

    @staticmethod
    def get_integer(data: Union[str, int, float]) -> int:
        result = float(data) * 10

        return int(result - 1 if result == 10 else result)

    @staticmethod
    def draw_html_table(data: str) -> str:

        return f'<table class="w-97 p-4 tooltip">{data}</table>'

    @staticmethod
    def draw_row_html_table(name: str, data: str) -> str:

        return f'<tr><th class="p-2 h7 ">{name}:</th><th>{data}</td></th></tr>'

    @staticmethod
    def draw_cell_html_table(color: str, data: str) -> str:

        return f'<td style="color: {color}" class="h7 w-auto text-center">{data}</td>'

    @staticmethod
    def get_shown_character(cell: str) -> str:

        return cell[-2] if cell.endswith(')') else cell

    @classmethod
    def set_root_by_minimum(cls, newick_data: str, rooting_method: str) -> str:
        """
        Args:
            newick_data (str): A Newick formatted string representing the tree structure.
            rooting_method (str): `mad` (Minimal Ancestor Deviation), `mvr`(Minimum Variance Rooting)

        Returns:
            str: A Newick formatted string representing the tree structure.
        """
        phylo_tree = Phylo.read(StringIO(newick_data), 'newick')
        all_leaves = phylo_tree.get_terminals()
        dists = {leaf1: {leaf2: phylo_tree.distance(leaf1, leaf2) for leaf2 in all_leaves} for leaf1 in all_leaves}

        best_score = float('inf')
        best_clade = None
        best_x = 0

        for clade in phylo_tree.find_clades():
            if clade != phylo_tree.root:
                if rooting_method == 'mad':
                    x_opt = cls.calculate_mad_x(clade, dists)
                    score = cls.calculate_mad_score(clade, x_opt, dists)
                else:
                    x_opt = cls.calculate_mvr_x(phylo_tree, clade, dists)
                    score = cls.calculate_mvr_score(phylo_tree, clade, x_opt, dists)

                if score < best_score:
                    best_score = score
                    best_clade = clade
                    best_x = x_opt

        if best_clade and best_clade != phylo_tree.root:
            remaining_dist = max(0.0, (best_clade.branch_length or 0.0) - best_x)
            outgroup_leaves = best_clade.get_terminals()
            og = outgroup_leaves if len(outgroup_leaves) > 1 else outgroup_leaves[0]
            phylo_tree.root_with_outgroup(og, outgroup_branch_length=best_x)
            for clade in phylo_tree.root.clades:
                if clade != best_clade:
                    clade.branch_length = remaining_dist

        return ''.join(Writer((phylo_tree, )).to_strings(format_branch_length='%1.10f'))

    @classmethod
    def get_alphabet_from_dict(cls, msa_dict: Dict[str, str]) -> Tuple[str, ...]:
        character_list = []
        for sequence in msa_dict.values():
            character_list += [i for i in sequence]

        return cls.get_alphabet(set(character_list) - set(cls.get_unknown_characters()))

    @classmethod
    def check_tree(cls, phylo_tree: Union[str, 'Tree']) -> 'Tree':
        if isinstance(phylo_tree, str):
            phylo_tree = cls(phylo_tree)

        return phylo_tree

    @classmethod
    def rename_nodes(cls, phylo_tree: Union[str, 'Tree'], node_name: str = 'N', fill_character: str = '0',
                     number_length: int = 0) -> 'Tree':
        phylo_tree = cls.check_tree(phylo_tree)
        nodes_list = phylo_tree.nodes_objects
        num = phylo_tree.__counter()
        for current_node in nodes_list:
            if re.fullmatch(r'^nd\d{4}$', current_node.name):
                current_node.name = f'{node_name}{str(num()).rjust(number_length, fill_character)}'
            # if re.search(r'^[\'\"](.*)[\'\"]$', current_node.name):
            #     current_node.name = re.sub(r'^[\'\"](.*)[\'\"]$', r'[\1]', current_node.name)
            if re.search(r'[\'\"]', current_node.name):
                current_node.name = re.sub(r'[\'\"]', '', current_node.name)

        phylo_tree.all_nodes = {current_node.name: current_node for current_node in phylo_tree.all_nodes_objects}

        return phylo_tree

    @classmethod
    def compute_correlation(cls, num_taxa: int = 8,
                            sites_quantity: int = 100,
                            categories_quantity: int = 4,
                            alpha: float = 0.5,
                            pi_1: Union[float, np.float64, int] = 0.5,
                            branch_lengths: Union[float, np.float64, int] = 0.5,
                            seed: Optional[int] = None,
                            newick_text: Optional[str] = None,
                            fasta_text: Optional[str] = None) -> Tuple[float, float, np.asarray, np.ndarray]:
        if not newick_text:
            newick_text = cls.build_symmetric_newick(num_taxa, branch_lengths)

        newick_tree = cls(newick_text)

        if fasta_text:
            msa = newick_tree.get_msa_dict(fasta_text)
            sites_quantity = len(next(iter(msa.values())))

        tree_data = {'pi_1': pi_1,
                     'alpha': alpha,
                     'categories_quantity': categories_quantity,
                     'seed': seed}
        newick_tree.set_tree_data(**tree_data)

        true_rates = newick_tree.generate_site_rates(sites_quantity)

        print(f'\ttrue_rates: {[round(float(r), 4) for r in true_rates]}')

        if not fasta_text:
            fasta_text = newick_tree.generate_msa(msa_type=str, site_rate=true_rates)

        gloome_tree = cls(newick_text, msa=fasta_text, **tree_data)

        print(f'\trate_vector (4 Gamma categories): {[round(float(r), 4) for r in gloome_tree.rate_vector]}')

        gloome_tree.set_posterior_rates_vector()
        print(true_rates, gloome_tree.posterior_rates, sep='\n')
        r_val, p_val = pearsonr(true_rates, gloome_tree.posterior_rates)

        print(f'\tPearson r = {r_val:.4f}  (p = {p_val:.3e})\n')

        return r_val, p_val, true_rates, gloome_tree.posterior_rates

    @classmethod
    def generate_scatter_plot(cls, taxa_list: Union[List[int], Tuple[int, ...], np.ndarray],
                              sites_quantity: int = 100,
                              categories_quantity: int = 4,
                              alpha: float = 0.5,
                              pi_1: Union[float, np.float64, int] = 0.5,
                              branch_lengths: Union[float, np.float64, int] = 0.5,
                              seed: Optional[int] = None,
                              newick_text: Optional[str] = None,
                              fasta_text: Optional[str] = None,
                              out_path: Optional[str] = None) -> None:

        print(f'Correlation estimation. Args: alpha={alpha}, sites={sites_quantity}, branch lengths={branch_lengths}, '
              f'π1={pi_1}, Gamma categories={categories_quantity}, seed for randomizer={seed} \n')

        results = {}
        for num_taxa in taxa_list:
            print(f"\nN = {num_taxa} taxa")
            results[num_taxa] = cls.compute_correlation(num_taxa, sites_quantity, categories_quantity, alpha, pi_1,
                                                        branch_lengths, seed, newick_text, fasta_text)

        fig, axes = plt.subplots(2, 3, figsize=(10, 8))
        fig.suptitle(f'True vs. estimated site rates  (alpha={alpha}, sites={sites_quantity}, '
                     f'branch lengths={branch_lengths})', fontsize=13)

        for ax, num_taxa in zip(axes.flat, taxa_list):
            r_val, p_val, true_rates, est_rates = results[num_taxa]
            ax.scatter(true_rates, est_rates, alpha=0.55, s=25, color='steelblue', edgecolors='none')
            lim = [0, max(true_rates.max(), est_rates.max()) * 1.05]
            ax.plot(lim, lim, 'r--', lw=1, alpha=0.6, label='y = x')
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            ax.set_xlabel('True rate', fontsize=10)
            ax.set_ylabel('E(r|D)  posterior mean', fontsize=10)
            ax.set_title(f'Taxa = {num_taxa};  Pearson r = {r_val:.3f}  (p = {p_val:.2e})', fontsize=9)
            ax.legend(fontsize=8)

        plt.tight_layout()
        if out_path:
            plt.savefig(out_path, dpi=150, bbox_inches='tight')
            print(f"Plot saved → {out_path}")
        else:
            plt.show()

    @classmethod
    def set_root(cls, newick_data: str, rooting_method: str = 'midpoint', leaf: Optional[Union[str, Node]] = None
                 ) -> str:
        """
        Args:
            newick_data (str): A Newick formatted string representing the tree structure.
            rooting_method (str, optional): `mad` (Minimal Ancestor Deviation), `mvr`(Minimum Variance Rooting),
            `midpoint` (Midpoint Rooting, default), `outgroup` (Outgroup Rooting)
            leaf (str, Node, optional): `None` (default)

        Returns:
            str: A Newick formatted string representing the tree structure.
        """
        phylo_tree = cls(newick_data)
        if len(phylo_tree.root.children) > 2:
            rooting_method = rooting_method.strip().lower()
            if leaf and rooting_method == 'outgroup':
                newick_data = cls.set_root_by_outgroup(newick_data, leaf.name if isinstance(leaf, Node) else leaf)
            else:
                if rooting_method in ('mad', 'mvr'):
                    newick_data = cls.set_root_by_minimum(newick_data, rooting_method)
                else:
                    newick_data = cls.set_root_by_midpoint(newick_data)
            phylo_tree = cls(newick_data)
            for current_node in phylo_tree.all_nodes_objects[1:]:
                if current_node.distance_to_father == 0:
                    current_node.distance_to_father = eps2

            return phylo_tree.get_newick()

        return newick_data

    @classmethod
    def __get_html_tree(cls, structure: dict, status: str) -> str:
        tags = (f'<details {status}>', '</details>', '<summary>', '</summary>') if structure['children'] else ('', '',
                                                                                                               '', '')
        str_html = (f'<li> {tags[0]}{tags[2]}{structure["name"].strip()} \t ({structure["distance_to_father"]}) '
                    f'{tags[3]}')
        for child in structure['children']:
            str_html += f'<ul>{cls.__get_html_tree(child, status)}</ul>\n' if child[
                'children'] else f'{cls.__get_html_tree(child, status)}'
        str_html += f'{tags[1]}</li>'

        return str_html

    @classmethod
    def get_robinson_foulds_distance(cls, tree1: Union['Tree', str], tree2: Union['Tree', str]) -> float:
        tree1 = cls(tree1) if type(tree1) is str else tree1
        tree2 = cls(tree2) if type(tree2) is str else tree2

        edges_list1 = sorted(tree1.get_edges_list(), key=lambda item: item[1])
        edges_list2 = sorted(tree2.get_edges_list(), key=lambda item: item[1])

        distance = 0
        for newick_node in edges_list1:
            distance += 0 if edges_list2.count(newick_node) else 1
        for newick_node in edges_list2:
            distance += 0 if edges_list1.count(newick_node) else 1

        return distance / 2

    @classmethod
    def structure_to_html_tree(cls, structure: dict, styleclass: str = '', status: str = '') -> str:

        return (f'<ul {f" class = {chr(34)}{styleclass}{chr(34)}" if styleclass else ""}>'
                f'{cls.__get_html_tree(structure, status)}</ul>')

    @classmethod
    def subtree_to_structure(cls, newick_node: Node) -> Dict[str, str]:
        dict_node = {'name': newick_node.name.strip(), 'distance_to_father': newick_node.distance_to_father}
        list_children = []
        if newick_node.children:
            for child in newick_node.children:
                list_children.append(cls.subtree_to_structure(child))
        dict_node.update({'children': list_children})

        return dict_node
