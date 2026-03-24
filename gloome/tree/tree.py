import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import re

from shutil import rmtree
from json import loads, dumps
from pathlib import Path
from d3blocks import D3Blocks
from typing import Optional, List, Union, Dict, Tuple, Set, Any, Callable
from Bio import Phylo
from scipy.stats import gamma
from scipy.special import gammainc
from scipy.optimize import minimize_scalar

from gloome.tree.node import Node
from gloome.jsonNpEncoder.npencoder import NpEncoder as npEncode


class Tree:
    root: Optional[Node] = None
    alphabet: Optional[Tuple[str, ...]] = None
    msa: Optional[Dict[str, str]] = None
    rate_vector: Tuple[Union[float, np.ndarray, int], ...] = (1.0, )
    alpha: Optional[Union[float, np.ndarray, int]] = None
    categories_quantity: Optional[int] = None
    pi_1: Optional[Union[float, np.ndarray, int]] = None
    coefficient_bl: Optional[Union[float, np.ndarray, int]] = 1,
    log_likelihood_vector: Optional[List[Union[float, np.ndarray]]] = None
    log_likelihood: Optional[Union[float, np.ndarray]] = None
    likelihood: Optional[Union[float, np.ndarray]] = None
    calculated_ancestor_sequence: bool = False
    calculated_tree: bool = False
    calculated_likelihood: bool = False

    def __init__(self, data: Optional[Union[str, Node]] = None, node_name: Optional[str] = None, **kwargs) -> None:
        """
        Parameters
        ----------
            data: Optional[Union[str, Node]] = None
            node_name: Optional[str] = None
            msa: Optional[Union[Dict[str, str], str]] = None
            categories_quantity: Optional[float] = None
            alpha: Optional[float] = None
            beta: Optional[float] = None
            pi_0: Optional[Union[float, np.ndarray, int]] = None
            pi_1: Optional[Union[float, np.ndarray, int]] = None
            coefficient_bl: Optional[Union[float, np.ndarray, int]] = None
            is_optimize_pi: Optional[bool] = None,
            is_optimize_pi_average: Optional[bool] = None
            is_optimize_alpha: Optional[bool] = None
            is_optimize_bl: Optional[bool] = None
    """
        available_parameters = {'data', 'node_name', 'msa', 'categories_quantity', 'alpha', 'beta', 'pi_0', 'pi_1',
                                'coefficient_bl', 'is_optimize_pi', 'is_optimize_pi_average', 'is_optimize_alpha',
                                'is_optimize_bl'}
        invalid_parameters = set(kwargs.keys()) - available_parameters
        for key in invalid_parameters:
            del kwargs[key]

        if isinstance(data, str):
            self.newick_to_tree(data)
            if node_name and isinstance(node_name, str):
                Tree.rename_nodes(self, node_name)
        elif isinstance(data, Node):
            self.root = data
        else:
            self.root = Node('root')

        self.msa, self.alphabet, self.categories_quantity, self.alpha = None, None, None, None
        self.rate_vector = (1.0,)
        self.pi_1, self.coefficient_bl = None, 1
        self.log_likelihood_vector, self.log_likelihood, self.likelihood = None, None, None
        self.calculated_ancestor_sequence = self.calculated_tree = self.calculated_likelihood = False

        if any(kwargs.values()):
            self.set_tree_data(**kwargs)
        elif invalid_parameters:
            print(f'There are invalid parameters: {", ".join(invalid_parameters)}')

    def __str__(self) -> str:
        return self.get_newick()

    def __dir__(self) -> List[str]:
        return ['root', 'alphabet', 'msa', 'rate_vector', 'alpha', 'categories_quantity', 'pi_1', 'coefficient_bl',
                'log_likelihood_vector', 'log_likelihood', 'likelihood', 'calculated_ancestor_sequence',
                'calculated_tree', 'calculated_likelihood']

    def __dict__(self) -> [Optional[Node], Optional[Tuple[str, ...]], Optional[Dict[str, str]],
                           Tuple[Union[float, np.ndarray, int], ...], Optional[Union[float, np.ndarray, int]],
                           Optional[int], Optional[Union[float, np.ndarray, int]],
                           Optional[Union[float, np.ndarray, int]], Optional[List[Union[float, np.ndarray]]],
                           Optional[Union[float, np.ndarray]], Optional[Union[float, np.ndarray]], bool, bool, bool]:
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
                'likelihood': self.likelihood,
                'calculated_ancestor_sequence': self.calculated_ancestor_sequence,
                'calculated_tree': self.calculated_tree,
                'calculated_likelihood': self.calculated_likelihood}

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

    def print_args(self, prefix_name: str = '', prefix: str = '', sort: bool = False) -> None:
        if all((prefix_name, prefix)):
            print(f'{prefix_name}\t\t>\t>\t>\t\t{prefix}')
        items = dict(sorted(self.__dict__().items())).items() if sort else self.__dict__().items()
        for key, value in items:
            print(f'{key}:\t{value}')

    def set_tree_data(self, msa: Optional[Union[Dict[str, str], str]] = None,
                      categories_quantity: Optional[int] = None,
                      alpha: Optional[float] = None,
                      beta: Optional[float] = None,
                      pi_0: Optional[Union[float, np.ndarray, int]] = None,
                      pi_1: Optional[Union[float, np.ndarray, int]] = None,
                      coefficient_bl: Optional[Union[float, np.ndarray, int]] = None,
                      is_optimize_pi: Optional[bool] = None,
                      is_optimize_pi_average: Optional[bool] = None,
                      is_optimize_alpha: Optional[bool] = None,
                      is_optimize_bl: Optional[bool] = None) -> None:
        if isinstance(msa, str):
            self.msa = self.get_msa_dict(msa)
        elif isinstance(msa, dict):
            self.msa = msa
        if isinstance(self.msa, dict) and self.msa:
            self.alphabet = Tree.get_alphabet_from_dict(self.msa)

        self.set_all(categories_quantity, alpha, beta, pi_0, pi_1, coefficient_bl)

        self.optimize_coefficient_bl(is_optimize_bl)
        self.optimize_pi(is_optimize_pi, is_optimize_pi_average)
        self.optimize_alpha(is_optimize_alpha)

        if (is_optimize_alpha or is_optimize_pi or is_optimize_pi_average) and is_optimize_bl:
            self.optimize_coefficient_bl(is_optimize_bl)

        self.set_distance_taking_into_coefficient()

    def set_distance_taking_into_coefficient(self):
        for current_node in self.root.get_list_nodes_info(only_node_list=True):
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
            mode (str, optional): `None` (default), 'pre-order', 'in-order', 'post-order', 'level-order'.
            filters (Dict, optional): `None` (default)
        Returns:
            None: This function does not return any value; it only prints the nodes to the standard output.
        """
        data_structure = self.root.get_list_nodes_info(with_additional_details, mode, filters)

        str_result = ''
        for i in data_structure:
            str_result = f'{str_result}\n{i}'
        print(str_result, '\n')

    def get_tree_info(self, filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None
                      ) -> pd.Series:
        nodes_info = self.get_list_nodes_info(True, 'pre-order', filters)

        return pd.Series([pd.Series(i) for i in nodes_info], index=[i.get('node') for i in nodes_info])

    def get_list_nodes_info(self, with_additional_details: bool = False, mode: Optional[str] = None, filters:
                            Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None, only_node_list:
                            bool = False) -> List[Union[Dict[str, Union[float, np.ndarray, bool, str, List[float],
                                                  List[np.ndarray]]], 'Node']]:
        """
        Args:
            with_additional_details (bool, optional): `False` (default).
            mode (Optional[str]): `None` (default), 'pre-order', 'in-order', 'post-order', 'level-order'.
            filters (Dict, optional): `None` (default).
            only_node_list (Dict, optional): `False` (default).
        """
        return self.root.get_list_nodes_info(with_additional_details, mode, filters, only_node_list)

    def get_node_count(self, filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None) -> int:
        """
        Args:
            filters (Dict, optional):
        """
        return len(self.get_list_nodes_info(True, None, filters))

    def get_node_by_name(self, name: str) -> Optional[Node]:

        return self.root.get_node_by_name(name)

    def get_newick(self, with_internal_nodes: bool = False, decimal_length: int = 0,
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

        return name in self.root.get_list_nodes_info()

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
                if current_node.node_type in ('node', ) and self.is_bootstrap_value(current_node.name):
                    current_node.name = 'nd' + str(num()).rjust(4, '0')

            return self

    def get_html_tree(self, style: str = '', status: str = '') -> str:
        """This method is for internal use only."""
        return self.structure_to_html_tree(self.tree_to_structure(), style, status)

    def tree_to_structure(self) -> Dict[str, str]:
        """This method is for internal use only."""
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
        """This method is for internal use only."""
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
        nodes_info = self.get_list_nodes_info(True, None, filters)

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
                                     'probabilities_sequence_characters', 'log_likelihood_vector', 'sequence',
                                     'ancestral_sequence', 'probability_vector_gain', 'probability_vector_loss')

        for node_info in nodes_info:
            for i in set(node_info.keys()) - set(columns.keys()):
                node_info.pop(i)
            if not node_info.get('father_name'):
                node_info.update({'father_name': 'root'})
            if columns.get(distance_name):
                distance_value = node_info.pop(distance_name)
                if distance_type is str:
                    distance_value = f'{distance_value:.10f}'.ljust(decimal_length, "0"
                                                                    ) if distance_value else ' ' * decimal_length
                else:
                    distance_value = distance_type(distance_value)
                node_info.update({distance_name: distance_value})
            for i in lists:
                if columns.get(i):
                    node_info.update({i: self.get_list_decimals(node_info.get(i), list_type, decimals)})

        tree_table = pd.DataFrame([i for i in nodes_info], index=None)
        tree_table = tree_table.rename(columns=columns)
        tree_table = tree_table.reindex(columns=columns.values())
        # tree_table = tree_table.fillna(0)
        if isinstance(list_type, (list, tuple, set)):
            lists_names = [v for k, v in columns.items() if k in lists]
            sort_values_by = tuple([i for i in sort_values_by if i not in lists_names])

        return tree_table.sort_values(by=list(sort_values_by)) if sort_values_by else tree_table

    @staticmethod
    def get_round(obj: Union[int, float, np.ndarray], decimals: int = 4) -> float:
        return float(np.round(obj, decimals))

    @staticmethod
    def get_list_decimals(obj: Union[int, float, np.ndarray], list_type: type = str, decimals: int = 4) -> Any:
        if list_type in (list, tuple, set):
            if isinstance(obj, (list, tuple, set)):
                return list_type(map(lambda x: Tree.get_round(x, decimals) if (isinstance(x, (int, float, np.ndarray))
                                                                               ) else Tree.get_list_decimals(x,
                                                                                                             list_type,
                                                                                                             decimals),
                                     obj))
            else:
                return obj
        else:
            return ' '.join(map(str, obj))

    def calculate_ancestral_sequence(self, newick_node: Optional[Union[Node, str]] = None) -> str:
        if self.alphabet and not self.calculated_ancestor_sequence:
            node_list = []
            if not newick_node:
                node_list = self.root.get_list_nodes_info(filters={'node_type': ['node', 'leaf']}, only_node_list=True)
            else:
                node_list.append(newick_node)

            ancestral_alphabet = Tree.get_ancestral_alphabet()
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

    def calculate_gl_probability(self) -> None:
        node_list = self.root.get_list_nodes_info(filters={'node_type': ['node', 'leaf']}, only_node_list=True)

        for current_node in node_list:
            current_node.calculate_gl_probability()

    def calculate_marginal(self, newick_node: Optional[Union[Node, str]] = None) -> None:
        if not newick_node:
            node_list = self.root.get_list_nodes_info(filters={'node_type': ['node', 'root']}, only_node_list=True)
        else:
            node_list = []
            if isinstance(newick_node, str):
                node_list.append(self.get_node_by_name(newick_node))
            elif isinstance(newick_node, Node):
                node_list.append(newick_node)

        for current_node in node_list:
            current_node.calculate_marginal()

    def calculate_up(self, msa: str) -> Union[Tuple[Union[List[np.ndarray], List[float]], float], float]:

        return self.root.calculate_up(self.get_msa_dict(msa, self.alphabet))

    def calculate_down(self) -> None:

        self.root.calculate_down(self.get_tree_info())

    def get_msa_dict(self, msa: str, alphabet: Optional[Union[Tuple[str, ...], str]] = None, only_leaves: bool = True
                     ) -> Dict[str, Union[Tuple[int, ...], str]]:
        node_types = ['leaf'] if only_leaves else ['leaf', 'node', 'root']
        nodes_info = self.get_list_nodes_info(True, 'pre-order', {'node_type': node_types})
        msa_list = msa.strip().split()
        msa_list_size, msa_dict = len(msa_list), dict()
        if msa_list_size == 1:
            for i, node_info in enumerate(nodes_info):
                if alphabet:
                    value = [0] * len(alphabet)
                    value[alphabet.index(msa[i])] = 1
                    value = tuple(value)
                else:
                    value = msa[i]
                msa_dict.update({node_info.get('node'): value})
        else:
            for j in range(msa_list_size // 2):
                node_name = msa_list[j + j][1::]
                if self.find_dict_in_iterable(nodes_info, 'node', node_name):
                    value = msa_list[j + j + 1]
                    value = ''.join(value)
                    msa_dict.update({node_name: value})

        return msa_dict

    def calculate_tree(self) -> Dict[str, Union[float, np.ndarray, int]]:
        if self.msa and not self.calculated_tree:
            self.clean_all()

            leaves = self.get_list_nodes_info(filters={'node_type': ['leaf']}, only_node_list=True)
            len_seq = len(min(list(self.msa.values())))
            for i in range(len_seq):
                self.likelihood *= self.calculate_up(msa=''.join([self.msa.get(leaf.name)[i] for leaf in leaves]))
                self.calculate_down()
                self.calculate_marginal()
                self.calculate_gl_probability()
            self.log_likelihood, self.log_likelihood_vector = self.root.log_likelihood, self.root.log_likelihood_vector
            self.calculated_tree = True
            self.calculated_likelihood = True

        return {'likelihood': self.likelihood, 'log_likelihood': self.log_likelihood,
                'log_likelihood_vector': self.log_likelihood_vector}

    def calculate_likelihood(self) -> None:
        if self.msa and self.alphabet and not self.calculated_likelihood:
            self.clean_all()
            self.log_likelihood_vector, self.log_likelihood, self.likelihood = (
                self.root.calculate_likelihood(self.msa))
            self.calculated_likelihood = True

    def get_fasta_text(self) -> str:
        fasta_text = ''
        for k, v in self.msa.items():
            fasta_text += f'>{k}\n{v}\n'

        return fasta_text[:-1]

    def get_json_structure(self, return_table: bool = False, columns: Optional[Dict[str, str]] = None,
                           mode: str = 'node', taking_into_coefficient: bool = True
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

        return loads(dumps(dict_json, cls=npEncode).replace(f'\'', r'"'))

    def tree_to_fasta_file(self, file_name: str = 'file.fasta') -> str:

        fasta_text = self.get_fasta_text()

        return self.write_file(file_name, fasta_text)

    def probability_to_tsv(self, file_name: str = 'ProbabilityPerPositionsPerBranches.tsv', sep: str = '\t',
                           taking_into_coefficient: bool = True) -> str:
        ancestral_comparison = ['absence', 'loss', 'gain', 'presence']
        probability_limit = 0.05
        rows = []

        suffix = '_taking_into_coefficient' if taking_into_coefficient else ''
        distance_to_father = f'distance_to_father{suffix}'
        distance_to_root = f'distance_to_root{suffix}'
        distance_to_nearest = f'distance_to_nearest{suffix}'

        list_nodes = self.get_list_nodes_info(only_node_list=True)
        for current_node in list_nodes:
            branch_probability_vector = current_node.branch_probability_vector
            for pos, value in enumerate(branch_probability_vector, start=1):
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
        df.to_csv(file_name, sep=sep, index=False)

        return file_name

    def attributes_to_tsv(self, file_name: str = 'TreeAttributes.tsv', sep: str = '\t') -> str:

        self.make_dir(file_name)
        data = loads(dumps({'π1 value': self.pi_1,
                            'Γ distribution α value': self.alpha,
                            'number of rate categories': self.categories_quantity,
                            'coefficient of branch lengths': self.coefficient_bl,
                            'rate vector': self.rate_vector,
                            'alphabet': self.alphabet,
                            'log_likelihood': self.log_likelihood}, cls=npEncode))
        df = pd.DataFrame({k: ((v, ) if isinstance(v, (set, tuple, list)) else v) for k, v in data.items()
                           if v is not None})
        df.to_csv(file_name, sep=sep, index=False)

        return file_name

    def likelihood_to_tsv(self, file_name: str = 'LogLikelihood.tsv', sep: str = '\t') -> str:

        self.make_dir(file_name)
        self.calculate_likelihood()
        df = pd.DataFrame({'POS': range(len(self.log_likelihood_vector)),
                           'log-likelihood': self.log_likelihood_vector})
        df.to_csv(file_name, sep=sep, index=False)

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
        table.to_csv(file_name, index=False, sep=sep)

        return file_name

    def tree_to_newick_file(self, file_name: str = 'tree_file.tree', with_internal_nodes: bool = False,
                            decimal_length: int = 0, taking_into_coefficient: bool = True) -> str:

        newick_text = self.get_newick(with_internal_nodes, decimal_length, taking_into_coefficient)

        return self.write_file(file_name, newick_text)

    def tree_to_visual_format(self, file_name: str = 'VisualTree.svg', with_internal_nodes: bool = False,
                              file_extensions: Optional[Union[str, Tuple[str, ...]]] = None, show_axes: bool = False,
                              taking_into_coefficient: bool = True) -> Dict[str, str]:
        file_extensions = Tree.check_file_extensions_tuple(file_extensions, 'svg')

        self.make_dir(file_name)
        tmp_dir = Path(file_name).parent.joinpath('tmp')
        tmp_file = f'{tmp_dir.joinpath(f"{Tree.get_random_name()}.tree")}'
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

        self.calculate_tree()
        self.calculate_ancestral_sequence()
        size_factor = min(1 + self.get_node_count({'node_type': ['leaf']}) // 9, 6)
        columns, lists, decimals = self.get_columns(mode='tree_html', taking_into_coefficient=taking_into_coefficient)
        df = self.tree_to_table(columns=columns, distance_type=float, filters={'node_type': ['leaf', 'node', 'root']},
                                list_type=list, taking_into_coefficient=taking_into_coefficient, lists=lists,
                                decimals=decimals)
        df_copy = df.copy()
        del df['sequence'], df['node_type'], df['prob_characters']
        df = df.iloc[1:]

        d3 = D3Blocks(verbose=100, chart='tree', frame=False)
        d3.set_node_properties(df)

        d3.font = {'size': 12}
        d3.hierarchy = [i for i in range(1, len(df_copy.T.count()) + 1)]
        d3.title = 'Phylogenetic tree'
        d3.filepath = file_name
        d3.figsize = (500, 500)
        d3.showfig, d3.overwrite, d3.reset_properties, d3.save_button = True, True, True, True
        d3.notebook = False
        d3.config = d3.chart.set_config(config=d3.config, filepath=d3.filepath, font=d3.font, title=d3.title,
                                        margin={"top": 20, "right": 40, "bottom": 20, "left": 40},
                                        showfig=d3.showfig, overwrite=d3.overwrite, figsize=d3.figsize,
                                        reset_properties=d3.reset_properties, notebook=d3.notebook,
                                        hierarchy=d3.hierarchy, save_button=d3.save_button)

        colors = ['crimson', 'orangered', 'darkorange', 'gold', 'yellowgreen', 'forestgreen', 'mediumturquoise',
                  'dodgerblue', 'slateblue', 'darkviolet']
        colors_as = {'A': 'crimson', 'L': 'darkorange', 'G': 'forestgreen', 'P': 'slateblue'}
        for i in df_copy.T:
            probability_coefficient = ancestral_sequence = ''
            sequence = ''.join([Node.draw_cell_html_table(colors[Node.get_integer(j)], j)
                                for j in df_copy['sequence'][i]])
            sequence = Node.draw_row_html_table('Sequence', sequence)
            if df_copy["node_type"][i] != 'root':
                ancestral_sequence = ''.join([Node.draw_cell_html_table(colors_as[j], j)
                                              for j in df_copy['ancestral_sequence'][i]])
                ancestral_sequence = Node.draw_row_html_table('Ancestral Comparison', ancestral_sequence)
            if df_copy["node_type"][i] != 'leaf':
                probability_coefficient = ''.join([Node.draw_cell_html_table(colors[Node.get_integer(j)], f'{j:.3f}')
                                                  for j in df_copy['prob_characters'][i]])
                probability_coefficient = Node.draw_row_html_table('Probability coefficient', probability_coefficient)
                if df_copy["node_type"][i] == 'node':
                    d3.node_properties.get(df_copy['target'][i])['color'] = 'darkorange'
                    d3.node_properties.get(df_copy['target'][i])['size'] = 15 / size_factor
                if df_copy["node_type"][i] == 'root':
                    d3.node_properties.get(df_copy['target'][i])['color'] = 'firebrick'
                    d3.node_properties.get(df_copy['target'][i])['size'] = 20 / size_factor
            else:
                d3.node_properties.get(df_copy['target'][i])['color'] = 'forestgreen'
                d3.node_properties.get(df_copy['target'][i])['size'] = 10 / size_factor
            distance = f'<td class="h7 w-auto text-center">{df_copy["weight"][i]}</td>'
            info = (f'{Node.draw_row_html_table("Distance", distance)}{sequence}{probability_coefficient}'
                    f'{ancestral_sequence}')
            d3.node_properties.get(df_copy['target'][i])['tooltip'] = Node.draw_html_table(info)
            d3.font.update({'type': 'Anonymous Pro'})

        d3.set_edge_properties(df)
        d3.show()

        return file_name

    def tree_to_graph(self, file_name: str = 'graph.svg', file_extensions: Optional[Union[str, Tuple[str, ...]]] = None
                      ) -> Union[str, Dict[str, str]]:
        file_extensions = Tree.check_file_extensions_tuple(file_extensions, 'png')

        size_factor = min(1 + self.get_node_count({'node_type': ['leaf']}) // 9, 6)
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

    def optimize(self, func: Union[Callable, str], bracket: Tuple[Union[float, np.ndarray], ...] = (0.5,),
                 bounds: Tuple[Union[float, np.ndarray], ...] = (0.001, 0.999), args: Optional[Tuple[Any, ...]] = None,
                 result_fild: Optional[str] = None):
        """
            result_fild: `str` (default), message, success, status, fun, x, nit, nfev
        """
        func = self.__getattribute__(func) if isinstance(func, str) else func
        min_scalar = minimize_scalar(func, bracket=bracket, bounds=bounds) if args is None else (
            minimize_scalar(func, args=args, bracket=bracket, bounds=bounds))

        return min_scalar[result_fild] if result_fild else min_scalar

    def pi_optimization(self, pi: Union[float, np.ndarray], mode: int = 1) -> Union[float, np.ndarray]:
        current_pi = (pi, None)
        self.clean_all()
        self.set_pi(current_pi[mode], current_pi[::-1][mode])
        self.set_vars()

        return -self.root.calculate_likelihood(self.msa)[1]

    def alpha_optimization(self, alpha: Union[int, float, np.ndarray]) -> Union[float, np.ndarray]:
        self.clean_all()
        self.set_gamma_distribution_categories_vector(alpha)
        self.set_vars()

        return -self.root.calculate_likelihood(self.msa)[1]

    def coefficient_bl_optimization(self, coefficient_bl: Union[int, float, np.ndarray]) -> Union[float, np.ndarray]:
        self.clean_all()
        self.set_coefficient_bl(coefficient_bl)
        self.set_vars()

        return -self.root.calculate_likelihood(self.msa)[1]

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
        elif is_optimize_pi_average:
            all_lines_list = list(self.msa.values())
            all_lines = ''.join(all_lines_list)
            self.pi_1 = all_lines.count(self.alphabet[mode]) / len(all_lines)

        self.set_vars()

    def clean_all(self):
        self.root.clean_all()
        self.likelihood, self.log_likelihood, self.log_likelihood_vector = 1, 0, []

    def set_all(self, categories_quantity: Optional[int] = None, alpha: Optional[float] = None,
                beta: Optional[float] = None, pi_0: Optional[Union[float, np.ndarray, int]] = None,
                pi_1: Optional[Union[float, np.ndarray, int]] = None,
                coefficient_bl: Optional[Union[float, np.ndarray, int]] = None) -> None:

        self.categories_quantity = categories_quantity
        self.set_alpha(alpha, beta)
        self.set_pi(pi_0, pi_1)
        self.set_coefficient_bl(coefficient_bl)
        self.set_gamma_distribution_categories_vector(self.alpha)
        self.set_vars()

    def set_gamma_distribution_categories_vector(self, alpha: Union[int, float, np.ndarray]) -> None:
        self.set_alpha(alpha)
        categories_vector = []
        gamma_percent_point = self.get_gamma_distribution_percent_point()
        for i in range(self.categories_quantity):
            lower_gamma_inc_1 = gammainc(self.alpha + 1, gamma_percent_point[i] * self.alpha)
            lower_gamma_inc_2 = gammainc(self.alpha + 1, gamma_percent_point[i + 1] * self.alpha)
            mean = (self.alpha / self.alpha) * (lower_gamma_inc_2 - lower_gamma_inc_1) / (1 / self.categories_quantity)
            categories_vector.append(mean)

        self.rate_vector = tuple(categories_vector)

    def set_coefficient_bl(self, coefficient_bl: Optional[Union[float, np.ndarray, int]] = None) -> None:
        self.coefficient_bl = 1.0 if coefficient_bl is None else coefficient_bl

    def set_alpha(self, alpha: Optional[float] = None, beta: Optional[float] = None) -> None:
        self.alpha = alpha if alpha else (beta if beta else 0.5)

    def set_pi(self, pi_0: Optional[Union[float, np.ndarray, int]] = None,
               pi_1: Optional[Union[float, np.ndarray, int]] = None) -> None:
        alphabet_size = len(self.alphabet)
        if pi_0:
            self.pi_1 = 1 - pi_0
        elif pi_1:
            self.pi_1 = pi_1
        else:
            self.pi_1 = 1 / alphabet_size

    def set_vars(self) -> None:
        alphabet_size = len(self.alphabet)
        rate_vector_size = len(self.rate_vector)
        if self.pi_1:
            frequency = (1 - self.pi_1, self.pi_1)
        else:
            frequency = (1 / alphabet_size, 1 / alphabet_size)
        for current_node in self.root.get_list_nodes_info(only_node_list=True):
            current_node.alphabet = self.alphabet
            current_node.alphabet_size = alphabet_size
            current_node.rate_vector_size = rate_vector_size
            current_node.frequency = frequency
            current_node.pi_1 = self.pi_1
            current_node.coefficient_bl = self.coefficient_bl
            current_node.pmatrix = tuple([current_node.get_pmatrix(r) for r in self.rate_vector])

    def get_gamma_distribution_percent_point(self) -> List[float]:
        probability_vector = np.linspace(0, 1, self.categories_quantity + 1)

        return gamma.ppf(probability_vector, a=self.alpha, scale=1/self.alpha)

    @staticmethod
    def is_bootstrap_value(number_str: str, lower: Union[float, np.ndarray, int] = 0,
                           upper: Union[float, np.ndarray, int] = 100) -> bool:
        re_result = bool(re.fullmatch(r'^-?\d+(\.\d+)?$', number_str)) if len(number_str.strip()) else False

        return lower <= float(number_str) <= upper if re_result else False

    @staticmethod
    def get_columns(mode: str = 'node', columns: Optional[Dict[str, str]] = None,
                    taking_into_coefficient: bool = True) -> Tuple[Dict[str, str], Tuple[str, ...], int]:

        suffix = '_taking_into_coefficient' if taking_into_coefficient else ''
        distance_name = f'distance{suffix}'

        lists = ('children', 'full_distance', 'full_distance_taking_into_coefficient', 'up_vector', 'down_vector',
                 'marginal_vector', 'marginal_bl_vector', 'probability_vector', 'probabilities_sequence_characters',
                 'log_likelihood_vector', 'probability_vector_gain', 'probability_vector_loss', 'ancestral_sequence',
                 'sequence')
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
                                               'probabilities_sequence_characters': 'Probability coefficient',
                                               'ancestral_sequence': 'Ancestral comparison', 'sequence_likelihood':
                                               'Likelihood of sequence', 'log_likelihood': 'Log-likelihood',
                                               'log_likelihood_vector': 'Vector of log-likelihood'}
            lists = ('children', 'probabilities_sequence_characters')
            decimals = 8
        elif mode == 'branch_tsv':
            columns = columns if columns else {'father_name': 'Parent node', 'node': 'Child node', distance_name:
                                               'Branch length', 'branch_probability_vector':
                                               'Branch probability vector'}
            lists = ('branch_probability_vector', )
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
    def get_alphabet_from_dict(msa_dict: Dict[str, str]) -> Tuple[str, ...]:
        character_list = []
        for sequence in msa_dict.values():
            character_list += [i for i in sequence]

        return Tree.get_alphabet(set(character_list))

    @staticmethod
    def get_columns_list_for_sorting(mode: str = 'node') -> Dict[str, List[str]]:
        if mode == 'node':
            result = {'List for sorting': ['Name', 'Node type', 'Distance to parent', 'Sequence',
                                           'Probability coefficient', 'Ancestral Comparison']}
        else:
            result = {'List for sorting': ['Parent node', 'Child node', 'Branch length', 'Gain probability',
                                           'Loss probability']}

        return loads(dumps(result, cls=npEncode).replace(f'\'', r'"'))

    @staticmethod
    def get_random_name(lenght: int = 24) -> str:
        abc_list = [_ for _ in 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz1234567890']

        return ''.join(np.random.choice(abc_list, lenght))

    @staticmethod
    def get_ancestral_alphabet() -> Tuple[str, ...]:

        return 'A', 'L', 'G', 'P'

    @staticmethod
    def get_alphabet(search_argument: Union[Set[str], int, str]) -> Tuple[str, ...]:
        alphabets = (('0', '1'), ('A', 'C', 'G', 'T'),
                     ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                      'V'))
        if isinstance(search_argument, int):
            return tuple(alphabets[search_argument])
        if isinstance(search_argument, str):
            search_argument = set([i for i in search_argument])
        if isinstance(search_argument, set):
            for alphabet in alphabets:
                if not search_argument - set(alphabet):
                    return alphabet

    @staticmethod
    def find_dict_in_iterable(iterable: Union[List[Union[Dict[str, Union[float, np.ndarray, bool, str, List[float],
                              List[np.ndarray]]], 'Node']], Tuple[Dict[str, Union[float, bool, str, List[float], Tuple[
                                   int, ...]]]]], key: str, value: Optional[Union[float, bool, str, List[float]]] = None
                              ) -> Dict[str, Union[float, bool, str, List[float], List[int], Tuple[int, ...]]]:
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
    def check_tree(newick_tree: Union[str, 'Tree']) -> 'Tree':
        if isinstance(newick_tree, str):
            newick_tree = Tree(newick_tree)

        return newick_tree

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
    def __set_node(node_str: str, num) -> Node:
        """This method is for internal use only."""
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
    def rename_nodes(newick_tree: Union[str, 'Tree'], node_name: str = 'N', fill_character: str = '0', number_length:
                     int = 0) -> 'Tree':
        newick_tree = Tree.check_tree(newick_tree)
        nodes_list = newick_tree.get_list_nodes_info(filters={'node_type': ['root', 'node']}, only_node_list=True,
                                                     mode='pre-order')
        num = newick_tree.__counter()
        for current_node in nodes_list:
            if re.fullmatch(r'^nd\d{4}$', current_node.name):
                current_node.name = f'{node_name}{str(num()).rjust(number_length, fill_character)}'

        return newick_tree

    @staticmethod
    def __counter():
        """This method is for internal use only."""
        count = 0

        def sub_function():
            nonlocal count
            count += 1
            return count

        return sub_function

    @classmethod
    def __get_html_tree(cls, structure: dict, status: str) -> str:
        """This method is for internal use only."""
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
        """This method is for internal use only."""
        tree1 = Tree(tree1) if type(tree1) is str else tree1
        tree2 = Tree(tree2) if type(tree2) is str else tree2

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
        """This method is for internal use only."""
        return (f'<ul {f" class = {chr(34)}{styleclass}{chr(34)}" if styleclass else ""}>'
                f'{cls.__get_html_tree(structure, status)}</ul>')

    @classmethod
    def subtree_to_structure(cls, newick_node: Node) -> Dict[str, str]:
        """This method is for internal use only."""
        dict_node = {'name': newick_node.name.strip(), 'distance_to_father': newick_node.distance_to_father}
        list_children = []
        if newick_node.children:
            for child in newick_node.children:
                list_children.append(cls.subtree_to_structure(child))
        dict_node.update({'children': list_children})

        return dict_node
