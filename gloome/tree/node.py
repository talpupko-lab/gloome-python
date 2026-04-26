import pandas as pd
import numpy as np

from typing import Optional, Dict, Union, List, Tuple, Any
from scipy.linalg import expm
from math import log, prod
from json import loads, dumps

from .npencoder import NpEncoder

eps = 5e-324


class Node:
    father: Optional['Node']
    children: List['Node']
    name: str
    node_type: str
    distance_to_father: Union[float, np.ndarray]
    distance_to_root: Union[float, np.ndarray]
    distance_to_root_vector: List[Union[float, np.ndarray]]
    distance_to_nearest: Union[float, np.ndarray]
    distance_to_father_taking_into_coefficient: Union[float, np.ndarray]
    distance_to_root_taking_into_coefficient: Union[float, np.ndarray]
    distance_to_root_vector_taking_into_coefficient: List[Union[float, np.ndarray]]
    distance_to_nearest_taking_into_coefficient: Union[float, np.ndarray]
    level: int
    levels_to_nearest: int
    alphabet: Tuple[str, ...]
    alphabet_size: int
    rate_vector_size: int
    pi_1: Union[float, np.ndarray]
    frequency: Tuple[Union[float, np.ndarray], ...]
    coefficient_bl: Union[float, np.ndarray, int]
    pmatrix: Optional[Tuple[np.ndarray, ...]]
    log_likelihood_vector: List[Union[float, np.ndarray]]
    log_likelihood: Union[float, np.ndarray]
    sequence_likelihood: Union[float, np.ndarray]
    likelihood: Union[float, np.ndarray]
    up_vector: List[List[Union[float, np.ndarray]]]
    down_vector: List[List[Union[float, np.ndarray]]]
    marginal_vector: List[List[Union[float, np.ndarray]]]
    marginal_bl_vector: List[List[Union[float, np.ndarray]]]
    probability_vector: List[List[Union[float, np.ndarray]]]
    branch_probability_vector: List[List[Union[float, np.ndarray]]]
    probability_vector_gain: List[Union[float, np.ndarray]]
    probability_vector_loss: List[Union[float, np.ndarray]]
    sequence: str
    probabilities_sequence_characters: List[Union[float, np.ndarray]]
    ancestral_sequence: str

    def __init__(self, name: Optional[str]) -> None:
        self.father = None
        self.children = []
        self.name = name
        self.node_type = ''
        self.distance_to_father = 0.0
        self.distance_to_root = 0.0
        self.distance_to_root_vector = []
        self.distance_to_nearest = 0.0
        self.distance_to_father_taking_into_coefficient = 0.0
        self.distance_to_root_taking_into_coefficient = 0.0
        self.distance_to_root_vector_taking_into_coefficient = []
        self.distance_to_nearest_taking_into_coefficient = 0.0
        self.level = 0
        self.levels_to_nearest = 0
        self.alphabet = ('0', '1')
        self.alphabet_size = 2
        self.rate_vector_size = 1
        self.pi_1 = 0.5
        self.frequency = (0.5, 0.5)
        self.coefficient_bl = 1.0
        self.pmatrix = None
        self.log_likelihood_vector = []
        self.log_likelihood = 0.0
        self.sequence_likelihood = 1.0
        self.likelihood = 0.0
        self.up_vector = []
        self.down_vector = []
        self.marginal_vector = []
        self.marginal_bl_vector = []
        self.probability_vector = []
        self.branch_probability_vector = []
        self.probability_vector_gain = []
        self.probability_vector_loss = []
        self.sequence = ''
        self.probabilities_sequence_characters = []
        self.ancestral_sequence = ''

    def __str__(self) -> str:
        return self.get_name(True)

    def __dir__(self) -> list:
        return ['father', 'children', 'name', 'distance_to_father', 'distance_to_root', 'distance_to_root_vector',
                'distance_to_nearest', 'distance_to_father_taking_into_coefficient',
                'distance_to_root_taking_into_coefficient', 'distance_to_root_vector_taking_into_coefficient',
                'distance_to_nearest_taking_into_coefficient', 'level', 'levels_to_nearest', 'alphabet',
                'alphabet_size', 'rate_vector_size', 'pi_1', 'frequency', 'coefficient_bl', 'pmatrix',
                'log_likelihood_vector', 'log_likelihood', 'sequence_likelihood', 'likelihood', 'up_vector',
                'down_vector', 'marginal_vector', 'marginal_bl_vector' 'probability_vector',
                'branch_probability_vector', 'probability_vector_gain', 'probability_vector_loss', 'sequence',
                'probabilities_sequence_characters', 'ancestral_sequence']

    def get_list_nodes_info(self, with_additional_details: bool = False,
                            mode: Optional[str] = None,
                            filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None,
                            only_node_list: bool = False
                            ) -> List[Union[Dict[str, Union[float, np.ndarray, bool, str, List[float],
                                      List[np.ndarray]]], 'Node']]:
        """
        Retrieve a list of descendant nodes from a given node, including the node itself.

        This function collects all child nodes of the specified `node`, including the `node` itself. The function
        returns a list of nodes or a list of dictionaries with information about these nodes.

        Args:
            with_additional_details (bool, optional): `False` (default).
            mode (str, optional): None (default), 'pre-order', 'in-order', 'post-order', 'level-order'.
            filters (Dict, optional):
            only_node_list (Dict, optional): `False` (default).
        Returns:
            list: A list of descendant nodes from a given node, including the node itself or a list of dictionaries
            with information about these nodes.
        """
        list_result = []
        mode = 'pre-order' if mode is None or mode.lower() not in ('pre-order', 'in-order', 'post-order', 'level-order'
                                                                   ) else mode.lower()
        condition = with_additional_details or only_node_list

        def get_list(trees_node: Node) -> None:
            nonlocal list_result, filters, mode, condition

            nodes_info = trees_node.get_node_info()
            list_item = trees_node if only_node_list else nodes_info
            if trees_node.check_filter_compliance(filters, nodes_info):
                if mode == 'pre-order':
                    list_result.append(list_item if condition else trees_node.name)

                for i, child in enumerate(trees_node.children):
                    get_list(child)
                    if mode == 'in-order' and not i:
                        list_result.append(list_item if condition else trees_node.name)

                if not trees_node.children:
                    if mode == 'in-order':
                        list_result.append(list_item if condition else trees_node.name)

                if mode == 'post-order':
                    list_result.append(list_item if condition else trees_node.name)
            else:
                for child in trees_node.children:
                    get_list(child)

        if mode == 'level-order':
            nodes_list = [self]
            while nodes_list:
                newick_node = nodes_list.pop(0)
                if newick_node.check_filter_compliance(filters, newick_node.get_node_info()):
                    level_order_item = newick_node if only_node_list else newick_node.get_node_info()
                    list_result.append(level_order_item if condition else newick_node.name)

                for nodes_child in newick_node.children:
                    nodes_list.append(nodes_child)
        else:
            get_list(self)

        return list_result

    def get_node_info(self) -> Dict[str, Union[float, np.ndarray, bool, str, List[float], List[np.ndarray]]]:

        result = {'node': self.name,
                  'distance': self.distance_to_father,
                  'distance_taking_into_coefficient': self.distance_to_father_taking_into_coefficient,
                  'distance_to_root': self.distance_to_root,
                  'distance_to_root_taking_into_coefficient': self.distance_to_root_taking_into_coefficient,
                  'distance_to_nearest': self.distance_to_nearest,
                  'distance_to_nearest_taking_into_coefficient': self.distance_to_nearest_taking_into_coefficient,
                  'level': self.level,
                  'levels_to_nearest': self.levels_to_nearest,
                  'node_type': self.node_type,
                  'father_name': self.father.name if self.father else '',
                  'full_distance': self.distance_to_root_vector,
                  'full_distance_taking_into_coefficient': self.distance_to_root_vector_taking_into_coefficient,
                  'children': [i.name for i in self.children],
                  'up_vector': self.up_vector,
                  'down_vector': self.down_vector,
                  'likelihood': self.likelihood,
                  'sequence_likelihood': self.sequence_likelihood,
                  'log_likelihood': self.log_likelihood,
                  'log_likelihood_vector': self.log_likelihood_vector,
                  'marginal_vector': self.marginal_vector,
                  'marginal_bl_vector': self.marginal_bl_vector,
                  'probability_vector': self.probability_vector,
                  'sequence': self.sequence,
                  'probabilities_sequence_characters': self.probabilities_sequence_characters,
                  'ancestral_sequence': self.ancestral_sequence,
                  'branch_probability_vector': self.branch_probability_vector,
                  'probability_vector_gain': self.probability_vector_gain,
                  'probability_vector_loss': self.probability_vector_loss,
                  'alphabet': self.alphabet,
                  'alphabet_size': self.alphabet_size,
                  'rate_vector_size': self.rate_vector_size,
                  'pi_1': self.pi_1,
                  'frequency': self.frequency,
                  'coefficient_bl': self.coefficient_bl,
                  'pmatrix': self.pmatrix}

        return loads(dumps(result, cls=NpEncoder))

    def get_node_by_name(self, node_name: str) -> Optional['Node']:
        if node_name == self.name:
            return self
        else:
            for child in self.children:
                newick_node = child.get_node_by_name(node_name)
                if newick_node:
                    return newick_node
        return None

    def get_pmatrix(self, rate: Union[float, np.ndarray] = 1.0):
        return self.get_one_parameter_pmatrix(rate)

    def calculate_sequence_likelihood(self) -> None:
        self.sequence_likelihood *= self.likelihood
        self.log_likelihood += log(max(self.likelihood, eps))
        self.log_likelihood_vector.append(log(max(self.likelihood, eps)))

    def calculate_gl_probability(self) -> None:
        self.marginal_bl_vector = []

        for r in range(self.rate_vector_size):
            current_marginal_bl_vector = []
            for j in range(self.alphabet_size):
                for i in range(self.alphabet_size):
                    current_marginal_bl_vector.append(self.frequency[i] * self.up_vector[r][j] *
                                                      self.pmatrix[r][i, j] * self.down_vector[r][i])
            self.marginal_bl_vector.append(current_marginal_bl_vector)

        likelihood = (np.sum([np.sum(self.marginal_bl_vector[r]) for r in range(self.rate_vector_size)]) /
                      self.rate_vector_size)
        likelihood = max(likelihood, eps)

        branch_probability_vector = []
        for i in range(self.alphabet_size * self.alphabet_size):
            branch_probability_vector.append(np.sum([self.marginal_bl_vector[r][i] for r in
                                             range(self.rate_vector_size)]) / self.rate_vector_size / likelihood)
        self.branch_probability_vector.append(branch_probability_vector)
        self.probability_vector_loss.append(branch_probability_vector[1])
        self.probability_vector_gain.append(branch_probability_vector[2])

    def calculate_marginal(self) -> Tuple[Union[Union[List[List[np.ndarray]], List[List[float]]], float],
                                          Union[np.ndarray, float]]:
        self.marginal_vector = []

        for r in range(self.rate_vector_size):
            current_marginal_vector = []
            for j in range(self.alphabet_size):
                marg = 0
                for i in range(self.alphabet_size):
                    marg += self.frequency[i] * self.pmatrix[r][i, j] * self.down_vector[r][i]
                current_marginal_vector.append(self.up_vector[r][j] * marg)
            self.marginal_vector.append(current_marginal_vector)

        likelihood = (np.sum([np.sum(self.marginal_vector[r]) for r in range(self.rate_vector_size)]) /
                      self.rate_vector_size)
        likelihood = max(likelihood, eps)

        probability_vector = []
        for i in range(self.alphabet_size):
            probability_vector.append(np.sum([self.marginal_vector[r][i] for r in range(self.rate_vector_size)]) /
                                      self.rate_vector_size / likelihood)
        self.probability_vector.append(probability_vector)
        probability = max(self.probability_vector[-1])
        self.sequence = f'{self.sequence}{self.alphabet[self.probability_vector[-1].index(probability)]}'
        self.probabilities_sequence_characters.append(probability)

        return self.marginal_vector, likelihood

    def calculate_up(self, nodes_dict: Dict[str, Tuple[int, ...]]
                     ) -> Union[Union[List[List[np.ndarray]], List[List[float]]], float]:
        self.up_vector = []
        self.likelihood = 0

        if not self.children:
            up_vector = list(nodes_dict.get(self.name))
            max_up_vector = max(up_vector)
            self.likelihood = (np.sum([self.frequency[s] * up_vector[s] for s in range(self.alphabet_size)]) /
                               self.rate_vector_size)
            probable_character = self.alphabet[up_vector.index(max_up_vector)]
            self.sequence = f'{self.sequence}{probable_character}'
            self.probabilities_sequence_characters.append(max_up_vector)
            self.up_vector = [up_vector for _ in range(self.rate_vector_size)]

            self.calculate_sequence_likelihood()

            return self.up_vector

        for child in self.children:
            child.calculate_up(nodes_dict)

        for r in range(self.rate_vector_size):
            current_up_vector = []
            for j in range(self.alphabet_size):
                probabilities = {}
                for i in range(self.alphabet_size):
                    for child in self.children:
                        p1 = child.pmatrix[r][j, i] * child.up_vector[r][i]
                        probabilities.update({child.name: probabilities.get(child.name, 0.0) + p1})

                current_up_vector.append(prod(probabilities.values()))
            self.up_vector.append(current_up_vector)
            self.likelihood += np.sum([self.frequency[i] * 1 / self.rate_vector_size * v for i, v in
                                       enumerate(current_up_vector)])

        self.calculate_sequence_likelihood()

        if self.father:
            return self.up_vector
        else:
            return self.likelihood

    def calculate_down(self, tree_info: pd.Series) -> None:
        self.down_vector = []

        father = self.father
        if father:
            brothers = tuple([father.get_node_by_name(i) for i in tree_info.get(father.name).get('children') if i !=
                              self.name])

            for r in range(self.rate_vector_size):
                current_down_vector = []
                for j in range(self.alphabet_size):
                    probabilities = {}
                    for i in range(self.alphabet_size):
                        for brother in brothers:
                            probabilities.update(
                                {brother.name:
                                 probabilities.get(brother.name, 0) + (brother.pmatrix[r][j, i] *
                                                                       brother.up_vector[r][i])})
                        if father.father:
                            probabilities.update(
                                {father.name: probabilities.get(father.name, 0) + (father.pmatrix[r][j, i] *
                                                                                   father.down_vector[r][i])})

                    current_down_vector.append(prod(probabilities.values()))
                self.down_vector.append(current_down_vector)

            for child in self.children:
                child.calculate_down(tree_info)
        else:
            self.down_vector = [[1] * self.alphabet_size for _ in range(self.rate_vector_size)]
            for child in self.children:
                child.calculate_down(tree_info)

    def clean_all(self):
        for current_node in self.get_list_nodes_info(only_node_list=True):
            current_node.log_likelihood_vector = []
            current_node.log_likelihood = 0.0
            current_node.sequence_likelihood = 1.0
            current_node.likelihood = 0.0
            current_node.up_vector = []
            current_node.down_vector = []
            current_node.marginal_vector = []
            current_node.marginal_bl_vector = []
            current_node.probability_vector = []
            current_node.branch_probability_vector = []
            current_node.probability_vector_gain = []
            current_node.probability_vector_loss = []
            current_node.sequence = ''
            current_node.probabilities_sequence_characters = []
            current_node.ancestral_sequence = ''

    def calculate_likelihood(self, msa_dict: Dict[str, str]) -> Tuple[List[float], float, float]:

        leaves_info = self.get_list_nodes_info(True, 'pre-order', {'node_type': ['leaf']})

        len_seq = len(min(list(msa_dict.values())))
        likelihood, log_likelihood, log_likelihood_list = 1, 0, []
        for i_char in range(len_seq):
            nodes_dict = dict()
            for i in range(len(leaves_info)):
                node_name = leaves_info[i].get('node')
                character = msa_dict.get(node_name)[i_char]
                nodes_dict.update({node_name: tuple([int(j == character) for j in self.alphabet])})

            char_likelihood = self.calculate_up(nodes_dict)
            likelihood *= char_likelihood
            log_likelihood += log(max(char_likelihood, eps))
            log_likelihood_list.append(log(max(char_likelihood, eps)))

        return log_likelihood_list, log_likelihood, likelihood

    def get_one_parameter_pmatrix(self, rate: Union[float, np.ndarray] = 1) -> np.ndarray:
        qmatrix = np.zeros((2, 2), dtype='float32')
        qmatrix[0, 0] = - 1 / (2 * (1 - self.pi_1))
        qmatrix[0, 1] = 1 / (2 * (1 - self.pi_1))
        qmatrix[1, 0] = 1 / (2 * self.pi_1)
        qmatrix[1, 1] = - 1 / (2 * self.pi_1)

        return expm(qmatrix * (self.distance_to_father * self.coefficient_bl * rate))

    def get_jukes_cantor_pmatrix(self, rate: Union[float, np.ndarray] = 1) -> np.ndarray:
        qmatrix = np.ones((self.alphabet_size, self.alphabet_size))
        np.fill_diagonal(qmatrix, 1 - self.alphabet_size)
        qmatrix = qmatrix * 1 / (self.alphabet_size - 1)

        return expm(qmatrix * (self.distance_to_father * self.coefficient_bl * rate))

    def node_to_json(self) -> Dict[str, Union[str, List[Any], float, np.ndarray]]:
        dict_json = dict()
        dict_json.update({'name': self.name})
        dict_json.update({'distance': f'{float(self.distance_to_father)}'})

        if self.children:
            dict_json.update({'children': []})
            for child in self.children:
                dict_json['children'].append(child.node_to_json())

        return dict_json

    def get_distance_to_father(self, taking_into_coefficient: bool) -> Union[float, np.ndarray]:
        return self.distance_to_father * self.coefficient_bl if taking_into_coefficient else self.distance_to_father

    def subtree_to_newick(self, with_internal_nodes: bool = False,
                          decimal_length: int = 0,
                          taking_into_coefficient: bool = False) -> str:
        node_list = self.children
        if node_list:
            result = '('
            for child in node_list:
                distance = child.get_distance_to_father(taking_into_coefficient)
                if child.children:
                    child_name = child.subtree_to_newick(with_internal_nodes, decimal_length,
                                                         taking_into_coefficient)
                else:
                    child_name = child.name
                result += f'{child_name}:' + f'{distance:.10f}'.ljust(decimal_length, '0') + ','
            result = f'{result[:-1]}){self.name if with_internal_nodes else ""}'
        else:
            distance = self.get_distance_to_father(taking_into_coefficient)
            result = f'{self.name}:' + f'{distance}'.ljust(decimal_length, '0')
        return result

    def get_name(self, is_full_name: bool = False) -> str:
        return (f'{self.subtree_to_newick() if self.children and is_full_name else self.name}:'
                f'{self.distance_to_father:.6f}')

    def add_child(self, child: 'Node', distance_to_father: Optional[float] = None) -> None:
        self.children.append(child)
        child.father = self
        if distance_to_father is not None:
            child.distance_to_father = distance_to_father

    def get_full_distance_to_father(self, return_list: bool = False) -> Union[List[float], float]:
        list_result = []
        father = self
        while father:
            list_result.append({'node': father, 'distance': father.distance_to_father})
            father = father.father
        result = [i['distance'] for i in list_result]
        return result if return_list else sum(result)

    def set_levels_and_distance_to_nearest(self) -> None:
        nodes_info_list = self.get_list_nodes_info(filters={'node_type': ['leaf']}, only_node_list=True)
        levels_list = []
        distance_list = []
        for newick_node in nodes_info_list:
            levels_list.append(round(newick_node.level - self.level))
            distance_list.append(round(newick_node.distance_to_root - self.distance_to_root, 14))
        self.levels_to_nearest = min(levels_list)
        self.distance_to_nearest = min(distance_list)

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
    def check_filter_compliance(filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]], info: Dict[str,
                                Union[float, bool, str, list[float]]]) -> bool:
        permission = 0
        if filters:
            for key in filters.keys():
                for value in filters.get(key):
                    permission += sum(k == key and info[k] == value for k in info)
        else:
            permission = 1

        return bool(permission)
