import numpy as np

from typing import Optional, Dict, Union, List, Tuple, Any, Set
from scipy.linalg import expm
from json import loads, dumps

from .npencoder import NpEncoder

eps = 5e-324


class Node:
    father: Optional['Node']
    children: List['Node']
    name: str
    node_type: str
    distance_to_father: Union[float, np.float64]
    distance_to_root: Union[float, np.float64]
    distance_to_root_vector: List[Union[float, np.float64]]
    distance_to_nearest: Union[float, np.float64]
    distance_to_father_taking_into_coefficient: Union[float, np.float64]
    distance_to_root_taking_into_coefficient: Union[float, np.float64]
    distance_to_root_vector_taking_into_coefficient: List[Union[float, np.float64]]
    distance_to_nearest_taking_into_coefficient: Union[float, np.float64]
    level: int
    levels_to_nearest: int
    alphabet: Tuple[str, ...]
    pi_1: Union[float, np.float64]
    frequency: Optional[np.ndarray]
    coefficient_bl: Union[float, np.float64, int]
    pmatrix: Optional[np.ndarray]
    log_likelihood_vector: Optional[np.ndarray]
    log_likelihood: Union[float, np.float64, None]
    likelihood_vector: Optional[np.ndarray]
    likelihood: Union[float, np.float64, None]
    sequence_likelihood: Union[float, np.float64]
    up_vector: Optional[np.ndarray]
    down_vector: Optional[np.ndarray]
    marginal_vector: Optional[np.ndarray]
    marginal_bl_vector: Optional[np.ndarray]
    probability_vector: Optional[np.ndarray]
    branch_probability_vector: Optional[np.ndarray]
    probability_vector_gain: Optional[np.ndarray]
    probability_vector_loss: Optional[np.ndarray]
    sequence: str
    probabilities_sequence_characters: Optional[np.ndarray]
    ancestral_sequence: str
    aliases = Dict[str, str]

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
        self.pi_1 = 0.5
        self.frequency = np.asarray((0.5, 0.5))
        self.coefficient_bl = 1.0
        self.pmatrix = None
        self.log_likelihood_vector = None
        self.log_likelihood = None
        self.likelihood_vector = None
        self.likelihood = None
        self.sequence_likelihood = 1.0
        self.up_vector = None
        self.down_vector = None
        self.marginal_vector = None
        self.marginal_bl_vector = None
        self.probability_vector = None
        self.branch_probability_vector = None
        self.probability_vector_gain = None
        self.probability_vector_loss = None
        self.sequence = ''
        self.probabilities_sequence_characters = None
        self.ancestral_sequence = ''
        self.aliases = dict()

    def __str__(self) -> str:
        return self.get_name(True)

    def __dir__(self) -> list:
        return ['father', 'children', 'name', 'distance_to_father', 'distance_to_root', 'distance_to_root_vector',
                'distance_to_nearest', 'distance_to_father_taking_into_coefficient',
                'distance_to_root_taking_into_coefficient', 'distance_to_root_vector_taking_into_coefficient',
                'distance_to_nearest_taking_into_coefficient', 'level', 'levels_to_nearest', 'alphabet', 'pi_1',
                'frequency', 'coefficient_bl', 'pmatrix', 'log_likelihood_vector', 'log_likelihood',
                'sequence_likelihood', 'likelihood', 'up_vector', 'down_vector', 'marginal_vector',
                'marginal_bl_vector', 'probability_vector', 'branch_probability_vector', 'probability_vector_gain',
                'probability_vector_loss', 'sequence', 'probabilities_sequence_characters', 'ancestral_sequence',
                'aliases']

    def get_list_nodes_info(self, with_additional_details: bool = False,
                            mode: Optional[str] = None,
                            filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]] = None,
                            only_node_list: bool = False,
                            fields: Optional[Set[str]] = None
                            ) -> List[Union[Dict[str, Union[float, np.float64, bool, str, np.ndarray, List[float],
                                      List[np.float64]]], 'Node']]:
        """
        Retrieve a list of descendant nodes from a given node, including the node itself.

        This function collects all child nodes of the specified `node`, including the `node` itself. The function
        returns a list of nodes or a list of dictionaries with information about these nodes.

        Args:
            with_additional_details (bool, optional): `False` (default).
            mode (str, optional): None (default), 'pre-order', 'in-order', 'post-order', 'level-order'.
            filters (Dict, optional):
            only_node_list (Dict, optional): `False` (default).
            fields (Set[str], optional): `None` (default) builds every field, as before. When given, only these
                keys are built (see ``get_node_info``) -- restricting this away from the big per-node arrays
                (up_vector/down_vector/marginal_vector/marginal_bl_vector/pmatrix) when a caller doesn't need them
                is the difference between a cheap call and one that JSON round-trips gigabytes of numpy arrays.
        Returns:
            list: A list of descendant nodes from a given node, including the node itself or a list of dictionaries
            with information about these nodes.
        """
        list_result = []
        mode = 'pre-order' if mode is None or mode.lower() not in ('pre-order', 'in-order', 'post-order', 'level-order'
                                                                   ) else mode.lower()

        def resolve_item(trees_node: 'Node') -> Union[str, 'Node', Dict[str, Any]]:
            if only_node_list:
                return trees_node
            if with_additional_details:
                return trees_node.get_node_info(fields)
            return trees_node.name

        def get_list(trees_node: Node) -> None:
            nonlocal list_result, filters, mode

            if trees_node.check_filter_compliance(filters):
                list_item = resolve_item(trees_node)
                if mode == 'pre-order':
                    list_result.append(list_item)

                for i, child in enumerate(trees_node.children):
                    get_list(child)
                    if mode == 'in-order' and not i:
                        list_result.append(list_item)

                if not trees_node.children:
                    if mode == 'in-order':
                        list_result.append(list_item)

                if mode == 'post-order':
                    list_result.append(list_item)
            else:
                for child in trees_node.children:
                    get_list(child)

        if mode == 'level-order':
            nodes_list = [self]
            while nodes_list:
                newick_node = nodes_list.pop(0)
                if newick_node.check_filter_compliance(filters):
                    list_result.append(resolve_item(newick_node))

                for nodes_child in newick_node.children:
                    nodes_list.append(nodes_child)
        else:
            get_list(self)

        return list_result

    _HEAVY_FIELDS = frozenset({'up_vector', 'down_vector', 'marginal_vector', 'marginal_bl_vector', 'pmatrix'})

    def get_node_info(self, fields: Optional[Set[str]] = None
                      ) -> Dict[str, Union[float, np.float64, bool, str, np.ndarray, List[float],
                                List[np.float64]]]:
        """
        Args:
            fields (Set[str], optional): `None` (default) returns every field, as before. When given, only
                these keys are included -- in particular, skipping the heavy per-node tensors
                (``_HEAVY_FIELDS``) here means they're never handed to the ``loads(dumps(...))`` round-trip
                below, which is what actually costs memory: that round-trip turns each packed numpy array
                into nested Python lists of boxed floats, at several times the array's own footprint, and
                CPython doesn't return that memory to the OS as calls repeat. A caller that doesn't need
                these fields (most don't -- see Tree.get_columns) should always pass `fields`.
        """
        want = (lambda k: fields is None or k in fields)

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
                  'likelihood': self.likelihood,
                  'sequence_likelihood': self.sequence_likelihood,
                  'log_likelihood': self.log_likelihood,
                  'log_likelihood_vector': self.log_likelihood_vector,
                  'probability_vector': self.probability_vector,
                  'sequence': self.sequence,
                  'probabilities_sequence_characters': self.probabilities_sequence_characters,
                  'ancestral_sequence': self.ancestral_sequence,
                  'branch_probability_vector': self.branch_probability_vector,
                  'probability_vector_gain': self.probability_vector_gain,
                  'probability_vector_loss': self.probability_vector_loss,
                  'alphabet': self.alphabet,
                  'pi_1': self.pi_1,
                  'frequency': self.frequency,
                  'coefficient_bl': self.coefficient_bl}
        result = {k: v for k, v in result.items() if want(k)}
        for key in self._HEAVY_FIELDS:
            if want(key):
                result[key] = getattr(self, key)

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

    def get_pmatrix(self, rate: Union[float, np.float64, np.ndarray] = 1.0):

        return self.get_one_parameter_pmatrix(rate)

    def calculate_up(self, rate_vector_length: int, alphabet_length: int, msa_length: int) -> None:
        total_up = np.ones((rate_vector_length, alphabet_length, msa_length))

        for child in self.children:
            child_contribution = np.einsum('rji,ril->rjl', child.pmatrix, child.up_vector)

            total_up *= child_contribution

        self.up_vector = total_up

        weighted_vector = self.up_vector * self.frequency[:, np.newaxis]

        likelihood_per_site = np.sum(np.mean(weighted_vector, axis=0), axis=0)

        invalid_mask = (likelihood_per_site <= 0.0) | np.isnan(likelihood_per_site)
        likelihood_per_site = np.where(invalid_mask, eps, likelihood_per_site)

        self.likelihood_vector = likelihood_per_site
        self.likelihood = np.prod(likelihood_per_site)
        self.log_likelihood_vector = np.log(likelihood_per_site)
        self.log_likelihood = np.sum(self.log_likelihood_vector)

    def calculate_down(self, rate_vector_length: int, alphabet_length: int, msa_length: int) -> None:
        total_down = np.ones((rate_vector_length, alphabet_length, msa_length))
        if not self.father:
            self.down_vector = total_down

        else:
            brothers = [b for b in self.father.children if b != self]
            for brother in brothers:
                brother_contrib = np.einsum('rji,ril->rjl', brother.pmatrix, brother.up_vector)
                total_down *= brother_contrib

            if self.father.father:
                father_contrib = np.einsum('rji,ril->rjl', self.father.pmatrix, self.father.down_vector)
                total_down *= father_contrib

        self.down_vector = total_down

    def calculate_marginal(self, rate_vector_length: int, msa_length: int) -> None:
        marg = np.einsum('i,rij,ril->rjl', self.frequency, self.pmatrix, self.down_vector)

        self.marginal_vector = self.up_vector * marg
        self.marginal_bl_vector = np.einsum('i,ril,rij,rjl->rjil', self.frequency, self.down_vector,
                                            self.pmatrix, self.up_vector)

        likelihoods = np.sum(self.marginal_vector, axis=(0, 1)) / rate_vector_length
        invalid_mask = (likelihoods == 0.0) | np.isnan(likelihoods) | np.isinf(likelihoods)
        likelihoods = np.where(invalid_mask, eps, likelihoods)
        summed_marginal = np.sum(self.marginal_vector, axis=0)
        summed_marginal_bl = np.sum(self.marginal_bl_vector, axis=0)
        current_branch_prob = summed_marginal_bl / (rate_vector_length * likelihoods[None, None, :])

        self.probability_vector = (summed_marginal / (rate_vector_length * likelihoods)).T
        self.branch_probability_vector = current_branch_prob.transpose(2, 0, 1).reshape(msa_length, -1)
        self.probability_vector_loss = self.branch_probability_vector[:, 1]
        self.probability_vector_gain = self.branch_probability_vector[:, 2]

        max_indices = np.argmax(self.probability_vector, axis=1)
        reconstructed_chars = np.array(self.alphabet)[max_indices]

        self.probabilities_sequence_characters = self.probability_vector[np.arange(msa_length), max_indices]
        self.sequence = ''.join(reconstructed_chars)

    def clean_all(self):
        self.log_likelihood_vector = None
        self.log_likelihood = None
        self.likelihood_vector = None
        self.likelihood = None
        self.sequence_likelihood = 1.0
        self.up_vector = None
        self.down_vector = None
        self.marginal_vector = None
        self.marginal_bl_vector = None
        self.probability_vector = None
        self.branch_probability_vector = None
        self.probability_vector_gain = None
        self.probability_vector_loss = None
        self.sequence = ''
        self.probabilities_sequence_characters = None
        self.ancestral_sequence = ''

    def get_one_parameter_pmatrix(self, rate: Union[float, np.float64, np.ndarray] = 1.0) -> np.ndarray:
        qmatrix = np.zeros((2, 2), dtype=np.float64)
        qmatrix[0, 0] = - 1 / (2 * (1 - self.pi_1))
        qmatrix[0, 1] = 1 / (2 * (1 - self.pi_1))
        qmatrix[1, 0] = 1 / (2 * self.pi_1)
        qmatrix[1, 1] = - 1 / (2 * self.pi_1)

        return expm(qmatrix * (self.distance_to_father * self.coefficient_bl * rate))

    def get_jukes_cantor_pmatrix(self, alphabet_length: int, rate: Union[float, np.float64, np.ndarray] = 1
                                 ) -> np.ndarray:
        qmatrix = np.ones((alphabet_length, alphabet_length))
        np.fill_diagonal(qmatrix, 1 - alphabet_length)
        qmatrix = qmatrix * 1 / (alphabet_length - 1)

        return expm(qmatrix * (self.distance_to_father * self.coefficient_bl * rate))

    def get_jukes_cantor_transition_probs(self, alphabet_length: int,
                                          rate: Union[float, np.float64, np.ndarray, Any] = 1.0
                                          ) -> Union[Tuple[np.ndarray, np.ndarray], Tuple[np.float64, np.float64],
                                                     Tuple[float, float]]:
        other_states = alphabet_length - 1

        branch_length = self.distance_to_father * self.coefficient_bl * rate
        exponent = np.exp((-alphabet_length / other_states) * branch_length)

        p_identity = (1 / alphabet_length) + (other_states / alphabet_length) * exponent
        p_mutation = (1 / alphabet_length) - (1 / alphabet_length) * exponent

        return p_identity, p_mutation

    def node_to_json(self) -> Dict[str, Union[str, List[Any], float, np.float64, np.ndarray]]:
        dict_json = {'name': self.name, 'distance': f'{float(self.distance_to_father)}'}

        if self.children:
            dict_json.update({'children': []})
            for child in self.children:
                dict_json['children'].append(child.node_to_json())

        return dict_json

    def get_distance_to_father(self, taking_into_coefficient: bool) -> Union[float, np.float64, np.ndarray]:

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

    def get_filter_value(self, key: str) -> Any:
        if key == 'father_name':
            return self.father.name if self.father else ''
        if key == 'children':
            return [i.name for i in self.children]

        return getattr(self, self.aliases.get(key, key), None)

    def check_filter_compliance(self, filters: Optional[Dict[str, List[Union[float, int, str, List[float]]]]]) -> bool:
        if not filters:
            return True

        return any(self.get_filter_value(key) == value for key, values in filters.items() for value in values)
