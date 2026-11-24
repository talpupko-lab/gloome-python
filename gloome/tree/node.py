import numpy as np

from typing import Optional, Dict, Union, List, Tuple, Any, Set
from scipy.linalg import expm
from scipy.special import logsumexp
from json import loads, dumps

from .npencoder import NpEncoder

eps = 5e-324
# Conditional vectors whose largest entry (per site pattern) drops below this are rescaled to 1, and the removed factor
# is kept as a log (up_log_scale / down_log_scale), so large trees do not underflow float64.
scaling_threshold = 1e-100


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
    site_pattern_index: Optional[np.ndarray]
    up_log_scale: Optional[np.ndarray]
    down_log_scale: Optional[np.ndarray]
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
        self.site_pattern_index = None
        self.up_log_scale = None
        self.down_log_scale = None
        self.aliases = dict()

    def __str__(self) -> str:
        return self.get_name(True)

    def __dir__(self) -> list:
        return ['father', 'children', 'name', 'distance_to_father', 'distance_to_root', 'distance_to_root_vector',
                'distance_to_nearest', 'distance_to_father_taking_into_coefficient',
                'distance_to_root_taking_into_coefficient', 'distance_to_root_vector_taking_into_coefficient',
                'distance_to_nearest_taking_into_coefficient', 'level', 'levels_to_nearest', 'alphabet', 'pi_1',
                'frequency', 'coefficient_bl', 'pmatrix', 'log_likelihood_vector', 'log_likelihood',
                'likelihood_vector', 'likelihood', 'up_vector', 'down_vector', 'marginal_vector', 'marginal_bl_vector',
                'probability_vector', 'branch_probability_vector', 'probability_vector_gain', 'probability_vector_loss',
                'sequence', 'probabilities_sequence_characters', 'ancestral_sequence', 'site_pattern_index',
                'up_log_scale', 'down_log_scale', 'aliases']

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
            filters (Dict, optional): `None` (default).
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
    # Engine tensors kept per unique site pattern (last axis) -- expanded back to sites for output.
    _PATTERN_FIELDS = frozenset({'up_vector', 'down_vector', 'marginal_vector', 'marginal_bl_vector'})
    _LIKELIHOOD_FIELDS = frozenset({'likelihood', 'likelihood_vector', 'log_likelihood', 'log_likelihood_vector'})

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

        # The up pass computes the subtree likelihood only at the root; other nodes compute it on request.
        if (self.likelihood_vector is None and self.up_vector is not None and self.children
                and any(want(k) for k in self._LIKELIHOOD_FIELDS)):
            self.calculate_likelihood()

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
                  'likelihood_vector': self.likelihood_vector,
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
                value = getattr(self, key)
                if key in self._PATTERN_FIELDS and value is not None:
                    value = self.expand_sites(self.unscale(key, value))
                result[key] = value

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

    def expand_sites(self, pattern_values: np.ndarray, axis: int = -1) -> np.ndarray:
        """Map values computed per unique site pattern back to every MSA site along `axis`."""
        if self.site_pattern_index is None:
            return pattern_values

        return np.take(pattern_values, self.site_pattern_index, axis=axis)

    @staticmethod
    def rescale(conditional_vector: np.ndarray) -> np.ndarray:
        """
        In place: for every (rate category, site pattern) whose largest entry over states is below
        `scaling_threshold`, divide by that entry. Rate categories are scaled separately because one subtree can fit
        a site only under a slow rate and its sibling only under a fast one: a shared factor would let every
        category underflow at their parent.

        Returns:
            np.ndarray: (rate, pattern) ln of the removed factor (0.0 where nothing was rescaled).
        """
        state_max = conditional_vector.max(axis=1)
        needs_scaling = (state_max > 0.0) & (state_max < scaling_threshold)
        log_scale = np.zeros(state_max.shape, dtype=np.float64)
        if needs_scaling.any():
            factor = np.where(needs_scaling, state_max, 1.0)
            conditional_vector /= factor[:, np.newaxis, :]
            log_scale[needs_scaling] = np.log(state_max[needs_scaling])

        return log_scale

    @staticmethod
    def get_rate_weights(log_scale: np.ndarray) -> Optional[np.ndarray]:
        """
        (rate, pattern) weights exp(log_scale - max over rates) that put rescaled rate categories back on a common
        scale before summing overrates; the per-pattern constant left out cancels in every per-site ratio.
        `None` when nothing was rescaled.
        """
        if not log_scale.any():
            return None

        return np.exp(log_scale - log_scale.max(axis=0))

    def unscale(self, key: str, value: np.ndarray) -> np.ndarray:
        """True (unscaled) value of a pattern-space engine tensor, for output."""
        log_scales = [scale for scale, used in ((self.up_log_scale, key != 'down_vector'),
                                                (self.down_log_scale, key != 'up_vector')
                                                ) if used and scale is not None]
        log_scale = sum(log_scales) if log_scales else 0.0
        if not np.any(log_scale):
            return value

        shape = (value.shape[0], ) + (1, ) * (value.ndim - 2) + (value.shape[-1], )

        return value * np.exp(log_scale).reshape(shape)

    def calculate_up(self, rate_vector_length: int, alphabet_length: int, patterns_length: int) -> None:
        total_up = np.ones((rate_vector_length, alphabet_length, patterns_length))
        log_scale = np.zeros((rate_vector_length, patterns_length), dtype=np.float64)

        for index, child in enumerate(self.children, start=1):
            child_contribution = np.matmul(child.pmatrix, child.up_vector)

            total_up *= child_contribution
            log_scale = log_scale + child.up_log_scale
            if index % 2 == 0 or index == len(self.children):
                log_scale = log_scale + self.rescale(total_up)

        self.up_log_scale = log_scale
        self.up_vector = total_up
        self.likelihood_vector = self.likelihood = self.log_likelihood_vector = self.log_likelihood = None

    def calculate_likelihood(self) -> None:
        weighted_vector = self.up_vector * self.frequency[:, np.newaxis]

        if not self.up_log_scale.any():
            likelihood_per_pattern = np.sum(np.mean(weighted_vector, axis=0), axis=0)

            invalid_mask = (likelihood_per_pattern <= 0.0) | np.isnan(likelihood_per_pattern)
            likelihood_per_pattern = np.where(invalid_mask, eps, likelihood_per_pattern)
            log_likelihood_per_pattern = np.log(likelihood_per_pattern)
        else:
            likelihood_per_rate = np.sum(weighted_vector, axis=1)
            with np.errstate(divide='ignore'):
                log_likelihood_per_rate = np.log(likelihood_per_rate) + self.up_log_scale
            log_likelihood_per_pattern = logsumexp(log_likelihood_per_rate, axis=0) - np.log(len(likelihood_per_rate))
            log_likelihood_per_pattern = np.where(np.isfinite(log_likelihood_per_pattern), log_likelihood_per_pattern,
                                                  np.log(eps))
            likelihood_per_pattern = np.maximum(np.exp(log_likelihood_per_pattern), eps)

        self.likelihood_vector = self.expand_sites(likelihood_per_pattern)
        self.likelihood = np.prod(self.likelihood_vector)
        self.log_likelihood_vector = self.expand_sites(log_likelihood_per_pattern)
        self.log_likelihood = np.sum(self.log_likelihood_vector)

    def calculate_down(self, rate_vector_length: int, alphabet_length: int, patterns_length: int) -> None:
        total_down = np.ones((rate_vector_length, alphabet_length, patterns_length))
        log_scale = np.zeros((rate_vector_length, patterns_length), dtype=np.float64)
        if self.father:
            brothers = [b for b in self.father.children if b != self]
            factors = [(brother.pmatrix, brother.up_vector, brother.up_log_scale) for brother in brothers]
            if self.father.father:
                factors.append((self.father.pmatrix, self.father.down_vector, self.father.down_log_scale))

            for index, (pmatrix, conditional_vector, factor_log_scale) in enumerate(factors, start=1):
                total_down *= np.matmul(pmatrix, conditional_vector)
                log_scale = log_scale + factor_log_scale
                if index % 2 == 0 or index == len(factors):
                    log_scale = log_scale + self.rescale(total_down)

        self.down_log_scale = log_scale
        self.down_vector = total_down

    def calculate_marginal(self, rate_vector_length: int, patterns_length: int) -> None:
        marg = np.einsum('i,rij,ril->rjl', self.frequency, self.pmatrix, self.down_vector)

        self.marginal_vector = self.up_vector * marg
        self.marginal_bl_vector = np.einsum('i,ril,rij,rjl->rjil', self.frequency, self.down_vector,
                                            self.pmatrix, self.up_vector)

        marginal_vector, marginal_bl_vector = self.marginal_vector, self.marginal_bl_vector
        rate_weights = self.get_rate_weights(self.up_log_scale + self.down_log_scale)
        if rate_weights is not None:
            marginal_vector = marginal_vector * rate_weights[:, np.newaxis, :]
            marginal_bl_vector = marginal_bl_vector * rate_weights[:, np.newaxis, np.newaxis, :]

        likelihoods = np.sum(marginal_vector, axis=(0, 1)) / rate_vector_length
        invalid_mask = (likelihoods == 0.0) | np.isnan(likelihoods) | np.isinf(likelihoods)
        likelihoods = np.where(invalid_mask, eps, likelihoods)
        summed_marginal = np.sum(marginal_vector, axis=0)
        summed_marginal_bl = np.sum(marginal_bl_vector, axis=0)
        current_branch_prob = summed_marginal_bl / (rate_vector_length * likelihoods[None, None, :])

        probability_per_pattern = (summed_marginal / (rate_vector_length * likelihoods)).T
        branch_probability_per_pattern = current_branch_prob.transpose(2, 0, 1).reshape(patterns_length, -1)
        max_indices = np.argmax(probability_per_pattern, axis=1)
        reconstructed_chars = np.array(self.alphabet)[max_indices]

        self.probability_vector = self.expand_sites(probability_per_pattern, axis=0)
        self.branch_probability_vector = self.expand_sites(branch_probability_per_pattern, axis=0)
        self.probability_vector_loss = self.branch_probability_vector[:, 1]
        self.probability_vector_gain = self.branch_probability_vector[:, 2]
        self.probabilities_sequence_characters = self.expand_sites(
            probability_per_pattern[np.arange(patterns_length), max_indices])
        self.sequence = ''.join(self.expand_sites(reconstructed_chars))

    def clean_all(self):
        self.up_log_scale = None
        self.down_log_scale = None
        self.log_likelihood_vector = None
        self.log_likelihood = None
        self.likelihood_vector = None
        self.likelihood = None
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

        return self.get_two_state_pmatrix(self.pi_1, self.distance_to_father * self.coefficient_bl * rate)

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

    @staticmethod
    def get_two_state_pmatrix(pi_1: Union[float, np.float64], branch_length: Union[float, np.float64, np.ndarray]
                              ) -> np.ndarray:
        """
        Exact expm(Q * t) for the 2-state generator Q = [[-a, a], [b, -b]], a = 1/(2(1-pi_1)), b = 1/(2 pi_1),
        whose mean rate is 1:  P(t) = Pi + exp(-(a + b) t) (I - Pi),  where every row of Pi is (1 - pi_1, pi_1).

        Args:
            pi_1 (float, np.float64, optional).
            branch_length: scalar or array of t values (branch length * coefficient * rate).

        Returns:
            np.ndarray: shape branch_length.shape + (2, 2), rows = parent state, columns = child state.
        """
        stationary = np.array([[1 - pi_1, pi_1], [1 - pi_1, pi_1]], dtype=np.float64)
        decay = np.exp(-np.asarray(branch_length, dtype=np.float64) / (2 * pi_1 * (1 - pi_1)))

        return stationary + decay[..., np.newaxis, np.newaxis] * (np.eye(2) - stationary)

    @staticmethod
    def get_two_state_pmatrix_derivative(pi_1: Union[float, np.float64],
                                         branch_length: Union[float, np.float64, np.ndarray]) -> np.ndarray:
        """
        d/dt of `get_two_state_pmatrix`:  P'(t) = -(a + b) exp(-(a + b) t) (I - Pi).

        Returns:
            np.ndarray: shape branch_length.shape + (2, 2), rows = parent state, columns = child state.
        """
        stationary = np.array([[1 - pi_1, pi_1], [1 - pi_1, pi_1]], dtype=np.float64)
        total_rate = 1 / (2 * pi_1 * (1 - pi_1))
        decay = total_rate * np.exp(-total_rate * np.asarray(branch_length, dtype=np.float64))

        return -decay[..., np.newaxis, np.newaxis] * (np.eye(2) - stationary)
