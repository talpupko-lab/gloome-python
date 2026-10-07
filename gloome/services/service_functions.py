import inspect
import json

from pathlib import Path
from typing import Callable, Any, Union, Tuple, Optional, Dict, List
from datetime import timedelta
from shutil import make_archive, move
from numpy import ndarray

from gloome.tree.tree import Tree

FILE_LIST = [
    'Log-likelihood (tsv)',
    'Table of nodes (tsv)',
    'Table of branches (tsv)',
    'Branch position probabilities (tsv)',
    'Table of coevolution (tsv)',
    'Parsimony and homoplasy scores (tsv)',
    'Table of posterior rates (tsv)',
    'Table of pearson correlation (tsv)',
    'Tree attributes (tsv)',
    'Phylogenetic tree (nwk)',
    'Simulated datasets (fastas)',
    'Log-File (log)',
    'JSON response file (json)',
    'Barplot of correlation (svg)',
    'Plot of correlation by rate-bin (svg)',
    'Plot of distribution of correlation (svg)',
    'Interactive tree (html)',
    'Newick tree (svg)',
    'Newick tree (png)',
    'Newick tree (txt)',
    'Archive (zip)'
]
SELECTED_FILES = {'file_interactive_tree_html': True,
                  'file_newick_tree_png': True,
                  'file_table_of_coevolution_tsv': True,
                  'file_simulated_datasets_fastas': True,
                  'file_table_of_posterior_rates_tsv': True,
                  'file_barplot_of_correlation_svg': True,
                  'file_plot_distribution_of_correlation_svg': True,
                  'file_plot_correlation_by_rate_bin_svg': True,
                  'file_table_of_pearson_correlation_tsv': True,
                  'file_table_of_nodes_tsv': True,
                  'file_branch_position_probabilities_tsv': True,
                  'file_table_of_branches_tsv': True,
                  'file_log_likelihood_tsv': True,
                  'file_table_of_attributes_tsv': True,
                  'file_table_of_parsimony_and_homoplasy_scores_tsv': True,
                  'file_phylogenetic_tree_nwk': True}


def get_digit(data: str) -> Union[int, float, str]:
    try:
        return int(data)
    except ValueError:
        try:
            return float(data)
        except ValueError:
            return str(data)


def get_variables(request_form: Dict[str, str]) -> Tuple[Union[str, int, float], ...]:
    result = [bool(int(v)) if k[:2] == 'is' or k[:4] == 'file' else get_digit(v) for k, v in request_form.items()]

    return tuple(result)


def get_dict(request_form: Dict[str, str]) -> Tuple[Union[str, int, float], ...]:
    result = {k: bool(int(v)) if k[:2] == 'is' or k[:4] == 'file' else get_digit(v) for k, v in request_form.items()}

    return tuple(result)


def get_path(path: Union[str, Path]) -> Path:
    if path and isinstance(path, str):
        return Path(path)

    return path


def create_file(file_path: Union[str, Path], data: Union[str, Any], file_name: Optional[str] = None) -> Path:
    file_path = get_path(file_path)
    if file_name and isinstance(file_name, str):
        file_path.joinpath(file_name)
    save_file(file_path, data)

    return file_path


def del_file(file_path: Union[str, Path]) -> None:
    file_path = get_path(file_path)
    if file_path.is_file():
        file_path.unlink(missing_ok=True)


def read_file(file_path: Union[str, Path], mode: str = 'r') -> str:
    file_path = get_path(file_path)
    if file_path.is_file():
        with open(file_path, mode) as f:
            return f.read()

    return ''


def save_file(file_path: Union[str, Path], data: Union[str, Any], mode: str = 'w') -> None:
    with open(file_path, mode) as f:
        if isinstance(data, str):
            f.write(data)
        else:
            f.write(dumps_json(data))


def loads_json(data: str) -> Any:

    return json.loads(data)


def dumps_json(data: Any) -> str:

    return json.dumps(data)


def del_files(file_list: Union[Union[str, Path], Tuple[Union[str, Path], ...]]) -> None:
    if isinstance(file_list, (str, Path)):
        del_file(file_list)
    else:
        for file in file_list:
            del_file(file)


def get_result_data(data: Union[Dict[str, Union[str, int, float, ndarray, List[Union[float, ndarray]]]],
                                List[Union[float, ndarray, Any]]],
                    action_name: str, form_data: Optional[Dict[str, Union[str, int, float, ndarray]]] = None
                    ) -> Dict[str, Union[str, int, float, ndarray, Dict[str, Union[str, int, float, ndarray]],
                                         List[Union[float, ndarray, Any]]]]:
    result = {'action_name': action_name, 'data': data}
    if form_data is not None:
        result.update({'form_data': form_data})

    return result


def check_tree_data(newick_tree: Union[str, Tree], msa: Union[Dict[str, str], str],
                    alphabet: Optional[Tuple[str, ...]]):
    if isinstance(newick_tree, str):
        newick_tree = Tree.rename_nodes(newick_tree)
    if isinstance(msa, str):
        msa = newick_tree.get_msa_dict(msa)
    if alphabet is None:
        alphabet = Tree.get_alphabet_from_dict(msa)

    return newick_tree, msa, alphabet


def execute_all_actions(newick_tree: Union[str, Tree], file_path: Union[str, Path],
                        log_file: Optional[str] = None, with_internal_nodes: bool = True,
                        actions: Optional[Dict[str, bool]] = None, selected_files: Optional[Dict[str, bool]] = None,
                        use_copap: Optional[bool] = None,
                        probability_lg: Union[float, ndarray] = 0.5,
                        number_lg: Union[float, ndarray, int] = 1,
                        number_datasets: int = 100
                        ) -> Union[Dict[str, str], Path]:
    result_data = {}
    if actions is None or actions.get('draw_tree', False):
        result_data.update({'draw_tree': draw_tree(newick_tree)})
    if actions is None or actions.get('compute_likelihood_of_tree', False):
        result_data.update({'compute_likelihood_of_tree': compute_likelihood_of_tree(newick_tree)})
    if actions is None or actions.get('create_all_file_types', False):
        result_data.update({'create_all_file_types': create_all_file_types(newick_tree, file_path, log_file,
                                                                           with_internal_nodes, selected_files,
                                                                           use_copap, probability_lg, number_lg,
                                                                           number_datasets)})

    return result_data


def compute_likelihood_of_tree(newick_tree: Union[str, Tree]) -> Union[List[Union[float, ndarray]], str]:
    newick_tree.calculate_likelihood()
    result = [newick_tree.log_likelihood]

    return result


def create_all_file_types(newick_tree: Union[str, Tree], file_path: Union[str, Path],
                          log_file: Optional[Union[str, Path]] = None,
                          with_internal_nodes: Optional[bool] = True,
                          selected_files: Optional[Dict[str, bool]] = None,
                          use_copap: Optional[bool] = None,
                          probability_lg: Union[float, ndarray] = 0.5,
                          number_lg: Union[float, ndarray, int] = 1,
                          number_datasets: int = 100
                          ) -> Union[Dict[str, str], str]:
    selected_files = (SELECTED_FILES if selected_files is None else selected_files)
    result = {}
    newick_tree = Tree.check_tree(newick_tree)
    taking_into_coefficient = newick_tree.coefficient_bl != 1
    use_correlation = newick_tree.msa_length > 1
    use_rates = newick_tree.rate_vector_length > 1
    use_simulated_datasets_file = selected_files.get('file_simulated_datasets_fastas', False)
    if use_correlation:
        use_coevolution_file = selected_files.get('file_table_of_coevolution_tsv', False)
        use_barplot_of_correlation_file = selected_files.get('file_barplot_of_correlation_svg', False)
        use_plot_distribution_of_correlation_file = selected_files.get('file_plot_distribution_of_correlation_svg',
                                                                       False)
        use_plot_correlation_by_rate_bin_file = selected_files.get('file_plot_correlation_by_rate_bin_svg', False)
    else:
        use_coevolution_file = False
        use_barplot_of_correlation_file = False
        use_plot_distribution_of_correlation_file = False
        use_plot_correlation_by_rate_bin_file = False

    if selected_files.get('file_interactive_tree_html', False):
        result.update({'Interactive tree (html)':
                       newick_tree.tree_to_interactive_html(file_name=f'{file_path}/InteractiveTree.html',
                                                            taking_into_coefficient=taking_into_coefficient)})
    if selected_files.get('file_newick_tree_png', False):
        result.update(newick_tree.tree_to_visual_format(file_name=f'{file_path}/VisualTree.svg',
                                                        with_internal_nodes=with_internal_nodes,
                                                        taking_into_coefficient=taking_into_coefficient,
                                                        file_extensions=('png', )))
    if selected_files.get('file_table_of_posterior_rates_tsv', False) and use_copap and use_rates:
        result.update({'Table of posterior rates (tsv)':
                       newick_tree.posterior_rates_to_tsv(file_name=f'{file_path}/PosteriorRates.tsv')})
    if selected_files.get('file_table_of_pearson_correlation_tsv', False) and use_copap and use_correlation:
        result.update({'Table of pearson correlation (tsv)':
                       newick_tree.pearson_correlation_to_tsv(file_name=f'{file_path}/PearsonCorrelation.tsv',
                                                              probability_lg=probability_lg,
                                                              number_lg=number_lg)})
    if selected_files.get('file_table_of_nodes_tsv', False):
        result.update({'Table of nodes (tsv)':
                       newick_tree.tree_to_tsv(file_name=f'{file_path}/Nodes.tsv',
                                               taking_into_coefficient=taking_into_coefficient,
                                               mode='node_tsv')})
    if selected_files.get('file_branch_position_probabilities_tsv', False):
        result.update({'Branch position probabilities (tsv)':
                       newick_tree.probability_to_tsv(file_name=f'{file_path}/BranchPositionProbabilities.tsv',
                                                      taking_into_coefficient=taking_into_coefficient)})
    if selected_files.get('file_table_of_branches_tsv', False):
        result.update({'Table of branches (tsv)':
                       newick_tree.tree_to_tsv(file_name=f'{file_path}/Branches.tsv',
                                               taking_into_coefficient=taking_into_coefficient,
                                               mode='branch_tsv')})
    if selected_files.get('file_log_likelihood_tsv', False):
        result.update({'Log-likelihood (tsv)':
                       newick_tree.likelihood_to_tsv(file_name=f'{file_path}/LogLikelihood.tsv')})
    if selected_files.get('file_table_of_attributes_tsv', False):
        result.update({'Tree attributes (tsv)':
                       newick_tree.attributes_to_tsv(file_name=f'{file_path}/TreeAttributes.tsv')})
    if selected_files.get('file_table_of_parsimony_and_homoplasy_scores_tsv', False):
        result.update({'Parsimony and homoplasy scores (tsv)':
                       newick_tree.parsimony_score_to_tsv(file_name=f'{file_path}/ParsimonyAndHomoplasyScores.tsv')})
    if selected_files.get('file_phylogenetic_tree_nwk', False):
        result.update({'Phylogenetic tree (nwk)':
                       newick_tree.tree_to_newick_file(file_name=f'{file_path}/PhylogeneticTree.nwk',
                                                       taking_into_coefficient=taking_into_coefficient,
                                                       with_internal_nodes=with_internal_nodes,
                                                       decimal_length=0)})
    if any((use_simulated_datasets_file, use_coevolution_file, use_barplot_of_correlation_file,
            use_plot_distribution_of_correlation_file, use_plot_correlation_by_rate_bin_file)) and use_copap:
        result.update(newick_tree.simulate_datasets(file_path=f'{file_path}',
                                                    number_datasets=number_datasets,
                                                    use_simulated_datasets_file=use_simulated_datasets_file,
                                                    use_coevolution_file=use_coevolution_file,
                                                    use_barplot_of_correlation_file=use_barplot_of_correlation_file,
                                                    use_plot_distribution_of_correlation_file=
                                                    use_plot_distribution_of_correlation_file,
                                                    use_plot_correlation_by_rate_bin_file=
                                                    use_plot_correlation_by_rate_bin_file))

    if result:
        file_path = get_path(file_path)
        archive_name = Path(make_archive(f'{file_path}', 'zip', f'{file_path}', '.'))
        new_archive_name = file_path.joinpath(archive_name.name)
        move(archive_name, new_archive_name)
        result.update({'Archive (zip)': f'{new_archive_name}'})

    if log_file is not None:
        result.update({'Log-File (log)': f'{log_file}'})

    return result


def draw_tree(newick_tree: Tree) -> Union[List[Any], str]:
    result = [newick_tree.get_json_structure(),
              newick_tree.get_json_structure(return_table=True),
              newick_tree.get_columns_list_for_sorting(),
              {'Size factor': min(1 + newick_tree.get_leaves_count() // 9, 6)},
              newick_tree.get_json_structure(return_table=True, mode='branch'),
              newick_tree.get_columns_list_for_sorting(mode='branch'),
              {'Sequence length': len(tuple(newick_tree.msa.values())[0])}]

    return result


def convert_seconds(seconds: float) -> str:

    return str(timedelta(seconds=seconds))


def del_bootstrap_values(newick_text: str) -> str:

    return Tree.del_bootstrap_values(newick_text)


def get_leaves(data) -> List[str]:

    return Tree(data).get_leaves(only_node_list=False)


def check_data(*args) -> List[Tuple[str, str]]:
    err_list = []
    newick_text = args[0].strip()
    msa = args[1].strip()
    categories_quantity = int(args[2])
    alpha = float(args[3])
    pi_1 = float(args[4])
    coefficient_bl = float(args[5])
    probability_lg = float(args[6])
    number_lg = int(args[7])
    number_datasets = int(args[8])
    is_optimize_pi = bool(args[9])
    is_optimize_pi_average = bool(args[10])
    is_optimize_alpha = bool(args[11])
    is_optimize_bl = bool(args[12])
    is_do_not_use_copap = bool(args[13])
    file_interactive_tree_html = bool(args[14])
    file_newick_tree_png = bool(args[15])
    file_table_of_coevolution_tsv = bool(args[16])
    file_simulated_datasets_fastas = bool(args[17])
    file_table_of_posterior_rates_tsv = bool(args[18])
    file_barplot_of_correlation_svg = bool(args[19])
    file_plot_distribution_of_correlation_svg = bool(args[20])
    file_plot_correlation_by_rate_bin_svg = bool(args[21])
    file_table_of_pearson_correlation_tsv = bool(args[22])
    file_table_of_nodes_tsv = bool(args[23])
    file_branch_position_probabilities_tsv = bool(args[24])
    file_table_of_branches_tsv = bool(args[25])
    file_log_likelihood_tsv = bool(args[26])
    file_table_of_attributes_tsv = bool(args[27])
    file_table_of_parsimony_and_homoplasy_scores_tsv = bool(args[28])
    file_phylogenetic_tree_nwk = bool(args[29])
    rooting_method = args[30].strip()
    leaf = args[31].strip()

    if not isinstance(categories_quantity, int) or not 1 <= categories_quantity <= 16:
        err_list.append((f'Number of rate categories value error [ {categories_quantity} ]',
                         f'The value must be between 1 and 16.'))

    if not isinstance(alpha, float) or not 0.1 <= alpha <= 20:
        err_list.append((f'Alpha value error [ {alpha} ]', f'The value must be between 0.1 and 20.'))

    if not isinstance(pi_1, float) or not 0.001 <= pi_1 <= 0.999:
        err_list.append((f'π1 value error [ {pi_1} ]', f'The value must be between 0.001 and 0.999.'))

    if not isinstance(coefficient_bl, float) or not 0.1 <= coefficient_bl <= 10:
        err_list.append((f'Branch lengths (BL) coefficient value error [ {coefficient_bl} ]',
                         f'The value must be between 0.1 and 10.'))

    if not isinstance(probability_lg, float) or not 0.01 <= probability_lg <= 0.99:
        err_list.append((f'Probability of loss/gain event value error [ {probability_lg} ]',
                         f'The value must be between 0.01 and 0.99.'))

    if not isinstance(number_lg, int) or not 1 <= number_lg <= 20:
        err_list.append((f'Number of loss/gain events value error [ {number_lg} ]',
                         f'The value must be between 1 and 20.'))

    if not isinstance(number_datasets, int) or not 1 <= number_datasets <= 1000:
        err_list.append((f'Number of simulation events value error [ {number_datasets} ]',
                         f'The value must be between 1 and 1000.'))

    if not isinstance(is_optimize_pi, bool):
        err_list.append((f'Optimize π1 value (algorithmic) error [ {is_optimize_pi} ]',
                         f'The value must be boolean type.'))

    if not isinstance(is_optimize_pi_average, bool):
        err_list.append((f'Optimize π1 value (empirical) error [ {is_optimize_pi_average} ]',
                         f'The value must be boolean type.'))

    if not isinstance(is_optimize_alpha, bool):
        err_list.append((f'Optimize α value error [ {is_optimize_alpha} ]', f'The value must be boolean type.'))

    if not isinstance(is_optimize_bl, bool):
        err_list.append((f'Optimize branch lengths coefficient value error [ {is_optimize_bl} ]',
                         f'The value must be boolean type.'))

    if not isinstance(is_do_not_use_copap, bool):
        err_list.append((f'Do not use CoPAP value error [ {is_do_not_use_copap} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_interactive_tree_html, bool):
        err_list.append((f'Interactive tree (html) value error [ {file_interactive_tree_html} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_newick_tree_png, bool):
        err_list.append((f'Newick tree (png) value error [ {file_newick_tree_png} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_table_of_coevolution_tsv, bool):
        err_list.append((f'Coevolution (tsv) value error '
                         f'[ {file_table_of_coevolution_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_simulated_datasets_fastas, bool):
        err_list.append((f'Simulated datasets (fastas) value error '
                         f'[ {file_simulated_datasets_fastas} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_table_of_posterior_rates_tsv, bool):
        err_list.append((f'Table of posterior rates (tsv) value error [ {file_table_of_posterior_rates_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_barplot_of_correlation_svg, bool):
        err_list.append((f'Barplot of correlation (svg) value error '
                         f'[ {file_barplot_of_correlation_svg} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_plot_distribution_of_correlation_svg, bool):
        err_list.append((f'Plot of distribution of correlation (svg) value error '
                         f'[ {file_plot_distribution_of_correlation_svg} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_plot_correlation_by_rate_bin_svg, bool):
        err_list.append((f'Plot of correlation by rate-bin (svg) value error '
                         f'[ {file_plot_correlation_by_rate_bin_svg} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_table_of_pearson_correlation_tsv, bool):
        err_list.append((f'Table of pearson correlation (tsv) value error [ {file_table_of_pearson_correlation_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_table_of_nodes_tsv, bool):
        err_list.append((f'Table of nodes (tsv) value error [ {file_table_of_nodes_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_branch_position_probabilities_tsv, bool):
        err_list.append((f'Branch position probabilities (tsv) value error [ '
                         f'{file_branch_position_probabilities_tsv} ]', f'The value must be boolean type.'))

    if not isinstance(file_table_of_branches_tsv, bool):
        err_list.append((f'Table of branches (tsv) value error [ {file_table_of_branches_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_log_likelihood_tsv, bool):
        err_list.append((f'Log-likelihood (tsv) value error [ {file_log_likelihood_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_table_of_attributes_tsv, bool):
        err_list.append((f'Tree attributes (tsv) value error [ {file_table_of_attributes_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_table_of_parsimony_and_homoplasy_scores_tsv, bool):
        err_list.append((f'Parsimony and homoplasy scores (tsv) value error [ '
                         f'{file_table_of_parsimony_and_homoplasy_scores_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_phylogenetic_tree_nwk, bool):
        err_list.append((f'Phylogenetic tree (nwk) value error [ {file_phylogenetic_tree_nwk} ]',
                         f'The value must be boolean type.'))

    if not isinstance(rooting_method, str):
        err_list.append((f'Rooting method value error [ {rooting_method} ]', f'The value must be string type.'))

    if (not isinstance(leaf, str) or not leaf) and rooting_method == 'outgroup':
        err_list.append((f'Leaf value error [ {leaf} ]', f'The value must be a non-empty string.'))

    if not msa:
        err_list.append(('MSA error', 'No MSA was provided.'))
    elif not msa.startswith('>'):
        err_list.append(('MSA error', 'Wrong MSA format. Please provide MSA in FASTA format.'))
    else:
        msa_list = msa.split()
        msa_list_size = len(msa_list)
        if msa_list_size / 2 < 2:
            err_list.append(('MSA error', 'There should be at least two sequences in the MSA.'))
        else:
            allowed = set('01?')
            all_chars = set()
            msa_taxa_set = set()

            first_len = None
            is_different_lengths = False
            has_duplicate_taxa = False

            for i in range(0, msa_list_size, 2):
                taxa = msa_list[i][1:]
                if taxa in msa_taxa_set:
                    has_duplicate_taxa = True
                msa_taxa_set.add(taxa)

                if i + 1 < msa_list_size:
                    line = msa_list[i + 1]
                    current_len = len(line)

                    if first_len is None:
                        first_len = current_len
                    elif current_len != first_len:
                        is_different_lengths = True

                    all_chars.update(line)

            unique_incorrect = all_chars - allowed
            incorrect_characters = ' '.join(unique_incorrect)

            if is_different_lengths:
                err_list.append(('MSA error', 'The MSA contains sequences of different lengths.'))

            if incorrect_characters:
                err_list.append(('MSA error',
                                 f'MSA file contains an illegal character(s) [ {incorrect_characters} ]. '
                                 f'Please note that “0”, “1” and “?” (missing data) are the only allowed characters '
                                 f'in the phyletic MSAs.'))

            if has_duplicate_taxa:
                err_list.append(('MSA error', 'Duplicate taxa names found.'))

            if not newick_text:
                err_list.append((f'TREE error', f'No Phylogenetic tree was provided.'))
            elif (not (newick_text.startswith('(') and newick_text.endswith(';')) or
                  (newick_text.count('(') != newick_text.count(')'))):
                err_list.append((f'TREE error',
                                 'Wrong Phylogenetic tree format. Please provide a tree in Newick format.'))
            else:
                try:
                    current_tree = Tree(newick_text)
                    Tree.rename_nodes(current_tree)
                except ValueError:
                    current_tree = None

                if current_tree:
                    for current_node in current_tree.get_list_nodes_info(with_additional_details=True,
                                                                         filters={'distance': [0.0, ]},
                                                                         only_node_list=True):
                        current_node.distance_to_father = float(f'{current_node.distance_to_father:.4f}1')
                    edges_distances_list = current_tree.tree_to_table(filters={'node_type': ['leaf', 'node']},
                                                                      columns={'distance': 'distance'},
                                                                      distance_type=float,
                                                                      taking_into_coefficient=False
                                                                      ).T.values[0].tolist()
                    if not all(edges_distances_list):
                        err_list.append((f'TREE error',
                                         f'One or more branches in the tree have zero length.\n'
                                         f'{edges_distances_list}'))
                    if not (current_tree.get_leaves_count() == len(msa.split('\n')) / 2 == msa.count('>')):
                        err_list.append((f'MSA error',
                                         f'A discrepancy exists between the number of leaves in the phylogenetic tree '
                                         f'and the number of sequences present in the MSA data.'))

                    tree_taxa_info = current_tree.get_leaves(only_node_list=False)

                    tree_taxa_set = set(tree_taxa_info)
                    if len(tree_taxa_info) != len(tree_taxa_set):
                        err_list.append((f'TREE error', f'Duplicate taxa names found.'))

                    if tree_taxa_set.difference(msa_taxa_set):
                        err_list.append((f'DATA MISMATCH error',
                                         f'Taxa names in the MSA and phylogenetic tree do not match.'))
                    if not current_tree.all_nodes.get(leaf) and rooting_method == 'outgroup':
                        err_list.append((f'TREE error', f'Leaf {leaf} not found.'))
                else:
                    err_list.append((f'TREE error',
                                     f'Wrong Phylogenetic tree format. Please provide a tree in Newick format.'))

    return err_list


def get_function_parameters(func: Callable) -> Tuple[str, ...]:

    return tuple(inspect.signature(func).parameters.keys())
