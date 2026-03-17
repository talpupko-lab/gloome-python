import inspect
import json

from pathlib import Path
from typing import Callable, Any
from datetime import timedelta
from shutil import make_archive, move
from numpy import ndarray
from validate_email import validate_email
from flask import url_for

from gloome.tree.tree import Tree
from gloome.services.design_functions import *

SELECTED_FILES = {'file_interactive_tree_html': True,
                  'file_newick_tree_png': True,
                  'file_table_of_nodes_tsv': True,
                  'file_probability_per_pos_per_branches_tsv': True,
                  'file_table_of_branches_tsv': True,
                  'file_log_likelihood_tsv': True,
                  'file_table_of_attributes_tsv': True,
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


def execute_all_actions(newick_tree: Union[str, Tree], file_path: Union[str, Path], create_new_file: bool = False,
                        form_data: Optional[Dict[str, Union[str, int, float, ndarray]]] = None,
                        log_file: Optional[str] = None, with_internal_nodes: bool = True,
                        actions: Optional[Dict[str, bool]] = None, selected_files: Optional[Dict[str, bool]] = None
                        ) -> Union[Dict[str, str], Path]:
    result_data = {}
    if actions is None or actions.get('draw_tree', False):
        result_data.update({'draw_tree': draw_tree(newick_tree)})
    if actions is None or actions.get('compute_likelihood_of_tree', False):
        result_data.update({'compute_likelihood_of_tree': compute_likelihood_of_tree(newick_tree)})
    if actions is None or actions.get('create_all_file_types', False):
        result_data.update({'create_all_file_types': create_all_file_types(newick_tree, file_path, log_file,
                                                                           with_internal_nodes, selected_files)})
    if create_new_file:
        file_path = file_path.joinpath('result.json') if isinstance(file_path, Path) else f'{file_path}/result.json'
        return create_file(file_path, get_result_data(result_data, 'execute_all_actions', form_data), 'result.json')

    return result_data


def compute_likelihood_of_tree(newick_tree: Union[str, Tree]) -> Union[List[Union[float, ndarray]], str]:

    newick_tree.calculate_likelihood()
    result = [newick_tree.log_likelihood]

    return result


def create_all_file_types(newick_tree: Union[str, Tree], file_path: Union[str, Path],
                          log_file: Optional[Union[str, Path]] = None,
                          with_internal_nodes: Optional[bool] = True, selected_files: Optional[Dict[str, bool]] = None
                          ) -> Union[Dict[str, str], str]:
    selected_files = (SELECTED_FILES if selected_files is None else selected_files)
    result = {}
    newick_tree = Tree.check_tree(newick_tree)
    taking_into_coefficient = newick_tree.coefficient_bl != 1
    # result.update(newick_tree.tree_to_graph(f'{file_path}/graph.txt', ('dot', 'png', 'svg')))
    # result.update(newick_tree.tree_to_visual_format(f'{file_path}/visual_tree.svg', True, ('txt', 'png', 'svg')))
    # result.update({'Newick text (tree)': newick_tree.tree_to_newick_file(f'{file_path}/newick_tree.tree', True)})
    # table = newick_tree.tree_to_table(columns=columns, list_type=list, lists=lists, distance_type=float)
    # result.update({'Fasta (fasta)': newick_tree.tree_to_fasta_file(f'{file_path}/fasta_file.fasta')})
    if selected_files.get('file_interactive_tree_html', False):
        result.update({'Interactive tree (html)':
                       newick_tree.tree_to_interactive_html(f'{file_path}/InteractiveTree.html',
                                                            taking_into_coefficient=taking_into_coefficient)})
    if selected_files.get('file_newick_tree_png', False):
        result.update(newick_tree.tree_to_visual_format(f'{file_path}/VisualTree.svg', with_internal_nodes,
                                                        ('png', ),
                                                        taking_into_coefficient=taking_into_coefficient))
    if selected_files.get('file_table_of_nodes_tsv', False):
        result.update({'Table of nodes (tsv)':
                       newick_tree.tree_to_tsv(f'{file_path}/Nodes.tsv', mode='node_tsv',
                                               taking_into_coefficient=taking_into_coefficient)})
    if selected_files.get('file_probability_per_pos_per_branches_tsv', False):
        result.update({'Probability per positions per branches (tsv)':
                       newick_tree.probability_to_tsv(f'{file_path}/ProbabilityPerPositionsPerBranches.tsv',
                                                      taking_into_coefficient=taking_into_coefficient)})
    if selected_files.get('file_table_of_branches_tsv', False):
        result.update({'Table of branches (tsv)':
                       newick_tree.tree_to_tsv(f'{file_path}/Branches.tsv', mode='branch_tsv',
                                               taking_into_coefficient=taking_into_coefficient)})
    if selected_files.get('file_log_likelihood_tsv', False):
        result.update({'Log-likelihood (tsv)':
                       newick_tree.likelihood_to_tsv(f'{file_path}/LogLikelihood.tsv')})
    if selected_files.get('file_table_of_attributes_tsv', False):
        result.update({'Tree attributes (tsv)':
                       newick_tree.attributes_to_tsv(f'{file_path}/TreeAttributes.tsv')})
    if selected_files.get('file_phylogenetic_tree_nwk', False):
        result.update({'Phylogenetic tree (nwk)':
                       newick_tree.tree_to_newick_file(f'{file_path}/PhylogeneticTree.nwk', True, 0,
                                                       taking_into_coefficient=taking_into_coefficient)})

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
              {'Size factor': min(1 + newick_tree.get_node_count({'node_type': ['leaf']}) // 9, 6)},
              newick_tree.get_json_structure(return_table=True, mode='branch'),
              newick_tree.get_columns_list_for_sorting(mode='branch'),
              {'Sequence length': len(tuple(newick_tree.msa.values())[0])}]

    return result


def convert_seconds(seconds: float) -> str:
    return str(timedelta(seconds=seconds))


def get_error(err_list: List[Tuple[str, str]]) -> str:
    return ''.join([f'{key_design(error_type, True, 14)}'
                    f'{value_design(error, True, 14)}\n' for error_type, error in err_list])


def check_data(*args) -> List[Tuple[str, str]]:
    err_list = []
    newick_text = args[0].strip()
    msa = args[1].strip()
    categories_quantity = int(args[2])
    alpha = float(args[3])
    pi_1 = float(args[4])
    coefficient_bl = float(args[5])
    e_mail = args[6]
    is_optimize_pi = bool(args[7])
    is_optimize_pi_average = bool(args[8])
    is_optimize_alpha = bool(args[9])
    is_optimize_bl = bool(args[10])
    is_do_not_use_e_mail = bool(args[11])
    file_interactive_tree_html = bool(args[12])
    file_newick_tree_png = bool(args[13])
    file_table_of_nodes_tsv = bool(args[14])
    file_probability_per_pos_per_branches_tsv = bool(args[15])
    file_table_of_branches_tsv = bool(args[16])
    file_log_likelihood_tsv = bool(args[17])
    file_table_of_attributes_tsv = bool(args[18])
    file_phylogenetic_tree_nwk = bool(args[19])

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

    if (not isinstance(e_mail, str) or not validate_email(e_mail)) and not is_do_not_use_e_mail:
        err_list.append((f'Invalid email address [ {e_mail} ]', f'Must be valid email address.'))

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

    if not isinstance(is_do_not_use_e_mail, bool):
        err_list.append((f'Do not use e-mail value error [ {is_do_not_use_e_mail} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_interactive_tree_html, bool):
        err_list.append((f'Interactive tree (html) value error [ {file_interactive_tree_html} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_newick_tree_png, bool):
        err_list.append((f'Newick tree (png) value error [ {file_newick_tree_png} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_table_of_nodes_tsv, bool):
        err_list.append((f'Table of nodes (tsv) value error [ {file_table_of_nodes_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_probability_per_pos_per_branches_tsv, bool):
        err_list.append((f'Probability per positions per branches (tsv) value error [ '
                         f'{file_probability_per_pos_per_branches_tsv} ]', f'The value must be boolean type.'))

    if not isinstance(file_table_of_branches_tsv, bool):
        err_list.append((f'Table of branches (tsv) value error [ {file_table_of_branches_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_log_likelihood_tsv, bool):
        err_list.append((f'Log-likelihood (tsv) value error [ {file_log_likelihood_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_table_of_attributes_tsv, bool):
        err_list.append((f'Table of attributes (tsv) value error [ {file_table_of_attributes_tsv} ]',
                         f'The value must be boolean type.'))

    if not isinstance(file_phylogenetic_tree_nwk, bool):
        err_list.append((f'Phylogenetic tree (nwk) value error [ {file_phylogenetic_tree_nwk} ]',
                         f'The value must be boolean type.'))

    if not msa:
        err_list.append(('MSA error', 'No MSA was provided.'))
    elif not msa.startswith('>'):
        err_list.append(('MSA error', 'Wrong MSA format. Please provide MSA in FASTA format.'))
    elif len(msa.split('\n')) / 2 < 2:
        err_list.append(('MSA error', 'There should be at least two sequences in the MSA.'))
    else:
        len_list = []
        incorrect_characters = ''
        for i, current_line in enumerate(msa.split()):
            if i % 2:
                current_line = current_line.strip()
                len_list.append(len(current_line))
                for j in current_line:
                    if j not in '01':
                        incorrect_characters += f'{j} '

        if min(len_list) != max(len_list):
            err_list.append((f'MSA error', f'The MSA contains sequences of different lengths.'))
        if incorrect_characters:
            err_list.append(('MSA error',
                             f'MSA file contains an illegal character(s) [ {incorrect_characters.strip()} ]. '
                             f'Please note that “0” and “1” are the only allowed characters in the phyletic MSAs.'))

        msa_list = msa.strip().split()
        msa_taxa_info = [msa_list[j + j][1::] for j in range(len(msa_list) // 2)]

        if len(msa_taxa_info) != len(msa_taxa_info):
            err_list.append((f'MSA error', f'Duplicate taxa names found.'))

        if not newick_text:
            err_list.append((f'TREE error', f'No Phylogenetic tree was provided.'))
        elif (not (newick_text.startswith('(') and newick_text.endswith(';')) or
              (newick_text.count('(') != newick_text.count(')'))):
            err_list.append((f'TREE error', f'Wrong Phylogenetic tree format. Please provide a tree in Newick format.'))
        else:
            try:
                current_tree = Tree(newick_text)
            except ValueError:
                current_tree = None

            if current_tree:
                for current_node in current_tree.get_list_nodes_info(with_additional_details=True, filters={'distance':
                                                                     [0.0, ]}, only_node_list=True):
                    current_node.distance_to_father = float(f'{current_node.distance_to_father:.4f}1')
                edges_distances_list = current_tree.tree_to_table(filters={'node_type': ['leaf', 'node']},
                                                                  columns={'distance': 'distance'},
                                                                  distance_type=float,
                                                                  taking_into_coefficient=False).T.values[0].tolist()
                if not all(edges_distances_list):
                    err_list.append((f'TREE error', f'One or more branches in the tree have zero length.\n'
                                                    f'{edges_distances_list}'))
                if not (current_tree.get_node_count({'node_type': ['leaf']}) == len(msa.split('\n')) / 2 ==
                        msa.count('>')):
                    err_list.append((f'MSA error',
                                     f'A discrepancy exists between the number of leaves in the phylogenetic tree and '
                                     f'the number of sequences present in the MSA data.'))

                tree_taxa_info = current_tree.tree_to_table(filters={'node_type': ['leaf']}, columns={'node': 'node'},
                                                            taking_into_coefficient=False).T.values[0].tolist()

                if len(tree_taxa_info) != len(set(tree_taxa_info)):
                    err_list.append((f'TREE error', f'Duplicate taxa names found.'))

                if set(tree_taxa_info).difference(set(msa_taxa_info)):
                    err_list.append((f'DATA MISMATCH error',
                                     f'Taxa names in the MSA and phylogenetic tree do not match.'))
            else:
                err_list.append((f'TREE error',
                                 f'Wrong Phylogenetic tree format. Please provide a tree in Newick format.'))

    return err_list


def get_function_parameters(func: Callable) -> Tuple[str, ...]:
    return tuple(inspect.signature(func).parameters.keys())


def recompile_json(output_file: Union[str, Path], process_id: int, create_link: bool) -> str:
    file_contents = read_file(file_path=output_file)
    json_object = loads_json(file_contents)
    action_name = json_object.pop('action_name')
    data = json_object.pop('data') if 'data' in json_object.keys() else json_object.copy()

    if 'execute_all_actions' in action_name:
        for key, value in data.items():
            data.update({key: get_response_design(value, key, create_link, output_file)})
    else:
        data = get_response_design(data, action_name, create_link, output_file)

    data.update({'title': process_id})
    data.update({'form_data': json_object.pop('form_data')})
    data.update({'action_name': action_name})
    create_file(file_path=output_file, data=data)

    return '; '.join(data.keys())


def get_response_design(json_object: Optional[Any], action_name: str, create_link: bool,
                        output_file:  Union[str, Path] = '') -> Optional[Any]:
    if 'create_all_file_types' in action_name and create_link:
        if output_file:
            json_object.update({'json response file (json)': output_file})
        json_object = link_design(json_object)
        json_object = result_design(json_object, change_value='compute_likelihood_of_tree' in action_name,
                                    change_value_style=False, change_key=True, change_key_style=False)
    return json_object


def link_design(json_object: Any) -> Any:
    for key, value in json_object.items():
        if key == 'execution_time':
            continue
        json_object.update(
            {f'{key}': [f'<a class="w-auto mw-auto form-control btn btn-outline-link rounded-pill" href="'
                        f'{url_for("get_file", file_path=value, mode="download")}" '
                        f'target="_blank">download</a>',
                        f'<a class="w-auto mw-auto form-control btn btn-outline-link rounded-pill" href="'
                        f'{url_for("get_file", file_path=value, mode="view")}" target="_blank">view</a>']})
    return json_object
