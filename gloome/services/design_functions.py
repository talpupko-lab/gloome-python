from typing import Union, Tuple, Optional, Dict, List, Set

STYLE_TAG = (('', ''),                                                  # 0
             ('<sub class="text-primary-emphasis">', '</sub>'),         # 1
             ('<sub class="text-secondary-emphasis">', '</sub>'),       # 2
             ('<sub class="text-success-emphasis">', '</sub>'),         # 3
             ('<sub class="text-info-emphasis">', '</sub>'),            # 4
             ('<sub class="text-warning-emphasis">', '</sub>'),         # 5
             ('<sub class="text-danger-emphasis">', '</sub>'),          # 6
             ('<sub class="text-light-emphasis">', '</sub>'),           # 7
             ('<sub class="text-dark-emphasis">', '</sub>'),            # 8

             ('<span class="text-primary">', '</span>'),                # 9
             ('<span class="text-secondary">', '</span>'),              # 10
             ('<span class="text-success">', '</span>'),                # 11
             ('<span class="text-info">', '</span>'),                   # 12
             ('<span class="text-warning">', '</span>'),                # 13
             ('<span class="text-danger">', '</span>'),                 # 14
             ('<span class="text-light">', '</span>'),                  # 15
             ('<span class="text-dark">', '</span>'),                   # 16

             ('<span class="text-primary-emphasis">', '</span>'),       # 17
             ('<span class="text-secondary-emphasis">', '</span>'),     # 18
             ('<span class="text-success-emphasis">', '</span>'),       # 19
             ('<span class="text-info-emphasis">', '</span>'),          # 20
             ('<span class="text-warning-emphasis">', '</span>'),       # 21
             ('<span class="text-danger-emphasis">', '</span>'),        # 22
             ('<span class="text-light-emphasis">', '</span>'),         # 23
             ('<span class="text-dark-emphasis">', '</span>')           # 24
             )


def key_design(key: str, change_style: bool = True, style: int = 22) -> str:
    style_tags = STYLE_TAG[style] if change_style else STYLE_TAG[0]
    return f'{style_tags[0]}{key.replace("_", " ")}:\t{style_tags[1]}'


def value_design(value: Optional[Union[str, Tuple[str, ...], List[str], Set[str]]], change_style: bool = True,
                 style: int = 20) -> str:
    style_tag = STYLE_TAG[style] if change_style else STYLE_TAG[0]
    if isinstance(value, (tuple, list, set)):
        return f'{style_tag[0]}{" ".join(map(str, value))}{style_tag[1]}'
    return f'{style_tag[0]}{value}{style_tag[1]}'


def dna_design(dna: str, different_color: Optional[Tuple[int, int]] = None, styles: Tuple[int, int] = (14, 13)) -> str:
    style_tag = STYLE_TAG[styles[0]] + STYLE_TAG[styles[1]]
    start = False
    start_dc = False
    str_result = ''
    for i in enumerate(dna):
        if different_color and different_color[0] == i[0]:
            j = f'{style_tag[2]}{i[1]}'
            start_dc = True
        elif i[1].upper() != i[1] and not start and not start_dc:
            j = f'{style_tag[0]}{i[1]}'
            start = True
        elif i[1].upper() == i[1] and start and not start_dc:
            j = f'{style_tag[1]}{i[1]}'
            start = False
        else:
            j = i[1]

        if different_color and different_color[1] == i[0]+1:
            j += f'{style_tag[3]}'
            start_dc = False

        str_result += j

    return f'<b>{str_result}</b>'


def result_design(data: Dict[str, Union[str, int, float]], change_key: bool = True, change_value: bool = True,
                  change_key_style: bool = True, change_value_style: bool = True) -> Dict[str, Union[str, int, float]]:
    result_data = dict()
    for key, value in data.items():
        if change_key:
            key = key_design(key, change_key_style)
        if change_value:
            value = value_design(value, change_value_style)
        result_data.update({key: value})
    return result_data
