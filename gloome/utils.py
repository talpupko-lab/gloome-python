import datetime
import logging

from time import time
from random import randint

from gloome.consts import *

formatter = logging.Formatter('%(asctime)s[%(levelname)s][%(filename)s][%(funcName)s]: %(message)s')


# init_dir_path()
# logging_file_name = LOGS_DIR.joinpath(datetime.datetime.now().strftime('%Y_%m_%d_%H_%M.log'))
# FORMAT = '%(asctime)s[%(levelname)s][%(filename)s][%(funcName)s]: %(message)s'
# formatter = logging.Formatter('%(asctime)s[%(levelname)s][%(filename)s][%(funcName)s]: %(message)s')
# logging.basicConfig(filename=logging_file_name, level=logging.WARNING, format=FORMAT)
# logger = logging.getLogger('main')
#
#
def init_dir_path():
    chdir(LOGS_DIR)


def current_time() -> str:
    now = datetime.datetime.now()
    return now.strftime("%H:%M:%S %Y-%m-%d")


def check_dir(file_path: Path, **kwargs) -> None:
    if not file_path.exists():
        file_path.mkdir(mode=kwargs.get('mode', 0o777), parents=kwargs.get('parents', True),
                        exist_ok=kwargs.get('exist_ok', True))


def get_job_logger(job_id: str, logs_dir: Path) -> logging.Logger:
    check_dir(logs_dir)
    current_logger = logging.getLogger(job_id)
    current_logger.setLevel(logging.INFO)
    if len(current_logger.handlers) == 0:
        log_file = logs_dir.joinpath(f'{job_id}.log')
        handler = logging.FileHandler(log_file, 'a', 'utf-8')
        handler.setFormatter(formatter)
        current_logger.addHandler(handler)
        current_logger.propagate = False
    return current_logger


def get_new_process_id() -> str:
    return f'{round(time())}{randint(1000, 9999)}'
