from gloome.config import Config

if __name__ == '__main__':
    conf = Config()
    conf.check_and_set_input_and_output_variables()
    conf.execute_calculation()
