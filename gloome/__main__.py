from gloome.config import Config


def main():
    conf = Config()
    conf.check_and_set_input_and_output_variables()
    conf.execute_calculation()


if __name__ == '__main__':
    main()
