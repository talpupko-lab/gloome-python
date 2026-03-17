### GLOOME: Gain Loss Mapping Engine
A bioinformatics tool for analyzing gene gain and loss events during evolution

GLOOME provides evolutionary analysis of presence and absence profiles (phyletic patterns). These patterns are assumed to result from gain and loss dynamics along a phylogenetic tree. Examples of characters represented by phyletic patterns include:

Restriction sites
Gene families
Introns
Indels
The primary purpose of the GLOOME server is to accurately infer branch-specific and site-specific gain and loss events. This inference is based on a stochastic mapping approach, using evolutionary models that accurately reflect the underlying biological processes.

Features
Support for various evolutionary models
Gain and loss inference using stochastic mapping or maximum parsimony
Estimation of gain/loss rates per character
Advanced optimization options
Likelihood and parsimony-based output
A full user manual is included in the release:
GLOOME.CoPAP.gainLoss.Manual.pdf

### Project structure
```
gloome
├── gloome
│   ├── data
│   │   └── initial_data
│   │       ├── tree
│   │       │   ├── newickTree0.tree
│   │       │   └── newickTree1.tree
│   │       └── msa
│   │           ├── patternMSA0.msa
│   │           └── patternMSA1.msa
│   ├── services
│   │   ├── __init__.py
│   │   ├── design_functions.py
│   │   └── service_functions.py
│   ├── tree
│   │   ├── __init__.py
│   │   ├── node.py
│   │   └── tree.py
│   ├── __init__.py
│   ├── __main__.py
│   ├── config.py
│   ├── consts.py
│   └── utils.py
├── logs
├── results
│   ├── out
│   └── in
├── pyproject.toml
└── README.md
```

This manual provides comprehensive Gloome on installation, input formats, command-line usage, interpretation of results, and troubleshooting.

### Installation
Requirements
UNIX/Linux system
GNU Compiler (gcc/g++)

To install the project, you need to run the command "pip install gloome" in the terminal.

### Program Execution
To get the project working, you need to run the command "python -m gloome" in the terminal with the following parameters:

Required parameters:

```
    --msa_file <type=str>
        Specify the msa filepath.

    --tree_file <type=str>
        Specify the newick filepath.
```

Optional parameters:

```
    --out_dir <type=str>
        Specify the outdir path.

    --process_id <type=str>
        Specify a process ID, otherwise it will be generated automatically.

    --mode <type=str>
        Execution mode style. 
        Possible options: ('draw_tree', 'compute_likelihood_of_tree', 'create_all_file_types', 'execute_all_actions'). 
        Default is 'execute_all_actions'.

    --with_internal_nodes <type=int> 
        Specify the tree has internal nodes. Default is 1.

    --categories_quantity <type=int>
        Specify categories quantity. Default is 4.

    --alpha <type=float>
        Specify alpha. Default is 0.5.

    --pi_1 <type=float> 
        Specify pi_1. Default is 0.5.

    --coefficient_bl <type=float> 
        Specify coefficient_bl. Default is 1.0.

    --is_optimize_pi <type=int> 
        Specify is_optimize_pi. Default is 1.

    --is_optimize_pi_average <type=int> 
        Specify is_optimize_pi_average. Default is 0.

    --is_optimize_alpha <type=int> 
        Specify is_optimize_alpha. Default is 1.

    --is_optimize_bl <type=int> 
        Specify is_optimize_bl. Default is 1.

    --file_interactive_tree_html <type=int> 
        Specify file_interactive_tree_html. Default is 0.

    --file_newick_tree_png <type=int> 
        Specify file_newick_tree_png. Default is 0.

    --file_table_of_nodes_tsv <type=int> 
        Specify file_table_of_nodes_tsv. Default is 1.

    --file_table_of_branches_tsv <type=int> 
        Specify file_table_of_branches_tsv. Default is 1.

    --file_log_likelihood_tsv <type=int> 
        Specify file_log_likelihood_tsv. Default is 1.

    --file_table_of_attributes_tsv <type=int> 
        Specify file_table_of_attributes_tsv. Default is 1.

    --file_phylogenetic_tree_nwk <type=int> 
        Specify file_phylogenetic_tree_nwk. Default is 1.

    --e_mail <type=str> 
        Specify e_mail (technical parameter, do not change).

    --is_do_not_use_e_mail <type=int> 
        Specify is_do_not_use_e_mail (technical parameter, do not change).
```

### Citing
If you use the GLOOME web server for your research, please make sure to cite the following publication:

Cohen, O., Ashkenazy, H., Belinky, F., Huchon, D., and Pupko, T. 2010. GLOOME: gain loss mapping engine. 
Bioinformatics 26(22):2914-2915. [[pdf](https://www.tau.ac.il/~talp/publications/GLOOME.pdf)] [[abs](https://academic.oup.com/bioinformatics/article/26/22/2914/228050)]

### BibTeX
Proper citation helps support the development and maintenance of this tool.

### Contact
For questions, issues, or contributions, please open an issue on the repository or contact the maintainer directly.

