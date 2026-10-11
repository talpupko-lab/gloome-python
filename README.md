# GLOOME: Gain Loss Mapping Engine
A high-performance bioinformatics tool for analyzing gene gain and loss events, and the co-evolution of
presence-absence patterns, along a phylogeny.

GLOOME provides evolutionary analysis of presence and absence profiles (phyletic patterns). These patterns are assumed
to result from gain and loss dynamics along a phylogenetic tree. Examples of characters represented by phyletic patterns
include:

* Restriction sites
* Gene families
* Introns
* Indels

The primary purpose of GLOOME is to infer branch-specific and site-specific gain and loss events, and to analyze the
co-evolution of presence-absence patterns (CoPAP). For every branch and every site, GLOOME computes the exact posterior
probability of a gain (0 → 1) and of a loss (1 → 0) under a continuous-time two-state model with gamma-distributed
rate variation among sites.

### Features
* **Gain and loss inference**: exact posterior probabilities of gain and loss events per branch and per site
* **Ancestral reconstruction** of presence/absence at every internal node
* **Posterior rate estimation** per character (site)
* **CoPAP analysis** (co-evolution of presence-absence patterns), with p-values from a simulation-based null
  distribution and Benjamini-Hochberg q-values
* **Parameter optimization** of π1 (frequency of presence), the gamma α parameter and a branch-length scaling coefficient
* **Parsimony and homoplasy scores** per character (parsimony score, minimum steps, homoplasy, consistency index)
* **Tree rooting** by MAD, MVR, midpoint or outgroup
* **High-performance computing** powered by NumPy vectorization

### Project structure

```text
gloome
├── gloome
│   ├── services
│   │   ├── __init__.py
│   │   └── service_functions.py
│   ├── tree
│   │   ├── __init__.py
│   │   ├── node.py
│   │   ├── npencoder.py
│   │   └── tree.py
│   ├── __init__.py
│   ├── __main__.py
│   ├── config.py
│   ├── consts.py
│   └── utils.py
├── pyproject.toml
└── README.md
```

This manual covers input formats, command-line usage, the webserver, and how to interpret the results.

## Input

* **MSA** (`--msa_file`): phyletic patterns in FASTA format. Each sequence is a 0/1 string
  (`1` = presence, `0` = absence, `?` = missing data). All sequences must have the same length, and sequence names
  must match the leaf names of the tree.
* **Tree** (`--tree_file`): a phylogenetic tree in Newick format with branch lengths. Bootstrap values are removed
  automatically, and the tree is rooted with the selected `--rooting_method`.

Example:

```text
>A
11101010111000101011111011111
>B1
...
```

```text
(A:0.1,((E:0.1,(F1:0.1,F2:0.1):0.2):0.2,((D1:1.1,D2:0.1):0.12,((B1:0.1,B2:0.1):0.2,C:0.1):0.2):0.2):0.2);
```

## Installation

GLOOME is available on [PyPI](https://pypi.org/project/gloome/) and requires Python 3.9 or later:

```bash
pip install gloome
```

This installs all required dependencies and the `gloome` command.

## Command-line usage

```bash
gloome --msa_file patterns.fasta --tree_file tree.nwk --out_dir results
```

`python -m gloome` works as well.

Results are written to `<out_dir>/<process_id>/`. If `--process_id` is not given, one is generated.

#### Required parameters

```text
    --msa_file <type=str>
        Path to the MSA in FASTA format. Input must be a 0/1 binary matrix
        ('1' for presence, '0' for absence, and '?' for missing data).

    --tree_file <type=str>
        Path to the tree in Newick format.
```

#### Optional parameters

```text
    --out_dir <type=str>
        Output directory.

    --process_id <type=str>
        Process ID (name of the result subdirectory). Generated automatically if omitted.

    --with_internal_nodes <type=int>
        Write internal node names into the output Newick tree (0/1). Default is 1.

    --rooting_method <type=str>
        Tree rooting method. Possible options: ('mad', 'mvr', 'midpoint', 'outgroup').
        mad - Minimal Ancestor Deviation
        mvr - Minimum Variance Rooting
        midpoint - Midpoint Rooting
        outgroup - Outgroup Rooting
        Default is 'midpoint'.

    --leaf <type=str>
        Leaf to use as the outgroup when --rooting_method is 'outgroup'. Default is ''.
```

Model parameters:

```text
    --pi_1 <type=float>
        π1, the stationary frequency of presence (state '1'), between 0.001 and 0.999. Default is 0.5.

    --alpha <type=float>
        Shape parameter α of the gamma distribution of rates among sites, between 0.1 and 20. Default is 0.5.

    --categories_quantity <type=int>
        Number of discrete gamma rate categories, between 1 and 16. Default is 4.

    --coefficient_bl <type=float>
        Scaling coefficient applied to all branch lengths, between 0.1 and 10. Default is 1.0.

    --is_optimize_pi <type=int>
        Optimize π1 by maximum likelihood (0/1). Default is 1.

    --is_optimize_pi_average <type=int>
        Set π1 to the observed frequency of '1' in the MSA (0/1). Used only when --is_optimize_pi is 0.
        Default is 0.

    --is_optimize_alpha <type=int>
        Optimize α by maximum likelihood (0/1). Default is 1.

    --is_optimize_bl <type=int>
        Optimize the branch-length coefficient by maximum likelihood (0/1). Default is 1.
```

CoPAP parameters:

```text
    --is_do_not_use_copap <type=int>
        Skip the CoPAP analysis (0/1). Default is 0.

    --probability_lg <type=float>
        Minimum gain/loss probability for an event to count on a branch, between 0.01 and 0.99.
        Used to select the sites for the Pearson correlation table. Default is 0.5.

    --number_lg <type=int>
        Minimum number of such gain/loss events a site needs to be included in the Pearson
        correlation table, between 1 and 20. Default is 1.

    --number_datasets <type=int>
        Number of datasets simulated under the fitted model to build the null distribution
        for the coevolution p-values, between 1 and 1000. Default is 100.
```

Output files (each flag is 0/1 and turns one output file on or off; see [Output files](#output-files)):

```text
    --file_interactive_tree_html <type=int>                          Default is 0.
    --file_newick_tree_png <type=int>                                Default is 0.
    --file_table_of_coevolution_tsv <type=int>                       Default is 1.
    --file_simulated_datasets_fastas <type=int>                      Default is 1.
    --file_barplot_of_correlation_svg <type=int>                     Default is 1.
    --file_plot_distribution_of_correlation_svg <type=int>           Default is 1.
    --file_plot_correlation_by_rate_bin_svg <type=int>               Default is 1.
    --file_table_of_posterior_rates_tsv <type=int>                   Default is 1.
    --file_table_of_pearson_correlation_tsv <type=int>               Default is 1.
    --file_table_of_nodes_tsv <type=int>                             Default is 1.
    --file_branch_position_probabilities_tsv <type=int>              Default is 1.
    --file_table_of_branches_tsv <type=int>                          Default is 1.
    --file_log_likelihood_tsv <type=int>                             Default is 1.
    --file_table_of_attributes_tsv <type=int>                        Default is 1.
    --file_table_of_parsimony_and_homoplasy_scores_tsv <type=int>    Default is 1.
    --file_phylogenetic_tree_nwk <type=int>                          Default is 1.
```

## Output files

All tables are tab-separated. **Positions (`POS`) are 0-based**: `POS 0` is the first column of the MSA.
Branches are named after their child node, so "branch `N3`" means the branch leading from the parent of `N3` to `N3`.
Internal nodes are named `N1`, `N2`, … in the output tree.

The CoPAP files (coevolution, Pearson correlation, posterior rates, simulated datasets and correlation plots) are
created only when CoPAP is enabled and the MSA has more than one position. The posterior rates table needs more than
one rate category.

### Model and likelihood

| File                   | Contents                                                                                                                                                                                  |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `TreeAttributes.tsv`   | Final model parameters after optimization: π1, gamma α, number of rate categories, branch-length coefficient, the rate of each gamma category, the alphabet and the total log-likelihood. |
| `LogLikelihood.tsv`    | Log-likelihood (natural log) of each position, plus a `TOTAL` row with the log-likelihood of the whole dataset. Positions with very low values fit the model poorly.                      |
| `PosteriorRates.tsv`   | `POS`, `rate`: posterior mean evolutionary rate of each position. Rates above 1 mean the character changes faster than average, rates near 0 mean it is nearly conserved.                 |
| `PhylogeneticTree.nwk` | The rooted tree used in the analysis, with internal node names and branch lengths after scaling by the coefficient.                                                                       |

### Gain and loss

| File                              | Contents                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|-----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Branches.tsv`                    | One row per branch: `Parent node`, `Child node`, `Branch length`, then per-position lists: `Gain probability` P(0→1), `Loss probability` P(1→0), and `Branch probability vector`, which holds four probabilities per position in the order [absence maintained (0→0), loss (1→0), gain (0→1), presence maintained (1→1)]. The four values add up to 1.                                                                                                                                                                         |
| `BranchPositionProbabilities.tsv` | The same gain/loss probabilities in long format, one row per (position, branch, event): `G/L` (`gain` or `loss`), `POS`, `branch`, `branchLength`, `distance2root`, `distance2NearestOTU`, `numOfNodes2NearestOTU`, `probability`. Only rows with probability > 0.05 are written. Use this file to filter events, e.g. all gains with probability > 0.9, or to count the expected number of gains/losses per site or per branch by summing `probability`.                                                                      |
| `Nodes.tsv`                       | One row per node: `Name`, `Parent`, `Distance to parent`, `Children`, `Sequence` (reconstructed presence/absence pattern, most probable state at each position), `Probability coefficient` (posterior probability of that reconstructed state), `Ancestral comparison` (change from parent, see below), plus the technical likelihood vectors (`P matrix`, `Up vector`, `Down vector`, `Marginal vector`, `Vector of current conditional likelihood (subtree)`, `Vector of current conditional log-likelihood (ln, subtree)`). |

`Ancestral comparison` codes, per position, compare a node's reconstructed state with its parent's:

| Code | Meaning                     |
|------|-----------------------------|
| `A`  | absence maintained (0 → 0)  |
| `L`  | loss (1 → 0)                |
| `G`  | gain (0 → 1)                |
| `P`  | presence maintained (1 → 1) |

The codes come from the most probable states, so they give a single best reconstruction. For uncertain events, use
the probabilities in `Branches.tsv` or `BranchPositionProbabilities.tsv`.

### Co-evolution (CoPAP)

| File                            | Contents                                                                                                                                                                                                                                                                                                                               |
|---------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Coevolution.tsv`               | Main CoPAP result, one row per pair of positions: `POS1`, `POS2`, `r` (Pearson correlation of the two positions' gain/loss probability profiles across all branches), `rate-bin` (pair of gamma rate categories of the two positions), `p-value`, `q-value`, `direction`. Sorted by q-value, so the most significant pairs come first. |
| `PearsonCorrelation.tsv`        | `POS1`, `POS2`, `correlation`: raw Pearson correlation of gain/loss profiles for pairs of positions that have at least `--number_lg` events with probability ≥ `--probability_lg`. No significance test; use `Coevolution.tsv` for that.                                                                                               |
| `SimulatedDatasets.fastas`      | The `--number_datasets` MSAs simulated on the same tree under the fitted model (each position keeps its posterior rate), separated by `iterations = i` headers. They form the null distribution for the p-values.                                                                                                                      |
| `DistributionOfCorrelation.svg` | Density (KDE) of `r` over all tested pairs, with the α and β of a Beta distribution fitted to it.                                                                                                                                                                                                                                      |
| `BarplotOfCorrelation.svg`      | Histogram of `r` over all tested pairs, with the same Beta fit.                                                                                                                                                                                                                                                                        |
| `CorrelationByRateBin.svg`      | `r` of every pair, grouped by rate-bin. Shows whether strong correlations concentrate in particular rate classes.                                                                                                                                                                                                                      |

How to read `Coevolution.tsv`:

* Only positions that vary and have at least one gain or loss with probability > 0.5 are tested.
* `r` > 0, `direction = co-occurrence`: the two characters tend to be gained and lost on the same branches.
  `r` < 0, `direction = avoidance`: when one is gained or lost, the other tends not to be.
* `p-value` is empirical: the fraction of simulated pairs **in the same rate-bin** with |r| at least as large as the
  observed one. Comparing within a rate-bin corrects for the fact that fast-evolving characters correlate more by
  chance. The smallest possible p-value is 1 / (number of null pairs in that bin + 1), so more simulated datasets
  give finer p-values.
* `q-value` is the Benjamini-Hochberg false discovery rate. A common cutoff is q < 0.05.

### Parsimony

| File                              | Contents                                                                                                                                                                                                                                                                                                                                                                                       |
|-----------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `ParsimonyAndHomoplasyScores.tsv` | Per position: `Minimum Possible Number Of Steps (M)` (number of states observed minus 1), `Parsimony Score (S)` (minimum number of changes on this tree, Fitch algorithm), `Homoplasy Score (H = S - M)` (extra changes caused by convergence or reversal), `Consistency Index (CI = M/S)` (1 means no homoplasy). A `TOTAL` row sums M, S and H over all positions and gives the ensemble CI. |

### Other files

| File                   | Contents                                                                                                       |
|------------------------|----------------------------------------------------------------------------------------------------------------|
| `InteractiveTree.html` | Standalone interactive tree that opens in a browser (off by default).                                          |
| `VisualTree.png`       | Static picture of the tree (off by default).                                                                   |
| `*.log`                | Log of the run, with parameters and timings.                                                                   |
| `<process_id>.zip`     | Archive of all output files.                                                                                   |

## Webserver

GLOOME is available as a webserver at [gloome.tau.ac.il](https://gloome.tau.ac.il).

### Running an analysis

1. Paste a phyletic pattern MSA (FASTA) and a Newick tree into the text boxes, or upload files. **Load an example**
   fills in an example dataset, **Clear** resets the form.
2. Optionally adjust the model parameters, the CoPAP parameters, the rooting method (for outgroup rooting, choose a
   leaf from the list, which fills in from the tree) and the output files. These options match the
   [command-line parameters](#command-line-usage).
3. Optionally enter an e-mail address to get a notification when the job finishes.
4. Click **SUBMIT**. The page checks the job status every few seconds and shows the results when it finishes.
   Each run has its own results page at `/results/<process_id>`, which you can bookmark.

### Results page

The results page has four parts.

**Log-likelihood bar.** Shows the total log-likelihood of the tree. Click it to copy the value to the clipboard.

**Interactive tree.**

* Node colors: red = root, orange = internal node, green = leaf.
* **Hover over a node** to see a short summary of the node in a tooltip. The node grows and changes color.
* **Hover over a branch** to see a short summary of the branch in a tooltip. The branch turns dark red.
* **Click a node** to open the **Node information** table below the tree.
* **Click a branch** to open the **Branch information** table below the tree.
* Zoom with the mouse wheel or a two-finger trackpad pinch, and drag to move the tree.
* **show distance to parent** writes each branch length next to its node.
* **show radial tree** switches between a rectangular and a circular layout.

**Node information table** (one column per position; the last row, `Index`, gives the 0-based position):

| Row                       | Meaning                                                                                                                                             |
|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
| `Name`, `Node type`       | Node name and type (root, node or leaf).                                                                                                            |
| `Distance to parent`      | Length of the branch above the node.                                                                                                                |
| `Sequence`                | Presence/absence at each position: observed for leaves, reconstructed for internal nodes.                                                           |
| `Probability coefficient` | Posterior probability of the shown state, colored on a scale from red (low) to violet (high). Values near 0.5 mean the reconstruction is uncertain. |
| `Ancestral Comparison`    | Change from the parent at each position: `A` absence maintained, `L` loss, `G` gain, `P` presence maintained (see the legend next to the tree).     |

Missing data (`?`) in a leaf is shown as `? (x)`, where `x` is the state GLOOME reconstructed for that position.

**Branch information table** (one column per position):

| Row                         | Meaning                                                                  |
|-----------------------------|--------------------------------------------------------------------------|
| `Parent node`, `Child node` | The two ends of the branch.                                              |
| `Branch length`             | Length of the branch.                                                    |
| `Gain probability`          | Posterior probability of a gain (0 → 1) on this branch at each position. |
| `Loss probability`          | Posterior probability of a loss (1 → 0) on this branch at each position. |

Probabilities are colored from red (near 0) to violet (near 1), so likely events stand out when you scan a row.
These rows contain the same values as `Branches.tsv`.

**View/Download Results.** Lists every output file described in [Output files](#output-files).
**view** opens a file in the browser, **download** saves it. The archive (zip) contains all files.

### Loading previous results

Use **read JSON response file** to upload the `result.json` from an earlier run. The page then shows that run's tree,
tables and log-likelihood again without recomputing anything.

## Citing
If you use GLOOME for your research, please cite the following publication:

Cohen, O., Ashkenazy, H., Belinky, F., Huchon, D., and Pupko, T. 2010. GLOOME: gain loss mapping engine.
Bioinformatics 26(22):2914-2915. [pdf](https://www.tau.ac.il/~talp/publications/GLOOME.pdf) [abs](https://academic.oup.com/bioinformatics/article/26/22/2914/228050)

### BibTeX
Proper citation helps support the development and maintenance of this tool.

```bibtex
@article{cohen2010gloome,
  author = {Cohen, Ofir and Ashkenazy, Haim and Belinky, Frida and Huchon, Doroth{\'e}e and Pupko, Tal},
  title = {GLOOME: gain loss mapping engine},
  journal = {Bioinformatics},
  volume = {26},
  number = {22},
  pages = {2914--2915},
  year = {2010},
  publisher = {Oxford University Press}
}
```

## Contact
For questions, issues, or contributions, please open an issue on the repository or contact the maintainer directly.
