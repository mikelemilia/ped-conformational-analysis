# Conformational analysis of protein structural ensembles

### Authors

[Milia Mikele](https://github.com/mikelemilia) -
[Pezzutti Giulia](https://github.com/giuliapezzutti) -
[Vettor Federica](https://github.com/FeVe98)

Structural Bioinformatics Project

Department of Mathematics - University of Padua, Italy

## Table of contents

- [Project specification](#project-specification)
    + [Task 1](#task-1)
    + [Task 2](#task-2)
- [Project structure](#project-structure)
    + [Features](#features)
        * [Extraction policies](#extraction-policies)
    + [Classes](#classes)
        * [Menu](#menu-class)
        * [Model features](#model-features-class)
        * [PED features](#ped-features-class)
- [Getting started](#getting-started)
    + [Run with custom dataset](#run-with-custom-dataset)

## Project specification

Intrinsically disordered proteins (IDPs) lack a fixed 3D structure and instead exhibit extreme conformational dynamics
in the free state. Similar to the unfolded state, IDPs and intrinsically disordered regions (IDRs) must be described as
ensembles of heterogeneous, rapidly interconverting conformations. Structural comparisons rely on quantitative
similarity measures. The most common measure is the root-mean-square deviation (RMSD) of the atomic positions between
two structures, which is minimized upon rigid-body superimposition of these structures. But, the RMSD is often not very
informative because it averages out differences across regions of the structures with varying similarity levels.
Characterizing and comparing IDP ensembles is therefore particularly challenging. First, their extreme conformational
heterogeneity makes it difficult to evaluate the degree of global similarity between two ensembles by any measure, let
alone by RMSD-based metrics. Second, the function of disordered proteins is often mediated by short, sequentially
contiguous binding motifs adopting locally relevant conformations. The latter are interconnected through more
structurally variable linkers that determine the relative overall configuration of these important motifs. Therefore,
the similarity of the IDP and IDR ensembles must be evaluated at both the local and global levels in a statistically
meaningful approach.

### Task 1

Relationships within an ensemble can be identified considering the structural features of single conformations: for this
purpose, from each conformation, a set of features must be extracted. The conformations should then be clusterized to
find the most significant ones, and these will be represented in a graph and in a Pymol image.

Features analyzed for a single conformation:

1. Radius of gyration of the structure.
2. Relative accessible surface area (ASA) for each residue.
3. Secondary structure (SS) for each residue, mapped into four classes based on phi/psi angles:
   α-helix, left-handed helix, β-strand, polyproline I and II.
4. Residue distance matrix considering Cα atoms.

### Task 2

Relationships between different ensembles of the same protein will be identified considering ensemble features
calculated from the output of the first software component. It is necessary to identify a measure (global score) to
quantify global differences between ensembles pairs - evaluated thanks to dendrogram and heatmap - and another measure (
local score) to identify low/high variance positions along the sequence - evaluated with a residues-scores plot.

Features analyzed for ensembles (multiple conformations)

1. Radius of gyration for each conformation in the ensemble.
2. Secondary structure entropy for each position across ensemble conformations.
3. Median solvent accessibility for each position across ensemble conformations.
4. Median RMSD for each position across ensemble conformations.
5. Median distance of each pair of equivalent positions across ensemble conformations.
6. Standard deviation of the distance of each pair of equivalent positions across ensemble conformations.

## Project structure

```
.
├── data/                 
│   ├── model-features/           # first task features data
│   ├── ped-features/             # second task features data
│   └── ...                       # .pdb / .ent files
├── output/               
│   ├── model-features/           # first task output
│   └── ped-features/             # second task output
├── src/                  
│   └── main.py                   
│   └── *.py                         
├── README.md                     # project overview
├── Makefile                      # makefile
└── requirements.txt              # required packages
```

The project is structured in three main folders: `data`, `output` and `src`.

Specifically, `data` contains all the input files related to the PED that we intend to analyze in .pdb (or .ent) format.
In addition, there are two sub-folders, `model-features` and `ped-features` which are used to better organize the
partial results (features file) produced by both project task.

The `output`, for coherence, is organized in the same way, in fact inside it there are two sub-folders  `model-features`
and `ped-features`, in which are respectively saved the graphical results produced by first and second task of the
project.

The `src` folder, instead, contains all the implemented code. A detailed analysis of the code structure is covered in
the following sections.

In addition to the main folders, the files `Makefile` and `requirements.txt` are made available, both of which can be
useful to install the minimum requirements to allow proper execution.

### Features extraction policies

As stated in the subsections of [Project specification](#project-specification), both tasks require the extraction of a
set of features. In order to make the code execution smoother, if the features files for the PED under analysis are already present 
inside `data/model-features`
or `data/ped-features` folders (respectively for the first or second task, calculated during a previous code run), 
every feature will be loaded from those files, instead of being recomputed from scratch.

To make the output files more understandable, we decided to map each set of features extracted during the analysis in a 
single row, applying a one to one policy. 
Clearly, this decision required some additions to each row:

- in task 1, it was necessary to insert both the model ID and the number of residues
- in task 2, it was necessary to insert the PED ID

The extraction functions `extract_vectors_model_feature(...)` and `extract_vectors_ped_feature(...)`, allow to obtain
different parts of features file depending on the needs. Specifically, once passed the data frame containing the
feature matrix, it is possible to obtain:

- all rows or a specific subset of them, containing a certain feature (i.e., RG, ASA, SS, etc ...)
- the interval extremes for a certain features (i.e., RG, ASA, SS, etc ...)
- all the feature intervals as slices

### Classes

To make it easier to read, the code has been divided into classes since the structure of the project, from a logical
point of view, suggested a modular organization. Specifically, three distinct classes have been implemented: `Menu`,
`ModelFeatures` and `PEDFeatures`. In subsequent sections, each class will be analyzed in detail.

#### Menu class

This class takes care of generating a useful menu for the user with the allowed options, which are:

1. Analyze features of models inside a PED
2. Compare inter and intra features of a PED
3. Exit

In order to select an option it is sufficient to report its number. To undo any wrong selection, it is enough to
digit `Q` (QUIT) or `B` (BACK) or `U` (UNDO).

Choosing to execute the first task or the second task, you are asked to insert the path to the folder containing the PED of interest. 
Subsequently, the PED ID list of the files found inside the provided folder is shown and you are asked to insert the 
index or the ID of the one of interest. As soon as the PED ID is correctly supplied, an
instance of the `ModelFeatures` class, if the first task has been chosen, or instance of the `PEDFeatures` class, if the 
case of the second one, is initialized. 

> **N.B. The folder path that can be inserted when requested can be absolute or relative. Note that whenever a relative
> path is provided, it must be referred to the working directory that has been set in the program as the project folder.**

An example of Menu usage for the `data` folder and PED00020 is here provided: 

```
== == == == == == == == == == == == == == == == == == ==
1.    Task 1 : Analyze features of models inside a PED
2.    Task 2 : Compare inter and intra features of a PED
3.    Exit
== == == == == == == == == == == == == == == == == == ==
You can reverse any wrong selection with
  - Q (QUIT)
  - B (BACK)
  - U (UNDO)
== == == == == == == == == == == == == == == == == == ==

Select what do you want to do: 1

Insert the folder path containing the PED you want to analyze: data

Which PED do you want to analyze: 

0 - PED00020

Your choice: 0
Your choice is PED00020

Analyzing PED00020e001...
```

#### Model features class

This class, named `ModelFeatures`, takes care of the whole first task workflow. In detail, once an instance is
initialized, the workflow is the following:

- firstly, class method `choice_maker` is called and the presence of feature files is checked, if they're found the
  class method `extract` is subsequently called, otherwise the pair of class method `compute` and `save` take care of
  computing each requested feature of the model and storing it inside the `data/model-features` folder in a specific
  file called `PEDxxxxxexxx_features.csv`. Once `choice_maker` concludes its execution, we're sure that features of the
  PED of interest are computed or at least loaded.

- then class method `compute_clustering` is called and **K-medoids** clustering is performed on the PED features. Since
  computed features were not comparable, an ad-hoc metric taking care of each one of them was implemented, and can be
  found inside the `metrics` class method. Thanks to the centroids (or better medoids) generated, through the class
  method `generate_graph` a fully connected graph of the representative conformations is given in output and can be
  found inside the `output/model-features` folder. In order to draw the graph, library `networkx` was used.

- finally, class method `generate_pymol_image` is called. This method accepts the computed graph and after computing
  each residue variability, through `pymol_metric` method, returns a single PYMOL image in which each residue is colored
  according to its variability.

#### PED features class

This class, named `PedFeatures`, takes care of the whole second task workflow. In detail, once an instance is
initialized, the workflow is the following:

- firstly, class method `choice_maker` is called. Initially it loads (or if not present, compute) the models features of
  corresponding to the first task with `load_models_files`: as a matter of fact, if the features files have not been generated 
  within the first task, the program will automatically generate them (but will not perform the first task analysis steps). 
  The function then checks whether ped features files should be generated or
  loaded (as happens in the first task). The only difference is the folder within each feature file is saved that is `data/ped-features`
  and the name of the file that now is simply `PEDxxxxx_features.csv`

- class method `global_dendrogram` is called in order to plot the weighted distance between a pair of ensembles with
  respect of a global metric that can be found inside `global_metric` method

- then the same metric is applied inside the class method `global_heatmap` in order to obtain the pairwise difference of
  ensembles

- finally, class method `local_metric` is called in order to plot how variable a residue is with respect to all the
  other ensambles

## Getting started

Before all, make sure to have [DSSP](https://ssbio.readthedocs.io/en/latest/instructions/dssp.html) installed on your
machine:

```shell
sudo apt-get install dssp

# symlink its name (dssp) from mkdssp, since installs itself as mkdssp
sudo ln -s /usr/bin/mkdssp /usr/bin/dssp
```

In order to get access to the code and replicate our results, follow these steps:

```shell
# Clone project repository
git clone https://github.com/mikelemilia/ped-conformational-analysis.git project-name

cd project-name                   # move inside project folder

make init                         # install project requirements (only if needed)

cd src                            # move inside src folder
python main.py                    # execute w/ default command line parameters
```

[comment]: <> (### Run with custom dataset)
[comment]: <> (If the path to folder containing PED files is not specified, `data` folder of this project is considered the default)
[comment]: <> (input folder. So, in order to test our code with different PED, you must provide a custom path &#40;pointing to their)
[comment]: <> (location&#41; as parameter or simply insert the PED .pdb &#40;or .ent&#41; file inside the `data` folder.)
[comment]: <> (To provide a custom path, you can use `-p` flag, as in the following example:)
[comment]: <> (```shell)
[comment]: <> (python main.py -p custom-path     # execute w/ custom path)
[comment]: <> (```)
