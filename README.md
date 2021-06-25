# Conformational analysis of protein structural ensembles

## Authors

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
        * [Menu](#menu)
        * [Model features](#model-features)
        * [PED features](#model-features)
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
quantify global differencesbetween ensembles pairs - evaluated thanks to dendrograms and heatmap - and another measure (
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
partial results (files of extracted features) produced by both project task.

The `output` for coherence is organized in the same way, in fact inside it there are two sub-folders  `model-features`
and `ped-features`, in which are respectively saved the graphical results produced by produced by first and second task
of the project.

As for the `src` folder, instead, it contains all the implemented code. A detailed analysis of the code structure is
covered in the following sections.

In addition to the main folders, the files `Makefile` and `requirements.txt` are made available, both of which can be
useful to install the minimum requirements to allow proper execution.

### Features

As stated in the subsections of [Project specification](#project-specification), both tasks require the extraction of a
set of features. To make the output files more understandable, regardless of the task, we decided to map them with this
policy: one line one feature set.

Clearly, this decision required some additions to each feature set:

- regarding the first task, it was necessary to insert both model ID and number of residues
- regarding the second task, it was necessary to insert only PED ID

#### Extraction policies

In order to make the code execution smoother, we thought to implement specialized functions for feature extraction.
Since we didn't think it was useful to recalculate every time the features produced and stored inside
data `model-features`
or `ped-features` folder (depending on the task) during previous code run.

The extraction functions `extract_vectors_model_feature(...)` and `extract_vectors_ped_feature(...)`, allow to obtain
different parts of the features file depending on the need. Specifically, once passed the data frame containing the
feature matrix, it is possible to obtain:

- all rows, or a specific rows subset from a certain column 'till the end
- all rows, or a specific rows subset between a certain interval
- all the rows containing a certain feature (i.e., RG, ASA, SS, etc ...)
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

Choosing to execute the first task, you are asked for the PED of interest for which you wish to generate the features,
which must be provided in the `PEDxxxxx` format (where `x` is a digit). As soon as the PED is correctly supplied, an
instance of the `ModelFeatures` class which will take care of the correct operation of the first task in all its aspects
is initialized.

Choosing to execute the second task, you are asked again for the PED of interest for which you wish to generate the
features, which must be provided in the same format as before. As soon as the PED is correctly supplied, an instance of
the `PEDFeatures` class which will take care of the correct operation of the second task in all its aspects is
initialized.

#### Model features class

Task 1 is implemented in ModelFeatures class. An example of its usage can be found in 'first_task' function inside Menu
class. Initially, if the features files are not already present in the correct folder, they are computed and saved,
otherwise they are loaded. Subsequently, functions to perform clustering and the correspondent graph and Pymol image are
implemented.

**- Task 1 : Abbiamo fatto clustering con K-medoids e una metrica ad hoc - generazione grafo con networkx e pymol image
con la variabilità dei residui calcolata con un'altra metrica fatta da noi.**

#### Ped features class

Task 2 is implemented in PedFeatures class. An example of its usage can be found in 'second_task' function inside Menu
class. Even in this case, if the features file is already saved in the folder, it is loaded, otherwise the entire
computation is performed and saved. The functions to build heatmaps and dendrograms (according to the defined global
metric) and the one to build the plot of local variability
(according to the local metric) are then implemented.

**- Task 2: Choice maker e feature extractor: implementazione di dendrogram, heatmap e variabilità locale con metriche
ad hoc.**

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

### Run with custom dataset

If the path to folder containing PED files is not specified, `data` folder of this project is considered the default
input folder. So, in order to test our code with different PED, you must provide a custom path (pointing to their
location) as parameter or simply insert the PED .pdb (or .ent) file inside the `data` folder.

To provide a custom path, you can use `-p` flag, as in the following example:

```shell
python main.py -p custom-path     # execute w/ custom path
```