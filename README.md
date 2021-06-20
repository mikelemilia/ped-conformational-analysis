# Conformational analysis of protein structural ensembles

## Authors
[Milia Mikele](https://github.com/mikelemilia) - 
[Pezzutti Giulia](https://github.com/giuliapezzutti) -
[Vettor Federica](https://github.com/FeVe98)

Structural Bioinformatics project

Department of Mathematics - University of Padua, Italy

## Project specification

Intrinsically disordered proteins (IDPs) lack a fixed 3D structure and instead exhibit extreme
conformational dynamics in the free state. Similar to the unfolded state, IDPs and intrinsically
disordered regions (IDRs) must be described as ensembles of heterogeneous, rapidly interconverting
conformations. Structural comparisons rely on quantitative similarity measures. The most common measure is the root
mean square deviation (RMSD) of the atomic positions between two structures, which is minimized
upon rigid-body superimposition of these structures. But, the RMSD is often not very informative
because it averages out differences across regions of the structures with varying similarity levels.
Characterizing and comparing IDP ensembles is therefore particularly challenging. First, their extreme
conformational heterogeneity makes it difficult to evaluate the degree of global similarity between
two ensembles by any measure, let alone by RMSD-based metrics. Second, the function of disordered
proteins is often mediated by short, sequentially contiguous binding motifs adopting locally
relevant conformations. The latter are interconnected through more structurally variable linkers that
determine the relative overall configuration of these important motifs. Therefore, the similarity of the IDP
and IDR ensembles must be evaluated at both the local and global levels in a statistically meaningful
approach.

#### Task 1

Relationships within an ensemble can be identified considering the structural features of single
conformations: for this purpose, from each conformation, a set of features must be extracted. The 
conformations should then be clusterized to find the most significant ones and these will be rappresented 
in a graph and in a Pymol image.

Features analyzed for a single conformation:
1. Radius of gyration of the structure.
2. Relative accessible surface area (ASA) for each residue. 
3. Secondary structure (SS) for each residue, mapped into four classes based on phi/psi angles:
α-helix, left-handed helix, β-strand, polyproline I and II. 
4. Residue distance matrix considering Cα atoms.

#### Task 2

Relationships between different ensembles of the same protein will be identified considering
ensemble features calculated from the output of the first software component. It is necessary to 
identify a measure (global score) to quantify global differencesbetween ensembles pairs - evaluated 
thanks to dendrograms and heatmap - and another measure (local score) to identify low/high variance positions
along the sequence - evaluated with a residues-scores plot.

Features analyzed for ensembles (multiple conformations)
1. Radius of gyration for each conformation in the ensemble. 
2. Secondary structure entropy for each position across ensemble conformations.
3. Median solvent accessibility for each position across ensemble conformations.
4. Median RMSD for each position across ensemble conformations.
5. Median distance of each pair of equivalent positions across ensemble conformations.
6. Standard deviation of the distance of each pair of equivalent positions across ensemble
conformations.

## Project Structure

All the script and classes implemented for this project can be found in 'scripts' folder. 
'Output' folder will contains all the output and plots generated during the execution.  

Menu class implements the different choices that can be taken from the inital menù. 

Task 1 is implemented in ModelFeatures class. An example of its usage can be found in 'first_task' function
inside Menu class. Initially, if the features files are not already present in the correct folder, they 
are computed and saved, otherwise they are loaded. Subsequently, functions to perform clustering and the 
correspondent graph and Pymol image are implemented.

Task 2 is implemented in PedFeatures class. An example of its usage can be found in 'second_task' function
inside Menu class. Even in this case, if the features file is already saved in the folder, it is loaded, 
otherwise the entire computation is performed and saved. The functions to build heatmaps and 
dendrograms (according to the defined global metric) and the one to build the plot of local variability 
(according to the local metric) are then implemented.

## How to run

*python main.py -p \<folder path>*

To execute it from command line, the path to folder already containg PED files of interest (.pdb or .ent formats are accepted) 
can be passed as input as reported above. If it is not specified, 'data' folder of this project is considered as input folder. 
Notice that inside the folder containing the data, for each task, a new folder will be created: 'model_features' will contain
files generated with the first task while 'ped_features' the ones of the second task. 

When the program starts, the user need to select the task to be performed from the menù (it is sufficient to report the task 
number: 1 for task 1, 2 for task 2 and 3 to exit). To subsequently undo the selection of the task, it is enough to digit 'Q'. 
The user is then asked to insert the PED ID of interest: it should be of the form PEDxxxxxx, where x corresponds 
to a digit. Not valid ID or ID of PEDs not present in the folder will be rejected!

