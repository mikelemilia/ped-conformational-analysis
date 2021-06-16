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

#### Task 2

Relationships between different ensembles of the same protein will be identified considering
ensemble features calculated from the output of the first software component. It is necessary to 
identify a measure (global score) to quantify global differencesbetween ensembles pairs - evaluated 
thanks to dendrograms and heatmap - and another measure (local score) to identify low/high variance positions
along the sequence - evaluated with a residues-scores plot.

## How to run

All the script and classes implemented for this project can be found in 'scripts' folder. 
'Output' folder will contains all the output and plots generated during the execution.  

Task 1 is implemented in ModelFeatures class. An example of its usage can be found in 'first_task' function
inside Menu class. 
Task 2 is implemented in PedFeatures class. An example of its usage can be found in 'second_task' function
inside Menu class. 
Menu class, instead, implements the different choices that can be taken from the inital menù. 

To execute the program from command line, the folder already containg the PED files of interest (.pdb or .ent formats are accepted) 
can be passed as input as follows. If the folder path is not specified, data folder of this project is considered as input folder. 
Notice that inside the folder containing the data, for each task, a new folder will be created: 'model features' will contain
files generated with the first task while 'ped features' files of the second one. 

python main.py -p \<folder path>

Once the user selects the task to be performed, he is asked to insert the PED ID of interest. Notice that it should
be of the form PEDxxxxxx, where x corresponds to a digit. Not valid ID will be rejected. 

