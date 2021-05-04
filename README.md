# Conformational analysis of protein structural ensembles

## Introduction
Intrinsically disordered proteins (IDPs) lack a fixed 3D structure and instead exhibit extreme
conformational dynamics in the free state. Similar to the unfolded state, IDPs and intrinsically
disordered regions (IDRs) must be described as ensembles of heterogeneous, rapidly interconverting
conformations. Conformational ensembles are representative sets of conformers reflecting on the
structural dynamics of IDPs sampling the space. Ensemble modeling usually relies on experimental data.
These measurements are then used to define local or nonlocal structural constraints for the
computational modeling of the conformational ensemble. Solving structural ensembles, however, is
fraught with uncertainties, because the number of degrees of freedom is inherently much larger than the
number of experimentally determined structural restraints. We don’t yet know how to select the ‘best’
ensemble from multiple alternatives, neither can we be sure if an actual ensemble is a faithful
representation of the real physical state of the IDP/IDR, nor is only a reasonable fit to experiment
observations. To help address these issues, solved IDP/IDR ensembles are collected and made
available in the dedicated Protein Ensemble Database (PED, [1][https://proteinensemble.org/]).

## Comparison of alternative ensembles
Structural comparisons rely on quantitative similarity measures. The most common measure is the root
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