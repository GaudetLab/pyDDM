# pyDDM
Code to build and analyze distance difference matrices (DDMs) to compare protein conformations. Also contains all code associated with the paper "They all rock: A systematic comparison of conformational movements in LeuT-fold transporters" by Jacob A. Licht, Samuel P. Berry, Michael A. Gutierrez and Rachelle Gaudet.

To quickly install this code, you can use pip:

`pip install git+https://github.com/GaudetLab/pyDDM.git`

and then import as

`import pyDDM as ddm`

For example usage as well as reproduction code for the manuscript, see the attached jupyter notebooks:

`Building DDMs for the APC superfamily.ipynb` shows how to construct a series of DDMs, ho-DDMs and hb-DDMs for a database of proteins in the LeuT-fold superfamily

`Clustermaps.ipynb` shows how to cluster the helices in these hb-DDMs to look for rigid bodies

`PCA on the DDMs.ipynb` shows how to run sparse PCA on the DDMs, along with cross-validation and reproduction of the plots included in the paper. This notebook would also be the place to start for anyone interested in exploring how varying the parameters in this analysis change the downstream conclusions (some of this is explored in the notebook as well!)

As a note, if you would like to make binned DDMs, you will need to determine and input a list of starts and ends for the segments of interest for each protein. In our paper those segments are helices, and we anticipate this will be the case for others studying membrane proteins, but the `helix_bin()` function (despite its name) could equally well bin the DDM based on other attributes, such as by domain for very large multi-domain proteins.

If you use this code, please consider citing:

Licht, J. A., Berry, S. P., Gutierrez, M. A., & Gaudet, R. (2024). They all rock: A systematic comparison of conformational movements in LeuT-fold transporters. bioRxiv, 2024-01.


