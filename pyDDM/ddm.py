import numpy as np
import pandas as pd
from bioviper import msa, pdb
from scipy.cluster import hierarchy
import seaborn as sns

def DDM(pdb1, chain1, pdb2, chain2, fasta_path,
        helix_starts=None, helix_ends=None,
        helix_only=False, helix_binned=False, pdb_dir = './',
        first_index=1):

    '''Calculate a distance-difference matrix (DDM) to compare the conformations of two structures.

        Parameters
        ----------------
        pdb1: the first pdb structure. Should match the filename exactly/*
        chain1: the chain number to use from the first PDB file
        pdb2: the second pdb structure
        chain2: the chain number to use from the second PDB file
        fasta_path: a path to the raw fasta sequence

        To calculate a *helix-only* DDM:
            helix_only: whether to filter out non-helical residues to focus on the structural core and ignore loops.
                        defaults to FALSE; if TRUE, you must also pass:
            helix_starts: the start value of each helix, using the numbering in the fasta file
            helix_ends: the end value of each helix, using the numbering in the fasta file

        To calculate a helix-binned DDM, pass the above and also set helix_binned=True.

        Other optional arguments:
            pdb_dir: the directory where your PDB files are. if you don't anything, assumes the current directory
            first_index: when you pass numberings for helix_starts, helix_ends, etc., whether the sequence is 1-indexed or 0-indexed
                        defaults to 1, as is conventional in biochemistry; since python is 0-indexed, we therefore have to subtract
                        1 from from the positions you put in.

        Returns
        ------------
        ddm: a distance-difference matrix as a numpy array
        (if not helix_binned==True): mapping_df, a pandas dataframe for mapping the new DDM back onto the positions of the input structure.
            The index of this matches the index of rows/columns in the DDM, and the columns of the dataframe are:
                pdb1_pos: the residue # in the first PDB
                pdb2_pos: the residue # in the second PDB (this should be the same for structures from the same protein!)
                pdb1_index: of the residues with coordinates, which residue it was in the PDB
                    (e.g. what position it would be in a distance matrix calculated from that PDB)
                pdb2_index: same for the 2nd pdb
                pdb1_aa: the amino acid corresponding to that residue in structure 1
                pdb2_aa: the amino acid corresponding to that residue in structure 2
        (if you set helix_only==True): helix_breaks, where on the new DDM the breaks between helices are (e.g. if you want to draw lines there)
    '''

    # Must be helix-only if helix binned
    if helix_binned:
        helix_only = True

    # Read in full fasta sequence
    seq = msa.readAlignment(fasta_path)

    # In case didn't pass an extension to either pdb code, add it
    if ".pdb" not in pdb1:
        pdb1 = pdb1+'.pdb'
    if ".pdb" not in pdb2:
        pdb2 = pdb2 + '.pdb'

    # Read in the two structures
    structA = msa.readPDB(pdb_dir + pdb1, name=pdb1.split('.')[0], chains=chain1)
    structB = msa.readPDB(pdb_dir + pdb2, name=pdb2.split('.')[0], chains=chain2)

    # FOR A HELIX-ONLY DDM: we want to redefine the sequence to include only helical positions
    if helix_only:

        # Initialize lists for helical positions and the new locations of the helix breaks
        helical_pos = []; helix_breaks = []

        # Loop through each helix
        n = 0  # keep track of the new indexing of how many positions we've added (for helix_breaks)
        for nhelix in range(len(helix_starts)):

            # Loop through each position in the helix
            for i in range(helix_starts[nhelix], helix_ends[nhelix]):

                # Add that position to the list
                helical_pos.append(i)
                n+=1

            # However many positions you've added at this point, you'll want to introduce a new "helix break"
            helix_breaks.append(n)

        helical_pos = np.array(helical_pos) - first_index   # subtract off first_index for zero-indexing reasons
        helix_breaks = np.array(helix_breaks)

        seq = seq[:,helical_pos]  # redefine the sequence as only

    # Align each sequence to the reference sequence
    seq.attach_structure(structA, 0)
    seq.attach_structure(structB, 0)

    # Identify the overlapping positions to use for the DDM
    s1, s2 = seq[0].structures.values()
    s1_overlap = s1._pdbpos[np.isin(s1._alipos, s2._alipos)]
    s2_overlap = s2._pdbpos[np.isin(s2._alipos, s1._alipos)]

    # calculate the ddm over the selected overlapping positions
    ddm = s2.distance_matrix()[s2_overlap,:][:,s2_overlap] - s1.distance_matrix()[s1_overlap,:][:,s1_overlap]
    
    # generate the mapping_df for interpretability
    mapping_df = pd.DataFrame(np.array([s1.residue_ids[s1_overlap], s2.residue_ids[s2_overlap],
                              s1_overlap, s2_overlap, 
                              np.array(list(s1.sequence))[s1_overlap], 
                              np.array(list(s2.sequence))[s2_overlap]]).T,
                              columns=("pdb1_pos", "pdb2_pos", "pdb1_index", "pdb2_index", "pdb1_aa", "pdb2_aa"))

    if helix_binned:

        # Need to add 0 to the start
        helix_breaks_0 = np.concatenate([[0], helix_breaks])
        hb_ddm = helix_bin(ho_ddm, np.concatenate([[0], helix_breaks[:-1]]), helix_breaks)

        return hb_ddm

    elif helix_only:
        return ddm, mapping_df, helix_breaks

    else:
        return ddm, mapping_df


def helix_bin(X, helix_starts, helix_ends,shift=0, norm=True):

    '''
    Bin a matrix (typically a DDM) based on helix assignments. By default takes the Frobenius
    norm (squared sum of all entries) but can also take the mean (set norm=False).
    '''

    num_helices = len(helix_starts)
    helix_binned_ddm = np.zeros((num_helices,num_helices))

    for n1 in range(num_helices):
        for n2 in range(num_helices):

            A = X[(helix_starts[n1]-shift):(helix_ends[n1]-shift),:][:,(helix_starts[n2]-shift):(helix_ends[n2]-shift)]

            # Whether to take the Frobenius norm
            if norm:
                helix_binned_ddm[(n1,n2)] = np.sqrt(np.mean(A**2))

            # Otherwise take a simple sum
            else:
                helix_binned_ddm[(n1,n2)] = np.mean(A)

    return helix_binned_ddm

def align_structures(pdb1, pdb2, fasta_path, chain1="A", chain2="A", pdb_dir="./"):

    '''
    This function is just the structure aligning part of the DDM() function above - it is otherwise identical
    '''

    # Read in full fasta sequence
    seq = msa.readAlignment(fasta_path)

    # In case didn't pass an extension to either pdb code, add it
    if ".pdb" not in pdb1:
        pdb1 = pdb1+'.pdb'
    if ".pdb" not in pdb2:
        pdb2 = pdb2 + '.pdb'

    # Read in the two structures
    structA = msa.readPDB(pdb_dir + pdb1, name=pdb1.split('.')[0], chains=chain1)
    structB = msa.readPDB(pdb_dir + pdb2, name=pdb2.split('.')[0], chains=chain2)

    # Align each sequence to the reference sequence
    seq.attach_structure(structA, 0)
    seq.attach_structure(structB, 0)

    return seq

def cluster_helices(D, helix_names=None, method='average', plot=True,
                       dendrogram_ratio=0.3, colormap='Greys', figsize=(5,5)):

    '''
    Treating a DDM as a distance matrix, perform hierarchical clustering on the helices.
    Optionally (if plot=True), visualize as a seaborn clustermap.
    '''

    L = D.shape[0]

    # Conver the DDM into a "proper" distance matrix
    for i in range(L):
        for j in range(i,L):

            # We need to make the diagonals into zeros
            if i==j:
                D[(i,j)] = 0

            # Fix numerical errors that make it not technically symmetric
            elif D[(i,j)] != D[(j,i)]:
                D[(i,j)] = D[(j,i)]

    # generates the "linkage" object for hierarchical clustering
    linkage = hierarchy.linkage(squareform(D), method=method, optimal_ordering=True)

    if plot:
        if type(helix_names) == type(None):
            print("Must pass helix names to plot!")

        else:
            Ddf = pd.DataFrame(D, columns=helix_names)
            Ddf.index = helix_names
            sns.clustermap(Ddf, row_linkage=linkage, col_linkage=linkage, cmap=colormap, figsize=figsize,
              dendrogram_ratio=dendrogram_ratio)

    return linkage

def runPCA(hbDDMs, n_components, sparse=True, alpha=2, nonnegative=True, nonneg_method="rezero", positive_components=True,
          spca_method="lars"):
    
    '''Run (sparse) principal component analysis on a set of hbDDMs.
    
        DDMs: a list or array of all DDMs
        n_components: how many components to use.
        Optional:
            sparse: whether to perform sparse pca (default TRUE, sparse PCA is better for this)
            nonnegative: whether to perform a correction to the PCA projections to make them all nonnegative
                and set the "zero" in PCA space to represent the identity DDM (zero matrix). Default TRUE.
            nonneg_method: if nonnegative==TRUE, two options for method:
                nonneg_method = rezero subtracts off the projection of the identity DDM from all other matrices
                nonneg_method = nnmf used nonnegative matrix factorization with a fixed H to perform nonnegative
                    linear regression on the components to generated nonnegative weights for each uncentered DDM
                    from the components. In practice, this should give nearly the same answer as rezero.
            '''
    
    hbDDMs = np.array(hbDDMs)
    N_samples = hbDDMs.shape[0]
    N_helices = hbDDMs.shape[1]
    
    all_hb_ddms_flattened = hbDDMs.reshape((N_samples, N_helices**2))
    
    if sparse:
        pca = SparsePCA(n_components = n_components, method=spca_method, alpha=alpha)
    else:
        pca = PCA(n_components=n_components)
        
    pca_proj = pca.fit_transform(all_hb_ddms_flattened)
    
    if positive_components:
        pos_arr = (np.mean(pca.components_, axis=1) > 0).astype('int') * 2 - 1
    else:
        pos_arr = np.ones((n_components))
        
    components = pca.components_ * pos_arr[:, None]
    
    if nonnegative:
        if nonneg_method == "rezero":
            projections = pca_proj * pos_arr - pca.transform([np.zeros((N_helices**2))])[0] * pos_arr
    else:
        projections = pca_proj * pos_arr
        
    return components, projections
