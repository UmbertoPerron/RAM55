# RAM55

This repository contains additional code and supplementary files from:

Perron et al., "Modelling structural constraints on protein evolution via side-chain conformational states"
https://doi.org/10.1101/530634

# MANUSCRIPT ABSTRACT
Few models of sequence evolution incorporate parameters describing protein structure, despite its high conservation, essential functional role and increasing availability. 
We present a structurally-aware empirical substitution model for amino acid sequence evolution in which proteins are expressed using an expanded alphabet that relays both amino acid identity 
and structural information. Each character specifies an amino acid as well information about the rotamer configuration of its side-chain: the discrete geometric pattern of permitted side-chain 
atomic positions, as defined by the dihedral angles between covalently linked atoms. By assigning rotamer states in 251,194 protein structures and identifying 4,508,390 substitutions between 
closely related sequences, we generate a 55-state ''Dayhoff-like'' model that shows that the evolutionary properties of amino acids depend strongly upon side-chain geometry. 
The model performs as well as or better than traditional 20-state models for divergence time estimation, tree inference and ancestral state reconstruction. 
We conclude that not only is rotamer configuration a valuable source of information for phylogenetic studies, but that modelling the concomitant evolution of sequence and structure may have 
important implications for understanding protein folding and function.

# Generated canonical sequence alignment
aligned_canonical_input.fasta

# Generated rotaequence alignment (gapped), uses 1-character symbols for rotamer states 
aligned_polypeptide_1-char.rotasequences.fasta

# Generated alignment of AA polypeptide sequences (one per chain model) from PDBe structures
aligned_polypeptide_input.fasta

# Ancestral rotamer or amino acid state reconstruction using the Joint reconstruction method,
# modified for our expanded state-set. 
# (Pupko, T., Pe’er, I., Shamir, R., and Graur, D. 2000.  A fast algorithm for joint reconstruction of ancestral amino acid sequences. Mol. Biol. Evol. , 17(6): 890–896)
ancestral_joint_reconstruct.py	

# Marginal rotamer or amino acid state reconstruction using the marginal reconstruction method, modified for our expanded state-set. 
# (Yang,  Z.,  Kumar,  S.,  and  Nei,  M.  1995.   A  new  method  of  inference  of  ancestral  nucleotide  and  amino  acid  sequences. Genetics, 141(4): 1641–1650.)
ancestral_marginal_reconstruct.py 

# Example canonical Uniprot sequences (un-aligned)
canonical_input.fasta

# A jupyter notebook example showing how to obtain a rotasequence alignment for a set of uniprot IDs
PDBe_API_rotamer-assignment_example.ipynb

# Maps Penultimate Library states (the ones PDBe provides) to Dunbrack Library states (the ones RAM55 uses)
Penultimate_2_Dunbrack.table

# Example AA polypepide sequences (un-gapped)
polypeptide_input.fasta

# Random tree generation script
tree_gen.py

# Substitution simulation along tree branches, analogous to Method 1 by Fletcher et al., modified for our expanded state-set. 
# (Fletcher, W. and Yang, Z. 2009.  Indelible: A flexible simulator of biological sequence evolution. Molecular Biology and Evolution, 26(8): 1879–1888)
tree_sim.py	

# All supplementary files from the Perron et al. manuscript
Supplementary_Files