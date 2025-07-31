# AlphaFold, AlphaFold-Multimer, AlphaFoldPulldown

Sources:

[EMBL-EBI Tutorial](https://www.ebi.ac.uk/training/online/courses/alphafold/])

[VIB Tutorial](https://elearning.vib.be/courses/alphafold/)

## Confidence scores - AlphaFold

### pLDDT: Predicted Local Distance Difference Test

Per-residue measure of local confidence. It measures confidence in the local structure, estimating how well the prediction would agree with an experimental structure.

Visualization: color in structures or line plot

Scale: 0-100, the higher the better. +70 is good.

Low confidence, two reasons:

(1) naturally highly flexible or intrinsically disordered, like linkers for example

(2) not enough information to predict with confidence

### PEA: Predited Aligned Error 

Is a measure of how confident AlphaFold2 is in the relative position of two residues within the predicted structure. PAE is defined as the expected positional error at residue X, measured in Ångströms (Å), if the predicted and actual structures were aligned on residue Y. Measure of how confident AlphaFold2 is that the domains are well packed and that the **relative placement of the domains** in the predicted structure is correct.

Visualization: heatmap of aligned residues (PEA plot)

Scale: 0 - 30+, the lower the better.

## Confidence scores - AlphaFold-Multimer

### pTM: Predicted Template Modelling

Is an integrated measure of how well AlphaFold-Multimer has predicted the overall structure of the complex. It is the predicted TM score for a superposition between the predicted structure and the hypothetical true structure.

**pTM score should be interpreted cautiously. For example, imagine a situation where one of the interacting proteins is larger, and its structure is predicted correctly, while a smaller partner protein is predicted incorrectly. The resulting pTM score of the complex could be dominated by the larger protein and show a pTM score above 0.5.**

Scale: 0 to 1, above 0.5 is good

### ipTM: Interface Predicted Template Modelling

Measures the accuracy of the predicted relative positions of the subunits forming the protein-protein complex.

Scale: 0 to 1, above 0.8 is good, between 0.5 and 0.8 is grey zone, below 0.5 is bad

These values assume modelling with multiple recycling steps, so the process of prediction reaches a degree of convergence. In large-scale screenings for protein-protein interactions, often settings optimised for the speed of prediction are used, e.g. very few or no recycling steps. In such cases ipTM thresholds as low as 0.3 have been used for initial screening; importantly though, all pairs of proteins with ipTM scores higher than 0.3 were then subjected to additional examination (e.g. Weeratunga et al., 2023). Disordered regions and regions with low pLDDT score may negatively impact the ipTM score even if the structure of the complex is predicted correctly.

## Confidence scores - AlphaFold-PullDown

### ipTM: see above

### PI-score: Protein Interface Score

A geometry-based measure that evaluates the quality of the predicted interface by assessing. It’s computed as a post-processing quality metric.

High PI-score = likely biologically relevant and well-packed interface.

Low PI-score = loose or unrealistic interface geometry, possible artifact.

### mpDockQ/pDockQ score

Original metric for comparing predicted and experimental protein-protein docking models. Combines measures like interface RMSD, fraction of native contacts, and ligand RMSD.

Scale: 0 to 1, above 0.5 likely correct docking arrangment, 0.23 to 0.5 possible but uncertain, less then 0.23: likely incorrect docking geometry.

# Issues/considerations for my project

- VCAM1 is natively in dimer (R&D protein) -> check this in literature if its true and how to deal with this 
- Interactors of VCAM1 have low pLDDT and/or PEA scores, however, they can be used as the folding is likely dependent on the interaction
- There are multiple DBs that predict the TM region of proteins: Membranone, AlphaFoldTm and TmAlphaFold. 
- Isoforms of VCAM1 
- Only use the ECD of VCAM1 for screening
- Use also the IgG proteins that don't bind to VCAM1 IP as a negative control for the AlpaFold Pulldown, as such you can set threshold on the false positives. (ROC curves, score histograms and percentile-based threshold). You can set these thresholds with the ELISA dataset from Nuno. Then apply this to your possible interactions.

# Considerations for running the software

- It uses random seeds to initialize the structure, re-run the prediction with several different seeds
- Increase the number of recycles (3 to 20) is an effective method for improving prediction quality
- Optimize the MSA: deeper MSA with many 1000s of sequences will generally lead to a better prediction
- AlphaBridge (post-processing toolkit designed to analyze and visualize protein–protein complex predictions made by AlphaFold 3 / AlphaFold-Multimer.)
- AlphaFold-Multimer is from AlphaFold package

# VSC considerations
- The database needs be on the VSC_SCRATCH, not the VSC_DATA. How do I connect to this db, where is it stored in the VSC? 
For VUB HPC, its stored in /databases/bio/alphafold-<version> with the latest version alphafold-2.3.1 & alphafold-3.0.1
For KUL HPC, its stored in /lustre1/project/res_00002/lp_alphafold/ with the latest AlphaFoldDBdir_20240326
For Ughent HPC, its stored in /arcanine/scratch/gent/apps/AlphaFold

# Roadmap
- First find out whether in at KUL, VUB or Ughent they have the most recent version of AlphaFold or AlphaPulldown.
VUB: alphafold-2.3.1 & alphafold-3.0.1 (based on db version)
KUL: AlphaFold/2.3.4 (with module avail)

Follow: https://github.com/hpcleuven/AlphaFold

1) install locally and make the docker image, push docker image to your docker repo

2) pull the docker image from the repo in a singularity image on the hpc

#### 31.07.2025

Tried building the image locally from the gitrepo but had issue with my M3 chip.
On HPC, the Apptainer container installation worked. Used this docker image (catgumag/alphafold).
Issue: do not have read permission for reading in the databases. Asked to join lp_alphafold group.


- Run a simple AlphaFold of VCAM1
- Run a simple AlphaMultimer of VCAM1-VCAM1
- Run the AlphaPulldown of VCAM1 with interactors




