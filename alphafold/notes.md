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

Visualization: heatmap of aligned residues (PEA plot).

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
- Only use the ECD of VCAM1 for screening, see (paper)[https://www.biorxiv.org/content/10.1101/2023.03.16.531341v2]
- Use also the IgG proteins that don't bind to VCAM1 IP as a negative control for the AlpaFold Pulldown, as such you can set threshold on the false positives. (ROC curves, score histograms and percentile-based threshold). You can set these thresholds with the ELISA dataset from Nuno. Then apply this to your possible interactions.

# Considerations for running the software

- It uses random seeds to initialize the structure, re-run the prediction with several different seeds
- Increase the number of recycles (3 to 20) is an effective method for improving prediction quality
- Optimize the MSA: deeper MSA with many 1000s of sequences will generally lead to a better prediction
- AlphaBridge (post-processing toolkit designed to analyze and visualize protein–protein complex predictions made by AlphaFold 3 / AlphaFold-Multimer.)
- AlphaFold-Multimer is from AlphaFold package
- monomer_ptm produces PAE erros, not monomer. Also multimer produces these PAE errors.

# HPC considerations
- The database needs be on the VSC_SCRATCH, not the VSC_DATA. How do I connect to this db, where is it stored in the VSC? 
- For VUB HPC, its stored in /databases/bio/alphafold-<version> with the latest version alphafold-2.3.1 & alphafold-3.0.1
- For KUL HPC, its stored in /lustre1/project/res_00002/lp_alphafold/ with the latest AlphaFoldDBdir_20240326
- For Ughent HPC, its stored in /arcanine/scratch/gent/apps/AlphaFold

# Running the AlphaFold on HPC
- First find out whether in at KUL, VUB or Ughent they have the most recent version of AlphaFold or AlphaPulldown.
VUB: alphafold-2.3.1 & alphafold-3.0.1 (module --show-hidden spider alphafold)
KUL: AlphaFold/2.3.4 (with module spider alphafold)
=> they both have the following most recent AlphaFold2 version: AlphaFold/2.3.4-foss-2022a-CUDA-11.7.0-ColabFold

For KUL, follow: https://github.com/hpcleuven/AlphaFold


#### 31.07.2025

Tried building the image locally from the gitrepo but had issue with my M3 chip.
On HPC, the Apptainer container installation worked. Used this docker image (catgumag/alphafold).
Issue: do not have read permission for reading in the databases. Asked to join lp_alphafold group, next day I joined the lp_alphafold group.

#### 04.08.2025 - 07.08.2025

I managed to build the apptainer container but now when I try a simple monomer run, it exits and warns about that a maximun number of residues are exceeded. Apparently, there is something wrong with the Uniref30 database that I'm using on the VSC KU Leuven, see [issue](https://arc.net/l/quote/hgezehdt), I cannot change this. Thus, I will try four things: 

(1) run AlphaFold without container using the module system. When running AlphaFold with the module system on KUL HPC, I get the same error of titin sequence.

(2) run it with reduced and full databases on KUL & VUB & UGent

    - KUL_reduced. I get an error, same as this [one](https://github.com/google-deepmind/alphafold/issues/743) Solution: Edit the Fasta file header with no spaces. KUL job nr: 64910853 (reduced). The output doesnt state why it failed, it probably got cancelled by the group account that I took, not due to the script. I requested more credits. I got access to lp_big_mem_gpu, to run GPU jobs here. I now have the error: "sbatch: error: Batch job submission failed: Invalid qos specification", I asked the HPC admin.

    - VUB. Error that it cant find the reduced db: "Could not find Jackhmmer database /databases/bio/alphafold-2.3.1/small_bfd/bfd-first_non_consensus_sequences.fasta". Solution: ran with a more newer version of AlphaFold (AlphaFold/2.3.4-foss-2022a-CUDA-11.7.0-ColabFold). Error: it could not find the database because there isnt any for 2.3.4 so I set the data-dir to the last version. VUB job nr: 10981208 (full) & 10981207 (reduced). Again the same error because in the previous version there isnt also anh "small_bfd". Seem that the reduced run isnt available in the VUB & it indeed is, no small_bfd in VUB. With the full run, I have again the issue with the titin sequence & again the uniref30 database used is 2021_03 (alphafold-2.3.1 db). Nextup: run it with the AlphaFold 2.2.0 verion and db as these do not seem to have the Uniref30 database (job nr 10983334). Issue: cannot load earlier AlphaFold modules (only v2.3.1 and newer) and these have all the Uniref30 db that is given errors. Next: run it for GPR37L1 with v2.3.4 & full db (job 10983466). This job ran for 4 hours but I got this (error)[https://github.com/google-deepmind/alphafold/issues/743]. This is related to version of AlphaFold, I need to run it with a more recent version (the header of the fasta file is short and wo spaces). Thus, I will now try to run alphafold inside an apptainer image from (this)[https://hub.docker.com/r/uvarc/alphafold/tags]. Next: run with version 2.3.1-foss-2022a-CUDA-11.7.0 on GPR37L1 (job nr 10983911). This run worked.

    GPR37L1 dimer. Tried running the dimer of GPR37L1 (job nr 10989886), this worked and took 9h.
    
    If this does not work: make apptainer container on KUL system and then transfer to this system.
    
    - from Ugent. Reduced_db failed. Two Ugent jobs: 40709080 (reduced) & 40709191 (full). With the full database job, I get the same titin error. With the reduced db run, I get just a truncated output log. Although I found the "small_bfd" db inside AlphaFold/20230310 and AlphaFold/2.3.2-foss-2023a-CUDA-12.1.1. Try to run it with GPR37L1 protein in reduced (40710071) and full db (40710070) mode for 16h. This failed because the system wasnt finding the GPU cores, so there this is a too long run. It's also the only version that is available on the Ugent system.

(3) run it with GPR37L1, is this titin problem VCAM1-specfic? Yes

(4) run with ColabFold/MMseqs2 databases

(5) ask HPC admins to use updated databases

- Run a simple AlphaFold of VCAM1
- Run a simple AlphaMultimer of VCAM1-VCAM1
- Run a simple AlphaMultimer of VCAM1-VLA4
- Run the AlphaPulldown of VCAM1 with interactors

Conclusions untill now:
VUB: Can't run the reduced_db option. Modules are outdated. Run (GPR37L1 - full_db - AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0) worked, took around 4h. Dimer of GPR37L1 also worked, took 9h.
KUL: 
UGent: only v2.3.2 available and it doesn't find the GPU.
General:
    - The titin error seems to be VCAM1-specific and not module specific

Can connect to login node of KUL & VUB with VSCode.

# Downstream analysis using ChimeraX, PyMOL and Python code for unpickling
There seems to be an issue with the compatibility with the Jaxlib version, see (issue)[https://github.com/jax-ml/jax/issues/18368]. I implement the suggested solution but then I got some more errors. Now I have a working downstream jupyter notebook for AlphaFold Multimer.



#### Good tutorials:
https://www.youtube.com/@Brown_Lab



## Questions for couse:

- ChimeraX: How do I select the chains that are only in contact with each other from outputs of AlphaFold Multimer



