**This repository contains the processing pipelines, analysis scripts, and figure-generation code associated with the paper:**

A Systematic Cross-Modal Approach Identifies Astrocytic VCAM1 and GPR37L1 as Regulators of Hippocampal Synapse Development
--------------------------------------------------------------------------------------------------------------------------

_Cydric Geyskens, Dan Dascenco, Elke Leysen, Efstathia Kotoula, Kjen Bogaert, Joris Vandenbempt, Francisco Pestana, Benjamin Pavia, Steffen Fieuws, Lars Lefever, Pedro Magalhaes and Joris de Wit._

The code spans multiple experimental modalities (spatial transcriptomics, imaging, and proteomics) and is organized by figure. Each sub-repository includes scripts for raw data processing, analysis, and figure reproduction.

### Repository Structure
```
phd/
├── spatial-transcriptomics/ # Spatial transcriptomics analysis scripts
├── image-analysis/
│ ├── synapse-counting/ # Automated image processing and synapse quantification
│ └── in-vitro-synapse-counting/ # Semi-automated Fiji pipeline for in-vitro data
└── ip-proteomics/ # DIA-NN proteomics analysis and downstream statistics
```

### Data & Code availability

| Figure | Experiment/Modality | Raw Data Repo | Code Repo | 
| --- | ---| --- | --- | 
| 1, S1     | Image-based Spatial Transcriptomics | | [phd/spatial-transcriptomics](https://github.com/cgeyskens/phd/tree/master/spatial-transcriptomics) |
| 3, S3, S4 | Synapse In-Vivo Imaging Data | | [phd/image-analysis/synapse-counting](https://github.com/cgeyskens/phd/tree/master/image-analysis/synapse-counting)
| 4,  S5 | co-IP Proteomics | ProteomeXchange: PXD070371 & PXD070422 | [phd/ip-proteomics](https://github.com/cgeyskens/phd/tree/master/ip-proteomics)
| 5 | Synapse In-Vitro Imaging Data | BioImage Archive: S-BIAD2430 | [phd/image-analysis/in-vitro-synapse-counting](https://github.com/cgeyskens/phd/tree/master/image-analysis/in-vitro-synapse-counting)        

### Dependencies
Each sub-repo contains its own dependencies, either Conda, Apptainer or Docker.

### Author Contributions

```markdown
C.G. developed and implemented all data processing pipelines and statistical analyses.  
F.P. advised on the Resolve Bioscience spatial transcriptomics workflow analysis.  
B.P. reviewed the synapse image processing workflow.  
S.F. advised on statistical modeling.  
P.M. and D.D. provided input on proteomics analysis.
```
### Contact

For questions or requests, please contact:
**Cydric Geyskens** (cydric.geyskens@kuleuven.be) or **Joris de Wit** (joris.dewit@kuleuven.be)

### Citation

If you use this code, please cite:

Geyskens, C., et al. (2026). *A Systematic Cross-Modal Approach Identifies Astrocytic VCAM1 and GPR37L1 as Regulators of Hippocampal Synapse Development.* [GitHub Repository](https://github.com/cgeyskens/phd)
