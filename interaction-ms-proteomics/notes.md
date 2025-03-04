## Notes for Co-Ip MS-based proteomics

### Input from Pedro

##### Instrument specific parameters:

ZenoTOF range: 350 to 1250 m/z range

##### DIA-NN:

Pedro used version 19.

Mass accuracy: 20 and MS1 accuracy: 12

Always use latest version

MBR: Match-between-runs. If a peptide is abundant in at least one sample, then it will be extracted with high sensitivity across the entire experiment. If you don't set this, it will take the first sample as "bleuprint". The more variable your samples, the more advantagous this MBR is. 

1 missed cleavage, try 2, this will be complexer spectra but could potentially increase the amount of protein identified.

Always change the output folder because it will overwrite it! 

##### Downstream analysis packages:

Peptide-level analysis: MS Stats, MS-DAP

Protein level analysis: Alphapepstats: doesnt use limma but t-tests & doesnt adjust the p values. But you can make it work with additional code.

Papers: Peng et al. 2024 Nature comms

### Good YT videos:

https://www.youtube.com/watch?v=cGKJzx9IIi4

### Development notes

TODO:

V- copy raw data to HPC

V - make apptainer container for DIA-NN (see notes in README)

- try to run DIA-NN with same settings as Pedro on the cluster
  
  - Write the outputfiles to another folder, be carefull for OVERWRITING!
  - 

- make the .config file with the parameters

- Finetune the running parameters (try different ones)

- I have got paired samples so I should do a paired analysis!!!!