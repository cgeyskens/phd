## Things to think about
- human or mouse interactors


## Next-up
- look into omnipath.R script again (make R markdown work in VSCode)

## Workflow
1. get data from PSICQUIC, STRING, OMNIPATH (and APID webscraping?)
2. remove duplicates between databases
3. try to have a score for the interactions from each database (STRING score, nr of publications)
3. only include one of the following interactors (using the Uniprot API): transmembrane signal receptor, cell adhesion molecule, cell junction protein, defense immune protein, extracellular matrix protein, intercellular signal molecule, membrane traffic protein, transmembrane signal receptor (for protein class); cell junction, extracellular region part, extracellular region, membrane part, membrane, synapse part, synapse (for cellular component); biological adhesion, developmental process, immune system process, response to stimulus, signaling (for biological process).
4. visualize the results in a venn diagram and heatmap


## Goal
- Make this automated with a nextflow script calling the R and python scripts with just one parameter: the uniprot ID of the protein
