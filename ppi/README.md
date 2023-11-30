# In silico protein-protein interactions 
With this script you can collect all the known interactions of a human or mouse protein from the following databases: [BioGRID](https://thebiogrid.org/), [STRING](https://string-db.org/), [IntAct](https://www.ebi.ac.uk/intact/home) & [APID](http://cicblade.dep.usal.es:8080/APID/init.action). You can find the script under the folder biogridStrinIntactApid/. 

<img src="image-readme.png" alt="drawing" width="500" height=200/>

### How to use
Requirments in conda environment: `ppi-env.yml`

##### Input
To run the script you must add the query protein right after script name. The format must be in Uniprot ID format and it can take both mouse and human uniprot IDs.  
If the query = `L1CAM_HUMAN`
eg: `python biogrid-string-intact-apid.py L1CAM_HUMAN`
eg: `python biogrid-string-intact-apid.py L1CAM_HUMAN`

##### Output
It will return three files inside a `/data` folder:
1. `interactors_of_query.csv`: Unfiltered dataset of all the know interactions of the query protein
2. `interactors_of_query_intersectionsPlot.png`: An intersection plot (upsetplot) from which database you gathered the interactions.
3. `interactors_of_query_filtered.csv` : Filtered dataset of all the known interaction of query protein. The dataset was filtered based on subcellular location (Uniprot). Interactors with the following subcellular locations were included: membrane, cell membrane, cell junction, plasma membrane, secreted, extracellular space, extracellular matrix and cell surface.

Both tables are in .csv format, with the following columns:

| Column                              | Description                                                                                                       |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------|
| interactor_of_query                 | interactor of query in uniprot ID format (eg NCHL1_HUMAN ) 
| interactor_of_query_Uniprot_AcNr    | interactor of query in uniprot Accession Number format (eg O00533)
| IntAct_method                       | detection method of the interaction, if in IntAct database                    
| IntAct_pubID                        | publication id of the interaction, if in IntAct database
| IntAct_publication                  | publication of the interaction, if in IntAct database



### Database selection  
Describe here your reasoning for database selection. Reference to other papers.
There is no real consensus which database to use. 

### What is does

### Background to script development
Which APIs that you use. Omnipath first, refer to issue.