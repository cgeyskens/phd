# In silico protein-protein interactions 
With this script, you can collect all the known **experimentally verified** interactions of a **human** or **mouse** protein from the following databases: [BioGRID](https://thebiogrid.org/), [STRING](https://string-db.org/), [IntAct](https://www.ebi.ac.uk/intact/home) & [APID](http://cicblade.dep.usal.es:8080/APID/init.action). You can find the script under the folder `biogridStrinIntactApid/`.  

<p align="center">
    <img src="image-readme.png" alt="drawing" width="500" height=200/>
</p>  

### How to use
Requirments in conda environment: `ppi-env.yml`

#### Input
To run the script, you must add the query protein as an arugment. The input format is Uniprot ID format and it can take both mouse and human uniprot IDs.  
If the query = `L1CAM_HUMAN`  of `L1CAM_MOUSE`

eg: `python biogrid-string-intact-apid.py L1CAM_HUMAN`   

eg: `python biogrid-string-intact-apid.py L1CAM_MOUSE`

#### Output
It will return three files inside a `/data` folder from your current directory:
1. `interactors_of_query.csv`: Unfiltered dataset of all the know interactions of the query protein
2. `interactors_of_query_intersectionsPlot.png`: An intersection plot (upsetplot) from which databases you gathered the interactions.
3. `interactors_of_query_filtered.csv` : Filtered dataset of all the known interaction of query protein. The dataset was filtered based on subcellular location (Uniprot). Interactors with the following subcellular locations were included: membrane, cell membrane, cell junction, plasma membrane, secreted, extracellular space, extracellular matrix and cell surface.

Both tables are in .csv format, with the following columns:

| Column                              | Description                                                                    |
| ----------------------------------- | -------------------------------------------------------------------------------|
| interactor_of_query                 | interactor of query in uniprot ID format (eg NCHL1_HUMAN ) 
| interactor_of_query_Uniprot_AcNr    | interactor of query in uniprot Accession Number format (eg O00533)
| IntAct_method                       | detection method of the interaction, if in IntAct database                    
| IntAct_pubID                        | publication ID of the interaction, if in IntAct database
| IntAct_publication                  | publication of the interaction, if in IntAct database
| apid_method                         | detection method of the interaction, if in APID database  
| apid_publication                    | publication of the interaction, if in APID database 
| apid_source                         | The primary ppi database source of APID, APID is a secondary database
| biogrid_method                      | detection method of the interaction, if in BioGRID database
| biogrid_pubID                       | publication ID of the interaction, if in BioGRID database
| biogrid_publication                 | publication of the interaction, if in BioGRID database 
| string_escore                       | experimental of the STRING database, if in STRING database
| string_score                        | general score of the STRING database, if in STRING database
| subcellularLocation                 | subcellular location of the interactor of the query, determined by Uniprot

### Database selection  
The goal of this project was to collect almost all the known experimentally verified protein protein interations. While researching which database I should use for collecting the known interactors of a certain protein, I saw that different papers use different databases i.e. there is no consensus. First, I tried the [omnipath](https://omnipathdb.org/) database, because it was suggested that it included multiple ppi resources. However, the main focus of omnipath is [causal interactions](https://github.com/saezlab/OmnipathR/issues/59) (pathway or activity flow) and it doesn't include all the interactions.  

This most recent benchmark [paper](https://pubmed.ncbi.nlm.nih.gov/32001390/) compared over 16 ppi databases. From the primary databases, I selected BioGRID and IntAct since they collectively have almost all of the exclusive experimentally verified protein interactions among all the primary databases (Figure 9). For the secondary databases, I selected STRING and APID since they are the top databases with the most experimentally verified interactions (Figure 6). Also, STRING includes data from BIND, DIP, GRID, HRPD, IntAct, MINT & PID. APID includes data from BioGRID, IntAct, MINT, HRPD, BioPlex and DIP. I also wanted to include HIPPIE, but it had problems with data retrieval during script development. These four databases, two primary (BioGRID & IntAct) and two secondary (STRING & APID), should have almost complete coverage of all the knwon experimentally verified ppi's.

### What is does in the background
This script will collect all the experimentally verified interactions of the input query protein form the STRING, IntAct, BioGRID and APID database using their APIs or through webscraping (APID). Since STRING also include non-experimentally verified interactions (such as from coexpression and textmining), I filtered those interactions based on whether the interaction has an experimental score (!= 0), and whether the combined score is larger then 0.4 (>0.4), which is medium to high confidence in the interaction based on all the different scores (such as the experimental score, coexpression score, gene neighborhood score, ...). The script combines all the data from each database and outputs a .csv file & an intersections plot, such that you can check the distributions from where each interaction comes and whether there is any overlap. Then, the script filters the interactors of the query based on their subcellular location (because our lab is mostly interested in cell-cell interactions) and will output this.  

### API usage
The following APIs were used: [Uniprot](https://www.uniprot.org/help/programmatic_access), [mygene](https://mygene.info/), [PANTHER](https://pantherdb.org/services/openAPISpec.jsp), [BioGRID](https://wiki.thebiogrid.org/doku.php/biogridrest), [PSICQUIC](https://psicquic.github.io/) & [STRING](https://www.string-db.org/help/api/). 