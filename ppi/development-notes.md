## Things to think about
- human or mouse interactors

## Next-up
- look into omnipath.R script again (make R markdown work in VSCode)

## Problems 
- HIPPIE website is out

## Done things
- get data from secondary: STRING, APID, HIPPIE & primary: BioGRID & IntAct
- plot the data using the upSetplot for intersections
- from the UniProt API, make a function that converts proteinID/Name
- try to have a score for the interactions from each database (STRING score, nr of publications)
- cleaning up the dataframes further seperatly before merge (merging rows based on one value in one column)
- make one large dataframe off all the databases
- filter using the Uniprot API: only include one of the following interactors (using the Uniprot API, or PANTHER?): transmembrane signal receptor, cell adhesion molecule, cell junction protein, defense immune protein, extracellular matrix protein, intercellular signal molecule, membrane traffic protein, transmembrane signal receptor (for protein class); cell junction, extracellular region part, extracellular region, membrane part, membrane, synapse part, synapse (for cellular component); biological adhesion, developmental process, immune system process, response to stimulus, signaling (for biological process).
- I'm losing a lot of data in the conversions (human_UniprotID -> human_UniprotAcNr -> human_ensembl -> mouse_ensembl -> mouseGeneSymbol, should be better to go from human_UniprotAcNr -> mouse_UniprotAcNr -> mouse_geneSymbol). Now much straightforward using the pantherDB API: from human_UniprotID -> human_UniprotAcNr -> MGIgene
- batch retrieval REST calling, especially Uniprot/Mygene/Ensemble. Use these APIs more effectively.
- properly document each function & each step
- pipe the dataframe operations with the pipe()
- if/print statements during each API request, 200 response or not
- make it into a .py script that is callable from the CLI (arggparse)
- debug last issues: converting from entrez gene IDs to uniprot IDs using the myGene API python package, a lot is lost in this conversion
- Also for mouse proteins? No for biogrid, string has taxID, made function for this at the beginning. Problem with mouse ID in string with the function: convert_stringID_to_uniprotName
- Check for the databases which uniprot format input is preferred. G37L1_HUMAN in APID for example & if they accept mouse Uniprot ID. IntAct prefers AcNR, Biogrid prefers AcNr, APID prefers AcNr -> do all input as AcNrs
- something wrong with the function -get_interactors_for_target- but only for string data cleaning, why? Dont know, just dropped it at the end.
- Also include statements if the get request from the API doesnt return anything, ie no results returned if there is actually no results and it handles it well.
- Format of .csv file should be the same for unfiltered and filtered, same order of column names. Reoder the last part of the script.
- check the filtering, because there is also mitochondrial membrane, endoplamic reticulum, ... in there.
- try all of your proteins, and debug the issues. Tested all 5 (mouse and human), and the isssue's are:

## TODO
- clean up the script:
    - write a documentation README file, with png and how to use the script

- make a heatmap & UpSetplot
    - add heatmap of GO terms (phagocytosis, synapse phagocytosis), now I just uploaded a .csv from MGI for these GO terms
    - add heatmap of the dangyang et al. 2022. Again need to reanalyze the data, now I just uploaded the xlsx file from the paper

## Example papers
Favuzzi et al. 2021 Cell
Gesuita et al. 2022 Cell Rep
Bernard et al. 2022 Science

## Goal
Make this automated e.g.  $ script.py -p "VCAM1_HUMAN"
outputs: 
- interactions.csv
- upsetplot.png
- interactions_filtered.csv (filtered on subcellular location)

Make separate script for plotting:
- heatmap 


## Webscraping tips from Simon (with APID db)

1. Go to the website
2. CNTRL + SHIFT + I (or go to settings -> more tools -> more developer tools), now you get a side window
3. put in the query (bv VCAM1_HUMAN)
4. click on the table icon to get the table
5. Now in the developer tools window, right-click on "interactions.actions?protein=P19320" -> copy -> copy as cURL (bash)
6. go to postman, as a GET request and paste the content in the bar 
7. this get you an HTML file, the content that copied are the headers of the request
8. from this content you can create a JSON object using BeautifulSoup


