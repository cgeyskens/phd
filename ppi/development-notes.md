## Things to think about
- human or mouse interactors

## Next-up
- look into omnipath.R script again (make R markdown work in VSCode)

## Workflow
- get data from PSICQUIC, STRING, OMNIPATH, APID and HIPPIE.
- from the UniProt API, make a function that converts proteinID/Name
- make JSON object, this can be used for other APIs to attach data to in the future
- remove duplicates between databases
- try to have a score for the interactions from each database (STRING score, nr of publications)
- filter: only include one of the following interactors (using the Uniprot API): transmembrane signal receptor, cell adhesion molecule, cell junction protein, defense immune protein, extracellular matrix protein, intercellular signal molecule, membrane traffic protein, transmembrane signal receptor (for protein class); cell junction, extracellular region part, extracellular region, membrane part, membrane, synapse part, synapse (for cellular component); biological adhesion, developmental process, immune system process, response to stimulus, signaling (for biological process).
- visualize the results in a venn diagram (upSetplot) and heatmap

## Example papers
Favuzzi et al. 2021 Cell
Gesuita et al. 2022 Cell Rep
Bernard et al. 2022 Science

## Goal
- Make this automated with a nextflow script calling the R and python scripts with just one parameter: the uniprot ID of the protein



## Webscraping tips from Simon (with APID db)

1. Go to the website
2. CNTRL + SHIFT + I (or go to settings -> more tools -> more developer tools), now you get a side window
3. put in the query (bv VCAM1_HUMAN)
4. click on the table icon to get the table
5. Now in the developer tools window, right-click on "interactions.actions?protein=P19320" -> copy -> copy as cURL (bash)
6. go to postman, as a GET request and paste the content in the bar 
7. this get you an HTML file, the content that copied are the headers of the request
8. from this content you can create a JSON object using BeautifulSoup


