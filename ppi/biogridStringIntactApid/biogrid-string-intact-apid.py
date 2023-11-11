# This script uses the BioGrid, STRING, IntAct and APID databases to look for protein-protein interactions for a given query

# The first section of the script will have the custom functions necessary for downstream purposes.

# importing the necessary libraries
import argparse
import requests
import pandas as pd
import json
from upsetplot import UpSet, from_contents
import matplotlib.pyplot as plt
import itertools
from functools import reduce
import mygene
from io import StringIO
import time
import numpy as np
from bs4 import BeautifulSoup
import re

# Parsing the argument form the command line
parser = argparse.ArgumentParser(description = "Retrieves PPI interactions")
parser.add_argument("protein_to_query", type = str, help = "enter your protein of interest")
args = parser.parse_args()

protein = args.protein_to_query

# the query protein
query = protein + "_HUMAN"

#########################################################################################################################################
############################################################### FUNCTIONS ###############################################################
#########################################################################################################################################

def get_interactors_for_target(df, column_a, column_b, target_protein):
    """
    Takes a ppi dataframe in the format column1:interactor-of-query & column2:query, 
    and outputs from those two columns only one column with the interactor of the query (target protein)
    
    Args:
        df: the input dataframe.
        column_a: a column that holds either the query of the interactor of the query.
        column_b: a column that holds either the query of the interactor of the query.
        target_protein: the query where you want to find the interactors of.
    
    Returns:
        df: dataframe with one column less then input dataframe but now with only the interactor
        of the query.
    """    
    def get_interactors(row):
        if row[column_a] == target_protein:
            return row[column_b]
        elif row[column_b] == target_protein:
            return row[column_a]
        else:
            return None
    df["interactor_of_" + target_protein] = df.apply(get_interactors, axis = 1)
    
    return df


def convert_uniprotID_uniprotAcNr(uniprotID_or_AcNr): # e.g. can be P19320 or VCAM1_HUMAN
    """
    Converts a uniprotID (eg VCAM1_HUMAN) to a uniprot Accession Number (eg P19320), or 
    the other way around. This is done one by one, thus not in batch retrieval.
    """

    uniprot_api_url = "https://rest.uniprot.org/uniprotkb" 
    format = "json"
    uniprot_request_url = f"{uniprot_api_url}/{uniprotID_or_AcNr}?format={format}"
    
    uniprot_response = requests.get(uniprot_request_url)
    
    if uniprot_response.ok:       
        uniprot_dict = uniprot_response.json() 

        if "_" in uniprotID_or_AcNr:
            value = uniprot_dict["primaryAccession"]
        else:
            value = uniprot_dict["uniProtkbId"]
        return (value)
    
    else: 
        print("Data retrieval through Uniprot API failed, for the function convert_uniprotID_uniprotAcNr")
        
        
def extract_main_protein(isoform_list):
    """
    Extracts only the main proteins from a isoform list of uniprot accession numbers.
    This function is used in the function "convert_uniprotID_uniprotAcNr_from_df_column".
    """
    convert_None_to_NaN = [np.nan if isoform is None else isoform for isoform in isoform_list] # convert None to NaN
    convert_float_to_str = [str(isoform) if isoform is np.nan else isoform for isoform in convert_None_to_NaN] # conversion of NaN floats to NaN strings
    protein_list = list(filter(lambda isoform: "-" not in isoform, convert_float_to_str)) # the iteration, filtering out the isoforms
    
    return protein_list


def convert_uniprotID_uniprotAcNr_from_df_column(df, column_to_convert, new_column_name):
    """     
    Converts a pandas dataframe column that contains uniprot accession numbers (e.g. Q13740) 
    to uniprot IDs (e.g. CD166_HUMAN) & vice versa. It uses the Uniprot-API ID mapping
    and does this in batch retrieval so it has one large payload for in POST request.
    
    Args:
        df: the dataframe that contains the column that needs to be converted.
        columnd_to_convert: the column in the dataframe that needs conversion, 
            this column holds either uniprot IDs or accession numbers.
        new_column_name: name of the new column where the converted items 
            will be stored.
    
    Returns:
        An exact copy of the original dataframe but it holds an extra column that 
        contains the converted items, either uniprot IDs or names. It will also print
        the job ID so that you can check with Postman for example if it is not returning
        the expected output.
    
    Raises:
        JOB status error: job status is not running & will not get results, thus something 
            went wrong with job
        ValueError: The returned list from UniProt-API is not equal in length to the nr of 
            rows in the original dataframe
    """
    
    # Conversions for input into Uniprot-API
    df[column_to_convert] = df[column_to_convert].astype(str) # making sure it is a string
    list_of_uniprotIDs_or_Names = df[column_to_convert].tolist() # column to list
    ls = ",".join(list_of_uniprotIDs_or_Names) # conversion to input as payload to uniprot
    
    # making the POST request
    r = requests.post("https://rest.uniprot.org/idmapping/run", data={
        "from": "UniProtKB_AC-ID",
        "to": "UniProtKB", 
        "ids": ls 
    })

    # getting the job ID nr
    job_id = r.json()["jobId"]
    print("job ID:", job_id)

    # for loop that checks whether the job is running and when job is finished it will return a list
    while True:
        response = requests.get(f"https://rest.uniprot.org/idmapping/status/{job_id}")
        data_json = json.loads(response.text)
    
        if "jobStatus" in data_json:
            job_status = data_json["jobStatus"]
            print(f"job status: {job_status}")
            
            if job_status != "RUNNING":
                print("JobError: job status is not running, error with posting the request")
                break
    
        if "results" in data_json:
            if "_" in ls:
                desired_id = []
                for entry in list_of_uniprotIDs_or_Names:
                    matched = False
                    for json_entry in data_json["results"]:
                        if entry == json_entry["from"]:
                            desired_id.append(json_entry["to"]["primaryAccession"])
                            matched = True
                            break
                    if not matched:
                        desired_id.append("nan")
                desired_id_proteins = extract_main_protein(desired_id)   
                break

            else:
                desired_id_proteins = []
                for entry in list_of_uniprotIDs_or_Names:
                    matched = False
                    for json_entry in data_json["results"]:
                        if entry == json_entry["from"]:  
                            desired_id_proteins.append(json_entry["to"]["uniProtkbId"])
                            matched = True
                            break
                    if not matched:
                        desired_id_proteins.append("nan")
                break
        
        time.sleep(5) # wait for 5sec before the looop begins again
    
    # check before appending list to original dataframe, if the list has the same nr of items as the original df has rows
    if len(desired_id_proteins) != len(list_of_uniprotIDs_or_Names):
        raise ValueError("Data appending conflict: the returned list from UniProt-API is not equal in length to the nr of rows in the original dataframe")

    # now append list to dataframe   
    df[new_column_name] = desired_id_proteins
    
    return df


def convert_geneID_uniprotID(geneID):

    mygene_api_url = "https://mygene.info/v3/gene"
    entrezGeneID = geneID
    mygene_request_url = f"{mygene_api_url}/{entrezGeneID}?fields=all&dotfield=false&size=10"

    mygene_response = requests.get(mygene_request_url)
    
    if mygene_response.ok:     
        mygene_json = mygene_response.json()

        if "uniprot" in mygene_json:
            uniprot_id = mygene_json["uniprot"]["Swiss-Prot"]
        elif "pantherdb" in mygene_json:
            uniprot_id = mygene_json["pantherdb"]["uniprot_kb"]
        else:
            uniprot_id = "NA"
            
        return uniprot_id
    
    else:
        print("Data retrieval through MyGene API failed, for the function convert_geneID_uniprotID")
        

def convert_geneIDs_uniprotAcNr(df, df_column_name, df_new_column_name):
    """
    Using the MyGene-API python package this function converts entrez gene IDs to uniprot ID.
    It takes as input a dataframe column of which it will retrieve the uniprot accession numbers
    in batch retrieval and then append this list to the original dataframe.
    
    Args:
        df: original dataframe
        df_column_name: name of the column of which the entrez gene ID data will be send as an endpoint.
        df_new_column_name: name of the new column with the Uniprot accession number data.
        
    Returns:
        df: original dataframe with a new column of uniprot accession numbers.
    """

    mg = mygene.MyGeneInfo()
    
    # convert column to list
    list_of_column = df[df_column_name].tolist()
    
    # get the data from MyGene-API package
    df_mg = mg.getgenes(list_of_column, fields = "uniprot", as_dataframe = True)
    
    # filter it, convert to list & add column to original df
    df_mg_filtered = df_mg[["uniprot.Swiss-Prot"]]
    list_filtered = df_mg_filtered["uniprot.Swiss-Prot"].tolist()
    
    # append the list to the dataframe
    df[df_new_column_name] = list_filtered
    
    return df


def convert_stringID_to_uniprotName(string_id):
    """
    Converts string gene ID to uniprot ID using the uniprot-API,
    not in batch retrieval.
    """
    uniprot_api_url = "https://rest.uniprot.org/uniprotkb/search?query=gene_exact:"
    uniprot_request_url = f"{uniprot_api_url}{string_id}+AND+organism_id:9606"

    uniprot_response = requests.get(uniprot_request_url)
    
    if uniprot_response.ok:        
        uniprot_json = uniprot_response.json()

        if "results" in uniprot_json and len(uniprot_json["results"]) > 0:
            uniprot_name = uniprot_json["results"][0]["uniProtkbId"]
        else:
            uniprot_name = None

        return uniprot_name
    
    else:
        print("Data retrieval through Uniprot API failed, for the function convert_stringID_to_uniprotName")


def removeDuplicateRow_butRetainInfo (df, columnWithDuplicate, columnRetain1, columnRetain2, columnRetain3, columnRetain4 = None):
    """
    This function will take a dataframe that has duplicate rows for a certain column,
    which in this case is the interactor of the query. It will take the data from 
    certain specified columns and merge them into one row and retaining the info from the columns.
    
    Args:
        df: input dataframe
        columnWithDuplicate: the column which can hold duplicate values
        columnRetain1: column from the data should be retained
        columnRetain2: column from the data should be retained
        columnRetain3: column from the data should be retained
        columnRetain4: column from the data should be retained, this column is not necessary
    
    Returns:
        A dataframe that has no duplicate rows and contains all the info from the original dataframe.
    """
    
    def join_columnValues(series):
        return ", ".join(str(value) for value in series)

    columns_to_aggregate = {
        columnRetain1: join_columnValues,
        columnRetain2: join_columnValues,
        columnRetain3: join_columnValues}
    
    if columnRetain4 is not None:
        columns_to_aggregate[columnRetain4] = join_columnValues
        
    df_final = df.groupby(columnWithDuplicate).agg(columns_to_aggregate).reset_index()
    
    return df_final


def convert_uniprotAcNr_to_mouseMGI(uniprotAcNr):
    """
    Converts a uniprot accession number to mouse MGI ID using the panther-API,
    not in batch retrieval.
    """
    panther_api_url = "https://pantherdb.org/services/oai/pantherdb/ortholog/matchortho?"
    originOrganism = "9606"
    targetOrganism = "10090"
    panther_request_url = f"{panther_api_url}geneInputList={uniprotAcNr}&organism={originOrganism}&targetOrganism={targetOrganism}&orthologType=all"
    
    panther_response = requests.post(panther_request_url)
    
    if panther_response.ok:      
        panther_json = panther_response.json()
        
        mgi_ids = []
        if "search" in panther_json and "mapping" in panther_json["search"]:
            mapping = panther_json["search"]["mapping"]
            
            if "mapped" in mapping:
                if isinstance(mapping["mapped"], list):
                    for mapped_gene in mapping["mapped"]:
                        try:
                            mgi_id = mapped_gene["target_gene"].split("|")[1].split("=")[2]
                        except IndexError:
                            mgi_id = "NA"
                        mgi_ids.append(mgi_id)
                elif isinstance(mapping["mapped"], dict):
                    try:
                        mgi_id = mapping["mapped"]["target_gene"].split("|")[1].split("=")[2]
                    except (KeyError, IndexError):
                        mgi_id = "NA"
                    mgi_ids.append(mgi_id)
        
        return mgi_ids
    
    else:
        print("Data retrieval through Panther API failed, for the function: convert_uniprotAcNr_to_mouseMGI")
        
        
def convert_column_to_list(df, column):
    """
    Converts a column to a flattened list, specifcally for the upsetplot.
    """
    column_list = df[[column]].values.tolist()
    
    def flatten_list(nested_list):
        return list(itertools.chain(*nested_list))

    interactors = flatten_list(column_list)

    return interactors


#######################################################################################################################################
########################################################## DATA RETRIEVAL #############################################################
#######################################################################################################################################


########################################################## BioGRID API ################################################################

# calling the API
biogrid_api_url = "https://webservice.thebiogrid.org/interactions"

geneList = [convert_uniprotID_uniprotAcNr(query)]

params = {
    "accesskey": "caf14dfd9a0b7be447d282c322b8362e", # need to request
    "additionalIdentifierTypes": "UNIPROT",
    "format": "json",
    "geneList": geneList,
    "taxId": 9606, # human
    "max": 100000
}

response = requests.get(biogrid_api_url, params=params)

if response.ok:
    biogrid_interactions = response.json()

    biogrid_data = {}
    for interaction_id, interaction in biogrid_interactions.items():
        biogrid_data[interaction_id] = interaction
        biogrid_data[interaction_id]["INTERACTION_ID"] = interaction_id
else:
    print("Access to BioGrid database failed")
    
# Loading into a dataframe
biogrid_df = pd.DataFrame.from_dict(biogrid_data, orient="index")
columns = [
    "INTERACTION_ID",
    "ENTREZ_GENE_A",
    "ENTREZ_GENE_B",
    "OFFICIAL_SYMBOL_A",
    "OFFICIAL_SYMBOL_B",
    "EXPERIMENTAL_SYSTEM",
    "PUBMED_ID",
    "PUBMED_AUTHOR",
    "THROUGHPUT",
    "QUALIFICATIONS"]
biogrid_df = biogrid_df[columns]

# Data cleaning
df_biogrid = (biogrid_df
              .pipe(convert_geneIDs_uniprotAcNr, "ENTREZ_GENE_A", "UniprotID_A") # converting entrez gene IDs to uniprot AcNr
              .pipe(convert_geneIDs_uniprotAcNr, "ENTREZ_GENE_B", "UniprotID_B") # again converting entrez gene IDs to uniprot AcNr but for other column
              .pipe(convert_uniprotID_uniprotAcNr_from_df_column, "UniprotID_A", "UniprotName_A") # converting uniprotID to uniprot AcNr
              .pipe(convert_uniprotID_uniprotAcNr_from_df_column, "UniprotID_B", "UniprotName_B") # converting uniprotID to uniprot AcNr
              .filter(items = { "UniprotName_A", 
                                "UniprotName_B", 
                                "EXPERIMENTAL_SYSTEM", 
                                "PUBMED_ID",
                                "PUBMED_AUTHOR"}) # filtering the columns
              .rename(columns = {       "UniprotName_A": "biogrid_interactor_a", 
                                        "UniprotName_B": "biogrid_interactor_b",
                                        "EXPERIMENTAL_SYSTEM": "biogrid_method",
                                        "PUBMED_ID": "biogrid_pubID",
                                        "PUBMED_AUTHOR": "biogrid_publication"}) # renaming the columns
              .pipe(get_interactors_for_target, "biogrid_interactor_a", "biogrid_interactor_b", query) # only retaining the interactors column
              .pipe(removeDuplicateRow_butRetainInfo, "interactor_of_" + query ,"biogrid_publication", "biogrid_method", "biogrid_pubID")) # removing duplicate rows but retainig info


####################################################### IntAct using the PSIQUIC API #####################################################

# calling the API
intact_api_url = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices"

version = "current"
method = "interactor"
protein = query
format = "tab25"

intact_request_url = f"{intact_api_url}/{version}/search/{method}/{protein}?format={format}"

intact_response = requests.get(intact_request_url)

if intact_response.ok:
    int_act_interactions = intact_response.text
else:
    print("Access to IntAct database failed")
    
# Loading into a dataframe

intact_columns = [  "Unique identifier for interactor A", 
                    "Unique identifier for interactor B", 
                    "Alternative identifier for interactor A", 
                    "Alternative identifier for interactor B", 
                    "Aliases for A", 
                    "Aliases for B", 
                    "Interaction detection methods", 
                    "First author", 
                    "Identifier of the publication", 
                    "NCBI Taxonomy identifier for interactor A", 
                    "NCBI Taxonomy identifier for interactor B", 
                    "Interaction types", 
                    "Source databases", 
                    "Interaction identifier(s)", 
                    "Confidence score"]

intact_rows = [row.split("\t") for row in int_act_interactions.split("\n")]
intact_df = pd.DataFrame(intact_rows, columns=intact_columns)

# Data cleaning
df_intact = (intact_df
                .filter(items = { "Unique identifier for interactor A", 
                                            "Unique identifier for interactor B", 
                                            "Interaction detection methods", 
                                            "First author", 
                                            "Identifier of the publication", 
                                            "Confidence score"}) # filter columns
                .iloc[:-1] # remove last row, it holds no data when retrieving it through the API
                .assign(**{"Unique identifier for interactor A": lambda x: x["Unique identifier for interactor A"].str.removeprefix("uniprotkb:"),
                            "Unique identifier for interactor B": lambda x: x["Unique identifier for interactor B"].str.removeprefix("uniprotkb:")}) # removing the uniprotkbID refix from the name
                .loc[~intact_df["Unique identifier for interactor B"].str.contains("intact:", na = False)] # filtering out rows with IntactID instead of UniprotID
                .pipe(convert_uniprotID_uniprotAcNr_from_df_column, "Unique identifier for interactor A", "UniprotName_A") # converting uniprot ID to uniprot Ac Nr
                .pipe(convert_uniprotID_uniprotAcNr_from_df_column, "Unique identifier for interactor B", "UniprotName_B") # converting uniprot ID to uniprot Ac Nr
                .drop_duplicates() # dropping duplicates
                .rename(columns= {  "UniprotName_A": "IntAct_interactor_a", 
                                    "UniprotName_B": "IntAct_interactor_b",
                                    "Interaction detection methods": "IntAct_method",
                                    "First author": "IntAct_publication",
                                    "Identifier of the publication": "IntAct_pubID",
                                    "Confidence score": "IntAct_score"}) # renaming columns
                .pipe(get_interactors_for_target, "IntAct_interactor_a", "IntAct_interactor_b", query) # get only the interactors column
                .assign(**{"IntAct_score": lambda x: x["IntAct_score"].str.removeprefix("intact-miscore:")}) # removing a prefix from a certain column value
                .pipe(removeDuplicateRow_butRetainInfo, "interactor_of_" + query, "IntAct_publication", "IntAct_method", "IntAct_pubID", "IntAct_score") ) 


################################################################## STRING API ###################################################################

# Calling the API
string_api_url = "https://string-db.org//api"

output_format = "json"
method = "interaction_partners"

string_request_url = "/".join([string_api_url, output_format, method])

params = {
    "identifiers": query,
    "species": 9606, # human
    "required_score": 400,
    "limit": 1000000}

response = requests.post(string_request_url, data = params)

if response.ok:
    string_raw_data = response.text
    
# Loading into dataframe
string_df = pd.read_json(string_raw_data)

# Converting the GeneID to UniprotID
string_df[["UniprotID_A", "UniprotID_B"]] = string_df[["preferredName_A", "preferredName_B"]].map(convert_stringID_to_uniprotName) # takes a long time

# Data cleaning
df_string = (string_df
                  .filter(items = {"UniprotID_A", "UniprotID_B", "score", "escore"}) # filter out the columns that are needed
                  .rename(columns = {"UniprotID_A": "string_interactor_a", 
                                     "UniprotID_B": "string_interactor_b",
                                     "score": "string_score",
                                     "escore": "string_escore"}) # renaming column headers
                  .pipe(get_interactors_for_target, "string_interactor_a", "string_interactor_b", query) # get one column, only the interactor of the query
                  .dropna(subset = ["interactor_of_" + query]) # removing NAs
                  .sort_values(by = "string_escore") # sore of experimental score
                  .drop_duplicates(subset  =["interactor_of_" + query]) # removing duplicates
                  .loc[lambda x: x["string_escore"] != 0]) # only getting interactors with an experimental score != 0


################################################################ APID webscraping ###################################################################

def extract_table(proteinid):
    url = "http://cicblade.dep.usal.es:8080/APID/InteractionsGrid.action?protein1="+str(proteinid)+"&protein2=NA"

    payload = {}
    headers = {
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
    "Accept-Language": "en-US,en;q=0.9",
    "Connection": "keep-alive",
    "Cookie": "JSESSIONID=086030D12C94A8248DA2B5B9A84C16FA; _ga=GA1.2.581619552.1696441740; _gid=GA1.2.818200239.1696441740; _ga_7JSDHY18SK=GS1.2.1696441740.1.1.1696442421.0.0.0",
    "Referer": "http://cicblade.dep.usal.es:8080/APID/searchProtein.action",
    "Upgrade-Insecure-Requests": "1",
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36"
    }

    response = requests.request("GET", url, headers=headers, data=payload)
    soup = BeautifulSoup(response.text, "html.parser")
    table = soup.find("table", id="interactions") #Find the table

    rows = table.find_all("tr") # Find all table rows
    interactions = [] #initialize empty list
    for idx, row in enumerate(rows): #loop over rows, keep index
        if idx == 0: # first row is the header, skip.
            pass
        else:
            try:
                data1 = row.find_all("td") # get column
                interaction = { # build intraction object
                "ProteinA": data1[0].get_text().strip(),
                "ProteinB": data1[1].get_text().strip(),
                "MethodType": data1[2].get_text().strip(),
                "Method": data1[3].get_text().strip(),
                "Publication": re.sub(" +", " ",data1[4].get_text().strip().replace("\n", "")),
                "Source": data1[5].get_text().strip()
                }
                interactions.append(interaction) #append interaction object to result list
            except:
                pass
    return interactions #return the result list


# To find the UniProtID from the protein name
proteinid = convert_uniprotID_uniprotAcNr(query)
apid_results = extract_table(proteinid)

# loading into dataframe
apid_df = pd.DataFrame(apid_results)

# Data cleaning
df_apid = (apid_df
            .filter(items = {"ProteinA", "ProteinB", "Method", "Publication", "Source"})
            .rename(columns = {"ProteinA": "apid_interactor_a", 
                                "ProteinB": "apid_interactor_b",
                                "Method": "apid_method",
                                "Publication": "apid_publication",
                                "Source": "apid_source"})
            .pipe(get_interactors_for_target, "apid_interactor_a", "apid_interactor_b", query)
            .pipe(removeDuplicateRow_butRetainInfo, "interactor_of_" + query, "apid_method", "apid_publication", "apid_source"))


#######################################################################################################################################
##################################################### Intersections plot ##############################################################
#######################################################################################################################################

# Make an intersecions diagram using the pyUpSet

# Put the data into the right format using a custom function
biogrid_interactors = convert_column_to_list(df_biogrid, "interactor_of_" + query)
IntAct_interactors = convert_column_to_list(df_intact, "interactor_of_" + query)
string_interactors = convert_column_to_list(df_string, "interactor_of_" + query)
apid_interactors = convert_column_to_list(df_apid, "interactor_of_" + query)

# Plot the data into an upSetplot
ppis = from_contents({"BioGrid": biogrid_interactors, "IntAct": IntAct_interactors, "STRING": string_interactors, "APID": apid_interactors})
ax_dict = UpSet(ppis, subset_size="count", show_counts=True).plot()

# saving the upsetplot
plt.savefig("data/interactors_of_" + query + "_intersectionsPlot.png", dpi = 300)


#######################################################################################################################################
######################################## Merging data & getting subcellular locations #################################################
#######################################################################################################################################

dfs = [df_biogrid, df_intact, df_string, df_apid]
final_df = reduce(lambda left, right: pd.merge(left,right, on=["interactor_of_" + query], how="outer"), dfs)

# get the subcellular location of the interactor of the query, using the Uniprot-API
def get_subcellular_location(protein):
    
    uniprot_api_url = "https://rest.uniprot.org/uniprotkb" 
    format = "json"
    uniprot_request_url = f"{uniprot_api_url}/{protein}?format={format}"
    uniprot_response = requests.get(uniprot_request_url)

    uniprot_json = uniprot_response.json()

    subcellular_locations = []
    for comment in uniprot_json.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for subcellular_location in comment.get("subcellularLocations", []):
                location_value = subcellular_location.get("location", {}).get("value", "")
                subcellular_locations.append(location_value)
    
    return subcellular_locations
    
# now go over the dataframe and make a new column with the subcellular location
final_df[["subcellularLocation"]] = final_df[["interactor_of_" + query]].map(get_subcellular_location) # takes a long time

# undo the list in the "subcellular locations" column
final_df["subcellularLocation"] = final_df["subcellularLocation"].apply(lambda x: ", ".join(x))

# exporting the filtered dataset to a csv
final_df.to_csv("data/interactors_of_" + query + ".csv", index = False)


#######################################################################################################################################
###################################################### Filtering interactions #########################################################
#######################################################################################################################################

# filter on the following subcellular locations
subcellular_locations = ["membrane",
                         "cell junction",
                         "cell projection",
                         "cell membrane",
                         "plasma membrane",
                         "secreted",
                         "extracellular space",
                         "extracellular matrix"
                         "extracellular exosome",
                         "cell surface"] 

final_df_filtered = final_df[final_df["subcellularLocation"].str.contains('|'.join(subcellular_locations), case=False)]

# Data cleaning
df_final = (final_df_filtered
            # Convert the human uniprotID to human uniprotAcNr, make new column "interactor_of_VCAM1_uniprotAcNr"
            .pipe(convert_uniprotID_uniprotAcNr_from_df_column, "interactor_of_" + query, "interactor_of_" + query + "_Uniprot_AcNr")
            # Convert from human uniprotAcNr to mouseGene MGI
            .assign(**{"interactor_of_" + query + "_mouseGene_MGI": lambda x: x["interactor_of_" + query + "_Uniprot_AcNr"].map(convert_uniprotAcNr_to_mouseMGI)})
            # undo the list in the mouseGene MGI column (MGI is necessary for further downstream gene ontology)
            .assign(**{"interactor_of_" + query + "_mouseGene_MGI": lambda x: x["interactor_of_" + query + "_mouseGene_MGI"].apply(lambda y: ", ".join(y))})
            # column order
            .loc[:, ["interactor_of_" + query, "interactor_of_" + query + "_Uniprot_AcNr", "interactor_of_" + query + "_mouseGene_MGI"] 
                + list(final_df_filtered.columns.difference(["interactor_of_" + query, "interactor_of_" + query + "_Uniprot_AcNr", "interactor_of_" + query + "_mouseGene_MGI"]))])

# exporting the filtered dataset to a csv
final_df_filtered.to_csv("data/interactors_of_" + query + "_filtered.csv", index = False)

