import requests
import json
import re
from bs4 import BeautifulSoup


#Function gets protein id from name value. (Example: VCAM1_HUMAN -> P19320). We will need this protein ID to generate our table.
def extract_protein_id(name):
    url = "http://cicblade.dep.usal.es:8080/APID/searchProtein.action" # Endpoint to search for protein by name

    payload = 'proteinName='+str(name)+'&taxon=0'
    headers = {
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
    'Accept-Language': 'en-US,en;q=0.9',
    'Cache-Control': 'max-age=0',
    'Connection': 'keep-alive',
    'Content-Type': 'application/x-www-form-urlencoded',
    'Cookie': 'JSESSIONID=086030D12C94A8248DA2B5B9A84C16FA; _ga=GA1.2.581619552.1696441740; _gid=GA1.2.818200239.1696441740; _gat=1; _ga_7JSDHY18SK=GS1.2.1696441740.1.1.1696443448.0.0.0',
    'Origin': 'http://cicblade.dep.usal.es:8080',
    'Referer': 'http://cicblade.dep.usal.es:8080/APID/init.action',
    'Upgrade-Insecure-Requests': '1',
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36'
    }

    response = requests.request("POST", url, headers=headers, data=payload) # Execute HTTP request, response is stored in 'response'
    soup = BeautifulSoup(response.text, "html.parser") #Parse response.text (which is HTML) with a HTML parser
    table = soup.find('table', id="proteins") #use soup.find function to find a table in the HTML with id "proteins" (You can get this value by checking the HTML in the response from your HTTP request above. Hardcoded it to "proteins" because this will not change.)
    row = table.find_next('td') #use the beautifulsoup library again to find a html element in the found table. (Give me the first column value.)
    
    return row.text # row.text gives us the protein ID, we will need this protein ID to generate our table. (Example: VCAM1_HUMAN -> P19320)


def extract_table(proteinid):
    url = "http://cicblade.dep.usal.es:8080/APID/InteractionsGrid.action?protein1="+str(proteinid)+"&protein2=NA"

    payload = {}
    headers = {
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
    'Accept-Language': 'en-US,en;q=0.9',
    'Connection': 'keep-alive',
    'Cookie': 'JSESSIONID=086030D12C94A8248DA2B5B9A84C16FA; _ga=GA1.2.581619552.1696441740; _gid=GA1.2.818200239.1696441740; _ga_7JSDHY18SK=GS1.2.1696441740.1.1.1696442421.0.0.0',
    'Referer': 'http://cicblade.dep.usal.es:8080/APID/searchProtein.action',
    'Upgrade-Insecure-Requests': '1',
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36'
    }

    response = requests.request("GET", url, headers=headers, data=payload)
    soup = BeautifulSoup(response.text, "html.parser")
    table = soup.find('table', id="interactions") #Find the table

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
                "Publication": re.sub(' +', ' ',data1[4].get_text().strip().replace("\n", "")),
                "Source": data1[5].get_text().strip()
                }
                interactions.append(interaction) #append interaction object to result list
            except:
                pass
    return interactions #return the result list


proteinid = extract_protein_id("VCAM1_HUMAN")
results = extract_table(proteinid)
print(json.dumps(results, indent=4))













