import requests
import pandas as pd

# the URL
url = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/cadm4_human?format=tab25"

# retrieving the data from the psicquic API
data = requests.get(url).text

# # Getting the data into a pandas df

# # Split the data into a list of rows
# rows = [row.strip() for row in data.split('\n')]

# # Split each row into a list of columns (according to PSI-MI TAB 2.5 Format)
# columns = ['Unique identifier for interactor A', 
#            'Unique identifier for interactor B', 
#            'Alternative identifier for interactor A', 
#            'Alternative identifier for interactor B', 
#            'Aliases for A', 
#            'Aliases for B', 
#            'Interaction detection methods', 
#            'First author', 
#            'Identifier of the publication', 
#            'NCBI Taxonomy identifier for interactor A', 
#            'NCBI Taxonomy identifier for interactor B', 
#            'Interaction types', 
#            'Source databases', 
#            'Interaction identifier(s)', 
#            'Confidence score']

# # Create a pandas DataFrame from the list of rows and columns
# df = pd.DataFrame(rows, columns=columns)

# Print the DataFrame
# print(df)

# Split each row into a list of columns
columns = ['Unique identifier for interactor A', 
           'Unique identifier for interactor B', 
           'Alternative identifier for interactor A', 
           'Alternative identifier for interactor B', 
           'Aliases for A', 
           'Aliases for B', 
           'Interaction detection methods', 
           'First author', 
           'Identifier of the publication', 
           'NCBI Taxonomy identifier for interactor A', 
           'NCBI Taxonomy identifier for interactor B', 
           'Interaction types', 
           'Source databases', 
           'Interaction identifier(s)', 
           'Confidence score']
rows = [row.split('\t') for row in data.split('\n')]

# Create a pandas DataFrame from the list of rows and columns
df = pd.DataFrame(rows, columns=columns)

print(df)