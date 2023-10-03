import psicquic

# Create a Psicquic client
client = psicquic.PsicquicClient()

# Get the interaction data for the protein with identifier P01234
protein_id = "P01234"
interactions = client.get_interactions_by_protein_id(protein_id)

# Print the interaction data
for interaction in interactions:
    print(interaction)
