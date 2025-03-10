import sys
from Bio import Entrez, SeqIO

# Set the email for Entrez
Entrez.email = "your_email@example.com"

def search_entrez(database, term):
    """
    Search Entrez database using a search term and return WebEnv and QueryKey for history API.
    """
    print(f"Searching Entrez database '{database}' for '{term}'...")  # Debugging output
    with Entrez.esearch(db=database, term=term, usehistory="y") as handle:
        record = Entrez.read(handle)

    print(f"Received record: {record}")  # Debugging output

    if not record["IdList"]:
        print("No IDs found.")
        return None, None, None

    return record["IdList"], record["WebEnv"], record["QueryKey"]

def fetch_sequences(database, webenv, query_key):
    """
    Fetch sequences using WebEnv and QueryKey.
    """
    print("Fetching sequences...")  # Debugging output
    print(f"WebEnv: {webenv}, QueryKey: {query_key}")  # Debugging output

    with Entrez.efetch(db=database, rettype="fasta", retmode="text", webenv=webenv, query_key=query_key) as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    print(f"Fetched {len(records)} sequences.")  # Debugging output
    return records


def main():
    """
    Main function to execute the search and fetch sequences.
    """
    if len(sys.argv) != 3:
        print("Usage: python script.py <database> <search_term>")
        sys.exit(1)

    database, term = sys.argv[1], sys.argv[2]
    ids, webenv, query_key = search_entrez(database, term)

    if not ids:
        print(f"No sequences found for the term '{term}' in the database '{database}'.")
        sys.exit(1)

    sequences = fetch_sequences(database, webenv, query_key)
    SeqIO.write(sequences, sys.stdout, "fasta")

if __name__ == "__main__":
    main()
