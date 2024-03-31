import sys
import urllib.request
from Bio import Entrez
from Bio import SeqIO

# Set the email for Entrez
Entrez.email = "your_email@example.com"


def search_entrez(database, term):
    '''
    Search Entrez database using a search term and return the list of IDs.

    Args:
    - database (str): The name of the Entrez database to query.
    - term (str): The search term to use for querying the database.

    Returns:
    - list: List of IDs retrieved from the Entrez database.
   ''' 
    handle = Entrez.esearch(db=database, term=term, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]


def fetch_sequences(database, ids):
    '''
    Fetch sequences from Entrez database using a list of IDs.

    Args:
    - database (str): The name of the Entrez database to query.
    - ids (list): List of IDs to fetch sequences for.

    Returns:
    - list: List of sequences retrieved from the Entrez database.
    '''
    handle = Entrez.efetch(db=database, id=ids, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return records


def main(database, term):
    '''
    Main function to search Entrez database and fetch sequences.

    Args:
    - database (str): The name of the Entrez database to query.
    - term (str): The search term to use for querying the database.
    '''

    ids = search_entrez(database, term)

    if not ids:
        print(f"No sequences found for the term '{term}' in the database '{database}'.")
        return

    sequences = fetch_sequences(database, ids)

    # Output sequences in FASTA format to STDOUT
    SeqIO.write(sequences, sys.stdout, "fasta")


if __name__ == "__main__":
    # Check command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <database> <search_term>")
        sys.exit(1)

    database = sys.argv[1]
    term = sys.argv[2]

    main(database, term)
