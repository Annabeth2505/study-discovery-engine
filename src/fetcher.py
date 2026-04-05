from Bio import Entrez 
import pandas as pd 
import io
import time 
import os 
from dotenv import load_dotenv


def configure_entrez():
    """
    Configure Entrez with credentials from the .env file. 
    Must be called before any Entrez API calls.
    """
    dotenv_path = os.path.join(os.path.dirname(__file__),'..', '.env')
    load_dotenv(dotenv_path)

    Entrez.email = os.getenv('NCBI_EMAIL')
    Entrez.api_key = os.getenv('NCBI_API_KEY')
    if not Entrez.email:
        raise ValueError("NCBI_EMAIL not set in .env file. Please provide an email for Entrez API access.")

    print("Entrez configured with email:", Entrez.email)


def fetch_sra_studies (search_term, max_results = 100, batch_size =100) :
    """ 
    Search SRA for studies matching a search term. 
    Handles pagination to retrieve results beyond the first batch. 

    Args:
        search_term: keyword or taxonomy query string
        max_results: maximum number of studies to retrieve
        batch_size: number of records to fetch per request (max 10000)
    
    Returns:
        List of SRA IDs 
    
    """
    try:
        handle = Entrez.esearch(
            db = 'sra',
            term = search_term, 
            retmax = 0
        )
        record = Entrez.read(handle)
        handle.close()

        total_count = int(record['Count'])

        if max_results is None:
            total_to_fetch = total_count
        else:
            total_to_fetch = min(max_results, total_count)
            
        print(f'Found {total_count} studies. Fetching {total_to_fetch}...')

        all_ids = []
        retstart =0 

        while retstart < total_to_fetch:
            batch = min(batch_size, total_to_fetch - retstart)

            handle = Entrez.esearch(
                db = 'sra',
                term = search_term,
                retstart = retstart,
                retmax = batch
            )
            record = Entrez.read(handle)
            handle.close()


            all_ids.extend(record['IdList'])
            retstart += batch

            print(f' Fetched {len(all_ids)} of {total_to_fetch}')
            time.sleep(0.5)  # Be respectful of NCBI servers

        return all_ids
    
    except Exception as e:
        print(f"Error fetching SRA studies for search term '{search_term}': {e}")
        return []           

def fetch_runs_for_bioproject(sra_id):
    """ 
    Fetch al runs for a given SRA study ID and return as DataFrame. 
    
    Args:
        sra_id : NCBI internal numeric ID for the study
        Returns:
            DataFrame of runs or None if fetch failed. 
    """

    try:
        fetch_handle = Entrez.efetch(
            db = 'sra',
            id = sra_id,
            rettype = 'runinfo',
            retmode = 'text'
        )
        
        data = fetch_handle.read()
        fetch_handle.close()

        decoded = data.decode('utf-8')
        df = pd.read_csv(io.StringIO(decoded))
        df = df.dropna(subset = ['Run'])
        time.sleep(0.3)  # Be respectful of NCBI servers
        return df 
    
    except Exception as e:
        print(f"Error fetching runs for SRA ID {sra_id}: {e}") 
        return None

