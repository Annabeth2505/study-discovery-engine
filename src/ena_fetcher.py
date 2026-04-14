import requests
import pandas as pd
import time 

ENA_PORTAL_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

GENERIC_SCIENTIFIC_NAMES = [
    "metagenome", "gut metagenome", "microbial mat metagenome",
    "soil metagenome", "food metagenome", "environmental metagenome",
    "marine metagenome", "freshwater metagenome", "human gut metagenome",
    "mouse gut metagenome", "bovine gut metagenome"
]

def search_ena_studies(query, max_results=100): 
    """ 
    Search ENA for studies matching a query. 

    Args:
        query: ENA search query string
        max_results: maximum number of studies to return 

    Returns:
        list of unique study accessions 
    """

    all_accessions = set()
    offset = 0 
    batch_size = 500

    while len(all_accessions) < max_results: 
        params = {
            "result": "read_run",
            "query": query,
            "fields": "study_accession",
            "limit": batch_size,
            "offset": offset,
            "format": "json"
        }

        try:
            response = requests.get(ENA_PORTAL_URL, params = params, timeout = 30)
            if response.status_code != 200:
                print(f'ENA search failed: {response.status_code}')
                break 

            data = response.json()
            if not data:
                break 

            for record in data:
                all_accessions.add(record['study_accession'])

            if len(data) < batch_size:
                break 

            offset += batch_size 
            time.sleep(0.3)

        except Exception as e:
            print(f'Error searching ENA: {e}')
            break 

def fetch_runs_for_study(study_accession):
    """  
    Fetch all runs for a study from ENA. 

    Args:
        study_accession: ENA/NCBI study accession 
    Returns:
        DataFrame of runs or none if fetch failed
    """
   
    params = {
        "result": "read_run",
        "query": f'study_accession="{study_accession}"',
        "fields": "run_accession,study_accession,sample_accession,scientific_name,host,library_strategy,fastq_ftp",
        "limit": 1000,
        "format": "json"
    }

    try:
        response = requests.get(ENA_PORTAL_URL, params = params, timeout = 30)
        if response.status_code != 200:
            print(f'ENA fetch failed for {study_accession}: {response.status_code}')
            return None 

        data = response.json()
        if not data:
            print(f'No runs found for {study_accession}')
            return None 

        df = pd.DataFrame(data)
        time.sleep(0.3)
        return df
    except Exception as e:
        print(f'Error fetching {study_accession}: {e}')
        return None


def resolve_host_species(runs_df):
    """ 
    Resolve the best available host species from a runs DataFrame.

    scientific_name is often generic like "gut metagenome" but host field contains
    the actual host species. This function falls back to host field whenscientific_name
    is uninformative. 

    Args:
        runs_df: DataFrame of runs for a single study 
    Returns:
        host_species string or None is not determinable 
    """

    sci_names = runs_df["scientific_name"].dropna()
    sci_name = sci_names.mode()[0] if len(sci_names) > 0 else None 

    if sci_name and sci_name.lower() not in [g.lower() for g in GENERIC_SCIENTIFIC_NAMES]:
        return sci_name 

    hosts = runs_df["host"].dropna()
    hosts = hosts[hosts != ""]
    if len(hosts) > 0:
        return hosts.mode()[0]

    return None
