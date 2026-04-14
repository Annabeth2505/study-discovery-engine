import pandas as pd

def has_16S(runs_df):
    """
    Check if a study has 16S amplicon sequencing runs.
    
    Args: 
        runs_df = DataFrame of runs for a single study_accession
    Returns:
        bool: True if any run has library_strategy == AMPLICON, False otherwise
    """
    return any(runs_df['library_strategy']=='AMPLICON')


def has_WGS(runs_df): 
    """
    Check if a study has WGS sequencing runs. 
    
    Args: 
        runs_df = DataFrame of runs for a single study_accession
    Returns: 
        bool: True if any runs has library_strategy == WGS or METAGENOMIC
    """
    wgs_strategies = ['WGS', 'METAGENOMIC', 'Metagenomic', 'WHOLE GENOME SHOTGUN', 'metagenomics', 'shotgun']
    return any(runs_df['library_strategy'].isin(wgs_strategies))

def is_paired_v2(runs_df):
    """
    Check if a study has both 16S and WGS data under the same study_accession. 
    version 2 handles the real world case where paired studies are submitted with separate sample accession IDs for each sequencing type, but under the same study_accession ID. 

    Limitation: this function assumes that if a study_accession has both 16S and WGS runs, then they are paired. This may not always be the case, but it's a reasonable heuristic given the data available in SRA.
    Will handle this limitation in the agent layer wherein the scientist will be asked to verify if the 16S and WGS runs are indeed from the same samples. 
    
    Args:
        runs_df = DataFrame of runs for a single study_accession
    Returns:
        bool: True if the study_accession has both 16S and WGS runs, False otherwise
    """

    has_amplicon = any(runs_df['library_strategy']=='AMPLICON')
    has_wgs = has_WGS(runs_df)
    same_project = runs_df['study_accession'].nunique() == 1

    return has_amplicon and has_wgs and same_project

def has_fastq(runs_df):
    """
    Check if the study has downloadable  FASTQ files. 
    
    Args:
        runs_df = DataFrame of runs for a single study_accession
    Returns:
        bool: True if the study has downloadable FASTQ files, False otherwise
    """

    return any(
        runs_df['fastq_ftp'].notna() & 
        (runs_df['fastq_ftp'] != "")
    )


def check_hard_filters(runs_df):
    """
    Check if a study passes all four hard filters: has 16S data, has WGS data, is paired, and has downloadable FASTQ files. 
    
    Args:
        runs_df = DataFrame of runs for a single study_accession
    Returns: 
        dict: A dictionary with pass/fail results foreach filter and an overall passes field 
    """

    results = {
        'has_16S': has_16S(runs_df),
        'has_WGS': has_WGS(runs_df),
        'is_paired': is_paired_v2(runs_df),
        'has_fastq': has_fastq(runs_df)
    }

    results['passes'] = all([
        results['has_16S'], 
        results['has_WGS'], 
        results['is_paired'], 
        results['has_fastq']
    ])

    return results