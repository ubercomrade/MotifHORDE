import os
import subprocess
import pandas as pd
from io import StringIO


def run_prosampler(foreground_path, background_path, output_dir, motif_length, number_of_motifs):
    """Execute ProSampler for de novo motif discovery.
    
    This function runs the ProSampler external tool to discover motifs from
    foreground and background FASTA sequences.
    
    Parameters
    ----------
    foreground_path : str
        Path to the foreground FASTA file containing sequences of interest
    background_path : str
        Path to the background FASTA file for comparison
    output_dir : str
        Directory where discovered motifs will be saved
    motif_length : int
        Length of motifs to discover
    number_of_motifs : int
        Number of motifs to discover
    
    Returns
    -------
    int
        Return code (0 indicates successful execution)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    args = ['ProSampler',
            '-i', foreground_path,
            '-b', background_path,
           '-k', f'{motif_length}', '-l', '0',
           '-m', f'{number_of_motifs}', '-z', '0', '-t', '4', '-w', '2',
           '-o', f'{output_dir}/motifs']
    
    result = subprocess.run(args, shell=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"ProSampler failed with return code {result.returncode}: {result.stderr}")
    
    return result.returncode


def run_tomtom(motifs_1, motifs_2):
    """Execute TomTom for motif comparison.
    
    This function runs the TomTom external tool to compare motifs from two
    different motif files and returns the comparison results as a DataFrame.
    
    Parameters
    ----------
    motifs_1 : str
        Path to the first motifs file (query motifs)
    motifs_2 : str
        Path to the second motifs file (target motifs)
    
    Returns
    -------
    pandas.DataFrame
        DataFrame containing the motif comparison results with columns
        as defined by TomTom output format
    """
    args = ['tomtom', motifs_1, motifs_2,
           '-thresh', '1',
            '-text' ]
    result = subprocess.run(args, shell=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"TomTom failed with return code {result.returncode}: {result.stderr}")
    table = pd.read_csv(StringIO(result.stdout), sep='\t', comment='#')
    return table


def run_motali_pwm_pwm(fasta_path, motif_1, motif_2, dist_1, dist_2,
                       overlap_path, all_path, sta_path):
    """Execute Motali for PWM-PWM motif comparison.
    
    This function runs the Motali external tool to perform pairwise comparison
    between Position Weight Matrices (PWMs) using specified distance metrics
    and outputs the similarity score.
    
    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file containing sequences for analysis
    motif_1 : str
        Path to the first PWM motif file
    motif_2 : str
        Path to the second PWM motif file
    dist_1 : str
        First distance metric/method to use for comparison
    dist_2 : str
        Second distance metric/method to use for comparison
    overlap_path : str
        Path where overlap results will be saved
    all_path : str
        Path where aggregate results will be saved
    sta_path : str
        Path where statistical results will be saved
    
    Returns
    -------
    float
        Similarity score calculated by Motali for the motif comparison
    """
    args = ['mco_prc.exe',
            fasta_path,
            'pwm', 'pwm',
            motif_1,
            motif_2,
            dist_1,
            dist_2,
            '50', '0.002',
            'hist_pass', '1',
            'prc_pass', '1',
            overlap_path,
            all_path,
            sta_path]
    result = subprocess.run(args, shell=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Motali failed with return code {result.returncode}: {result.stderr}")
    with open(all_path) as file:
       score = float(file.readline().strip())
    return score


def run_sitega(tmp_dir, motif_length, number_of_lpds, number_of_motifs):
    """Execute SiteGA for de novo motif discovery using genetic algorithm.
    
    This function runs the SiteGA external tool to discover motifs using
    a genetic algorithm approach with locally positioned dinucleotides (LPDs).
    
    Parameters
    ----------
    tmp_dir : str
        Temporary directory for SiteGA execution and output
    motif_length : int
        Length of motifs to discover
    number_of_lpds : int
        Number of locally positioned dinucleotides to consider
    number_of_motifs : int
        Number of motifs to discover
    
    Returns
    -------
    subprocess.CompletedProcess
        Result object from subprocess execution containing return code,
        stdout, stderr, and other execution details
    """
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    args = ['andy05cell.exe',  # Changed from hardcoded path to just executable name
            f'{tmp_dir}/',
            'train.fa',
            'background.fa',
            '6',
            f'{motif_length}',
            f'{number_of_lpds}',
            '6',
            '1',
            f'{tmp_dir}/',
            '5000',
            'sitega.log']
    result = subprocess.run(args, shell=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"SiteGA failed with return code {result.returncode}: {result.stderr}")
    return result
