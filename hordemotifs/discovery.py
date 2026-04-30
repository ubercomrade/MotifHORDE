"""
discovery
=========

Classes for de‑novo motif discovery.  Each discovery tool implementation
wraps an external program (e.g. Streme, ProSampler, BaMM) or a
native algorithm and returns a collection of motif models.  All
classes inherit from :class:`MotifDiscoveryTool` which defines a
common interface.
"""

from __future__ import annotations

import os
import subprocess
import glob
from abc import ABC, abstractmethod
from typing import List

from .models import MotifModel, PwmMotif, BammMotif, SitegaMotif
from .execute import run_sitega


class MotifDiscoveryTool(ABC):
    """Abstract base class for motif discovery tools.

    Concrete implementations must implement :meth:`discover` which runs
    the underlying tool and returns a list of motif models.  Tools may
    accept additional configuration parameters via their constructor.
    """

    def __init__(self, name: str) -> None:
        self.name = name

    @abstractmethod
    def discover(self, foreground: str, background: str, output_dir: str, number_of_motifs: int, *args, **kwargs) -> List[MotifModel]:
        """Abstract method to perform de novo motif discovery.

        Concrete implementations must implement this method to run
        the underlying motif discovery tool and return a list of motif models.

        Parameters
        ----------
        foreground : str
            Path to the FASTA file containing positive sequences.
        background : str
            Path to the FASTA file containing negative/background sequences.
        output_dir : str
            Directory where intermediate files may be stored.  The
            directory will be created if it does not exist.
        number_of_motifs : int
            Maximum number of motifs to discover.

        Returns
        -------
        List[MotifModel]
            List of discovered motif models.  If no motifs are
            discovered an empty list should be returned.
        """
        raise NotImplementedError


class StremeDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around the MEME suite STREME program.

    STREME discovers PWM motifs by optimising the enrichment of a motif
    in the foreground relative to the background.  The discovered
    motifs are parsed from the MEME output and converted to PWMs.
    """

    def __init__(self, nmotifs: int = 5) -> None:
        super().__init__(name="streme")
        self.nmotifs = nmotifs

    def discover(self, foreground: str, background: str, output_dir: str, number_of_motifs: int, *args, **kwargs) -> List[MotifModel]:
        """Perform de novo motif discovery using STREME (Short Tandem Repeat EM).

        STREME is part of the MEME suite and identifies overrepresented motifs in DNA sequences
        by optimizing the enrichment of motifs in the foreground sequences relative to the
        background sequences. The discovered motifs are represented as Position Weight Matrices (PWMs).

        Parameters
        ----------
        foreground : str
            Path to the FASTA file containing positive sequences.
        background : str
            Path to the FASTA file containing negative/background sequences.
        output_dir : str
            Directory where intermediate files may be stored. The
            directory will be created if it does not exist.
        number_of_motifs : int
            Maximum number of motifs to discover.

        Returns
        -------
        List[MotifModel]
            List of discovered PWM motif models. If no motifs are
            discovered an empty list should be returned.

        Raises
        ------
        ValueError
            If 'length' parameter is not provided in kwargs.
        """

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        length = kwargs.get('length')
        if length is None:
            raise ValueError("Parameter 'length' is required for discovery")

        tmp_meme = os.path.join(output_dir, "motifs.meme")
        args = ['streme', '--p', foreground,
                '--n', background,
               '--objfun', 'de',
               '--w', f'{length}',
               '-nmotifs', f'{number_of_motifs}',
               '--text']

        p = subprocess.run(args, shell=False, capture_output=True)

        # Ensure we write stdout even if empty, but handle error
        with open(tmp_meme, 'wb') as file:
            file.write(p.stdout)

        if not os.path.exists(tmp_meme) or os.path.getsize(tmp_meme) == 0:
             return []

        # parse motifs from the generated file using the from_file class method
        motifs: List[MotifModel] = []
        
        for index in range(number_of_motifs):
            motif = PwmMotif.from_file(tmp_meme, index=index)
            motif.name = f'Streme-{index + 1}'
            motifs.append(motif)

        return motifs


class BammDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around the BaMM (Bayesian Markov Model) motif discovery program.

    BaMM discovers higher-order motifs by extending position weight matrices (PWMs)
    to model dependencies between positions in DNA sequences. This implementation
    first runs STREME to generate initial PWMs, then converts these to BaMM motifs.
    """

    def __init__(self) -> None:
        super().__init__(name="bamm")

    def discover(self, foreground: str, background: str, output_dir: str, number_of_motifs: int, *args, **kwargs) -> List[MotifModel]:
        """Perform de novo motif discovery using BaMM (Bayesian Markov Model).

        This method executes the BaMM motif discovery pipeline which involves first running STREME
        to generate initial PWM motifs, then converting these to higher-order BaMM motifs using
        the BaMMmotif tool. BaMMs capture dependencies between positions in DNA sequences beyond
        simple position weight matrices.

        Parameters
        ----------
        foreground : str
            Path to the FASTA file containing positive sequences.
        background : str
            Path to the FASTA file containing negative/background sequences.
        output_dir : str
            Directory where intermediate files may be stored. The
            directory will be created if it does not exist.
        number_of_motifs : int
            Maximum number of motifs to discover.

        Returns
        -------
        List[MotifModel]
            List of discovered BaMM motif models. If no motifs are
            discovered an empty list should be returned.

        Raises
        ------
        ValueError
            If 'length' parameter is not provided in kwargs.
        """

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        length = kwargs.get('length')
        order = kwargs.get('order', 2)  # Default order=2 for BaMM

        tmp_meme = os.path.join(output_dir, "motifs.meme")
        # call the helper from the original pipeline
        args = ['streme', '--p', foreground,
                '--n', background,
               '--objfun', 'de',
               '--w', f'{length}',
               '-nmotifs', f'{number_of_motifs}',
               '--text']
        p = subprocess.run(args, shell=False, capture_output=True)
        with open(tmp_meme, 'wb') as file:
            file.write(p.stdout)

        # parse motifs from the generated file
        args = ['BaMMmotif', output_dir, foreground, '--PWMFile', tmp_meme,
                '--EM', '--order', str(order), '--Order', str(order),
               '--basename', 'bamm',
               '--negSeqFile', background]
        _ = subprocess.run(args, capture_output=True)

        motifs = []
        for index in range(1, number_of_motifs + 1):
            bamm_path = f'{output_dir}/bamm_motif_{index}.ihbcp'
            bg_path = f'{output_dir}/bamm.hbcp'
            motif = BammMotif.from_file(bamm_path, bg_path, order)
            motifs.append(motif)

        return motifs


class SitegaDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around the SiteGA motif discovery program.

    SiteGA discovers motifs using a linear probabilistic grammar approach.
    """

    def __init__(self, nmotifs: int = 5) -> None:
        super().__init__(name="sitega")
        self.nmotifs = nmotifs

    def discover(self, foreground: str, background: str, output_dir: str, number_of_motifs: int, *args, **kwargs) -> List[MotifModel]:
        """Perform de novo motif discovery using SiteGA (Site Genetic Algorithm).

        SiteGA employs a genetic algorithm approach to identify overrepresented motifs in
        DNA sequences. It uses a linear probabilistic grammar model to represent motifs
        and optimize their detection in the input sequences.

        Parameters
        ----------
        foreground : str
            Path to the FASTA file containing positive sequences.
        background : str
            Path to the FASTA file containing negative/background sequences.
        output_dir : str
            Directory where intermediate files may be stored.
        number_of_motifs : int
            Maximum number of motifs to discover.

        Returns
        -------
        List[MotifModel]
            List of discovered SiteGA motif models. If no motifs are
            discovered an empty list should be returned.

        Raises
        ------
        ValueError
            If 'length' parameter is not provided in kwargs.
        """

        length = kwargs.get('length')
        number_of_lpd = kwargs.get('lpd', 20)  # Default lpd=20

        if length is None:
            raise ValueError("Parameter 'length' is required for discovery")

        # Run SiteGA
        run_sitega(output_dir, length, number_of_lpd, number_of_motifs)
        
        output_files = glob.glob(f"{output_dir}/train.fa_mat*")
        if len(output_files) == 0:
            return []  # No motifs found

        motifs = []
        for file in output_files:
            motif = SitegaMotif.from_file(file)
            motifs.append(motif)

        return motifs
