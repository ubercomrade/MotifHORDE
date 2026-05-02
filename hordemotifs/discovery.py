"""Wrappers for external de novo motif discovery tools."""

from __future__ import annotations

import glob
import logging
import os
import shutil
from abc import ABC, abstractmethod
from typing import List

from .execute import run_sitega
from .external import (
    DEFAULT_BAMM_COMMAND,
    DEFAULT_DIMONT_JAR,
    DEFAULT_MEME_COMMAND,
    DEFAULT_SLIM_JAR,
    DEFAULT_STREME_COMMAND,
    resolve_command,
    resolve_existing_path,
    run_checked,
)
from .io import write_jstacs_fasta
from .models import GenericModel, read_model

logger = logging.getLogger(__name__)


class MotifDiscoveryTool(ABC):
    """Base interface for de novo discovery wrappers."""

    def __init__(self, name: str) -> None:
        self.name = name

    @abstractmethod
    def discover(
        self,
        foreground: str,
        background: str,
        output_dir: str,
        number_of_motifs: int,
        *args,
        **kwargs,
    ) -> List[GenericModel]:
        """Run the discovery tool and return discovered models."""
        raise NotImplementedError


def _require_length(kwargs: dict) -> int:
    length = kwargs.get("length")
    if length is None:
        raise ValueError("Parameter 'length' is required for discovery")
    return int(length)


def _read_indexed_models(path: str, model_type: str, prefix: str, number_of_motifs: int) -> List[GenericModel]:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []

    motifs: List[GenericModel] = []
    for index in range(number_of_motifs):
        try:
            motif = read_model(path, model_type, index=index)
        except IndexError:
            break
        except ValueError:
            if index == 0:
                return []
            break
        motif.name = f"{prefix}-{index + 1}"
        motifs.append(motif)
    return motifs


def _run_streme(
    foreground: str,
    background: str,
    output_dir: str,
    length: int,
    number_of_motifs: int,
    command: str | None = None,
) -> str:
    os.makedirs(output_dir, exist_ok=True)
    tmp_meme = os.path.join(output_dir, "motifs.meme")
    resolved_command = command or resolve_command(
        "streme",
        DEFAULT_STREME_COMMAND,
        "HORDEMOTIFS_STREME_COMMAND",
    )
    args = [
        resolved_command,
        "--p",
        foreground,
        "--n",
        background,
        "--objfun",
        "de",
        "--w",
        str(length),
        "-nmotifs",
        str(number_of_motifs),
        "--text",
    ]
    result = run_checked(args)
    with open(tmp_meme, "w") as handle:
        handle.write(result.stdout)
    return tmp_meme


def _existing_xml_paths(output_dir: str, patterns: list[str]) -> list[str]:
    paths: list[str] = []
    for pattern in patterns:
        paths.extend(sorted(glob.glob(os.path.join(output_dir, pattern))))
        if paths:
            return paths
    return paths


def _read_xml_models(paths: list[str], model_type: str, prefix: str, number_of_motifs: int) -> List[GenericModel]:
    motifs: List[GenericModel] = []
    for path in paths:
        try:
            motif = read_model(path, model_type)
        except ValueError as exc:
            logger.warning("Skipping invalid %s output %s: %s", model_type, path, exc)
            continue
        motif.name = f"{prefix}-{len(motifs) + 1}"
        motifs.append(motif)
        if len(motifs) >= number_of_motifs:
            break
    return motifs


class StremeDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around STREME PWM discovery."""

    def __init__(self, nmotifs: int = 5, command: str | None = None) -> None:
        super().__init__(name="streme")
        self.nmotifs = nmotifs
        self.command = command

    def discover(
        self,
        foreground: str,
        background: str,
        output_dir: str,
        number_of_motifs: int,
        *args,
        **kwargs,
    ) -> List[GenericModel]:
        length = _require_length(kwargs)
        tmp_meme = _run_streme(foreground, background, output_dir, length, number_of_motifs, self.command)
        return _read_indexed_models(tmp_meme, "pwm", "Streme", number_of_motifs)


class MemeDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around MEME PWM discovery."""

    def __init__(
        self,
        command: str | None = None,
        objfun: str = "classic",
        model: str = "zoops",
        minsites: int | None = None,
        maxsites: int | None = None,
        seed: int | None = None,
        threads: int | None = None,
    ) -> None:
        super().__init__(name="meme")
        self.command = command
        self.objfun = objfun
        self.model = model
        self.minsites = minsites
        self.maxsites = maxsites
        self.seed = seed
        self.threads = threads

    def discover(
        self,
        foreground: str,
        background: str,
        output_dir: str,
        number_of_motifs: int,
        *args,
        **kwargs,
    ) -> List[GenericModel]:
        length = _require_length(kwargs)
        os.makedirs(output_dir, exist_ok=True)
        tmp_meme = os.path.join(output_dir, "motifs.meme")
        command = self.command or resolve_command("meme", DEFAULT_MEME_COMMAND, "HORDEMOTIFS_MEME_COMMAND")
        meme_args = [
            command,
            foreground,
            "-dna",
            "-revcomp",
            "-nmotifs",
            str(number_of_motifs),
            "-w",
            str(length),
            "-text",
            "-objfun",
            self.objfun,
            "-mod",
            self.model,
        ]
        if self.objfun != "classic":
            meme_args[4:4] = ["-neg", background]
        optional_args = [
            ("-minsites", self.minsites),
            ("-maxsites", self.maxsites),
            ("-seed", self.seed),
            ("-p", self.threads),
        ]
        for flag, value in optional_args:
            if value is not None:
                meme_args.extend([flag, str(value)])

        result = run_checked(meme_args)
        if result.stdout:
            with open(tmp_meme, "w") as handle:
                handle.write(result.stdout)
        elif os.path.exists(os.path.join(output_dir, "meme.txt")):
            tmp_meme = os.path.join(output_dir, "meme.txt")
        return _read_indexed_models(tmp_meme, "pwm", "Meme", number_of_motifs)


class BammDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around BaMM motif discovery."""

    def __init__(self, bamm_command: str | None = None, streme_command: str | None = None) -> None:
        super().__init__(name="bamm")
        self.bamm_command = bamm_command
        self.streme_command = streme_command

    def discover(
        self,
        foreground: str,
        background: str,
        output_dir: str,
        number_of_motifs: int,
        *args,
        **kwargs,
    ) -> List[GenericModel]:
        length = _require_length(kwargs)
        order = int(kwargs.get("order", 2))

        tmp_meme = _run_streme(foreground, background, output_dir, length, number_of_motifs, self.streme_command)
        if not os.path.exists(tmp_meme) or os.path.getsize(tmp_meme) == 0:
            return []

        bamm_command = self.bamm_command or resolve_command(
            DEFAULT_BAMM_COMMAND,
            DEFAULT_BAMM_COMMAND,
            "HORDEMOTIFS_BAMM_COMMAND",
        )
        bamm_args = [
            bamm_command,
            output_dir,
            foreground,
            "--PWMFile",
            tmp_meme,
            "--EM",
            "--order",
            str(order),
            "--Order",
            str(order),
            "--basename",
            "bamm",
            "--negSeqFile",
            background,
        ]
        run_checked(bamm_args)

        motifs: List[GenericModel] = []
        for index in range(1, number_of_motifs + 1):
            bamm_path = os.path.join(output_dir, f"bamm_motif_{index}.ihbcp")
            if not os.path.exists(bamm_path):
                logger.warning("Skipping missing BaMM output: %s", bamm_path)
                continue
            motif = read_model(bamm_path, "bamm", order=order)
            motif.name = f"Bamm-{index}"
            motifs.append(motif)
        return motifs


def _build_dimont_args(
    java_command: str,
    java_xmx: str,
    jar_path: str,
    output_dir: str,
    data_name: str,
    length: int,
    position_tag: str,
    value_tag: str,
    bg_order: int,
    motif_order: int,
    ess: float,
    starts: int,
    threads: int | None,
) -> list[str]:
    args = [
        java_command,
        "-Djava.awt.headless=true",
        f"-Xmx{java_xmx}",
        "-jar",
        jar_path,
        f"home={output_dir}",
        f"data={data_name}",
        "infix=dimont",
        f"position={position_tag}",
        f"value={value_tag}",
        f"motifWidth={length}",
        f"motifOrder={motif_order}",
        f"bgOrder={bg_order}",
        f"ess={ess}",
        f"starts={starts}",
    ]
    if threads is not None:
        args.append(f"threads={threads}")
    return args


def _build_slim_args(
    java_command: str,
    java_xmx: str,
    jar_path: str,
    output_dir: str,
    data_name: str,
    length: int,
    position_tag: str,
    value_tag: str,
    bg_order: int,
    motif_order: int,
    modify: bool,
    starts: int,
    threads: int | None,
) -> list[str]:
    args = [
        java_command,
        "-Djava.awt.headless=true",
        f"-Xmx{java_xmx}",
        "-jar",
        jar_path,
        f"home={output_dir}",
        f"data={data_name}",
        "infix=slim",
        f"position={position_tag}",
        f"value={value_tag}",
        f"motifWidth={length}",
        f"motifOrder={motif_order}",
        f"bgOrder={bg_order}",
        f"starts={starts}",
        f"modify={str(modify).lower()}",
    ]
    if threads is not None:
        args.append(f"threads={threads}")
    return args


class DimontDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around Jstacs Dimont discovery."""

    def __init__(
        self,
        jar_path: str | None = None,
        java_command: str = "java",
        java_xmx: str = "4G",
        threads: int | None = None,
        position_tag: str = "position",
        value_tag: str = "value",
        bg_order: int = -1,
        motif_order: int = 0,
        ess: float = 4.0,
        starts: int = 20,
    ) -> None:
        super().__init__(name="dimont")
        self.jar_path = jar_path
        self.java_command = java_command
        self.java_xmx = java_xmx
        self.threads = threads
        self.position_tag = position_tag
        self.value_tag = value_tag
        self.bg_order = bg_order
        self.motif_order = motif_order
        self.ess = ess
        self.starts = starts

    def discover(
        self,
        foreground: str,
        background: str,
        output_dir: str,
        number_of_motifs: int,
        *args,
        **kwargs,
    ) -> List[GenericModel]:
        length = _require_length(kwargs)
        os.makedirs(output_dir, exist_ok=True)
        jar_path = resolve_existing_path(
            self.jar_path,
            "HORDEMOTIFS_DIMONT_JAR",
            DEFAULT_DIMONT_JAR,
            "Dimont jar",
        )
        java_command = resolve_command(self.java_command)
        data_name = "train.annot.fa"
        write_jstacs_fasta(
            foreground,
            os.path.join(output_dir, data_name),
            self.position_tag,
            self.value_tag,
        )
        cmd = _build_dimont_args(
            java_command,
            self.java_xmx,
            jar_path,
            output_dir,
            data_name,
            length,
            self.position_tag,
            self.value_tag,
            self.bg_order,
            self.motif_order,
            self.ess,
            self.starts,
            self.threads,
        )
        run_checked(cmd)
        paths = _existing_xml_paths(output_dir, ["*dimont*.xml", "*Dimont*.xml", "*.xml"])
        return _read_xml_models(paths, "dimont", "Dimont", number_of_motifs)


class SlimDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around Jstacs SlimDimont discovery."""

    def __init__(
        self,
        jar_path: str | None = None,
        java_command: str = "java",
        java_xmx: str = "4G",
        threads: int | None = None,
        position_tag: str = "position",
        value_tag: str = "value",
        bg_order: int = -1,
        motif_order: int = -5,
        modify: bool = True,
        starts: int = 20,
    ) -> None:
        super().__init__(name="slim")
        self.jar_path = jar_path
        self.java_command = java_command
        self.java_xmx = java_xmx
        self.threads = threads
        self.position_tag = position_tag
        self.value_tag = value_tag
        self.bg_order = bg_order
        self.motif_order = motif_order
        self.modify = modify
        self.starts = starts

    def discover(
        self,
        foreground: str,
        background: str,
        output_dir: str,
        number_of_motifs: int,
        *args,
        **kwargs,
    ) -> List[GenericModel]:
        length = _require_length(kwargs)
        os.makedirs(output_dir, exist_ok=True)
        jar_path = resolve_existing_path(
            self.jar_path,
            "HORDEMOTIFS_SLIM_JAR",
            DEFAULT_SLIM_JAR,
            "SlimDimont jar",
        )
        java_command = resolve_command(self.java_command)
        data_name = "train.annot.fa"
        write_jstacs_fasta(
            foreground,
            os.path.join(output_dir, data_name),
            self.position_tag,
            self.value_tag,
        )
        cmd = _build_slim_args(
            java_command,
            self.java_xmx,
            jar_path,
            output_dir,
            data_name,
            length,
            self.position_tag,
            self.value_tag,
            self.bg_order,
            self.motif_order,
            self.modify,
            self.starts,
            self.threads,
        )
        run_checked(cmd)
        paths = _existing_xml_paths(output_dir, ["*slim*.xml", "*Slim*.xml", "*.xml"])
        return _read_xml_models(paths, "slim", "Slim", number_of_motifs)


class SitegaDiscoveryTool(MotifDiscoveryTool):
    """Wrapper around SiteGA motif discovery."""

    def __init__(self, nmotifs: int = 5) -> None:
        super().__init__(name="sitega")
        self.nmotifs = nmotifs

    def discover(
        self,
        foreground: str,
        background: str,
        output_dir: str,
        number_of_motifs: int,
        *args,
        **kwargs,
    ) -> List[GenericModel]:
        length = _require_length(kwargs)
        number_of_lpd = int(kwargs.get("lpd", 20))
        os.makedirs(output_dir, exist_ok=True)
        shutil.copyfile(foreground, os.path.join(output_dir, "train.fa"))
        shutil.copyfile(background, os.path.join(output_dir, "background.fa"))
        run_sitega(output_dir, length, number_of_lpd, number_of_motifs)

        motifs: List[GenericModel] = []
        for index, path in enumerate(sorted(glob.glob(os.path.join(output_dir, "train.fa_mat*"))), start=1):
            motif = read_model(path, "sitega")
            motif.name = f"Sitega-{index}"
            motifs.append(motif)
        return motifs
