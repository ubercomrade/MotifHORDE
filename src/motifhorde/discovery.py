"""Subprocess adapters for external motif discovery tools."""

from __future__ import annotations

import glob
import json
import logging
import os
from abc import ABC, abstractmethod
from dataclasses import replace
from pathlib import Path

from mimosa import MotifModel, read_model
from mimosa.errors import ModelFormatError

from .external import (
    DEFAULT_BAMM_COMMAND,
    DEFAULT_MEME_COMMAND,
    DEFAULT_SITEGA_COMMAND,
    DEFAULT_STREME_COMMAND,
    resolve_command,
    resolve_existing_path,
    run_checked,
)
from .io import write_jstacs_fasta

logger = logging.getLogger(__name__)


class MotifDiscoveryTool(ABC):
    """Common interface for external discovery adapters."""

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
    ) -> list[MotifModel]:
        raise NotImplementedError


def _require_length(kwargs: dict) -> int:
    length = kwargs.get("length")
    if length is None:
        raise ValueError("Parameter 'length' is required for discovery")
    length = int(length)
    if length < 1:
        raise ValueError("motif length must be positive")
    return length


def _require_number_of_motifs(number_of_motifs: int) -> int:
    if (
        isinstance(number_of_motifs, bool)
        or not isinstance(number_of_motifs, int)
        or number_of_motifs < 1
    ):
        raise ValueError("number_of_motifs must be positive")
    return number_of_motifs


def _check_length(model: MotifModel, expected_length: int | None) -> None:
    if expected_length is not None and model.motif_length != expected_length:
        raise ValueError(
            f"{type(model).__name__} model {model.name!r} has length "
            f"{model.motif_length}, expected {expected_length}"
        )


def _rename(model: MotifModel, prefix: str, index: int) -> MotifModel:
    return replace(model, name=f"{prefix}-{index}")


def _read_indexed_models(
    path: str,
    model_format: str,
    prefix: str,
    number_of_motifs: int,
    expected_length: int | None = None,
    **kwargs,
) -> list[MotifModel]:
    _require_number_of_motifs(number_of_motifs)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"Discovery output is missing or empty: {path}")

    motifs: list[MotifModel] = []
    for index in range(number_of_motifs):
        try:
            model = read_model(path, format=model_format, index=index, **kwargs)
        except IndexError:
            break
        except ModelFormatError as exc:
            if index > 0 and "motif index" in str(exc) and "out of range" in str(exc):
                break
            raise
        _check_length(model, expected_length)
        motifs.append(_rename(model, prefix, index + 1))
    if not motifs:
        raise RuntimeError(f"No models found in discovery output: {path}")
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
    output = os.path.join(output_dir, "motifs.meme")
    executable = resolve_command(
        DEFAULT_STREME_COMMAND,
        explicit=command,
        env_var="HORDEMOTIFS_STREME_COMMAND",
    )
    result = run_checked(
        [
            executable,
            "--p",
            foreground,
            "--n",
            background,
            "--objfun",
            "de",
            "--minw",
            str(length),
            "--maxw",
            str(length),
            "-nmotifs",
            str(number_of_motifs),
            "--text",
        ]
    )
    if not result.stdout.strip():
        raise RuntimeError("STREME produced no MEME output")
    Path(output).write_text(result.stdout, encoding="utf-8")
    return output


def _existing_xml_paths(output_dir: str, patterns: list[str]) -> list[str]:
    for pattern in patterns:
        paths = sorted(glob.glob(os.path.join(output_dir, pattern)))
        if paths:
            return paths
    raise FileNotFoundError(f"No XML model output found in {output_dir}")


def _read_xml_models(
    paths: list[str],
    model_format: str,
    prefix: str,
    number_of_motifs: int,
    expected_length: int | None = None,
) -> list[MotifModel]:
    _require_number_of_motifs(number_of_motifs)
    motifs: list[MotifModel] = []
    for path in paths:
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise FileNotFoundError(f"Discovery output is missing or empty: {path}")
        model = read_model(path, format=model_format)
        _check_length(model, expected_length)
        motifs.append(_rename(model, prefix, len(motifs) + 1))
        if len(motifs) == number_of_motifs:
            break
    if not motifs:
        raise RuntimeError(f"No models found in XML discovery output: {paths}")
    return motifs


class StremeDiscoveryTool(MotifDiscoveryTool):
    def __init__(self, nmotifs: int = 5, command: str | None = None) -> None:
        super().__init__("streme")
        self.nmotifs = nmotifs
        self.command = command

    def discover(
        self, foreground, background, output_dir, number_of_motifs, *args, **kwargs
    ):
        _require_number_of_motifs(number_of_motifs)
        length = _require_length(kwargs)
        output = _run_streme(
            foreground, background, output_dir, length, number_of_motifs, self.command
        )
        return _read_indexed_models(
            output, "meme", "Streme", number_of_motifs, expected_length=length
        )


class MemeDiscoveryTool(MotifDiscoveryTool):
    def __init__(
        self,
        command: str | None = None,
        objfun: str = "classic",
        model: str = "zoops",
        minsites: int | None = None,
        maxsites: int | None = None,
        seed: int | None = None,
        threads: int | None = 1,
    ) -> None:
        super().__init__("meme")
        self.command = command
        self.objfun = objfun
        self.model = model
        self.minsites = minsites
        self.maxsites = maxsites
        self.seed = seed
        self.threads = threads

    def discover(
        self, foreground, background, output_dir, number_of_motifs, *args, **kwargs
    ):
        _require_number_of_motifs(number_of_motifs)
        length = _require_length(kwargs)
        os.makedirs(output_dir, exist_ok=True)
        output = os.path.join(output_dir, "motifs.meme")
        executable = resolve_command(
            DEFAULT_MEME_COMMAND,
            explicit=self.command,
            env_var="HORDEMOTIFS_MEME_COMMAND",
        )
        command = [
            executable,
            foreground,
            "-dna",
            "-revcomp",
            "-nmotifs",
            str(number_of_motifs),
            "-minw",
            str(length),
            "-maxw",
            str(length),
            "-nomatrim",
            "-oc",
            output_dir,
            "-objfun",
            self.objfun,
            "-mod",
            self.model,
        ]
        if self.objfun != "classic":
            command.extend(["-neg", background])
        for flag, value in (
            ("-minsites", self.minsites),
            ("-maxsites", self.maxsites),
            ("-seed", self.seed),
            ("-p", self.threads),
        ):
            if value is not None:
                command.extend([flag, str(value)])

        result = run_checked(command)
        report = os.path.join(output_dir, "meme.txt")
        if os.path.isfile(report) and os.path.getsize(report) > 0:
            output = report
        elif result.stdout.strip():
            Path(output).write_text(result.stdout, encoding="utf-8")
        else:
            raise RuntimeError("MEME produced no model output")
        return _read_indexed_models(
            output, "meme", "Meme", number_of_motifs, expected_length=length
        )


class BammDiscoveryTool(MotifDiscoveryTool):
    def __init__(
        self, bamm_command: str | None = None, streme_command: str | None = None
    ):
        super().__init__("bamm")
        self.bamm_command = bamm_command
        self.streme_command = streme_command

    def discover(
        self, foreground, background, output_dir, number_of_motifs, *args, **kwargs
    ):
        _require_number_of_motifs(number_of_motifs)
        length = _require_length(kwargs)
        order = int(kwargs.get("order", 2))
        pwm_output = _run_streme(
            foreground,
            background,
            output_dir,
            length,
            number_of_motifs,
            self.streme_command,
        )
        executable = resolve_command(
            DEFAULT_BAMM_COMMAND,
            explicit=self.bamm_command,
            env_var="HORDEMOTIFS_BAMM_COMMAND",
        )
        run_checked(
            [
                executable,
                output_dir,
                foreground,
                "--PWMFile",
                pwm_output,
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
        )

        motifs = []
        for index in range(1, number_of_motifs + 1):
            path = os.path.join(output_dir, f"bamm_motif_{index}.ihbcp")
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise FileNotFoundError(f"BaMM output is missing or empty: {path}")
            model = read_model(path, format="bamm", order=order)
            _check_length(model, length)
            motifs.append(_rename(model, "Bamm", index))
        return motifs


def _build_jstacs_args(
    java_command: str,
    java_xmx: str,
    jar_path: str,
    output_dir: str,
    data_name: str,
    tool_name: str,
    length: int,
    position_tag: str,
    value_tag: str,
    bg_order: int,
    motif_order: int,
    starts: int,
    threads: int | None,
    extra: list[str] | None = None,
) -> list[str]:
    args = [
        java_command,
        "-Djava.awt.headless=true",
        f"-Xmx{java_xmx}",
        "-jar",
        jar_path,
        f"home={output_dir}",
        f"data={data_name}",
        f"infix={tool_name}",
        f"position={position_tag}",
        f"value={value_tag}",
        f"motifWidth={length}",
        f"motifOrder={motif_order}",
        f"bgOrder={bg_order}",
        f"starts={starts}",
    ]
    if extra:
        args.extend(extra)
    if threads is not None:
        args.append(f"threads={threads}")
    return args


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
    return _build_jstacs_args(
        java_command,
        java_xmx,
        jar_path,
        output_dir,
        data_name,
        "dimont",
        length,
        position_tag,
        value_tag,
        bg_order,
        motif_order,
        starts,
        threads,
        extra=[f"ess={ess}"],
    )


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
    return _build_jstacs_args(
        java_command,
        java_xmx,
        jar_path,
        output_dir,
        data_name,
        "slim",
        length,
        position_tag,
        value_tag,
        bg_order,
        motif_order,
        starts,
        threads,
        extra=[f"modify={str(modify).lower()}"],
    )


class DimontDiscoveryTool(MotifDiscoveryTool):
    def __init__(
        self,
        jar_path: str | None = None,
        java_command: str = "java",
        java_xmx: str = "4G",
        threads: int | None = 1,
        position_tag: str = "position",
        value_tag: str = "value",
        bg_order: int = -1,
        motif_order: int = 0,
        ess: float = 4.0,
        starts: int = 20,
    ) -> None:
        super().__init__("dimont")
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
        self, foreground, background, output_dir, number_of_motifs, *args, **kwargs
    ):
        _require_number_of_motifs(number_of_motifs)
        length = _require_length(kwargs)
        os.makedirs(output_dir, exist_ok=True)
        jar_path = resolve_existing_path(
            self.jar_path, "HORDEMOTIFS_DIMONT_JAR", "Dimont jar"
        )
        java = resolve_command(self.java_command)
        data_name = "train.annot.fa"
        write_jstacs_fasta(
            foreground,
            os.path.join(output_dir, data_name),
            self.position_tag,
            self.value_tag,
        )
        run_checked(
            _build_jstacs_args(
                java,
                self.java_xmx,
                jar_path,
                output_dir,
                data_name,
                "dimont",
                length,
                self.position_tag,
                self.value_tag,
                self.bg_order,
                self.motif_order,
                self.starts,
                self.threads,
                extra=[f"ess={self.ess}"],
            )
        )
        return _read_xml_models(
            _existing_xml_paths(output_dir, ["*dimont*.xml", "*Dimont*.xml", "*.xml"]),
            "dimont",
            "Dimont",
            number_of_motifs,
            expected_length=length,
        )


class SlimDiscoveryTool(MotifDiscoveryTool):
    def __init__(
        self,
        jar_path: str | None = None,
        java_command: str = "java",
        java_xmx: str = "4G",
        threads: int | None = 1,
        position_tag: str = "position",
        value_tag: str = "value",
        bg_order: int = -1,
        motif_order: int = -5,
        modify: bool = True,
        starts: int = 20,
    ) -> None:
        super().__init__("slim")
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
        self, foreground, background, output_dir, number_of_motifs, *args, **kwargs
    ):
        _require_number_of_motifs(number_of_motifs)
        length = _require_length(kwargs)
        os.makedirs(output_dir, exist_ok=True)
        jar_path = resolve_existing_path(
            self.jar_path, "HORDEMOTIFS_SLIM_JAR", "SlimDimont jar"
        )
        java = resolve_command(self.java_command)
        data_name = "train.annot.fa"
        write_jstacs_fasta(
            foreground,
            os.path.join(output_dir, data_name),
            self.position_tag,
            self.value_tag,
        )
        run_checked(
            _build_slim_args(
                java,
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
        )
        return _read_xml_models(
            _existing_xml_paths(output_dir, ["*slim*.xml", "*Slim*.xml", "*.xml"]),
            "slim",
            "Slim",
            number_of_motifs,
            expected_length=length,
        )


def _manifest_model_path(root: Path, value: object) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError("SiteGA manifest model_file must be a non-empty string")
    if os.path.isabs(value) or "\\" in value or "\x00" in value:
        raise ValueError("SiteGA manifest model_file must be a relative POSIX path")
    parts = value.split("/")
    if any(part in ("", ".", "..") for part in parts):
        raise ValueError("SiteGA manifest model_file contains path traversal")
    for index, _ in enumerate(parts):
        component = root.joinpath(*parts[: index + 1])
        if component.is_symlink():
            raise ValueError("SiteGA manifest model_file contains a symlink")
    candidate = (root / Path(*parts)).resolve()
    try:
        if os.path.commonpath((str(root), str(candidate))) != str(root):
            raise ValueError("SiteGA manifest model_file escapes output directory")
    except ValueError as exc:
        raise ValueError("SiteGA manifest model_file escapes output directory") from exc
    return candidate


def _read_sitega_manifest(
    manifest_path: Path, output_dir: str, expected_length: int, number_of_motifs: int
) -> list[MotifModel]:
    _require_number_of_motifs(number_of_motifs)
    if not manifest_path.is_file():
        raise FileNotFoundError(f"SiteGA manifest is missing: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid SiteGA manifest: {manifest_path}") from exc
    if not isinstance(manifest, dict):
        raise ValueError("SiteGA manifest must be an object")
    required = ("schema_version", "status", "message", "models")
    missing = [key for key in required if key not in manifest]
    if missing:
        raise ValueError(
            "SiteGA manifest is missing required fields: " + ", ".join(missing)
        )
    version = manifest["schema_version"]
    if isinstance(version, bool) or not isinstance(version, int) or version != 1:
        raise ValueError("SiteGA manifest schema_version must be 1")
    status = manifest["status"]
    message = manifest["message"]
    if not isinstance(status, str) or not status:
        raise ValueError("SiteGA manifest status must be a non-empty string")
    if not isinstance(message, str):
        raise ValueError("SiteGA manifest message must be a string")
    entries = manifest["models"]
    if not isinstance(entries, list):
        raise ValueError("SiteGA manifest models must be a list")
    if status == "failure":
        raise RuntimeError(f"SiteGA failed: {message}")
    if status != "success":
        raise ValueError(f"unsupported SiteGA manifest status: {status!r}")
    if not entries:
        raise ValueError("SiteGA manifest models must be a non-empty list")

    root = Path(output_dir).resolve()
    models: list[MotifModel] = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("SiteGA manifest model entry must be an object")
        required = ("model_file", "model_type", "length", "name")
        missing = [key for key in required if key not in entry]
        if missing:
            raise ValueError(
                "SiteGA manifest model is missing required fields: "
                + ", ".join(missing)
            )
        if not isinstance(entry["model_type"], str):
            raise ValueError("SiteGA manifest model_type must be a string")
        if entry["model_type"] != "sitega":
            raise ValueError("SiteGA manifest model_type must be 'sitega'")
        length = entry["length"]
        if (
            isinstance(length, bool)
            or not isinstance(length, int)
            or length != expected_length
        ):
            raise ValueError("SiteGA manifest model length is incompatible")
        name = entry["name"]
        if not isinstance(name, str) or not name:
            raise ValueError("SiteGA manifest model name must be non-empty")
        path = _manifest_model_path(root, entry["model_file"])
        if not path.is_file():
            raise FileNotFoundError(f"SiteGA model is missing: {path}")
        model = read_model(path, format="sitega")
        if model.motif_length != expected_length:
            raise ValueError(
                f"SiteGA model {path} has length {model.motif_length}, expected {expected_length}"
            )
        models.append(replace(model, name=name))
    selected = models[:number_of_motifs]
    if not selected:
        raise RuntimeError("SiteGA manifest contains no requested models")
    return selected


class SitegaDiscoveryTool(MotifDiscoveryTool):
    def __init__(
        self,
        command: str | None = None,
        threads: int | None = 1,
        seed: int | None = None,
    ) -> None:
        super().__init__("sitega")
        self.command = command
        self.threads = threads
        self.seed = seed

    def discover(
        self, foreground, background, output_dir, number_of_motifs, *args, **kwargs
    ):
        _require_number_of_motifs(number_of_motifs)
        length = _require_length(kwargs)
        lpd = int(kwargs.get("lpd", 6))
        features = int(kwargs.get("features", 20))
        if not 1 <= lpd <= 6:
            raise ValueError("SiteGA lpd must be between 1 and 6")
        if features < 1:
            raise ValueError("SiteGA features must be positive")
        os.makedirs(output_dir, exist_ok=True)
        manifest = Path(output_dir) / "manifest.json"
        executable = resolve_command(
            DEFAULT_SITEGA_COMMAND,
            explicit=self.command,
            env_var="HORDEMOTIFS_SITEGA_COMMAND",
        )
        command = [
            executable,
            "--foreground",
            foreground,
            "--background",
            background,
            "--output-dir",
            output_dir,
            "--length",
            str(length),
            "--lpd",
            str(lpd),
            "--features",
            str(features),
            "--number-of-motifs",
            str(number_of_motifs),
            "--manifest",
            os.fspath(manifest),
        ]
        if self.threads is not None:
            command.extend(["--threads", str(self.threads)])
        if kwargs.get("seed", self.seed) is not None:
            command.extend(["--seed", str(kwargs.get("seed", self.seed))])
        run_checked(command)
        return _read_sitega_manifest(manifest, output_dir, length, number_of_motifs)
