from __future__ import annotations

import os
import stat
import subprocess

import pytest

from motifhorde.discovery import (
    BammDiscoveryTool,
    DimontDiscoveryTool,
    MemeDiscoveryTool,
    SlimDiscoveryTool,
    StremeDiscoveryTool,
    _build_jstacs_args,
)
from motifhorde.external import resolve_command, resolve_existing_path, run_checked
from motifhorde.io import write_jstacs_fasta


def _executable(path, source):
    path.write_text("#!/usr/bin/env python3\n" + source)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def _meme_text(*names, length=4):
    rows = "\n".join("0.25 0.25 0.25 0.25" for _ in range(length))
    return "\n".join(
        [
            "MEME version 4",
            "",
            "ALPHABET= ACGT",
            "",
            "strands: + -",
            "",
            *sum(
                (
                    [
                        f"MOTIF {name}",
                        f"letter-probability matrix: alength= 4 w= {length} nsites= 20 E= 0",
                        rows,
                        "",
                    ]
                    for name in names
                ),
                [],
            ),
        ]
    )


def _inputs(tmp_path):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    foreground.write_text(">f\nACGTACGT\n")
    background.write_text(">b\nTGCATGCA\n")
    return foreground, background


def _dimont_xml():
    return (
        "<root><MarkovModelDiffSM><bayesianNetworkSF><trees>"
        "<pos><parameterTree><root><treeElement><contextPos>0</contextPos>"
        "<pars><pos><parameter><symbol>0</symbol><value>0.25</value></parameter>"
        "</pos><pos><parameter><symbol>1</symbol><value>0.25</value></parameter>"
        "</pos><pos><parameter><symbol>2</symbol><value>0.25</value></parameter>"
        "</pos><pos><parameter><symbol>3</symbol><value>0.25</value></parameter>"
        "</pos></pars></treeElement></root></parameterTree></pos>"
        "</trees></bayesianNetworkSF></MarkovModelDiffSM></root>"
    )


def _slim_xml():
    return (
        "<root><SLIM><length>1</length>"
        "<componentMixtureParameters><pos><pos>1</pos></pos>"
        "</componentMixtureParameters>"
        "<ancestorMixtureParameters><pos><pos><pos>1</pos></pos></pos>"
        "</ancestorMixtureParameters>"
        "<dependencyParameters><pos><pos><pos><pos>0.25</pos><pos>0.25</pos>"
        "<pos>0.25</pos><pos>0.25</pos></pos></pos></pos>"
        "</dependencyParameters></SLIM></root>"
    )


def test_meme_reads_standard_output_in_manifest_order(tmp_path):
    foreground, background = _inputs(tmp_path)
    command = _executable(
        tmp_path / "meme-fake",
        "import pathlib, sys\n"
        "args = sys.argv[1:]\n"
        "root = pathlib.Path(args[args.index('-oc') + 1])\n"
        "(root / 'meme.txt').write_text(" + repr(_meme_text("Second", "First")) + ")\n",
    )

    models = MemeDiscoveryTool(command=os.fspath(command)).discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(tmp_path / "out"),
        2,
        length=4,
    )

    assert [(model.name, model.motif_length) for model in models] == [
        ("Meme-1", 4),
        ("Meme-2", 4),
    ]


def test_streme_reads_direct_meme_output(tmp_path):
    foreground, background = _inputs(tmp_path)
    command = _executable(
        tmp_path / "streme-fake",
        "print(" + repr(_meme_text("Only")) + ")\n",
    )

    models = StremeDiscoveryTool(command=os.fspath(command)).discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(tmp_path / "out"),
        1,
        length=4,
    )

    assert [(model.name, model.motif_length) for model in models] == [("Streme-1", 4)]


def test_bamm_runs_streme_then_reads_models(tmp_path):
    foreground, background = _inputs(tmp_path)
    streme = _executable(
        tmp_path / "streme-fake",
        "print(" + repr(_meme_text("PWM")) + ")\n",
    )
    bamm = _executable(
        tmp_path / "bamm-fake",
        "import pathlib, sys\n"
        "root = pathlib.Path(sys.argv[1])\n"
        "(root / 'bamm_motif_1.ihbcp').write_text('\\n\\n'.join(['0.25 0.25 0.25 0.25'] * 4))\n",
    )

    models = BammDiscoveryTool(
        bamm_command=os.fspath(bamm), streme_command=os.fspath(streme)
    ).discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(tmp_path / "out"),
        1,
        length=4,
        order=0,
    )

    assert [(model.name, model.motif_length) for model in models] == [("Bamm-1", 4)]


@pytest.mark.parametrize(
    ("tool_cls", "infix", "xml"),
    [
        (DimontDiscoveryTool, "dimont", _dimont_xml()),
        (SlimDiscoveryTool, "slim", _slim_xml()),
    ],
    ids=["dimont", "slim"],
)
def test_jstacs_adapters_use_fake_process_and_mimosa_reader(
    tmp_path, tool_cls, infix, xml
):
    foreground, background = _inputs(tmp_path)
    jar = tmp_path / f"{infix}.jar"
    jar.touch()
    java = _executable(
        tmp_path / "java-fake",
        "import pathlib, sys\n"
        "args = dict(item.split('=', 1) for item in sys.argv[1:] if '=' in item)\n"
        "home = pathlib.Path(args['home'])\n"
        "home.joinpath(args['infix'] + '.xml').write_text(" + repr(xml) + ")\n",
    )

    models = tool_cls(jar_path=os.fspath(jar), java_command=os.fspath(java)).discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(tmp_path / "out"),
        1,
        length=1,
    )

    assert len(models) == 1
    assert models[0].motif_length == 1


def test_discovery_rejects_wrong_meme_length(tmp_path):
    foreground, background = _inputs(tmp_path)
    command = _executable(
        tmp_path / "meme-fake",
        "import pathlib, sys\n"
        "root = pathlib.Path(sys.argv[sys.argv.index('-oc') + 1])\n"
        "(root / 'meme.txt').write_text(" + repr(_meme_text("Wrong", length=5)) + ")\n",
    )

    with pytest.raises(ValueError, match="expected 4"):
        MemeDiscoveryTool(command=os.fspath(command)).discover(
            os.fspath(foreground),
            os.fspath(background),
            os.fspath(tmp_path / "out"),
            1,
            length=4,
        )


def test_run_checked_reports_process_details(tmp_path):
    command = _executable(
        tmp_path / "fail",
        "import sys\nprint('out')\nprint('err', file=sys.stderr)\nsys.exit(3)\n",
    )

    with pytest.raises(RuntimeError, match=r"Return code: 3") as error:
        run_checked([os.fspath(command)])
    assert "out" in str(error.value)
    assert "err" in str(error.value)


def test_resolvers_use_cli_then_environment_then_path(tmp_path, monkeypatch):
    cli = tmp_path / "cli"
    env = tmp_path / "env"
    path = tmp_path / "path"
    for value in (cli, env, path):
        value.touch()
    monkeypatch.setenv("HORDEMOTIFS_TEST_COMMAND", os.fspath(env))
    monkeypatch.setattr("motifhorde.external.shutil.which", lambda _: os.fspath(path))

    assert resolve_command(
        "tool", os.fspath(cli), "HORDEMOTIFS_TEST_COMMAND"
    ) == os.fspath(cli)
    assert resolve_command("tool", env_var="HORDEMOTIFS_TEST_COMMAND") == os.fspath(env)
    monkeypatch.delenv("HORDEMOTIFS_TEST_COMMAND")
    assert resolve_command("tool", env_var="HORDEMOTIFS_TEST_COMMAND") == os.fspath(
        path
    )
    assert resolve_existing_path(
        os.fspath(cli), "HORDEMOTIFS_TEST_JAR", "jar"
    ) == os.fspath(cli)
    monkeypatch.setenv("HORDEMOTIFS_TEST_JAR", os.fspath(env))
    assert resolve_existing_path(None, "HORDEMOTIFS_TEST_JAR", "jar") == os.fspath(env)


def test_jstacs_fasta_and_single_internal_thread(tmp_path):
    source = tmp_path / "input.fa"
    output = tmp_path / "annotated.fa"
    source.write_text(">first description\nACGT\n>second\nTGCA\n")
    write_jstacs_fasta(os.fspath(source), os.fspath(output))
    lines = output.read_text().splitlines()
    assert lines[0] == "> position: 2; value: 1.0"
    assert lines[2] == "> position: 2; value: 2.0"

    args = _build_jstacs_args(
        "java",
        "4G",
        "tool.jar",
        "out",
        "train.annot.fa",
        "dimont",
        8,
        "position",
        "value",
        -1,
        0,
        20,
        1,
    )
    assert args[-1] == "threads=1"


def test_run_checked_does_not_use_shell(monkeypatch):
    seen = {}

    def fake_run(args, **kwargs):
        seen.update(kwargs)
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr("motifhorde.external.subprocess.run", fake_run)
    run_checked(["echo", "$(touch should-not-run)"])
    assert seen["shell"] is False
