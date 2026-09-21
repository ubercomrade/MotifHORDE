from __future__ import annotations

import json
import os
import stat

import pytest

from motifhorde.discovery import SitegaDiscoveryTool, _read_sitega_manifest


def _sitega_mat(path, length=4):
    path.write_text(
        "\n".join(
            [
                "SiteGATest",
                "1\tLPD count",
                f"{length}\tModel length",
                "0.0\tMinimum",
                "1.0\tRazmah",
                "0\t2\t1.5\t0\tac",
            ]
        )
    )


def _fake_sitega(path):
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import json, pathlib, sys\n"
        "args = dict(zip(sys.argv[1::2], sys.argv[2::2]))\n"
        "root = pathlib.Path(args['--output-dir'])\n"
        "root.mkdir(parents=True, exist_ok=True)\n"
        "for name in ('first.mat', 'second.mat'):\n"
        "    (root / name).write_text('SiteGATest\\n1\\tLPD count\\n4\\tModel length\\n0.0\\tMinimum\\n1.0\\tRazmah\\n0\\t2\\t1.5\\t0\\tac\\n')\n"
        "pathlib.Path(args['--manifest']).write_text(json.dumps({'schema_version': 1, 'status': 'success', 'message': '', 'models': [{'model_file': 'second.mat', 'model_type': 'sitega', 'length': 4, 'name': 'Second'}, {'model_file': 'first.mat', 'model_type': 'sitega', 'length': 4, 'name': 'First'}]}))\n"
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_sitega_external_manifest_order_and_command(tmp_path):
    command = tmp_path / "fake-sitega"
    _fake_sitega(command)
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    foreground.write_text(">f\nACGTACGT\n")
    background.write_text(">b\nACGTACGT\n")

    models = SitegaDiscoveryTool(command=os.fspath(command)).discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(tmp_path / "out"),
        1,
        length=4,
        lpd=6,
        features=10,
    )

    assert [(model.name, model.motif_length) for model in models] == [("Second", 4)]


def test_sitega_manifest_rejects_traversal(tmp_path):
    output = tmp_path / "out"
    output.mkdir()
    manifest = output / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": "success",
                "models": [
                    {
                        "model_file": "../escape.mat",
                        "model_type": "sitega",
                        "length": 4,
                        "name": "x",
                    }
                ],
            }
        )
    )

    with pytest.raises(ValueError, match="traversal"):
        _read_sitega_manifest(manifest, os.fspath(output), 4, 1)
