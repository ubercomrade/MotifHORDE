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
                "message": "",
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


def _write_manifest(path, models, *, status="success", message=""):
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": status,
                "message": message,
                "models": models,
            }
        )
    )


def _entry(model_file, *, name="Model", model_type="sitega", length=4):
    return {
        "model_file": model_file,
        "model_type": model_type,
        "length": length,
        "name": name,
    }


def test_sitega_manifest_validates_entries_before_limit(tmp_path):
    output = tmp_path / "out"
    output.mkdir()
    _sitega_mat(output / "first.mat")
    manifest = output / "manifest.json"
    _write_manifest(
        manifest,
        [_entry("first.mat"), _entry("../escape.mat")],
    )

    with pytest.raises(ValueError, match="traversal"):
        _read_sitega_manifest(manifest, os.fspath(output), 4, 1)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("schema_version", 2),
        ("status", "unknown"),
        ("message", 1),
        ("models", {}),
    ],
)
def test_sitega_manifest_rejects_invalid_root_fields(tmp_path, field, value):
    output = tmp_path / "out"
    output.mkdir()
    manifest = output / "manifest.json"
    data = {
        "schema_version": 1,
        "status": "success",
        "message": "",
        "models": [],
    }
    data[field] = value
    _write_manifest(
        manifest, data["models"], status=data["status"], message=data["message"]
    )
    raw = json.loads(manifest.read_text())
    raw["schema_version"] = data["schema_version"]
    raw["models"] = data["models"]
    raw[field] = value
    manifest.write_text(json.dumps(raw))

    with pytest.raises(ValueError):
        _read_sitega_manifest(manifest, os.fspath(output), 4, 1)


def test_sitega_manifest_failure_requires_message(tmp_path):
    output = tmp_path / "out"
    output.mkdir()
    manifest = output / "manifest.json"
    _write_manifest(manifest, [], status="failure", message="training failed")

    with pytest.raises(RuntimeError, match="training failed"):
        _read_sitega_manifest(manifest, os.fspath(output), 4, 1)


@pytest.mark.parametrize(
    "entry",
    [
        _entry("missing.mat"),
        _entry("missing.mat", model_type="pwm"),
        _entry("missing.mat", length=5),
    ],
)
def test_sitega_manifest_rejects_invalid_model_entries(tmp_path, entry):
    output = tmp_path / "out"
    output.mkdir()
    manifest = output / "manifest.json"
    _write_manifest(manifest, [entry])

    with pytest.raises((ValueError, FileNotFoundError)):
        _read_sitega_manifest(manifest, os.fspath(output), 4, 1)


def test_sitega_manifest_rejects_symlink(tmp_path):
    output = tmp_path / "out"
    output.mkdir()
    target = tmp_path / "target.mat"
    _sitega_mat(target)
    (output / "link.mat").symlink_to(target)
    manifest = output / "manifest.json"
    _write_manifest(manifest, [_entry("link.mat")])

    with pytest.raises(ValueError, match="symlink"):
        _read_sitega_manifest(manifest, os.fspath(output), 4, 1)
