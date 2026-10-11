"""Storage regressions use temporary synthetic PDBs only."""
import csv
from pathlib import Path
import shutil
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_esmfold_predictions as runner


def pdb_text():
    # The storage contract validates CA sequence, not geometry.
    return "ATOM      1  CA  GLY A   1       1.000   2.000   3.000  1.00 80.00           C  \nEND\n"


def sibling(path, suffix):
    return path.with_name(path.name + suffix)


def assert_clean(path):
    assert not sibling(path, ".lock").exists()
    assert not list(path.parent.glob(".esmfold_tmp_*"))


def test_success_completion_and_no_overwrite(tmp_path):
    path = tmp_path / "nested" / "prediction.pdb"
    runner.save_prediction(pdb_text(), path, "G")
    assert path.read_text() == pdb_text()
    assert sibling(path, ".complete").read_text() == "COMPLETE\n"
    assert runner.check_prediction(path, "G") == "PASS"
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        runner.save_prediction("replacement", path, "G")
    assert path.read_bytes() == before
    assert_clean(path)


def test_existing_lock_retained(tmp_path):
    path = tmp_path / "prediction.pdb"
    lock = sibling(path, ".lock")
    lock.write_text("another writer")
    with pytest.raises(FileExistsError):
        runner.save_prediction(pdb_text(), path, "G")
    assert lock.read_text() == "another writer"
    assert not path.exists()
    assert not list(tmp_path.glob(".esmfold_tmp_*"))


def test_reentrant_writer_excluded(tmp_path, monkeypatch):
    path = tmp_path / "prediction.pdb"
    extract = runner.extract_pdb_sequence
    attempts = []
    def during_validation(source):
        with pytest.raises(FileExistsError):
            runner.save_prediction(pdb_text(), path, "G")
        attempts.append(True)
        return extract(source)
    monkeypatch.setattr(runner, "extract_pdb_sequence", during_validation)
    runner.save_prediction(pdb_text(), path, "G")
    assert len(attempts) == 2
    assert_clean(path)


def test_validation_failure_cleans_temp_and_keeps_destination_absent(tmp_path):
    path = tmp_path / "prediction.pdb"
    with pytest.raises(ValueError, match="failed validation"):
        runner.save_prediction(pdb_text(), path, "A")
    assert not path.exists()
    assert not sibling(path, ".complete").exists()
    assert_clean(path)


def test_orphan_marker_never_overwritten(tmp_path):
    path = tmp_path / "prediction.pdb"
    marker = sibling(path, ".complete")
    marker.write_text("previous evidence")
    assert runner.check_prediction(path, "G") == "ORPHAN_COMPLETION_MARKER"
    with pytest.raises(FileExistsError):
        runner.save_prediction(pdb_text(), path, "G")
    assert marker.read_text() == "previous evidence"
    assert not path.exists()
    assert_clean(path)


def test_interrupted_copy_preserved_without_completion(tmp_path, monkeypatch):
    path = tmp_path / "prediction.pdb"
    def interrupted(source, destination):
        destination.write(b"partial evidence")
        raise OSError("injected copy interruption")
    monkeypatch.setattr(shutil, "copyfileobj", interrupted)
    with pytest.raises(OSError, match="injected"):
        runner.save_prediction(pdb_text(), path, "G")
    assert path.read_bytes() == b"partial evidence"
    assert not sibling(path, ".complete").exists()
    assert runner.check_prediction(path, "G") == "INCOMPLETE_NO_MARKER"
    with pytest.raises(FileExistsError):
        runner.save_prediction(pdb_text(), path, "G")
    assert_clean(path)


@pytest.mark.parametrize("contents,expected", [(None, "INCOMPLETE_NO_MARKER"), ("BAD\n", "INVALID_COMPLETION_MARKER"), ("COMPLETE\n", "PASS")])
def test_restart_classification(tmp_path, contents, expected):
    path = tmp_path / "prediction.pdb"
    assert runner.check_prediction(path, "G") == "MISSING"
    path.write_text(pdb_text())
    if contents is not None:
        sibling(path, ".complete").write_text(contents)
    assert runner.check_prediction(path, "G") == expected


def test_completed_sequence_mismatch_is_not_pass(tmp_path):
    path = tmp_path / "prediction.pdb"
    path.write_text(pdb_text())
    sibling(path, ".complete").write_text("COMPLETE\n")
    assert runner.check_prediction(path, "A") == "SEQUENCE_MISMATCH"


def test_lock_inspection_is_read_only(tmp_path):
    path = tmp_path / "prediction.pdb"
    lock = sibling(path, ".lock")
    lock.write_text("owner evidence")
    manifest = tmp_path / "manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb_path", "sequence"])
        writer.writeheader()
        writer.writerow(dict(pdb_path=path, sequence="G"))
    before = lock.stat().st_mtime_ns
    findings = runner.inspect_prediction_locks(manifest, minimum_age_minutes=0)
    assert findings[0]["review_recommended"] is True
    assert findings[0]["prediction_status"] == "MISSING"
    assert lock.read_text() == "owner evidence"
    assert lock.stat().st_mtime_ns == before


def test_destination_created_during_validation_is_preserved(tmp_path, monkeypatch):
    path = tmp_path / "prediction.pdb"
    extract = runner.extract_pdb_sequence
    def competing_writer(source):
        path.write_text("noncooperating writer evidence")
        return extract(source)
    monkeypatch.setattr(runner, "extract_pdb_sequence", competing_writer)
    with pytest.raises(FileExistsError):
        runner.save_prediction(pdb_text(), path, "G")
    assert path.read_text() == "noncooperating writer evidence"
    assert not sibling(path, ".complete").exists()
    assert_clean(path)


def test_copied_validation_failure_preserves_copy_without_marker(tmp_path, monkeypatch):
    path = tmp_path / "prediction.pdb"
    extract = runner.extract_pdb_sequence
    def inconsistent_copy(source):
        return "A" if Path(source) == path else extract(source)
    monkeypatch.setattr(runner, "extract_pdb_sequence", inconsistent_copy)
    with pytest.raises(ValueError, match="Copied PDB failed"):
        runner.save_prediction(pdb_text(), path, "G")
    assert path.read_text() == pdb_text()
    assert not sibling(path, ".complete").exists()
    assert_clean(path)


@pytest.mark.parametrize("completed", [True, False])
def test_cli_restart_skips_complete_and_blocks_incomplete_without_model(tmp_path, monkeypatch, capsys, completed):
    path = tmp_path / "prediction.pdb"
    path.write_text(pdb_text())
    if completed:
        sibling(path, ".complete").write_text("COMPLETE\n")
    manifest = tmp_path / "manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["target", "method", "index", "sequence", "pdb_path"])
        writer.writeheader()
        writer.writerow(dict(target="9NH7", method="baseline", index=0, sequence="G", pdb_path=path))
    monkeypatch.setattr(sys, "argv", ["runner", "--manifest", str(manifest)])
    def forbidden(*args, **kwargs):
        pytest.fail("Restart path reached model/device/prediction code")
    monkeypatch.setattr(runner, "load_model", forbidden)
    monkeypatch.setattr(runner, "predict", forbidden)
    monkeypatch.setattr(runner.torch.cuda, "is_available", forbidden)
    if completed:
        runner.main()
        assert "already completed" in capsys.readouterr().out
    else:
        with pytest.raises(RuntimeError, match="Investigate before continuing"):
            runner.main()
    assert path.read_text() == pdb_text()
