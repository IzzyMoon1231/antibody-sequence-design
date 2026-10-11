"""CPU-only contract regression; no model/tokenizer construction is needed.

UIRJ_SERIALIZER_BASELINE=1 selects the conversion block from the pinned Git
revision, not a rewritten approximation. The oracle is the official HF export
function extracted with AST to avoid importing/initializing the folding model.
"""
import ast
import hashlib
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
import torch
import transformers
from transformers.models.esm.openfold_utils import residue_constants as rc
from transformers.models.esm.openfold_utils.data_transforms import make_atom14_masks
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
from transformers.models.esm.openfold_utils.protein import Protein as OFProtein, to_pdb

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_esmfold_predictions as runner

BASE = "3be2136e0863b705aee42715a681600b64530ae4"
# Source hashes from official Transformers v5.18.0, commit a906d3c4b65095f2308b6a6a193e934d03b8eb5d.
HELPER_HASHES = {
    "feats.py": "40262ea6c548368e6326ec18937f137a335890f9469ba29f735629354c6923cb",
    "tensor_utils.py": "fcb02e726101ff8995258d680082f596cff9b893696d262ff49d0f34303a350f",
    "data_transforms.py": "43927f06927ff056b97d9f273fa90f079a29b3e638332752424c19fbfc8c0c10",
    "protein.py": "300ea24a17e5c24529651c0479f993e392e5e4a3d4f2527622e24fc87d8dfe62",
}


def fixture_output(batch=1):
    # Axes: structure-module stage, batch, residue, atom14, xyz.
    # The stage axis is not the number of outer recycling iterations.
    aatype = torch.tensor([[rc.restype_order[x] for x in "GAVW"]] * batch, device="cpu")
    output = make_atom14_masks({"aatype": aatype})
    output["positions"] = torch.arange(2 * batch * 4 * 14 * 3, device="cpu", dtype=torch.float32).reshape(2, batch, 4, 14, 3) / 10
    output["residue_index"] = torch.tensor([[0, 2, 5, 9]] * batch, device="cpu")
    output["plddt"] = torch.arange(batch * 4 * 37, device="cpu", dtype=torch.float32).reshape(batch, 4, 37) / 10 + 20
    return output


def official_oracle():
    esm = Path(transformers.__file__).parent / "models" / "esm"
    for name, digest in HELPER_HASHES.items():
        assert hashlib.sha256((esm / "openfold_utils" / name).read_bytes()).hexdigest() == digest, f"Unverified oracle helper: {name}"
    source = Path(os.environ.get("UIRJ_ESMFOLD_ORACLE", str(esm / "modeling_esmfold.py")))
    tree = ast.parse(source.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "EsmForProteinFolding")
    function = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "output_to_pdb")
    function.decorator_list = []
    # Pinned AST digest also permits the identical installed 5.8.1 export.
    assert hashlib.sha256(ast.dump(function).encode()).hexdigest() == ORACLE_AST_SHA256
    namespace = dict(atom14_to_atom37=atom14_to_atom37, OFProtein=OFProtein, to_pdb=to_pdb)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])), str(source), "exec"), namespace)
    return namespace["output_to_pdb"]


def original_serializer():
    source = subprocess.check_output(["git", "show", f"{BASE}:scripts/run_esmfold_predictions.py"], cwd=ROOT, text=True)
    tree = ast.parse(source)
    predict = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "predict")
    start = next(i for i, n in enumerate(predict.body) if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "atom37_positions" for t in n.targets))
    fn = ast.parse("def original(output):\n    pass\n").body[0]
    fn.body = predict.body[start:]
    namespace = dict(OFProtein=OFProtein, to_pdb=to_pdb)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), f"{BASE}:predict conversion", "exec"), namespace)
    return lambda output: namespace["original"](SimpleNamespace(**output))


def serializer(output, batch_index=0):
    if os.environ.get("UIRJ_SERIALIZER_BASELINE") == "1":
        assert batch_index == 0
        return original_serializer()(output)
    return runner.serialize_esmfold_output(output, batch_index=batch_index)


def atoms(pdb):
    return {(line[21], int(line[22:26]), line[17:20].strip(), line[12:16].strip()):
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[60:66])])
            for line in pdb.splitlines() if line.startswith("ATOM  ")}


def assert_atoms(actual, expected):
    left, right = atoms(actual), atoms(expected)
    assert left.keys() == right.keys(), "Atom identities/residue numbering differ"
    for identity in right:
        np.testing.assert_allclose(left[identity], right[identity], rtol=0, atol=0.005, err_msg=str(identity))


def test_realistic_atom14_export_matches_official_oracle(tmp_path):
    output = fixture_output()
    expected = official_oracle()(output)[0]
    actual = serializer(output)
    path = tmp_path / "fixture.pdb"
    path.write_text(actual)
    assert runner.validate_sequence("GAVW", runner.extract_pdb_sequence(path))[0] == "PASS"
    assert_atoms(actual, expected)
    assert actual == expected
    assert len(atoms(actual)) == 30  # GLY 4 + ALA 5 + VAL 7 + TRP 14.
    assert set(key[1] for key in atoms(actual)) == {1, 3, 6, 10}


def test_missing_atom_and_padded_batch_match_oracle():
    output = fixture_output(batch=2)
    output["atom37_atom_exists"][0, 3, rc.atom_order["CH2"]] = 0
    output["atom37_atom_exists"][1, 2:] = 0
    expected = official_oracle()(output)
    for i in range(2):
        actual = serializer(output, batch_index=i)
        assert_atoms(actual, expected[i])
    assert ("A", 10, "TRP", "CH2") not in atoms(serializer(output))
    assert len(atoms(serializer(output, batch_index=1))) == 9


def test_sequence_correct_atom_wrong_control_is_rejected(tmp_path):
    good = official_oracle()(fixture_output())[0]
    # Alter a side-chain coordinate while preserving CA sequence and all labels.
    bad = "\n".join(line[:30] + " 999.000" + line[38:] if line.startswith("ATOM  ") and line[12:16].strip() == "CH2" else line for line in good.splitlines()) + "\n"
    path = tmp_path / "scrambled.pdb"
    path.write_text(bad)
    assert runner.validate_sequence("GAVW", runner.extract_pdb_sequence(path))[0] == "PASS"
    with pytest.raises(AssertionError):
        assert_atoms(bad, good)


@pytest.mark.parametrize("field,mutation", [
    ("positions", lambda x: x[..., :13, :]),
    ("positions", lambda x: x[-1]),
    ("aatype", lambda x: x[:, :-1]),
    ("atom37_atom_exists", lambda x: x[..., :36]),
    ("residx_atom37_to_atom14", lambda x: x[..., :36]),
    ("residue_index", lambda x: x[:, :-1]),
    ("plddt", lambda x: x[..., :36]),
])
def test_mismatched_shapes_fail_clearly(field, mutation):
    output = fixture_output()
    output[field] = mutation(output[field])
    with pytest.raises(ValueError, match="shape"):
        serializer(output)


@pytest.mark.parametrize("field", ["positions", "plddt", "atom37_atom_exists"])
def test_nonfinite_fields_fail_clearly(field):
    output = fixture_output()
    output[field].reshape(-1)[-1] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        serializer(output)


@pytest.mark.parametrize("value", [-1, 14, 2])
def test_invalid_mapping_rejected_including_in_range_corruption(value):
    output = fixture_output()
    output["residx_atom37_to_atom14"][0, 0, rc.atom_order["N"]] = value
    with pytest.raises(ValueError, match="mapping"):
        serializer(output)


def test_missing_required_field():
    output = fixture_output()
    del output["residx_atom37_to_atom14"]
    with pytest.raises(ValueError, match="Missing"):
        serializer(output)


def test_invalid_batch_index():
    with pytest.raises(ValueError, match="batch_index"):
        serializer(fixture_output(), batch_index=1)

ORACLE_AST_SHA256 = '8d2f694046921d169c1a225add3b41d86f85e54f415c5269f386b3cc55c901ba'


@pytest.mark.parametrize("field", ["aatype", "residue_index", "residx_atom37_to_atom14"])
def test_floating_index_fields_rejected(field):
    output = fixture_output()
    output[field] = output[field].float()
    with pytest.raises(ValueError, match="integer dtype"):
        serializer(output)


@pytest.mark.parametrize("value", [-1, 21])
def test_out_of_range_aatype(value):
    output = fixture_output()
    output["aatype"][0, 0] = value
    with pytest.raises(ValueError, match="aatype"):
        serializer(output)


@pytest.mark.parametrize("value", [0.5, 1])
def test_invalid_atom_mask(value):
    output = fixture_output()
    output["atom37_atom_exists"][0, 0, rc.atom_order["CB"]] = value
    with pytest.raises(ValueError, match="binary subset"):
        serializer(output)


def test_ca_and_sidechain_coordinates_and_plddt_have_known_sentinels():
    output = fixture_output()
    actual = atoms(serializer(output))
    for letter, residue, residue_number in [("G", 0, 1), ("A", 1, 3), ("V", 2, 6), ("W", 3, 10)]:
        resname = rc.restype_1to3[letter]
        for name in ["CA", "O"] + (["CH2"] if letter == "W" else []):
            slot14 = rc.restype_name_to_atom14_names[resname].index(name)
            slot37 = rc.atom_order[name]
            expected_xyz = output["positions"][-1, 0, residue, slot14].numpy()
            expected_b = output["plddt"][0, residue, slot37].item()
            np.testing.assert_allclose(actual[("A", residue_number, resname, name)], np.r_[expected_xyz, expected_b], atol=0.005, rtol=0)
