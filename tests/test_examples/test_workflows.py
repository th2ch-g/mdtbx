import os
import re
import subprocess
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).parents[2]
EXAMPLE_ROOT = REPOSITORY_ROOT / "example"
WORKFLOW_ROOT = EXAMPLE_ROOT / "workflows"
WORKFLOW_SCRIPTS = (
    WORKFLOW_ROOT / "solution_setup.sh",
    WORKFLOW_ROOT / "membrane_setup.sh",
    WORKFLOW_ROOT / "run_slurm.sh",
    WORKFLOW_ROOT / "analyze.sh",
)
MDP_STAGES = ("mini", "eq1", "eq2", "eq3", "eq4", "eq5", "eq6", "prd")


def run_check(script: Path, working_directory: Path, environment: dict[str, str]):
    merged_environment = dict(os.environ)
    merged_environment.update(environment)
    return subprocess.run(
        ["bash", str(script), "--check"],
        cwd=working_directory,
        env=merged_environment,
        check=True,
        capture_output=True,
        text=True,
    )


def make_mdp_directory(path: Path):
    path.mkdir()
    for stage in MDP_STAGES:
        (path / f"{stage}.mdp").write_text("integrator = md\n")


def test_all_shell_examples_have_valid_bash_syntax():
    scripts = sorted(EXAMPLE_ROOT.rglob("*.sh"))
    assert scripts
    subprocess.run(
        ["bash", "-n", *map(str, scripts)],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("script", WORKFLOW_SCRIPTS)
def test_canonical_workflows_are_copy_ready(script: Path):
    text = script.read_text()

    assert script.stat().st_mode & 0o111
    assert "set -euo pipefail" in text
    assert "--check" in text
    assert re.search(r"^\s*(?:sbatch|pjsub)\b", text, re.MULTILINE) is None
    assert re.search(r"^\s*rm\b", text, re.MULTILINE) is None
    assert re.search(r"pixi\s+run\s+gmx.*\bmdrun\b", text) is None


def test_solution_setup_check_is_side_effect_free(tmp_path: Path):
    input_pdb = tmp_path / "prepared.pdb"
    input_pdb.write_text("END\n")
    mdp_directory = tmp_path / "mdp"
    make_mdp_directory(mdp_directory)
    output_directory = tmp_path / "run"

    result = run_check(
        WORKFLOW_ROOT / "solution_setup.sh",
        tmp_path,
        {
            "INPUT_PDB": str(input_pdb),
            "MDP_SOURCE_DIR": str(mdp_directory),
            "OUTPUT_DIR": str(output_directory),
        },
    )

    assert "build_solution" in result.stdout
    assert not output_directory.exists()


def test_membrane_setup_check_is_side_effect_free(tmp_path: Path):
    input_pdb = tmp_path / "oriented.pdb"
    input_pdb.write_text("END\n")
    mdp_directory = tmp_path / "mdp"
    make_mdp_directory(mdp_directory)
    output_directory = tmp_path / "run"

    result = run_check(
        WORKFLOW_ROOT / "membrane_setup.sh",
        tmp_path,
        {
            "INPUT_PDB": str(input_pdb),
            "MDP_SOURCE_DIR": str(mdp_directory),
            "OUTPUT_DIR": str(output_directory),
        },
    )

    assert "packmol-memgen" in result.stdout
    assert not output_directory.exists()


def test_membrane_restraint_macros_match_generated_prefixes():
    script = (WORKFLOW_ROOT / "membrane_setup.sh").read_text()
    assert '"${OUTPUT_DIR}/posres"' in script
    assert '"${OUTPUT_DIR}/posres_lipid"' in script

    for stage in ("mini", "eq1", "eq2", "eq3", "eq4", "eq5", "eq6"):
        mdp = (EXAMPLE_ROOT / "mdp" / "membrane" / f"{stage}.mdp").read_text()
        assert "-DPOSRES_LIPID" in mdp
        assert "-DPOSRES_LIPID_FC=" in mdp
        assert "-DPOSRES_FC_LIPID=" not in mdp


def test_slurm_check_requires_installer_gromacs_without_running_it(tmp_path: Path):
    for filename in (
        "gmx.gro",
        "gmx.top",
        "index.ndx",
        *(f"{stage}.mdp" for stage in MDP_STAGES),
    ):
        (tmp_path / filename).write_text("\n")
    fake_gromacs = tmp_path / "installer-gmx"
    fake_gromacs.write_text("#!/usr/bin/env bash\nexit 99\n")
    fake_gromacs.chmod(0o755)

    result = run_check(
        WORKFLOW_ROOT / "run_slurm.sh",
        tmp_path,
        {
            "GROMACS_BIN": str(fake_gromacs),
            "SYSTEM_DIR": str(tmp_path),
            "PRODUCTION_SEGMENTS": "2",
        },
    )

    assert "prd1..prd2" in result.stdout
    assert not (tmp_path / "mini.tpr").exists()


def test_slurm_starts_eq1_without_a_minimization_checkpoint():
    text = (WORKFLOW_ROOT / "run_slurm.sh").read_text()

    assert "run_stage eq1 eq1.mdp mini.gro" in text
    assert "run_stage eq1 eq1.mdp mini.gro mini.cpt" not in text


def test_slurm_rejects_a_pixi_gromacs_alias(tmp_path: Path):
    for filename in (
        "gmx.gro",
        "gmx.top",
        "index.ndx",
        *(f"{stage}.mdp" for stage in MDP_STAGES),
    ):
        (tmp_path / filename).write_text("\n")
    pixi_gromacs = tmp_path / ".pixi" / "envs" / "default" / "bin" / "gmx"
    pixi_gromacs.parent.mkdir(parents=True)
    pixi_gromacs.write_text("#!/usr/bin/env bash\nexit 99\n")
    pixi_gromacs.chmod(0o755)
    alias = tmp_path / "gmx-alias"
    alias.symlink_to(pixi_gromacs)
    environment = dict(os.environ)
    environment.update(
        {
            "GROMACS_BIN": str(alias),
            "SYSTEM_DIR": str(tmp_path),
        }
    )

    result = subprocess.run(
        ["bash", str(WORKFLOW_ROOT / "run_slurm.sh"), "--check"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Use installer GROMACS" in result.stderr


def test_analysis_check_is_side_effect_free(tmp_path: Path):
    run_directory = tmp_path / "run"
    run_directory.mkdir()
    (run_directory / "index.ndx").write_text("\n")
    for segment in (1, 2):
        (run_directory / f"prd{segment}.tpr").write_text("\n")
        (run_directory / f"prd{segment}.xtc").write_text("\n")

    result = run_check(
        WORKFLOW_ROOT / "analyze.sh",
        tmp_path,
        {
            "RUN_DIR": str(run_directory),
            "OUTPUT_DIR": "analysis",
            "PRODUCTION_SEGMENTS": "2",
        },
    )

    assert "rmsd, rmsf, contactmap" in result.stdout
    assert not (run_directory / "analysis").exists()


def test_example_catalog_covers_every_top_level_area():
    catalog = (EXAMPLE_ROOT / "README.md").read_text()
    areas = sorted(path.name for path in EXAMPLE_ROOT.iterdir() if path.is_dir())

    for area in areas:
        assert f"]({area}/)" in catalog
