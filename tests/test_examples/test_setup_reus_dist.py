import os
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).parents[2]
SETUP_SCRIPT = REPO_ROOT / "example" / "remd" / "setup_reus_dist.sh"
TEMPLATE_MDP = REPO_ROOT / "example" / "mdp" / "solution" / "re.mdp"
TEMPLATE_PLUMED = REPO_ROOT / "example" / "plumed" / "reus_dist.dat"


def test_setup_reus_dist_generates_plumed_biases(tmp_path):
    topology = tmp_path / "gmx.top"
    index = tmp_path / "index.ndx"
    targets = tmp_path / "targets.txt"
    submit_template = tmp_path / "submit-template.sh"
    fake_gmx = tmp_path / "fake-gmx"

    topology.write_text("[ system ]\nExample\n")
    index.write_text("[ TARGET1 ]\n1\n[ TARGET2 ]\n2\n")
    targets.write_text("init1.gro 1.5\ninit2.gro 2.0\n")
    submit_template.write_text(
        'REPLEX=1000\nN_REPLICA=16\nDEFFNM="grest" # or rest2, reus\n'
    )
    fake_gmx.write_text(
        "#!/bin/bash\n"
        'while [ "$#" -gt 0 ]; do\n'
        '    if [ "$1" = "-o" ]; then\n'
        '        touch "$2"\n'
        "        exit 0\n"
        "    fi\n"
        "    shift\n"
        "done\n"
    )
    fake_gmx.chmod(0o755)

    for name in ("init1.gro", "init2.gro"):
        (tmp_path / name).write_text("Example structure\n")

    env = os.environ.copy()
    env.update(
        {
            "TEMPLATE_MDP": str(TEMPLATE_MDP),
            "TEMPLATE_PLUMED": str(TEMPLATE_PLUMED),
            "SUBMIT_SCRIPT": str(submit_template),
            "TARGETS_FILE": str(targets),
            "TOPOLOGY": str(topology),
            "INDEX": str(index),
            "ITP_GLOB": "missing-*.itp",
            "GMX_CMD": str(fake_gmx),
            "REPLEX": "250",
        }
    )

    subprocess.run(
        ["bash", str(SETUP_SCRIPT)],
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "AT=1.5" in (tmp_path / "rep1" / "plumed.dat").read_text()
    assert "AT=2.0" in (tmp_path / "rep2" / "plumed.dat").read_text()
    assert "TARGET_DISTANCE" not in (tmp_path / "rep1" / "plumed.dat").read_text()
    assert (tmp_path / "rep1" / "reus.tpr").exists()
    assert (tmp_path / "rep2" / "reus.tpr").exists()
    mdp_lines = (tmp_path / "rep1" / "reus.mdp").read_text().splitlines()
    assert not any(line.lstrip().startswith("pull") for line in mdp_lines)

    submit_script = (tmp_path / "submit-template.sh").read_text()
    assert "REPLEX=250" in submit_script
    assert "N_REPLICA=2" in submit_script
    assert 'DEFFNM="reus"' in submit_script
