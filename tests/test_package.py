import os
import subprocess
import sys


def test_import_has_no_path_or_pymol_side_effects():
    code = """
import os
import sys
before = os.environ.get("PATH")
import mdtbx
assert os.environ.get("PATH") == before
assert "pymol_plugins" not in sys.modules
"""
    environment = dict(os.environ)
    environment["PYTHONPATH"] = "src"
    subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        cwd=os.getcwd(),
        env=environment,
        capture_output=True,
        text=True,
    )
