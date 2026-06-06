"""mdtbx package.

Runtime side effects that must run once before any subcommand executes live
here (they previously lived in config.py and were triggered indirectly via
`from ..config import *`):

- Prepend the pixi environment's bin directory to PATH so bundled external
  tools (gmx, tleap, antechamber, ...) resolve.
- Register the PyMOL plugins package (mocked out in tests via sys.modules).
"""

import os
from pathlib import Path

_current_path = os.environ.get("PATH", "")
_pixi_bin = Path(__file__).parent.parent / ".pixi/envs/default/bin"
os.environ["PATH"] = os.pathsep.join([str(_pixi_bin), _current_path])

import pymol_plugins  # noqa: E402,F401
