"""mdtbx package.

Runtime side effects that must run once before any subcommand executes live
here (they previously lived in config.py and were triggered indirectly via
`from ..config import *`):

- Prepend the pixi environment's bin directory to PATH so bundled external
  tools (gmx, tleap, antechamber, ...) resolve.
- Register the PyMOL plugins package (mocked out in tests via sys.modules).
"""

import importlib
import os
from pathlib import Path

_current_path = os.environ.get("PATH", "")
_pixi_bin = Path(__file__).parent.parent / ".pixi/envs/default/bin"
_path_entries = _current_path.split(os.pathsep) if _current_path else []
if _pixi_bin.is_dir() and str(_pixi_bin) not in _path_entries:
    os.environ["PATH"] = os.pathsep.join([str(_pixi_bin), *_path_entries])

importlib.import_module("pymol_plugins")
