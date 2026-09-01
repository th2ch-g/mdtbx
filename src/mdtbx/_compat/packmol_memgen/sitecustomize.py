"""Provide the Python 2 name still used by PACKMOL-Memgen 2025.1.29."""

import builtins


if not hasattr(builtins, "unicode"):
    setattr(builtins, "unicode", str)
