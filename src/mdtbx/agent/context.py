"""Execution context shared by the CLI and subprocess wrapper."""

from contextvars import ContextVar


JSON_MODE = ContextVar("mdtbx_json_mode", default=False)
