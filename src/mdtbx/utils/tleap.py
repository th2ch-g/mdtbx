"""Shared tleap invocation used by the AMBER build subcommands.

Replaces the duplicated "write tleap.in -> run tleap -> rm leap.log tleap.in"
sequence in build_solution / build_vacuum / gen_loop_aa / gen_resp / gen_am1bcc.
"""

import re
from pathlib import Path

from .proc import run_cmd


def run_tleap(input_text, *, keepfiles=False, extra_cleanup=(), cwd=None):
    """Write tleap.in, run ``tleap -f tleap.in``, then clean up.

    Parameters
    ----------
    input_text : str
        Full contents of the tleap input script.
    keepfiles : bool
        When True, leave tleap.in / leap.log (and extra_cleanup) on disk.
    extra_cleanup : Iterable[str]
        Additional files to remove unless keepfiles is set.
    cwd : str or Path, optional
        Isolated tleap working directory.
    """
    workdir = Path(cwd).expanduser().resolve() if cwd is not None else Path.cwd()
    input_path = workdir / "tleap.in"
    log_path = workdir / "leap.log"
    existing = [str(path) for path in (input_path, log_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing tleap working file(s): "
            + ", ".join(existing)
        )

    input_path.write_text(input_text)
    try:
        run_cmd(["tleap", "-f", str(input_path)], cwd=workdir, check=True)
        if log_path.is_file():
            log_text = log_path.read_text(errors="replace")
            summaries = re.findall(r"Exiting LEaP: Errors = (\d+)", log_text)
            if summaries and int(summaries[-1]) != 0:
                raise RuntimeError(
                    f"tleap reported {summaries[-1]} error(s); see {log_path}"
                )
    finally:
        if not keepfiles:
            cleanup = [log_path, input_path]
            for name in map(Path, extra_cleanup):
                cleanup.append(name if name.is_absolute() else workdir / name)
            for path in cleanup:
                path.unlink(missing_ok=True)
