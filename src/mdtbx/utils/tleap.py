"""Shared tleap invocation used by the AMBER build subcommands.

Replaces the duplicated "write tleap.in -> run tleap -> rm leap.log tleap.in"
sequence in build_solution / build_vacuum / gen_loop_aa / gen_resp / gen_am1bcc.
"""

from pathlib import Path

from .proc import run_cmd


def run_tleap(input_text, *, keepfiles=False, extra_cleanup=()):
    """Write tleap.in, run ``tleap -f tleap.in``, then clean up.

    Parameters
    ----------
    input_text : str
        Full contents of the tleap input script.
    keepfiles : bool
        When True, leave tleap.in / leap.log (and extra_cleanup) on disk.
    extra_cleanup : Iterable[str]
        Additional files to remove unless keepfiles is set.
    """
    input_path = Path("tleap.in")
    log_path = Path("leap.log")
    existing = [str(path) for path in (input_path, log_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing tleap working file(s): "
            + ", ".join(existing)
        )

    input_path.write_text(input_text)
    try:
        run_cmd(["tleap", "-f", str(input_path)], check=True)
    finally:
        if not keepfiles:
            for name in (log_path, input_path, *map(Path, extra_cleanup)):
                name.unlink(missing_ok=True)
