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
    Path("tleap.in").write_text(input_text)
    run_cmd(["tleap", "-f", "tleap.in"], check=True)
    if not keepfiles:
        for name in ("leap.log", "tleap.in", *extra_cleanup):
            Path(name).unlink(missing_ok=True)
