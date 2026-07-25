import shlex
from pathlib import Path

from ..utils.proc import run_cmd


def convert_to_gaussian_input(
    structure: str,
    filetype: str,
    output: str | Path,
) -> None:
    if not filetype:
        raise ValueError("Input structure must have a file extension")
    with Path(output).open("w") as output_file:
        run_cmd(
            ["obabel", "-i", filetype, structure, "-o", "gjf"],
            stdout=output_file,
        )


def configure_gaussian_input(
    path: str | Path,
    *,
    checkpoint: str,
    memory_gb: int,
    threads: int,
    route: str,
    charge: int,
    multiplicity: int,
) -> None:
    input_path = Path(path)
    lines = input_path.read_text().splitlines(keepends=True)
    if len(lines) < 2:
        raise ValueError(f"Gaussian input is incomplete: {input_path}")

    lines[0] = f"%chk={checkpoint}\n"
    lines[1] = f"%mem={memory_gb}GB\n"
    lines[2:2] = [f"%nprocshared={threads}\n", f"{route}\n"]

    for index, line in enumerate(lines):
        fields = line.split()
        if len(fields) != 2:
            continue
        try:
            int(fields[0])
            int(fields[1])
        except ValueError:
            continue
        lines[index] = f"{charge} {multiplicity}\n"
        break
    else:
        raise ValueError(f"Charge and multiplicity line not found in {input_path}")

    input_path.write_text("".join(lines))


def run_gaussian(
    command: str,
    input_path: str | Path,
    output_path: str | Path,
) -> None:
    command_args = shlex.split(command)
    if not command_args:
        raise ValueError("Gaussian command must not be empty")
    with (
        Path(input_path).open("rb") as input_file,
        Path(output_path).open("wb") as output_file,
    ):
        run_cmd(command_args, stdin=input_file, stdout=output_file)
