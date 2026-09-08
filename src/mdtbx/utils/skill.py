"""Display the bundled agent guide."""

from importlib.resources import files


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "skill",
        help="Print the agent guide",
        description="Print the bundled English agent guide as Markdown to stdout.",
    )
    parser.set_defaults(func=run)


def run(args):
    print(files("mdtbx").joinpath("MDTBX_SKILL.md").read_text(encoding="utf-8"), end="")
