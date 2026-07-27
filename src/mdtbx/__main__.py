"""Command-line entry point for humans and autonomous agents."""


def main():
    from . import prepare_runtime
    from .cli import cli

    prepare_runtime()
    cli()


if __name__ == "__main__":
    main()
