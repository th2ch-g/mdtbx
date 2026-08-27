"""Command-line entry point."""


def main():
    from . import prepare_runtime
    from .cli import cli

    prepare_runtime()
    cli()


if __name__ == "__main__":
    main()
