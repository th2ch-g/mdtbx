import warnings

warnings.simplefilter("ignore")

from .cli import cli  # noqa: E402


def main():
    cli()


if __name__ == "__main__":
    main()
