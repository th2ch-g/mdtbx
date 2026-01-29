import warnings

warnings.simplefilter('ignore')

from .cli import cli

def main():
    cli()


if __name__ == "__main__":
    main()
