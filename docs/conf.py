from importlib.metadata import version as package_version

project = "mdtbx"
author = "mdtbx contributors"
release = package_version("mdtbx")
version = release

extensions = [
    "sphinxcontrib.autoprogram",
]

root_doc = "index"
language = "en"
locale_dirs = ["locale/"]
gettext_compact = False

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

html_theme = "furo"
html_title = f"{project} {release}"
html_theme_options = {
    "announcement": (
        '<a href="https://th2ch-g.github.io/mdtbx/">English</a>'
        " / "
        '<a href="https://th2ch-g.github.io/mdtbx/ja/">日本語</a>'
    ),
    "navigation_with_keys": True,
    "source_repository": "https://github.com/th2ch-g/mdtbx/",
    "source_branch": "main",
    "source_directory": "docs/",
    "top_of_page_buttons": ["view", "edit"],
}
