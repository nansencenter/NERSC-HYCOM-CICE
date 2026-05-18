project = "NERSC-HYCOM-CICE"
author = "Nansen Environmental and Remote Sensing Center"
copyright = "2024, Nansen Environmental and Remote Sensing Center"
release = "develop"

extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_copybutton",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "rm.bak.md"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "navigation_depth": 3,
    "titles_only": False,
}
