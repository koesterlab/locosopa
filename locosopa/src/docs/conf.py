# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import json
import os
from pathlib import Path
from datetime import datetime
import yaml
from sphinxawesome_theme.postprocess import Icons

lsc = yaml.load(open(os.environ["LOCOSOPA_CONFIG"]), Loader=yaml.SafeLoader)
lsc_dir = Path(os.environ["LOCOSOPA_CONFIG"]).parent

date = datetime.now()

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
]

# templates_path = ['_templates']
exclude_patterns = []
# pygments_style = "sphinx"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinxawesome_theme'
html_static_path = ['_static']
html_css_files = ["custom.css"]
html_theme_options = {
    "logo_light": (lsc_dir / lsc["logo"]["light"]).name,
    "logo_dark": (lsc_dir / lsc["logo"]["dark"]).name,
    "awesome_external_links": True,
    "awesome_headerlinks": True,
    "show_prev_next": False,
    "main_nav_links": dict(),
}
links = lsc["docs"].get("links")
if links is not None:
    html_theme_options["main_nav_links"].update(links)

html_permalinks_icon = Icons.permalinks_icon

project = lsc["project"]["name"]
author = lsc["project"]["authors"]
release = lsc["project"]["version"]
copyright = lsc["project"]["copyright"]

globals().update(lsc["docs"].get("sphinx-config", dict()))