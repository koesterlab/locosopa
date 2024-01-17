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

html_theme = 'lutra'
html_static_path = ['_static']
html_css_files = [str(Path(__file__).parent / "theme.css")]
html_theme_options = {
    "primary_color": lsc["color"],
    "secondary_color": lsc["color"],
    "dark_logo": (lsc_dir / lsc["logo"]["dark"]).name,
    "light_logo": (lsc_dir / lsc["logo"]["dark"]).name,
    "navigation_style": "plain",
}
sidebar_links = lsc["docs"].get("sidebar-links")
if sidebar_links is not None:
    html_theme_options["sidebar_links"] = sidebar_links

globals().update(lsc["docs"]["sphinx-config"])