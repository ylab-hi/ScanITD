# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

# Make the src layout importable so autodoc can find the package
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

# -- Project information -----------------------------------------------------
project = "ScanITD"
copyright = "2024, Ting-You Wang"
author = "Ting-You Wang"
release = "0.9.1"
version = "0.9"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",        # Auto-generate API docs from docstrings
    "sphinx.ext.napoleon",       # Google/NumPy docstring support
    "sphinx.ext.viewcode",       # Add [source] links to API pages
    "sphinx.ext.autosummary",    # Generate summary tables
    "sphinx.ext.intersphinx",    # Cross-reference other projects
    "sphinx.ext.todo",           # Support .. todo:: directives
    "myst_parser",               # Parse Markdown files
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Napoleon settings (Google-style docstrings) ----------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True

# -- Autodoc settings -------------------------------------------------------
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": False,
    "exclude-members": "__weakref__, __dict__, __module__",
    "show-inheritance": True,
}
autodoc_typehints = "description"
autodoc_typehints_format = "short"
add_module_names = False

# -- Autosummary settings ---------------------------------------------------
autosummary_generate = True

# -- Intersphinx mapping ----------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pysam": ("https://pysam.readthedocs.io/en/latest", None),
}

# -- MyST (Markdown) settings -----------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
]

# -- Todo extension ---------------------------------------------------------
todo_include_todos = True

# -- Options for HTML output ------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "navigation_depth": 4,
    "titles_only": False,
    "logo_only": False,
    "prev_next_buttons_location": "both",
    "style_external_links": True,
    "collapse_navigation": False,
    "sticky_navigation": True,
}

html_context = {
    "display_github": True,
    "github_user": "ylab-hi",
    "github_repo": "ScanITD",
    "github_version": "main",
    "conf_py_path": "/docs/source/",
}
