# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cosmos'
copyright = '2025, Chul-Moon Yoo et al.'
author = 'Chul-Moon Yoo et al.'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['breathe', 'exhale']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'
# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

# -- Set up breathe extension ------------------------------------------------
breathe_projects = { "cosmos": "./_doxygen/xml" }
breathe_default_project = "cosmos"

# -- Set up the exhale extension ---------------------------------------------
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "cosmos-api-reference.rst",
    "doxygenStripFromPath":  "..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "COSMOS API Reference",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    "INPUT = ./../source/\nRECURSIVE = YES\nGENERATE_HTML = YES",
    "verboseBuild" : True
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
