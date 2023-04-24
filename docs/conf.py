# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import sys,os

sys.path.append(os.path.abspath("./_ext"))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pkdgrav3'
copyright = '2023, Douglas Potter'
author = 'Douglas Potter'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# If there is an error like WARNING: using "math" markup without a Sphinx math extension active
# switch to the second line
extensions = ['make_parameters']
# extensions = ['make_parameters','sphinx.ext.imgmath']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Uncomment if there is an error like ....contents.rst not found
# master_doc = 'index'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
