# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import sys,os,sphinx

sys.path.append(os.path.abspath("./_ext"))
sys.path.append(os.path.abspath("../modules"))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pkdgrav3'
copyright = '2023, Douglas Potter'
author = 'Douglas Potter'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions=['sphinx.ext.autodoc','sphinxcontrib.bibtex','make_parameters']
if sphinx.version_info < (1,8):
    extensions += ['sphinx.ext.imgmath']
bibtex_bibfiles = ['references.bib']
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
if sphinx.version_info >= (4,0):
    root_doc = 'index'
else:
    master_doc = 'index'
html_theme = 'alabaster'
html_static_path = ['_static']
