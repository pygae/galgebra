# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.append(os.path.abspath('_sphinxext'))


# -- Imports

import sphinx_rtd_theme
import sphinx
import releases_hack

# -- Project information -----------------------------------------------------

project = 'galgebra'
copyright = '2014-2019, Alan Bromborsky and GAlgebra team'
author = 'Alan Bromborsky'

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'nbsphinx',
    'nbsphinx_link',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'IPython.sphinxext.ipython_directive',
    'IPython.sphinxext.ipython_console_highlighting',
    'sphinx_markdown_tables',
    'releases',
    'sphinxcontrib.bibtex',
    'sphinx_rtd_theme',

    # local extensions
    'md_include',
]

# This is not actually used as we have a template overload to render this with
# latex. Leaving it here in case we change theme.
html_logo = "images/galgebra.svg"

bibtex_bibfiles = ['refs.bib']

# -- nbsphinx configuration ---------------------------------------------------

import galgebra
import nbsphinx
# nbsphinx_execute = 'always'
nbsphinx_execute = 'never'
nbsphinx_allow_errors=True
nbsphinx_kernel_name='python'
nbsphinx_timeout = 60

# This is processed by Jinja2 and inserted before each notebook
# Some change in dependencies made us need to replace `var` with
# `env.config.html_context['var']`.
nbsphinx_prolog = r"""
{% set docname = 'doc/' + env.doc2path(env.docname, base=None) %}
{% set git_ref = 'master' if not env.config.html_context['READTHEDOCS'] else
                 env.config.html_context['github_version']
                 if '.' not in env.config.html_context['current_version'] else
                 'v' + env.config.release %}
.. raw:: html

    <div class="admonition note">
      <p>This page was generated from
        <a class="reference external" href="https://github.com/pygae/galgebra/blob/{{ git_ref|e }}/{{ docname|e }}">{{ docname|e }}</a>.
        <!--
            This does not work yet due to nbsphinx-link
            Interactive online version:
            <a href="https://mybinder.org/v2/gh/pygae/galgebra/{{ git_ref|e }}?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>.
        -->
      </p>
      <script>
        if (document.location.host) {
          var p = document.currentScript.previousSibling.previousSibling;
          var a = document.createElement('a');
          a.innerHTML = 'View in <em>nbviewer</em>';
          a.href = `https://nbviewer.jupyter.org/url${
            (window.location.protocol == 'https:' ? 's/' : '/') +
            window.location.host +
            window.location.pathname.slice(0, -4) }ipynb`;
          a.classList.add('reference');
          a.classList.add('external');
          p.appendChild(a);
          p.appendChild(document.createTextNode('.'));
        }
      </script>
    </div>

.. raw:: latex

    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from
    \sphinxcode{\sphinxupquote{\strut {{ docname | escape_latex }}}} \dotfill}}
"""

# This is processed by Jinja2 and inserted after each notebook
nbsphinx_epilog = r"""
{% set docname = 'doc/' + env.doc2path(env.docname, base=None) %}
.. raw:: latex

    \nbsphinxstopnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{\dotfill\ \sphinxcode{\sphinxupquote{\strut
    {{ docname | escape_latex }}}} ends here.}}
"""

# -- latex configuration ---------------------------------------------------

galgebra_latex_macros = R"""
    \newcommand{\bm}[1]{\boldsymbol{#1}}
    \newcommand{\ubh}{\bm{\hat{u}}}
    \newcommand{\ebh}{\bm{\hat{e}}}
    \newcommand{\ebf}{\bm{e}}
    \newcommand{\mat}[1]{\left [ {#1} \right ]}
    \newcommand{\bra}[1]{{#1}_{\mathcal{G}}}
    \newcommand{\ket}[1]{{#1}_{\mathcal{D}}}
    \newcommand{\ds}{\displaystyle}
    \newcommand{\bfrac}[2]{\displaystyle\frac{#1}{#2}}
    \newcommand{\lp}{\left (}
    \newcommand{\rp}{\right )}
    \newcommand{\half}{\frac{1}{2}}
    \newcommand{\llt}{\left <}
    \newcommand{\rgt}{\right >}
    \newcommand{\abs}[1]{\left |{#1}\right |}
    \newcommand{\pdiff}[2]{\bfrac{\partial {#1}}{\partial {#2}}}
    \newcommand{\pdifftwo}[3]{\bfrac{\partial^{2} {#1}}{\partial {#2}\partial {#3}}}
    \newcommand{\lbrc}{\left \{}
    \newcommand{\rbrc}{\right \}}
    \newcommand{\set}[1]{\lbrc {#1} \rbrc}
    \newcommand{\W}{\wedge}
    \newcommand{\R}{\dagger}
    \newcommand{\lbrk}{\left [}
    \newcommand{\rbrk}{\right ]}
    \newcommand{\com}[1]{\lbrk {#1} \rbrk}
    \newcommand{\proj}[2]{\llt {#1} \rgt_{#2}}
    %\newcommand{\bm}{\boldsymbol}
    \newcommand{\braces}[1]{\left \{ {#1} \right \}}
    \newcommand{\grade}[1]{\left < {#1} \right >}
    \newcommand{\f}[2]{{#1}\lp {#2} \rp }
    \newcommand{\paren}[1]{\lp {#1} \rp }
    \newcommand{\eval}[2]{\left . {#1} \right |_{#2}}
    \newcommand{\prm}[1]{{#1}'}
    \newcommand{\ddt}[1]{\bfrac{d{#1}}{dt}}
    \newcommand{\deriv}[3]{\bfrac{d^{#3}#1}{d{#2}^{#3}}}
    \newcommand{\be}{\begin{equation}}
    \newcommand{\ee}{\end{equation}}
    \newcommand{\eb}{\bm{e}}
    \newcommand{\ehb}{\bm{\hat{e}}}
    \newcommand{\Tn}[2]{\f{\mathcal{T}_{#2}}{#1}}
    \newcommand{\tr}{\mbox{tr}}
    \newcommand{\T}[1]{\texttt{#1}}
    \newcommand{\grd}{\bm{\nabla}}
    \newcommand{\indices}[1]{#1}
    \newcommand{\xRightarrow}[1]{\overset{#1}{\Rightarrow}}
"""

# https://www.sphinx-doc.org/en/master/usage/extensions/math.html#module-sphinx.ext.mathjax
# mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

# # enable autonumbering
# mathjax_config = dict(TeX=dict(
#     equationNumbers=dict(autoNumber="AMS"),
#     extensions=["[Contrib]/preamble/preamble.js", "color.js"],
#     preamble=[galgebra_latex_macros]
# ))
# https://www.sphinx-doc.org/en/master/usage/extensions/math.html#module-sphinx.ext.mathjax
# https://docs.mathjax.org/en/latest/input/tex/extensions/autoload.html#autoload-options
# https://docs.mathjax.org/en/latest/input/tex/eqnumbers.html
# mathjax3_config = dict(
#     tex=dict(
#         tags='ams'
#     )
# )


# -- other extension configuration --------------------------------------------

napoleon_include_init_with_doc= False

autoclass_content = "both"  # include both class docstring and __init__
autodoc_default_options = {
        # Make sure that any autodoc declarations show the right members
        "members": True,
        # do not show inherited, as this pull in all the sympy members
        # "inherited-members": True,
        "member-order": "bysource",
        # "undoc-members": True,
        # "special-members": True,
        # "private-members": True,
        "show-inheritance": True,
}

#autodoc_default_flags='members'
# you have to list all files with automodule here due to bug in sphinx and nbsphinx
# https://github.com/spatialaudio/nbsphinx/issues/14
autosummary_generate=['api']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# from recommonmark.parser import CommonMarkParser

# source_parsers = {
#     '.md': CommonMarkParser,
# }

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ['.rst']

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    '_build', '**.ipynb_checkpoints', 'Thumbs.db', '.DS_Store',
    # these are here for users, not for Sphinx
    'books', 'old_installation.md', 'old_introduction.md',
    # this is converted into ipynb elsewhere
    'galgebra.md',
]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If your project is hosted on Github, set the releases_github_path setting instead, 
# to e.g. account/project. Releases will then use an appropriate Github URL for both
# releases and issues.
releases_github_path = "pygae/galgebra"
releases_release_uri = releases_hack.release_uri(releases_github_path)

# You may optionally set releases_debug = True to see debug output while building your docs.
releases_debug = True

# If your changelog includes “simple” pre-1.0 releases derived from a single branch
# (i.e. without stable release lines & semantic versioning) you may want to set 
# releases_unstable_prehistory = True.
releases_unstable_prehistory = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = dict(
    logo_only=True,
)

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# html_context = {
#     'css_files': [
#         '_static/css/theme.css',
#         '_static/theme_overrides.css',  # override wide tables in RTD theme
#         ],
#      }

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'galgebradoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    'preamble': galgebra_latex_macros,

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'galgebra.tex', 'galgebra Documentation',
     'Alan Bromborsky', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'galgebra', 'galgebra Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'galgebra', 'galgebra Documentation',
     author, 'galgebra', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'sympy': ('https://docs.sympy.org/latest', None),
}
