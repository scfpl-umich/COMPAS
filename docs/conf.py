# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

# Sphinx configuration for the COMPAS documentation.
# Build locally with:  pip install -r docs/requirements.txt
#                      sphinx-build -W --keep-going -b html docs docs/_build/html

project = "COMPAS"
author = "The Regents of the University of Michigan"
copyright = "2026, The Regents of the University of Michigan. GPL-3.0-or-later"

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
]

# Pages are Markdown (MyST). $...$ and $$...$$ are math, ::: fences hold directives such as
# notes, and headings get anchors so pages can link to a section with page.md#section.
source_suffix = {".md": "markdown"}
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist"]
myst_heading_anchors = 4

exclude_patterns = ["_build", "requirements.txt", "media"]

html_theme = "furo"
html_title = "COMPAS"
html_static_path = ["_static", "media/videos"]   # the videos are fetched into docs/media
html_css_files = ["custom.css"]
html_theme_options = {
    "source_repository": "https://github.com/scfpl-umich/COMPAS/",
    "source_branch": "main",
    "source_directory": "docs/",
}
