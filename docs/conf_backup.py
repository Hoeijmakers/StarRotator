import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'StarRotator'
author = 'Your Name'
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints'
]


