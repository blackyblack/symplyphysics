#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os

from setuptools import find_packages, setup  # type: ignore

# Package meta-data.
NAME = 'symplyphysics'
DESCRIPTION = 'Physics laws implemented as code.'
URL = 'https://github.com/blackyblack/symplyphysics'
EMAIL = 'sam.and.tetris@gmail.com'
AUTHOR = 'blackyblack'
REQUIRES_PYTHON = '>=3.10.0'
VERSION = '1.0.0'

# What packages are required for this module to be executed?
REQUIRED = ['sympy']

# What packages are optional?
EXTRAS = {'plots': ['matplotlib']}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
long_description: str = DESCRIPTION
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    pass

# Where the magic happens:
setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='MIT',
    classifiers=[
    'License :: OSI Approved :: MIT License', 'Programming Language :: Python',
    'Programming Language :: Python :: 3', 'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: Implementation :: CPython',
    'Programming Language :: Python :: Implementation :: PyPy'
    ],
)
