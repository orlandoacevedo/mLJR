#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
from setuptools import setup
import mljr
from mljr.mljr import __version__


def read(*names):
    values = {}
    extensions = ['', '.txt', '.rst', '.md',]
    for name in names:
        v = ''
        for ext in extensions:
            filename = name + ext
            if os.path.isfile(filename):
                with open(filename) as f:
                    v = f.read()
        values[name] = v
    return values

long_description = """%(README)s""" % read('README')

setup(
    name='mljr',
    version=__version__,
    description='Modified Lydersen-Joback-Reid method',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    author='Orlando Acevedo',
    author_email='orlando.acevedo@miami.edu',
    keywords='Modified Lydersen-Joback-Reid method',
    maintainer='Xiang Zhong',
    maintainer_email='xxz385@miami.edu',
    url='https://github.com/orlandoacevedo/mLJR',
    license='MIT',
    packages=['mljr'],
    entry_points={
        'console_scripts': [
            'mljr = mljr.mljr:main',
        ]
    },
    python_requires='>=3.5',
)



