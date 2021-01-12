#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

from setuptools import setup
import mljr
from mljr.mljr import __version__


setup(
    name='mljr',
    version=__version__,
    description='Modified Lydersen-Joback-Reid method',
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



