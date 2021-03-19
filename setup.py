#!/usr/bin/env python3
from setuptools import setup

with open("README.md") as fh:
    description = fh.read()

setup(
    name='jasper-vh',
    version='2.0',
    description="Just a simple virus's hosts predictor",
    long_description=description,
    long_description_content_type='text/markdown',
    author='Milosz Chodkowski',
    keywords='bioinformatics sequence DNA trna CRISPR blast virus host',
    license="GPLv3",
    author_email='milosz.chodkowski@student.put.poznan.pl',
    url="https://github.com/777moneymaker/jasper",
    download_url="https://github.com/777moneymaker/jasper/archive/v2.0.tar.gz",
    packages=['jasper'],
    install_requires=[
        "biopython",
        "numpy",
        "pandas",
        "argparse",
    ],
    scripts=['jasper-vh'],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        'Environment :: Console',
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Information Technology",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
