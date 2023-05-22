#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import setuptools
from glob import glob
setuptools.setup(name="benj",
                 version="0.0.3",
                 author="Benjamin James",
                 author_email="benjames@mit.edu",
                 url="https://github.com/KellisLab/benj",
                 license="GPL",
                 install_requires=[
                     "numpy",
                     "scipy",
                     "scikit-learn",
                     "numba",
                     "anndata>=0.8.0",
                     "pandas",
                     "umap-learn",
                     "pynndescent>=0.5.7",
                     "scanpy",
                     "matplotlib",
                     "seaborn",
                     # "python-igraph",
                     # "leidenalg",
                     "tqdm",
                     "pyranges",
                 ],
                 packages=setuptools.find_packages("."),
                 test_suite="test",
                 scripts=glob("scripts/*.py") + glob("scripts/*.sh")
                 )
