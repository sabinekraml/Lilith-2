# Lilith-2

Lilith (Light Likelihood Fit for the Higgs) is a light and easy-to-use Python tool for constraining new physics from signal strength measurements of the 125 GeV Higgs boson. Lilith is provided with the latest experimental measurements from the ATLAS and CMS collaborations at the LHC. The Higgs likelihood is based on experimental results stored in an easily extensible XML database; it is evaluated from the user input, given in XML format in terms of reduced couplings or signal strengths. 

New from version 2.0 onwards (more in changelog.txt):

- Use of variable Gaussian and generalized Poisson likelihoods for a __better treatment of assymetric uncertainties__. 
The generalized Poisson likelihood can be used for experimental results in 1 or 2 dimensions (the latter with correlation), while the variable Gaussian approximation is available for results of any dimension. 

- Use of __N-dim correlation matrices__ for ordinary and variable Gaussian likelihoods.

- Database 19.06 contains the __complete ATLAS and CMS Run 2 results for 36 fb-1__.

Prerequisites: Python 2.7 (not yet 3.X!), SciPy and NumPy; 
the example codes doing a likelihood profile analysis require iminuit.

Tested with Python 2.7.x (x>4), SciPy 0.13.0b1, NumPy 1.8.0rc1, iminuit 1.2.

Notes:

- The __master__ branch contains the latest official version. If you want to see the validation (plots and scripts, plus the various versions of xml files used), check out the __validation__ branch. 

- In case of problems running the code, check whether the `__init.py__` file exists in lilith/internal/ and is executable. If not, create it (as an empty file) and declare it as executable. If the code still does not work, check that all the Lilith Python (`.py`) files are executable.  

- If you get an error `ImportError: No module named lilith`, your path to Lilith is probably not correct. (e.g, when this happens with the Python example codes, check where `lilith_dir` points to)

Usage: see the original Lilith manual https://arxiv.org/abs/1502.04138 (a new manual is in preparation)

For usage in micrOMEGAs, see https://arxiv.org/abs/1606.03834

