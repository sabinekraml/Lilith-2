# Lilith-2

Lilith (Light Likelihood Fit for the Higgs) is a light and easy-to-use Python tool for constraining new physics from signal strength measurements of the 125 GeV Higgs boson. Lilith is provided with the latest experimental measurements from the ATLAS and CMS collaborations at the LHC. The Higgs likelihood is based on experimental results stored in an easily extensible XML database, and is evaluated from the user input, given in XML format in terms of reduced couplings or signal strengths. 

Lilith-2.0 features variable Gaussian and Poisson likelihoods as well as a complete database of 36 fb-1 Run 2 results.

Notes:

- The __master__ branch contains the latest official version. If you want to see the validation (plots and scripts, plus the various versions of xml files used), check out the __validation__ branch. 

- In case of problems running the code, check whether the `__init.py__` file exists in lilith/internal/ and is executable. If not, create it (as an empty file) and declare it as executable. If the code still does not work, check that all the Lilith Python (`.py`) files are executable.  
