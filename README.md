# emt
Eccentric Mass Transfer (emt)

Code to integrate the orbit-averaged equations of motion describing mass transfer in eccentric orbits.

Requires: Python, numpy, scipy, and included emtfunctions library (type `make` to install in current directory)

Testing emtfunctions: use nosetests on testemtlibrary.py: `nosetests testemtlibrary`

Usage: `python integrator.py`. Use command-line arguments to specify parameters,
or, for more flexibility, use the included applications.py' script to define parameters
in-script: `python applications.py --id x`, where `x` refers to the application.

Adrian Hamers, November 2018

Updated July 2019 (equilibrium tides)
