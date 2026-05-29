#!/usr/bin/env python3

## \file getConvergence.py
# \brief Extract convergence information from a Newton_GMRES verbose log file.
# \author Kevin K. Chen (Princeton University)

## \namespace getConvergence
# \brief The namespace for getConvergence.py.
# \author Kevin K. Chen (Princeton University)

"""Extract convergence information from a Newton_GMRES verbose log file.

When Newton is selected as the solver, this script searches the standard input
for lines

    Newton iteration <n>: ||fx|| = <fxNorm>; tolerance = <tol>.

and reprints all such lines as

    <n> <fxNorm> <tol>

When Armijo is given as an argument, this script searches the standard input for
lines

    Armijo iteration <n>: lambda = <lambda>; ||f(x + lambda * h)|| = <fCurrNorm>; tolerance = <tol>.

and reprints all such lines as

    <n> <lambda> <fCurrNorm> <tol>

When GMRES is given as an argument, this script searches the standard input for
lines

    GMRES iteration <n>: ||A * x - b|| = <rho>; tolerance = <errtol>.

and reprints all such lines as

    <n> <rho> <errtol>

The convergence information can also be plotted with matplotlib.

"""

import sys
import Convergence


## Read the command line arguments and data, and print and plot as requested.
def main(argv):
    convergence = Convergence.Convergence(argv, __doc__)
    convergence.run(sys.stdin.read())


if __name__ == '__main__':
    main(sys.argv[1:])
