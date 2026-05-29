## \file Convergence.py
# \brief Class for performing tasks for getConvergence.py.
# \author Kevin K. Chen (Princeton University)

## \namespace Convergence
# \brief The namespace for the Convergence class.
# \author Kevin K. Chen (Princeton University)

"""Class for performing tasks for getConvergence.py."""

import argparse
import re

import numpy
import matplotlib as mpl
from matplotlib import pyplot as plt


## Class for performing tasks for getConvergence.py.
class Convergence:
    ## Parse the input arguments and construct the object.
    #
    # \param[in] argv List of command-line arguments.
    # \param[in] doc Docstring describing the program.
    def __init__(self, argv, doc):
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

        parser.add_argument(
            '-s', '--suppress',
            action='store_true',
            help='suppress printing the convergence data to the standard output'
        )
        parser.add_argument(
            '-d', '--display',
            action='store_true',
            help='display the plot of the convergence data'
        )
        parser.add_argument(
            '-p', '--save',
            help='print the figure to the given file name'
        )
        parser.add_argument(
            'solver',
            choices=('Newton', 'Armijo', 'GMRES'),
            help='solver whose convergence information is printed or plotted'
        )

        parser.parse_args(argv, self)

    ## Process the text, and plot if necessary.
    #
    # \param[in] text The verbose output from Newton_GMRES.
    def run(self, text):
        self._readData(text)

        # Print the data.
        if not self.suppress:
            self._print()

        # Plot the data.
        if self.display or self.save:
            self._plot()

    ## Read the verbose output from Newton_GMRES and extract convergence data.
    #
    # \param[in] text The verbose output from Newton_GMRES.
    def _readData(self, text):
        num = r'([-+]?(\d+\.?\d*|\.\d+)([Ee][-+]?\d+)?)'

        if self.solver == 'Newton':
            self._data = numpy.array(
                re.findall(
                    r'{} iteration (\d+): \|\|fx\|\| = {}; tolerance = {}'
                    .format(self.solver, num, num),
                    text
                )
            )
            self._data = self._data[:, [0, 1, 4]]
        elif self.solver == 'Armijo':
            self._data = numpy.array(
                re.findall(
                    (r'{} iteration (\d+): lambda = {}; ' +
                     r'\|\|f\(x \+ lambda \* h\)\|\| = {}; tolerance = {}.')
                    .format(self.solver, num, num, num),
                    text
                )
            )
            self._data = self._data[:, [0, 1, 4, 7]]
        elif self.solver == 'GMRES':
            self._data = numpy.array(
                re.findall(
                    r'{} iteration (\d+): \|\|A \* x - b\|\| = {}; tolerance = {}.'
                    .format(self.solver, num, num),
                    text
                )
            )
            self._data = self._data[:, [0, 1, 4]]

    ## Print the convergence data to the standard output.
    def _print(self):
        for match in self._data:
            print(' '.join(match))

    ## Plot the convergence data, and save if requested.
    def _plot(self):
        ylabels = {'Newton': r'$\|\mathbf{f} (\mathbf{x})\|$',
                   'Armijo': r'$\|\mathbf{f} (\mathbf{x} + \lambda \mathbf{h})\|$',
                   'GMRES': (r'$\|(\partial \mathbf{f} / \partial \mathbf{x})' +
                             r'\mathbf{h} - \mathbf{f} (\mathbf{x})\|$')}

        mpl.rcParams['text.usetex'] = True

        plt.figure()
        ax1 = plt.axes()
        plt.hold(True)

        plt.semilogy(self._data[:, -2], 'b')
        plt.semilogy(self._data[:, -1], 'r')

        plt.xlabel('Cumulative iterations')
        plt.ylabel(ylabels[self.solver])
        plt.title(self.solver + ' convergence')
        plt.legend(('value', 'tolerance'), loc='best')
        plt.grid(True)

        if self.solver == 'Armijo': # Plot lambda.
            ax2 = ax1.twinx()

            plt.semilogy(self._data[:, 1], 'g')

            ax2.tick_params(axis='y', colors='g')
            plt.ylabel('$\lambda$', color='g')

        if self.save is not None:
            plt.savefig(self.save)
        if self.display:
            plt.show()
