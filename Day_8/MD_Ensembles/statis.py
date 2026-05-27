"""
File containing methods for loading statistics data from DL_POLY_4
Author: alin m elena, jacob wilkins, 2019-2021
"""

import numpy as np
import matplotlib.pyplot as plt

key = {'temperature': 4, 'total_energy':3}
class Statis():
    """Type to parse and interpret STATIS file

    :param source: STATIS to read
    :param control: Associated CONTROL
    :param config: Associated CONFIG

    """
    __version__ = "0"

    def __init__(self, source=None, control=None, config=None):
        self.rows = 0
        self.columns = 0
        self.data = None
        if source is not None:
            self.source = source
            self.read(source)

    _labelPos = property(lambda self: (len(self.labels)//5+1, len(self.labels) % 5+1))

    def add_label(self, arg):
        """Add a label to the list of labels

        :param arg: Label to add
        """
        self.labels.append("{0:d}-{1:d} {2:s}".format(*self._labelPos, arg))

    def read(self, filename="STATIS"):
        """Read and parse a STATIS file

        :param filename: File to read
        :returns: Parsed statis
        """
        with open(filename, 'r') as f:
            a = f.readline().split()[0]
            with open(filename, 'r') as f:
                h1, h2, s = f.read().split('\n', 2)
                self.data = np.array(s.split(), dtype=float)
                self.columns = int(self.data[2])
                self.rows = self.data.size//(self.columns + 3)
                self.data.shape = self.rows, self.columns + 3
                np.delete(self.data, 2, axis=1)
                self.columns += 2
        return self

    def gen_labels(self, control=None, config=None):
        """Generate labels for headers in STATIS file

        :param control: Control file relating to statis
        :param config: Config file relating to statis
        :returns: Set labels for further reference
        """
        self.labels = ["1-1 Total Extended System Energy",
                       "1-2 System Temperature",
                       "1-3 Configurational Energy",
                       "1-4 Short Range Potential Energy",
                       "1-5 Electrostatic Energy",
                       "2-1 Chemical Bond Energy",
                       "2-2 Valence Angle And 3-Body Potential Energy",
                       "2-3 Dihedral, Inversion, And 4-Body Potential Energy",
                       "2-4 Tethering Energy",
                       "2-5 Enthalpy (Total Energy + Pv)",
                       "3-1 Rotational Temperature",
                       "3-2 Total Virial",
                       "3-3 Short-Range Virial",
                       "3-4 Electrostatic Virial",
                       "3-5 Bond Virial",
                       "4-1 Valence Angle And 3-Body Virial",
                       "4-2 Constraint Bond Virial",
                       "4-3 Tethering Virial",
                       "4-4 Volume",
                       "4-5 Core-Shell Temperature",
                       "5-1 Core-Shell Potential Energy",
                       "5-2 Core-Shell Virial",
                       "5-3 Md Cell Angle Α",
                       "5-4 Md Cell Angle Β",
                       "5-5 Md Cell Angle Γ",
                       "6-1 Pmf Constraint Virial",
                       "6-2 Pressure",
                       "6-3 External Degree Of Freedom",
                       "6-4 stress xx",
                       "6-5 stress xy",
                       "7-1 stress xz",
                       "7-2 stress yx",
                       "7-3 stress yy",
                       "7-4 stress yz",
                       "7-5 stress zx",
                       "8-1 stress zy",
                       "8-2 stress zz"]

        # Catch Remainder
        for i in range(len(self.labels)+1, self.columns+1):
            self.add_label("col_{0:d}".format(i))
        self.labels = ['iter', 'time', 'vars'] + self.labels

    def flatten(self):
        """FIXME! briefly describe function"""
        for i in range(self.columns-3):
            with open(self.labels[i], 'w') as f:
                for j in range(self.rows):
                    f.write("{} {}\n".format(self.data[j, 1], self.data[j, i+3]))

    def get_total_energy(self):
        return self.data[:, key['total_energy']]

    def plot(self,label):
        if label in key.keys():
            plt.plot(self.data[:, 1],self.data[:, key[label]])
            plt.xlabel("t [ps])")
            plt.ylabel(label + "[a.u.])")
            plt.show()

