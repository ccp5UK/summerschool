"""
Module containing classes for loading rdf data from DL_POLY_4
"""

import numpy as np
import matplotlib.pyplot as plt

class rdf():
    """ class for reading RDFDAT

        :param source: Source RDF to read

        """
    __version__ = "0"

    def __init__(self, source=None):
        self.nRDF = 0
        self.nPoints = 0
        self.x = None
        self.data = None
        self.labels = None
        if source is not None:
            self.source = source
            self.read(source)

    def read(self, source="RDFDAT"):
        """ Read an RDF file into data

        :param source: File to read

        """
        with open(source, 'r') as f:
            a = f.readline().split()[0]
            with open(source, 'r') as fileIn:
                # Discard title
                _ = fileIn.readline()
                self.nRDF, self.nPoints = map(int, fileIn.readline().split())

                self.x = np.zeros(self.nPoints)
                self.data = np.zeros((self.nRDF, self.nPoints))
                self.labels = []
                s = True
                for sample in range(self.nRDF):
                    species = fileIn.readline().split()
                    if len(species) == 0:
                        break
                    self.labels.append(species)
                    for point in range(self.nPoints):
                        r, g_r, _ = map(float, fileIn.readline().split())
                        if s:
                            self.x[point] = r
                        self.data[sample, point] = g_r
                    s = False

    def plot(self):
        for i in range(len(self.labels)):
            plt.plot(self.x, self.data[i,:],label = "-".join(self.labels[i]))
        plt.xlabel("r [Å])")
        plt.ylabel("gofr [a.u.])")
        plt.legend()
        plt.show()
