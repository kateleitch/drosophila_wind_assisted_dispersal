from __future__ import print_function
import sys
import matplotlib.pyplot as plt
import kml_gen as kml

dir = sys.argv[1]
generator = kml.kmlGenerator(dir)
generator.run()
