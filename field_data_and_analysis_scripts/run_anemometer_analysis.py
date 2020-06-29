from __future__ import print_function
import sys
import anemometer_analysis as anem
import matplotlib.pyplot as plt

dir = sys.argv[1]
fig = plt.figure(figsize=(8,8))
ax = plt.subplot(1,1,1)
analyzer = anem.AnemometerAnalyzer(directory =dir, ax_handle=ax, time_shift = 0, annie_sim_data = False, bin_duration = 30.)
analyzer.run()
