from __future__ import print_function
import sys
import matplotlib.pyplot as plt
import trap_layout as t

dir = sys.argv[1]
planned_or_actual = sys.argv[2]

turn_off_text = False
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='polar')
analyzer = t.TrapLayoutVisualizer(dir, planned_or_actual, ax, turn_off_text)
if planned_or_actual == 'planned':
    analyzer.run_planned()
