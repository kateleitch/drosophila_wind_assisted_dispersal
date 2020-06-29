from __future__ import print_function
import cv2 # opencv
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import time
import json
from scipy.optimize import curve_fit
from scipy.stats import circmean
from scipy.stats import vonmises
from scipy.stats import circvar
import trap_layout as trap_layout
import trapcam_analysis_in_progress_2019_july as t
import anemometer_analysis as anemometer_analysis

# this 2nd version of this script is geared towards Michael's idea of downwindedness calculation -- declare 90 degrees of downwind-most traps and 90 degrees of upwind-most traps as the "downwind" and "upwind" bins. then calculate the ratio of downwind_bin_catches/downwind+upwind_bins_catches

# the exact choice of 90 degrees is somewhat arbitrary, but the overall idea is great because its level of precision is well-matched to our understanding of how to meaningfully summarize the real-world wind.

def adjust_spines(ax_handle, spines):
    for loc, spine in ax_handle.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            #spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine
    # turn off ticks where there is no spine
    if 'left' in spines:
        ax_handle.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax_handle.yaxis.set_ticks([])
    if 'bottom' in spines:
        ax_handle.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax_handle.xaxis.set_ticks([])

def giraldo_circvar(alpha,axis=None):
    N = len(alpha)
    R = np.sqrt(np.sum(np.sin(alpha),axis)**2 + np.sum(np.cos(alpha),axis)**2)/N
    V = 1-R
    return V

def sum_vectors_from_angles_and_mags(angle_list, magnitude_list, axis = None):
    angle_list = [x*np.pi/180. for x in angle_list]
    y_component = np.sum(np.sin(angle_list)*magnitude_list,axis)
    print ('y: ' +str(y_component))
    x_component = np.sum(np.cos(angle_list)*magnitude_list,axis)
    print ('x: ' +str(x_component))
    avg_vector_mag = np.sqrt(y_component**2 + x_component**2)/len(angle_list)
    avg_vector_dir = np.arctan2(y_component, x_component)
    return avg_vector_mag, avg_vector_dir

def calculate_downwindedness(alpha, wind_dir):
    N = len(alpha)
    list_of_scaling_factors = [np.abs(np.cos(a-wind_dir)) for a in alpha]
    y_components = [np.sin(a)*list_of_scaling_factors[index] for (index, a) in enumerate(alpha)]
    x_components = [np.cos(a)*list_of_scaling_factors[index] for (index, a) in enumerate(alpha)]
    downwindedness = np.sqrt(np.sum(y_components)**2 + np.sum(x_components)**2)/N
    return downwindedness

def calculate_coarse_binned_downwindedness(alpha, wind_dir, number_of_traps, traps_per_bin = 3):
    trap_angles = np.linspace(0,2*np.pi, number_of_traps, endpoint = False)
    rotated_traps = (trap_angles-wind_dir +np.pi/2.)%(2*np.pi)# rotates all trap angles into frame in which downwind points due east (lazy way to ignore wrapping)
    ang_distance_from_downwind = np.array([np.abs(x-np.pi/2) for x in rotated_traps])
    indices = ang_distance_from_downwind.argsort()
    downwind_traps = trap_angles[indices[0:traps_per_bin]]

    ang_distance_from_upwind = np.array([np.abs(x-3*np.pi/2) for x in rotated_traps])
    upwind_indices = ang_distance_from_upwind.argsort()
    upwind_traps = trap_angles[upwind_indices[0:traps_per_bin]]

    downwind_bin_counts=0
    upwind_bin_counts = 0
    for a in alpha:
        if a in downwind_traps:
            downwind_bin_counts +=1
        if a in upwind_traps:
            upwind_bin_counts +=1
    print (upwind_bin_counts)
    print (downwind_bin_counts)
    return downwind_bin_counts/float(upwind_bin_counts+downwind_bin_counts)

def color(normalized_trap_counts):
    cmap = plt.cm.magma
    v_min = 0
    v_max = 0.5
    norm = matplotlib.colors.Normalize(vmin=v_min ,vmax=v_max)
    dotcolors = plt.cm.ScalarMappable(norm, cmap)
    rgba_dotcolors = [dotcolors.to_rgba(x) for x in normalized_trap_counts]
    return rgba_dotcolors

d = {#'2017_04_15':              {'trap_radius_meters':250, '20min_wind_avg_mph':3.60, 'kappa':2.44, 'circvar':.24, 'scattermarker':'.', 'trap_number':8},
    # '2017_04_30_inner':         {'trap_radius_meters':250, '20min_wind_avg_mph':2.88, 'kappa':1.36, 'circvar':.44, 'scattermarker':'.', 'trap_number' :8},
    # '2017_04_30_outer':         {'trap_radius_meters':1000,'20min_wind_avg_mph':2.88, 'kappa':2.46, 'circvar':.24, 'scattermarker':'.', 'trap_number' :8},
    # '2017_06_16':               {'trap_radius_meters':1000,'20min_wind_avg_mph':4.47, 'kappa':1.90, 'circvar':.32, 'scattermarker':'.', 'trap_number' :8},
    '2017_10_26':               {'trap_radius_meters':1000,'20min_wind_avg_mph':2.80, 'kappa':1.82, 'circvar':.33, 'scattermarker':'.', 'trap_number':8},
    '2019_04_19_first_release': {'trap_radius_meters':1000,'20min_wind_avg_mph':6.20, 'kappa':4.49, 'circvar':.12, 'scattermarker':'.', 'trap_number':10},
    '2019_05_08':               {'trap_radius_meters':1000,'20min_wind_avg_mph':4.70, 'kappa':2.97, 'circvar':.19, 'scattermarker':'.', 'trap_number':10},
    '2019_06_11':               {'trap_radius_meters':1000,'20min_wind_avg_mph':1.39, 'kappa':0.86, 'circvar':.61, 'scattermarker':'.', 'trap_number':10},
    '2019_07_06':               {'trap_radius_meters':1000,'20min_wind_avg_mph':3.90, 'kappa':2.41, 'circvar':.25, 'scattermarker':'.', 'trap_number':10}}

# acurite_anemometer_experiments =
#     {'2017_04_15':              {'trap_radius_meters':250, '20min_wind_avg_mph':3.60, 'kappa':2.44, 'circvar':.24, 'scattermarker':'.', 'trap_number':8},
#     # '2017_04_30_inner':         {'trap_radius_meters':250, '20min_wind_avg_mph':2.88, 'kappa':1.36, 'circvar':.44, 'scattermarker':'.', 'trap_number' :8},
#     # '2017_04_30_outer':         {'trap_radius_meters':1000,'20min_wind_avg_mph':2.88, 'kappa':2.46, 'circvar':.24, 'scattermarker':'.', 'trap_number' :8},
#     '2017_06_16':               {'trap_radius_meters':1000,'20min_wind_avg_mph':4.47, 'kappa':1.90, 'circvar':.32, 'scattermarker':'.', 'trap_number' :8}}

# now adding to the dictionary a vector-summation metric for the wind
sum_vectors_for_x_min_post_release = 20

fig = plt.figure(figsize=(6,4))
ax = plt.subplot(1,1,1)

plot_counter = 0
for experiment in d:

    print ()
    print (experiment)
    anemometer_analyzer = anemometer_analysis.AnemometerAnalyzer(directory = experiment, ax_handle = 'blah', time_shift = 0)
    wind_file_list = anemometer_analyzer.get_filenames(experiment+'/weather_data/metone_anemometer_data/', contains ='.txt', does_not_contain=['~', '.pyc'])
    winddirection_list = []
    windspeed_list =[]
    timestamp_list = []
    print ("MAGNETIC DECLINATION: " +str(anemometer_analyzer.mag_declin_deg))
    with open(wind_file_list[0], 'r') as wind_file:
        for line in wind_file:
            winddirection_list.append(float(line.split(' ')[1])-anemometer_analyzer.mag_declin_deg) # in degrees from TRUE north
            timestamp_list.append(float(line.split(' ')[0])) # unix epoch
            windspeed_list.append(float(line.split(' ')[2])) # in meters per second
    print (len(winddirection_list))
    sec_since_release_list =[]
    for timestamp in timestamp_list:
        sec_since_release_list.append(anemometer_analyzer.get_time_since_release_from_epoch_timestamp(timestamp))
    indices_to_plot=[]
    for index, sec_since_release in enumerate(sec_since_release_list):
        if sec_since_release - anemometer_analyzer.time_shift< 0:
            continue
        if sec_since_release - anemometer_analyzer.time_shift > 60.*sum_vectors_for_x_min_post_release:
            print ('breaking')
            break
        indices_to_plot.append(int(index))

    bin_duration = 15. #5. #always specify in seconds
    binned_wind_dir = []
    binned_windspeeds = []
    # print ('length of winddirection list:' +str(len(winddirection_list)))
    direction_sublist = winddirection_list[indices_to_plot[0]:indices_to_plot[-1]] # in degrees
    speed_sublist = windspeed_list[indices_to_plot[0]:indices_to_plot[-1]] # in m/s
    sec_since_release_sublist = sec_since_release_list[indices_to_plot[0]:indices_to_plot[-1]] # in seconds

    start_time = sec_since_release_sublist[0]
    bin_count = 0

    vector_average_mag, vector_average_dir = sum_vectors_from_angles_and_mags(direction_sublist,speed_sublist)
    d[experiment]['wind_vector_avg_mph'] = vector_average_mag
    print ('vector average: ' + str(vector_average_mag))
    print ('vector average direction:' +str(vector_average_dir*180/np.pi))
    vector_average_downwind_dir = (vector_average_dir +np.pi)%(2*np.pi)
    print ('vector average downwind dir: '+str(vector_average_downwind_dir*180/np.pi))

    ####### now, retrieving the trap count dictionaries for the given experiment ########
    filename = './'+experiment+'/field_parameters.json'
    with open(filename) as f:
        field_param_dicts = json.load(f)
    trap_count_dictionary = field_param_dicts["trap counts"]
    trap_number = d[experiment]['trap_number']
    trap_lettering_list = ['A','B','C','D','E','F','G','H','I','J']

    trap_counts = []
    trap_angle_hist_list = []
    for trap in trap_count_dictionary:
        trap_index = trap_lettering_list.index(trap.split('_')[1])
        angle_to_trap = (2*np.pi/trap_number)*trap_index
        trap_counts.append(trap_count_dictionary[trap])
        trap_angle_hist_list.extend([float(angle_to_trap)]*(trap_count_dictionary[trap]))

    coarse_binning_nums = [1,2,3,4]
    for num in coarse_binning_nums:
        coarse_downwindedness = calculate_coarse_binned_downwindedness(alpha = trap_angle_hist_list, number_of_traps=trap_number,wind_dir = vector_average_downwind_dir, traps_per_bin = num)
        ax.scatter(vector_average_mag, coarse_downwindedness, s = 20, color = 'black', label = 'coarse-binned downwindedness' if plot_counter == 0 else "")
        plot_counter +=1
    downwindedness = calculate_downwindedness(alpha= trap_angle_hist_list, wind_dir=vector_average_downwind_dir)
    ax.scatter(vector_average_mag, downwindedness, s = 20, color = 'gray', label = 'downwindedness' if plot_counter == len(coarse_binning_nums) else '')



ax.set_ylim([0,1.1])
ax.set_xlim([0,3.0])
# plt.title('downwindedness')
ax.set_ylabel('measure of downwindedness')
ax.set_xlabel('vector-averaged windspeed '+str(sum_vectors_for_x_min_post_release)+' min post release')
adjust_spines(ax, spines=['bottom','left'])
plt.legend(loc = 2)
plt.tight_layout()

plt.savefig('./coarse_and_fine_binned_downwindedness.svg')
plt.show()

# plt.show()
