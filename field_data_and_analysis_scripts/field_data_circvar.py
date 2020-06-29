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
#### vector strength = 1 - circvar
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
    vector_avg = np.sqrt(y_component**2 + x_component**2)/len(angle_list)
    return vector_avg

def color_and_size(latency_dictionary, trap_catch_dictionary):
    cmap = plt.cm.viridis_r
    v_min = min([latency_dictionary[key] for key in latency_dictionary if latency_dictionary[key] is not None])
    print ('vmin: ' +str(v_min))
    v_max = max([latency_dictionary[key] for key in latency_dictionary])
    print ('vmax: ' +str(v_max))
    norm = matplotlib.colors.Normalize(vmin=v_min ,vmax=v_max)
    dotcolors = plt.cm.ScalarMappable(norm, cmap)

    return_dict = {}
    for trap_name in trap_catch_dictionary:
        if trap_name in latency_dictionary:
            print ('latency: ' + str(latency_dictionary[trap_name]))
            if latency_dictionary[trap_name] != None:
                dotcolor = dotcolors.to_rgba(latency_dictionary[trap_name])
            else:
                dotcolor = (0.75, 0.75, 0.75, 1)
        else:
            dotcolor = (0.75, 0.75, 0.75, 1)
        dotsize = 5* trap_catch_dictionary[trap_name]
        return_dict[trap_name.split('_')[-1]]={'dotcolor':dotcolor, 'dotsize':dotsize}
    return return_dict, cmap, norm

def calculate_circvar(dir):
    # analyzer = t.TrapcamAnalyzer(dir)
    with open(dir+'/field_parameters.json') as f:
        field_parameters = json.load(f)
    trap_catch_dictionary = field_parameters["trap counts"]  # dictionary, keyed by trap name (e.g. trap_A)

    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(1,1,1, polar = True)
    trap_layoutter = trap_layout.TrapLayoutVisualizer(directory = dir,
            planned_or_actual = 'actual',
            ax_handle = ax,
            turn_off_text = True,
            return_histogram = True)
    latency_dictionary = {'placeholder':20}
    dot_colors_and_sizes, cmap, norm = color_and_size(latency_dictionary, trap_catch_dictionary)
    (inner_hist_list, outer_hist_list) = trap_layoutter.run(dot_colors_and_sizes= dot_colors_and_sizes, trap_catch_dictionary = trap_catch_dictionary)

    if len(inner_hist_list)>0:
        circvar_to_report = giraldo_circvar(inner_hist_list)
        print ('inner giraldo circvar= '+str(giraldo_circvar(inner_hist_list)))
    if len(outer_hist_list)>0:
        circvar_to_report = giraldo_circvar(outer_hist_list)
        print ('outer giraldo circvar= '+str(giraldo_circvar(outer_hist_list)))
    return circvar_to_report

d = {'2017_10_26': {},
    '2019_04_19_first_release': {},
    '2019_05_08': {},
    '2019_06_11': {},
    '2019_07_06': {}}
sum_vectors_for_x_min_post_release = 20
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
for experiment in d:
    print ()
    print (experiment)
    circvar = calculate_circvar(experiment)
    anemometer_analyzer = anemometer_analysis.AnemometerAnalyzer(directory = experiment, ax_handle = 'blah', time_shift = 0, minutes_to_vector_average_wind = 20.)
    wind_file_list = anemometer_analyzer.get_filenames(experiment+'/weather_data/metone_anemometer_data/', contains ='.txt', does_not_contain=['~', '.pyc'])
    winddirection_list = []
    windspeed_list =[]
    timestamp_list = []
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

    direction_sublist = winddirection_list[indices_to_plot[0]:indices_to_plot[-1]] # in degrees
    speed_sublist = windspeed_list[indices_to_plot[0]:indices_to_plot[-1]] # in m/s
    sec_since_release_sublist = sec_since_release_list[indices_to_plot[0]:indices_to_plot[-1]] # in seconds

    start_time = sec_since_release_sublist[0]
    bin_count = 0

    vector_average = sum_vectors_from_angles_and_mags(direction_sublist,speed_sublist)
    #d[experiment]['wind_vector_avg_meters_per_sec'] = vector_average
    print ('vector average: ' + str(vector_average))
    ax.scatter(vector_average, circvar)
    d[experiment]['vector average'] = vector_average
    d[experiment]['circvar'] = circvar
dictionary_name = './circular_variance_as_func_of_wind_vector_averaged_%d_min.json' %(sum_vectors_for_x_min_post_release)
with open(dictionary_name, mode = 'w') as f:
    json.dump(d,f, indent = 1)
plt.show()
