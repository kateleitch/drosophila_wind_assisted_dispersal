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

# dir = raw_input('type the experiment directory to analyze: ')
# #dir = sys.argv[1]
# analyzer = t.TrapcamAnalyzer(dir)
# with open(dir+'/field_parameters.json') as f:
#     field_parameters = json.load(f)
# trap_catch_dictionary = field_parameters["trap counts"]  # dictionary, keyed by trap name (e.g. trap_A)
#
#
# fig = plt.figure(figsize=(8,8))
# ax = plt.subplot(1,1,1, polar = True)
# trap_layoutter = trap_layout.TrapLayoutVisualizer(directory = dir,
#         planned_or_actual = 'actual',
#         ax_handle = ax,
#         turn_off_text = True,
#         return_histogram = True)
# latency_dictionary = {'placeholder':20}
# dot_colors_and_sizes, cmap, norm = color_and_size(latency_dictionary, trap_catch_dictionary)
# (inner_hist_list, outer_hist_list) = trap_layoutter.run(dot_colors_and_sizes= dot_colors_and_sizes, trap_catch_dictionary = trap_catch_dictionary)
#
# x = np.linspace(0, 2*np.pi, 1000)
# fig = plt.figure(figsize=(8,8))
# ax = plt.subplot(1,1,1, polar = True)
# ax.set_theta_offset(np.pi/2)
# ax.set_theta_direction(-1)
# # ax = plt.subplot(1,1,1)
#
# if len(inner_hist_list)>0:
#     (inner_kappa, inner_loc, inner_scale) = vonmises.fit(inner_hist_list, fscale=1)
#     print ('inner_kappa= '+str(inner_kappa))
#     print ('inner_loc= ' +str(inner_loc))
#     print ('inner circvar= ' +str(circvar(inner_hist_list)))
#     print ('inner giraldo circvar= '+str(giraldo_circvar(inner_hist_list)))
#     ax.hist(inner_hist_list, normed=True,alpha=0.2)
# if len(outer_hist_list)>0:
#     (outer_kappa, outer_loc, outer_scale) = vonmises.fit(outer_hist_list, fscale =1)
#     print ('outer_kappa= '+str(outer_kappa))
#     print ('outer_loc= ' +str(outer_loc))
#     print ('outer_scale= '+str(outer_scale))
#     print ('outer circvar= ' +str(circvar(outer_hist_list)))
#     print ('outer giraldo circvar= '+str(giraldo_circvar(outer_hist_list)))
#     ax.hist(outer_hist_list, bins = 10, alpha=0.2)
#     ax.plot(x, 100*vonmises.pdf(x, outer_kappa, loc = outer_loc),'r-', lw=5, alpha=0.6, label='vonmises pdf')
# plt.show()

d = {'2017_04_15':              {'trap_radius_meters':250, '20min_wind_avg_mph':3.60, 'kappa':2.44, 'circvar':.24, 'scattermarker':'.'},
    '2017_04_30_inner':         {'trap_radius_meters':250, '20min_wind_avg_mph':2.88, 'kappa':1.36, 'circvar':.44, 'scattermarker':'v'},
    '2017_04_30_outer':         {'trap_radius_meters':1000,'20min_wind_avg_mph':2.88, 'kappa':2.46, 'circvar':.24, 'scattermarker':'v'},
    '2017_06_16':               {'trap_radius_meters':1000,'20min_wind_avg_mph':4.47, 'kappa':1.90, 'circvar':.32, 'scattermarker':'<'},
    '2017_10_26':               {'trap_radius_meters':1000,'20min_wind_avg_mph':2.80, 'kappa':1.82, 'circvar':.33, 'scattermarker':'p'},
    '2019_04_19_first_release': {'trap_radius_meters':1000,'20min_wind_avg_mph':6.20, 'kappa':4.49, 'circvar':.12, 'scattermarker':'d'},
    '2019_05_08':               {'trap_radius_meters':1000,'20min_wind_avg_mph':4.70, 'kappa':2.97, 'circvar':.19, 'scattermarker':'*'},
    '2019_06_11':               {'trap_radius_meters':1000,'20min_wind_avg_mph':1.39, 'kappa':0.86, 'circvar':.61, 'scattermarker':'s'},
    '2019_07_06':               {'trap_radius_meters':1000,'20min_wind_avg_mph':3.90, 'kappa':2.41, 'circvar':.25, 'scattermarker':'>'}}

d = {#'2017_04_15':              {'trap_radius_meters':250, '20min_wind_avg_mph':3.60, 'kappa':2.44, 'circvar':.24, 'scattermarker':'.'},
    #'2017_04_30_inner':         {'trap_radius_meters':250, '20min_wind_avg_mph':2.88, 'kappa':1.36, 'circvar':.44, 'scattermarker':'.'},
    #'2017_04_30_outer':         {'trap_radius_meters':1000,'20min_wind_avg_mph':2.88, 'kappa':2.46, 'circvar':.24, 'scattermarker':'.'},
    #'2017_06_16':               {'trap_radius_meters':1000,'20min_wind_avg_mph':4.47, 'kappa':1.90, 'circvar':.32, 'scattermarker':'.'},
    '2017_10_26':               {'trap_radius_meters':1000,'20min_wind_avg_mph':2.80, 'kappa':1.82, 'circvar':.33, 'scattermarker':'.'},
    '2019_04_19_first_release': {'trap_radius_meters':1000,'20min_wind_avg_mph':6.20, 'kappa':4.49, 'circvar':.12, 'scattermarker':'.'},
    '2019_05_08':               {'trap_radius_meters':1000,'20min_wind_avg_mph':4.70, 'kappa':2.97, 'circvar':.19, 'scattermarker':'.'},
    '2019_06_11':               {'trap_radius_meters':1000,'20min_wind_avg_mph':1.39, 'kappa':0.86, 'circvar':.61, 'scattermarker':'.'},
    '2019_07_06':               {'trap_radius_meters':1000,'20min_wind_avg_mph':3.90, 'kappa':2.41, 'circvar':.25, 'scattermarker':'.'}}

# now adding to the dictionary a vector-summation metric for the wind
sum_vectors_for_x_min_post_release = 30
for experiment in d:
    print ()
    print (experiment)
    anemometer_analyzer = anemometer_analysis.AnemometerAnalyzer(directory = experiment, ax_handle = 'blah', time_shift = 0)
    wind_file_list = anemometer_analyzer.get_filenames(experiment+'/weather_data/metone_anemometer_data/', contains ='.txt', does_not_contain=['~', '.pyc'])
    winddirection_list = []
    windspeed_list =[]
    timestamp_list = []
    # print ("MAGNETIC DECLINATION: " +str(anemometer_analyzer.mag_declin_deg))
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

    # bin_duration = 5. #always specify in seconds
    # binned_wind_dir = []
    # binned_windspeeds = []
    # print ('length of winddirection list:' +str(len(winddirection_list)))
    direction_sublist = winddirection_list[indices_to_plot[0]:indices_to_plot[-1]] # in degrees
    speed_sublist = windspeed_list[indices_to_plot[0]:indices_to_plot[-1]] # in m/s
    sec_since_release_sublist = sec_since_release_list[indices_to_plot[0]:indices_to_plot[-1]] # in seconds

    start_time = sec_since_release_sublist[0]
    bin_count = 0

    vector_average = sum_vectors_from_angles_and_mags(direction_sublist,speed_sublist)
    d[experiment]['wind_vector_avg_meters_per_sec'] = vector_average
    print ('vector average: ' + str(vector_average))



# now plotting von mises kappa as a function of average wind speed during first 20 minutes after release
fig = plt.figure(figsize=(8,8))
ax = plt.subplot(2,1,1)
for experiment in d:
    if d[experiment]['trap_radius_meters']==250:
        markercolor = 'gray'
    else:
        markercolor = 'black'
    ax.scatter(d[experiment]['wind_vector_avg_meters_per_sec'],d[experiment]['kappa'], s = 120, marker =d[experiment]['scattermarker'],color=markercolor, label = '_'.join(experiment.split('_')[0:3]))
#ax.legend(loc =2, prop={'size': 8},scatterpoints = 1)
# ax.set_xlabel('avg wind speed 20 min post release, mph')
ax.set_ylabel('trap distribution - kappa')
ax.set_ylim(0, 5)
ax.set_xlim(0, 3)
ax.tick_params(labelsize=12)
# ax.text(1.8,4.3,'250 m traps', color = 'gray', size = 10)
# ax.text(1.8,4.0, '1000 m traps', color = 'black', size = 10)
adjust_spines(ax_handle=ax, spines = ['left','bottom'])
plt.title('Field trap distributions vs vector avg windspeed', size = 12)

ax2 = plt.subplot(2,1,2)
for experiment in d:
    if d[experiment]['trap_radius_meters']==250:
        markercolor = 'gray'
    else:
        markercolor = 'black'
    ax2.scatter(d[experiment]['wind_vector_avg_meters_per_sec'],d[experiment]['circvar'], s = 120,marker =d[experiment]['scattermarker'],color=markercolor, label = '_'.join(experiment.split('_')[0:3]))
ax2.set_xlabel('vector-averaged windspeed '+str(sum_vectors_for_x_min_post_release)+' min post release, m/s')
ax2.set_ylabel('trap distribution - circular variance')
ax2.set_ylim(0, 1)
ax2.set_xlim(0, 3)
ax2.tick_params(labelsize=12)
adjust_spines(ax_handle=ax2, spines = ['left','bottom'])
# analyzer.format_matplotlib_ax_object(ax_handle = ax2)

plt.savefig('./trap_distributions_vs_vector_averaged_windspeed_'+str(sum_vectors_for_x_min_post_release)+'_min.png')
plt.savefig('./trap_distributions_vs_vector_averaged_windspeed_'+str(sum_vectors_for_x_min_post_release)+'_min.svg')
plt.show()
