from __future__ import print_function
import sys
import json
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib import gridspec
import trap_layout as trap_layout
import anemometer_analysis as anemometer_analysis
import acurite_anemometer_analysis as acurite_anemometer_analysis


#******* some figure defaults *******
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 9}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['xtick.labelsize']=9
plt.rcParams['ytick.labelsize']=9
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.5
plt.rcParams['axes.linewidth'] = 1.5
#************************************



# acurite_directory_list = ['./2017_04_15', '2017_04_30']

#directory_list  = ['./2019_07_06', './2019_06_11', './2017_10_26', './2019_05_08','./2019_04_19_first_release']
directory_list  = ['./2017_04_15+inner', './2017_04_30+inner', './2017_04_30+outer','./2019_07_06', './2019_06_11', './2017_10_26', './2019_04_19_first_release']

anemometer_dictionary = {'./2017_04_15':'acurite', './2017_04_30':'acurite','./2019_07_06':'metone', './2019_06_11':'metone', './2017_10_26':'metone', './2019_05_08':'metone','./2019_04_19_first_release':'metone'}

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
        ax_handle.yaxis.set_ticks_position('right')
        # ax_handle.yaxis.set_ticks([])
    if 'bottom' in spines:
        ax_handle.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax_handle.xaxis.set_ticks([])

def ask_user_which_traps_to_analyze():
    print ('')
    while True:
        analyze_trap_list = []
        letter = raw_input("Enter a trap letter to plot: ")
        analyze_trap_list.append('trap_'+letter)
        while True:
            letter = raw_input("Enter another trap letter to plot, or enter 'go': ")
            if letter == 'go':
                break
            else:
                analyze_trap_list.append('trap_'+letter)
        print ('')
        print ('you said you want to plot: ')
        for an_trap in analyze_trap_list:
            print (an_trap)
        user_go_ahead = raw_input("Are those the traps you'd like to plot? (y/n) ")
        if user_go_ahead == 'y':
            return analyze_trap_list
            break
        if user_go_ahead == 'n':
            continue
    print ('')

def calculate_latency_to_arrival_wave(flies_on_trap, flies_in_trap, sec_since_release, baseline_cutoff_seconds = 0, sigma_cutoff= 2, window_size = 10):
    #here I'm just calculating a coarse metric, but in a few days I'd like to revisit this and come up with something more principled
    sec_since_release = sec_since_release[window_size:]
    baseline_indices = np.where([x <baseline_cutoff_seconds for x in np.array(sec_since_release)])
    low_pass_flies_on_trap = np.zeros(len(flies_on_trap)-window_size)
    low_pass_flies_in_trap = np.zeros(len(flies_in_trap)-window_size)
    for i in range (window_size, len(flies_on_trap)):
        low_pass_flies_on_trap[i-window_size] = (np.mean(flies_on_trap[i-window_size:i]))
        low_pass_flies_in_trap[i-window_size] = (np.mean(flies_in_trap[i-window_size:i]))
    last_baseline_index = baseline_indices[0][-1]
    print ('last baseline index: ' + str(last_baseline_index))
    baseline_flies_on_trap = low_pass_flies_on_trap[0:last_baseline_index]
    baseline_flies_on_trap_mean = np.mean(baseline_flies_on_trap)
    baseline_flies_on_trap_std = np.std(baseline_flies_on_trap)
    test_flies_on_trap = low_pass_flies_on_trap[last_baseline_index:]
    threshold = baseline_flies_on_trap_mean+baseline_flies_on_trap_std*sigma_cutoff
    print ('threshold: ' +str(threshold))
    on_trap_arrival_wave_indices = np.where([x > threshold  for x in np.array(test_flies_on_trap)])[0]
    try:
        index_of_trap_arrival_wave = on_trap_arrival_wave_indices[0]# < --reports first index at which low-passed fly signal rises above sigma*std of baseline
        print ('on-trap wave calculated at index ' +str(index_of_trap_arrival_wave))
        return [sec_since_release[index_of_trap_arrival_wave], low_pass_flies_on_trap, low_pass_flies_in_trap]
    except:
        print ('no arrival wave calculated for this trap')
        return [None, low_pass_flies_on_trap, low_pass_flies_in_trap]

def find_first_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_first_above_a_threshold(array,threshold):
    above_array = []
    for value in array:
        if value > threshold:
            above_array.append(1)
        else:
            above_array.append(0)
    first_crossing_index = find_first_nearest(above_array,1)
    return first_crossing_index

def calculate_latency_to_some_percent_of_max(flies_on_trap, flies_in_trap, sec_since_release, percent_of_max = 37, window_size = 10):
    sec_since_release = sec_since_release[window_size:]
    low_pass_flies_on_trap = np.zeros(len(flies_on_trap)-window_size)
    low_pass_flies_in_trap = np.zeros(len(flies_in_trap)-window_size)
    for i in range (window_size, len(flies_on_trap)):
        low_pass_flies_on_trap[i-window_size] = (np.mean(flies_on_trap[i-window_size:i]))
        low_pass_flies_in_trap[i-window_size] = (np.mean(flies_in_trap[i-window_size:i]))
    print ('low pass flies on trap max: '+ str(low_pass_flies_on_trap.max()))
    on_trap_thresh = (percent_of_max/100.)*(low_pass_flies_on_trap.max())
    in_trap_thresh = (percent_of_max/100.)*(low_pass_flies_in_trap.max())
    print ('on trap thresh: '+str(on_trap_thresh))
    on_trap_index = find_first_above_a_threshold(low_pass_flies_on_trap, on_trap_thresh)
    in_trap_index = find_first_above_a_threshold(low_pass_flies_in_trap, in_trap_thresh)
    on_trap_time = sec_since_release[on_trap_index]
    in_trap_time = sec_since_release[in_trap_index]
    return [on_trap_time, in_trap_time, low_pass_flies_on_trap, low_pass_flies_in_trap]

def color_and_size(latency_dictionary, trap_catch_dictionary):
    cmap = plt.cm.viridis_r
    v_min = min([latency_dictionary[key] for key in latency_dictionary if latency_dictionary[key] is not None])
    print ('vmin: ' +str(v_min))
    v_max = max([latency_dictionary[key] for key in latency_dictionary])
    print ('vmax: ' +str(v_max))
    norm = mpl.colors.Normalize(vmin=v_min ,vmax=v_max)
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

def round_to_nearest_ten(x):
    return int(np.round(x / 10.0)) * 10

def downwind_vs_upwind_traps (wind_dir, trap_angles):
    rotated_traps = (trap_angles-wind_dir +np.pi/2.)%(2*np.pi)# rotates all trap angles into frame in which upwind points due east (lazy way to ignore wrapping)
    ang_distance_from_upwind = np.array([np.abs(x-np.pi/2) for x in rotated_traps])
    indices = ang_distance_from_upwind.argsort()
    upwind_traps = trap_angles[indices[0:len(trap_angles)/2.]]

    ang_distance_from_downwind = np.array([np.abs(x-3*np.pi/2) for x in rotated_traps])
    downwind_indices = ang_distance_from_downwind.argsort()
    downwind_traps = trap_angles[downwind_indices[0:len(trap_angles)/2.]]
    # downwind_traps = trap_angles[downwind_indices[0:2]]
    return upwind_traps, downwind_traps

#########------------- run
# fig = plt.figure(figsize=(2,2))
ax = plt.subplot(111)
# on_trap_full_hist_list =[]
# in_trap_full_hist_list = []
minutes_to_vector_average_wind = 20.
#
# downwind_on_trap_full_hist_list=[]
# downwind_in_trap_full_hist_list=[]
# upwind_on_trap_full_hist_list=[]
# upwind_in_trap_full_hist_list=[]

#
# wind_speed_cutoff_meters_per_second =

for dir in directory_list:
    print ()
    print ()
    print ()
    print ('............')
    downwind_on_trap_full_hist_list=[]
    downwind_in_trap_full_hist_list=[]
    upwind_on_trap_full_hist_list=[]
    upwind_in_trap_full_hist_list=[]

    directory = dir.split('+')[0]
    try:
        inner_or_outer  = dir.split('+')[1]
    except:
        inner_or_outer = 'single_ring'
    with open(directory+'/all_traps_final_analysis_output.json') as f:
        data = json.load(f)
    with open(directory+'/field_parameters.json') as f:
        field_params = json.load(f)
    trap_catch_dictionary = field_params["trap counts"]
    num_flies_released = field_params["estimated_number_of_flies_released"]
    release_time = field_params["time_of_fly_release"]

    # anemometer_time_shift = int(raw_input("If all went well, you'll be entering zero, but here's where you can declare the number of seconds by which the   anemometer timestamps were advanced: "))
    anemometer_time_shift = 0
    if anemometer_dictionary[directory] == 'metone':
        anemometer_analyzer = anemometer_analysis.AnemometerAnalyzer(directory = directory, ax_handle = ax, time_shift = anemometer_time_shift, minutes_to_vector_average_wind = minutes_to_vector_average_wind)
    if anemometer_dictionary[directory] == 'acurite':
        anemometer_analyzer = acurite_anemometer_analysis.AcuriteAnemometerAnalyzer(directory = directory,ax_handle = ax,time_shift = 0,turn_off_scalebar = False,orient_circmean_from_north = False,turn_off_title = True,turn_off_text = False,bin_duration = 60.,declare_fixed_scale = 0,plot_vectors_as_function_of_time = False,savefig = False)

    degrees_to_rotate, vector_averaged_dir_radians_cw_from_east, vector_averaged_mag_meters_per_second = anemometer_analyzer.run()
    print (vector_averaged_dir_radians_cw_from_east)
    print (inner_or_outer)
    if inner_or_outer == 'single_ring':
        trap_angles = np.linspace(5*np.pi/2.,np.pi/2., len(trap_catch_dictionary), endpoint = False) #
        trap_angles = trap_angles%(2*np.pi)
        trap_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if inner_or_outer == 'inner':
        trap_angles = np.linspace(5*np.pi/2.,np.pi/2., 8, endpoint = False) #
        trap_angles = trap_angles%(2*np.pi)
        trap_names = 'ABCDEFGH'
    if inner_or_outer == 'outer':
        trap_angles = np.linspace(17*np.pi/8.,1*np.pi/8., 4, endpoint = False) #
        trap_angles = trap_angles%(2*np.pi)
        trap_names = 'JLNP'
    trap_name_list = ['trap_'+trap_names[x] for x in range(0, len(trap_angles))]
    # if directory == './2017_04_30': # this experiment had inner and outer traps at 250 and 1000 meters respectively. outer traps offset by 22.5 degrees
    #     inner_trap_angles = np.linspace(5*np.pi/2.,np.pi/2., 8, endpoint = False)
    #     outer_trap_angles = np.linspace(19*np.pi/8.,3*np.pi/8., 8, endpoint = False)
    #     trap_angles = np.hstack((inner_trap_angles,outer_trap_angles))
    #     trap_name_list = ['trap_'+trap_names[x] for x in range(0, len(trap_angles))]

    upwind_trap_angles, downwind_trap_angles = downwind_vs_upwind_traps (wind_dir = vector_averaged_dir_radians_cw_from_east, trap_angles = trap_angles)
    # print ([u*180/np.pi for u in upwind_trap_angles])
    # print ([d*180/np.pi for d in downwind_trap_angles])
    upwind_trap_names = []
    downwind_trap_names = []
    for up_ang in upwind_trap_angles:
        idx = np.where(trap_angles == up_ang)[0]
        upwind_trap_names.append(trap_name_list[idx])
    for down_ang in downwind_trap_angles:
        idx = np.where(trap_angles == down_ang)[0]
        downwind_trap_names.append(trap_name_list[idx])
    print ('upwind traps: ')
    print (upwind_trap_names)
    print ('downwind traps: ')
    print (downwind_trap_names)
    exp_date = directory.split("/")[-1]

    for trap_name in data:

        sec_since_release = data[trap_name]['seconds since release:']
        print (trap_name + ' ' + str(sec_since_release[-1]/60.)+ ' minutes analyzed')
        flies_on_trap = data[trap_name]["flies on trap over time:"]
        flies_in_trap = data[trap_name]["flies in trap over time:"]
        if trap_name in upwind_trap_names:
            for index,sec in enumerate(sec_since_release):
                upwind_on_trap_full_hist_list.extend([sec]*int(flies_on_trap[index]))
                upwind_in_trap_full_hist_list.extend([sec]*int(flies_in_trap[index]))
        if trap_name in downwind_trap_names:
            for index,sec in enumerate(sec_since_release):
                downwind_on_trap_full_hist_list.extend([sec]*int(flies_on_trap[index]))
                downwind_in_trap_full_hist_list.extend([sec]*int(flies_in_trap[index]))
    fig = plt.figure(figsize=(3,5))
    ax = plt.subplot(211)
    axb = ax.twinx()
    ax2 = plt.subplot(212)
    ax2b = ax2.twinx()

    downwind_hist = ax.hist(downwind_on_trap_full_hist_list, color = 'gray', bins = 200, histtype = 'step')
    upwind_hist = axb.hist(upwind_on_trap_full_hist_list, color = 'black', bins = 200, histtype = 'step')
    ax.plot([0,1], [0,1], color = 'gray', label ='downwind traps', linewidth = 1.5)
    ax.plot([0,1], [0,1], color = 'black', label ='upwind traps', linewidth = 1.5)
    ax.set_ylabel('flies on downwind traps')
    axb.set_ylabel('flies on upwind traps')
    ax.legend(loc =1)

    ax2.hist(downwind_in_trap_full_hist_list, color = 'gray', bins = 200, histtype = 'step')
    ax2b.hist(upwind_in_trap_full_hist_list, color = 'black', bins = 200, histtype = 'step')
    ax2.set_xlabel('seconds since fly release')
    ax2.set_ylabel('flies in downwind traps')
    ax2b.set_ylabel('flies in upwind traps')
    adjust_spines(ax, spines = ['bottom', 'left'])
    adjust_spines(ax2, spines = ['bottom', 'left'])
    adjust_spines(axb, spines = ['bottom','right'])
    adjust_spines(ax2b, spines = ['bottom','right'])
    #plt.savefig(directory+'/lumped_time_series_data_' +inner_or_outer+ '.svg')
    plt.show()
