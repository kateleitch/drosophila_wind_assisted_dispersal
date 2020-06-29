from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit
from scipy.stats import circmean
from scipy.stats import vonmises
from scipy.stats import circvar
import trap_layout as trap_layout
import trapcam_analysis_in_progress_2019_july as t
import anemometer_analysis as anemometer_analysis
from operator import itemgetter
sys.path.append('/home/kate/src/flytrap_model')
import flytrap_model

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
def find_first_index_of_longest_string_of_ones(sequence):
    first_index_and_bout_length = {}
    bout_length = 0
    for index, value in enumerate(sequence):
        if value ==1:
            if bout_length == 0:
                first_index = index
            if index +1 == len(sequence): #if we're at the end of the sequence yet still in a supra-threshold bout:
                first_index_and_bout_length[first_index] = bout_length
            bout_length +=1
        else: #if we're not currently in a supra-threshold bout
            if bout_length >0: #... and if we just stepped out of a supra-thresh bout
                first_index_and_bout_length[first_index] = bout_length
            bout_length = 0

    sequence_index = max(first_index_and_bout_length.iteritems(), key=itemgetter(1))[0]
    return sequence_index

def adjust_spines(ax_handle, spines):
    for loc, spine in ax_handle.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)
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

def calculate_latency_to_arrival_wave(flies_on_trap, sec_since_release, baseline_cutoff_seconds = 0, sigma_cutoff= 2, window_size = 10):
    low_pass_flies_on_trap = np.zeros(len(flies_on_trap)-window_size)
    for i in range (window_size, len(flies_on_trap)):
        low_pass_flies_on_trap[i-window_size] = (np.mean(flies_on_trap[i-window_size:i]))

    baseline_indices = np.where([x <baseline_cutoff_seconds for x in np.array(sec_since_release)])
    last_baseline_index = baseline_indices[0][-1]
    print ('last baseline time: ' + str(sec_since_release[last_baseline_index]))
    baseline_flies_on_trap = flies_on_trap[0:last_baseline_index]
    baseline_flies_on_trap_mean = np.mean(baseline_flies_on_trap)
    print ('mean: '+ str(baseline_flies_on_trap_mean))
    baseline_flies_on_trap_std = np.std(baseline_flies_on_trap)
    print ('std: ' +str(baseline_flies_on_trap_std))

    test_flies_on_trap = low_pass_flies_on_trap

    threshold = np.max([baseline_flies_on_trap_mean + baseline_flies_on_trap_std*sigma_cutoff, 0.5])
    print ('threshold: ' +str(threshold))

    ###########################trying a schmitt trigger ############################################
    #boolean_lowpassed_flies_on_trap = np.where(np.array(test_flies_on_trap) > threshold, 1, 0)
    upper_thresh = threshold*1.0 #np.max([baseline_flies_on_trap_mean + baseline_flies_on_trap_std*sigma_cutoff*1.0, 0.5])
    lower_thresh = threshold*0.2 #np.max([baseline_flies_on_trap_mean + baseline_flies_on_trap_std*sigma_cutoff*0.2, 0.25])
    boolean_lowpassed_flies_on_trap=[]
    toggle = False
    for value in low_pass_flies_on_trap:
        if not toggle:
            if value > upper_thresh:
                boolean_lowpassed_flies_on_trap.append(1)
                toggle = True
            else:
                boolean_lowpassed_flies_on_trap.append(0)
        if toggle:
            if value < lower_thresh:
                boolean_lowpassed_flies_on_trap.append(0)
                toggle = False
            else:
                boolean_lowpassed_flies_on_trap.append(1)

    index_of_trap_arrival_wave = find_first_index_of_longest_string_of_ones(boolean_lowpassed_flies_on_trap)
    return [sec_since_release[index_of_trap_arrival_wave+window_size], low_pass_flies_on_trap, upper_thresh, lower_thresh, boolean_lowpassed_flies_on_trap]


d = {'2017_10_26':              {'trap_num':8, 'trap_radius_meters':1000., 'traps_with_decent_machine_vision': ['trap_A', 'trap_F', 'trap_G', 'trap_H']},
    '2019_04_19_first_release': {'trap_num':10,'trap_radius_meters':1000., 'traps_with_decent_machine_vision': []},
    '2019_05_08':               {'trap_num':10,'trap_radius_meters':1000., 'traps_with_decent_machine_vision': ['trap_C']},
    '2019_06_11':               {'trap_num':10,'trap_radius_meters':1000., 'traps_with_decent_machine_vision': ['trap_F', 'trap_G']},
    '2019_07_06':               {'trap_num':10,'trap_radius_meters':1000., 'traps_with_decent_machine_vision': ['trap_A', 'trap_B', 'trap_J']}}

d_scored_by_eye ={'2017_10_26': {'trap_num':8, 'trap_radius_meters':1000., '2nd_fly':
    {'trap_A':1137,'trap_B':1194,'trap_C':1190,'trap_D':1674,'trap_E':3279,'trap_F':547,'trap_G':749,'trap_H':403},
    '1st_fly':{'trap_A':1080,'trap_B':1080,'trap_C':1140,'trap_D':1442,'trap_E':655,'trap_F':544,'trap_G':665,'trap_H':391}},
    '2019_04_19_first_release': {'trap_num':10,'trap_radius_meters':1000.,'2nd_fly':{'trap_H':488,'trap_I':482,'trap_J':735}, '1st_fly':{'trap_H':472,'trap_I':468,'trap_J':520}},
    '2019_05_08': {'trap_num':10,'trap_radius_meters':1000.,'2nd_fly':{'trap_C':644,'trap_D':963,'trap_E':1174,'trap_F':3640,'trap_J':4451}, '1st_fly':{'trap_C':395,'trap_D':844,'trap_E':1142,'trap_F':3640,'trap_J':4414}},
    '2019_06_11':               {'trap_num':10,'trap_radius_meters':1000.,'2nd_fly':{'trap_A':933,'trap_B':912,'trap_C':1275,'trap_D':868,'trap_E':898,'trap_F':738,'trap_G':1122,'trap_H':1099,'trap_I':1735,'trap_J':955},'1st_fly':
     {'trap_A':674,'trap_B':856,'trap_C':968,'trap_D':844,'trap_E':689,'trap_F':587,'trap_G':839,'trap_H':968,'trap_I':426,'trap_J':929} },
    '2019_07_06':{'trap_num':10,'trap_radius_meters':1000.,'2nd_fly':{'trap_A':662,'trap_B':464,'trap_I':1860,'trap_J':2214},
    '1st_fly':{'trap_A':573,'trap_B':382,'trap_I':1641,'trap_J':1081} }}

# d_scored_by_eye ={
#     '2019_05_08': {'trap_num':1,'trap_radius_meters':1000.,'2nd_fly':{'trap_C':644}, '1st_fly':{'trap_C':395}},
#     '2019_06_11':               {'trap_num':1,'trap_radius_meters':1000.,'2nd_fly':{'trap_G':1122},'1st_fly':
#      {'trap_G':839}},
#      }

fig = plt.figure(figsize=(10,10))
ax2 = plt.subplot(1,1,1)
color_list = ['red','blue','black','cyan','green']
color_list = ['black']*10
plt.set_cmap('jet')

all_windspeeds_along_trajectory=[]
all_groundspeeds_along_trajectory=[]
all_trap_names =[]

for experiment_date in d:
    print (experiment_date)
    with open(experiment_date+'/field_parameters.json') as f:
        field_parameters = json.load(f)
    trap_catch_dictionary = field_parameters["trap counts"]  # dictionary, keyed by trap name (e.g. trap_A)
    with open(experiment_date+'/all_traps_final_analysis_output.json') as g:
        trap_arrival_dictionary = json.load(g)

# STEP ONE: for each trap, determine arrival time of first fraction of arrivers
#     using distance to trap, calculate the minimal ground speed (assuming straight flight to trap, and assuming minimal plume tracking of these lucky first few flies).

    moving_avg_window_size = 10
    for trap in d[experiment_date]['traps_with_decent_machine_vision']:
        print ()
        print (trap)
        seconds_since_release = trap_arrival_dictionary[trap]['seconds since release:']
        flies_on_trap_over_time = trap_arrival_dictionary[trap]['flies on trap over time:']
        latency_list = calculate_latency_to_arrival_wave(flies_on_trap = flies_on_trap_over_time, sec_since_release=seconds_since_release, baseline_cutoff_seconds = 500, sigma_cutoff= 2.5, window_size = moving_avg_window_size)

        print ('wave calculated at '+str(latency_list[0])+' seconds post release')

        # fig = plt.figure(figsize=(20,5))
        # ax = plt.subplot(1,1,1)
        # ax.set_title(experiment_date +' '+ trap)
        #
        # #above_thresh_indices = [x+moving_avg_window_size-1 for x in latency_list[4]]
        # #above_thresh_indices = [x for x in latency_list[4]]
        # thresholded_data = latency_list[4]
        # #ax.scatter([seconds_since_release[value] for value in above_thresh_indices],[flies_on_trap_over_time[value] for value in above_thresh_indices], color = 'yellow')
        #
        # for index,datapoint in enumerate(flies_on_trap_over_time):
        #     scattercolor = 'black'
        #     if thresholded_data[index-moving_avg_window_size] ==1:
        #         scattercolor = 'green'
        #     ax.scatter(seconds_since_release[index],datapoint, color = scattercolor)
        # #ax.scatter(seconds_since_release[x -moving_avg_window_size for x in latency_list[3]:], flies_on_trap_over_time[x -moving_avg_window_size for x in latency_list[3]:])
        # ax.plot(seconds_since_release[moving_avg_window_size:], latency_list[1], color ='red')
        # ax.set_xlabel('seconds since release')
        # ax.axvline(latency_list[0])
        # ax.axhline(latency_list[2])
        # ax.axhline(latency_list[3])
        # plt.show()

        arrival_wave_time = latency_list[0]
        min_ground_speed = d[experiment_date]['trap_radius_meters']/arrival_wave_time
        print ('minimal ground speed of first arrivers: %02f m/s' %(min_ground_speed))

        # plt.show()

        # STEP TWO: using the arrival time  calculated in the previous step, generate a summary of the anemometer data from time zero to arrival time:
        #     - in say 10-second bins, calculate magnitude of the component of the wind parallel to fly's trajectory to trap.
        #     (Will generate a list of magnitudes, in mph, whose variabity is related to the variability in windspeed and wind direction)
        anemometer_analyzer = anemometer_analysis.AnemometerAnalyzer(directory = experiment_date, ax_handle = 'blah', time_shift = 0, minutes_to_vector_average_wind = 40)
        wind_file_list = anemometer_analyzer.get_filenames(experiment_date+'/weather_data/metone_anemometer_data/', contains ='.txt', does_not_contain=['~', '.pyc'])
        winddirection_list = []
        windspeed_list =[]
        timestamp_list = []
        # print ("MAGNETIC DECLINATION: " +str(anemometer_analyzer.mag_declin_deg))
        with open(wind_file_list[0], 'r') as wind_file:
            for line in wind_file:
                winddirection_list.append(float(line.split(' ')[1])-anemometer_analyzer.mag_declin_deg) # in degrees from TRUE north
                timestamp_list.append(float(line.split(' ')[0])) # unix epoch
                windspeed_list.append(float(line.split(' ')[2])) # in meters per second

        sec_since_release_list =[]
        for timestamp in timestamp_list:
            sec_since_release_list.append(anemometer_analyzer.get_time_since_release_from_epoch_timestamp(timestamp))

        indices_to_plot=[]
        for index, sec_since_release in enumerate(sec_since_release_list):
            if sec_since_release - anemometer_analyzer.time_shift< 0:
                continue
            if sec_since_release - anemometer_analyzer.time_shift > arrival_wave_time:
                print ('breaking')
                break
            indices_to_plot.append(int(index))

        bin_duration = 5. #always specify in seconds
        binned_wind_dir = []
        binned_windspeeds = []
        # print ('length of winddirection list:' +str(len(winddirection_list)))
        direction_sublist = winddirection_list[indices_to_plot[0]:indices_to_plot[-1]] # in degrees
        speed_sublist = windspeed_list[indices_to_plot[0]:indices_to_plot[-1]] # in m/s
        sec_since_release_sublist = sec_since_release_list[indices_to_plot[0]:indices_to_plot[-1]] # in seconds

        start_time = sec_since_release_sublist[0]
        bin_count = 0
        #total_bin_number = np.ceil(len(speed_sublist)/(20*bin_duration))
        # print ('length of speed sublist: '+str(np.ceil(len(speed_sublist))))
        # print ('length of direction sublist: ' +str(len(direction_sublist)))
        #######################################3
        #total_seconds = int((arrival_wave_time/60.)*bin_duration)
        total_seconds = int(arrival_wave_time)
        print ('total seconds of data to analyze: ' +str(total_seconds))
        total_bin_number = arrival_wave_time/bin_duration
        print ('total bin number: ' +str(total_bin_number))
        while True:
            start_index = anemometer_analyzer.find_nearest(sec_since_release_sublist, start_time)
            end_time = start_time+bin_duration
            end_index = anemometer_analyzer.find_nearest(sec_since_release_sublist, end_time)

            if bin_count < total_bin_number:
                direction_slice = direction_sublist[start_index:end_index]
                speed_slice = speed_sublist[start_index:end_index] #still in m/s
            else:
                print ('this last bin is not a full %d seconds in length ' %(bin_duration))
                binned_wind_dir.append(circmean(direction_sublist[start_index:-1], high = 360, low =0))
                binned_windspeeds.append(np.mean(speed_sublist[start_index:-1]))
                break
            binned_wind_dir.append(circmean(direction_slice, high = 360, low =0))
            binned_windspeeds.append(np.mean(speed_slice)) # still in m/s
            start_time = end_time
            bin_count += 1

        # now calculating the angle of the path from the release site to the trap in question

        trap_number = d[experiment_date]['trap_num']
        trap_lettering_list = ['A','B','C','D','E','F','G','H','I','J']
        trap_index = trap_lettering_list.index(trap.split('_')[1])
        angle_to_trap = (360./trap_number)*trap_index
        print ('angle to trap: ' + str(angle_to_trap))

        list_of_headwind_magnitudes = [] # in meters/second
        for index, wind_angle in enumerate(binned_wind_dir[:-1]):
            wind_speed = binned_windspeeds[index]
            if np.isnan(angle_to_trap-wind_angle):
                # print (str(wind_angle))
                continue
            headwind_mag =  -1 * np.cos((angle_to_trap - wind_angle)*np.pi/180)*wind_speed #defined so that, if the wind is pushing a fly along its trajectory, we'll have a positive value
            list_of_headwind_magnitudes.append(headwind_mag)

        headwind_mag_mean = np.mean(list_of_headwind_magnitudes)
        headwind_mag_std = np.std(list_of_headwind_magnitudes)
# step three: plot the result from step one on y axis, and that from step two on x-axis (maybe with horiz bars to indicate variance in wind)
        ax2.errorbar(headwind_mag_mean, min_ground_speed, xerr=headwind_mag_std, fmt='o', label = experiment_date +'  ' +trap)
        ax2.legend(loc =2, prop={'size': 10},scatterpoints = 1)

        all_windspeeds_along_trajectory.append(headwind_mag_mean)
        all_groundspeeds_along_trajectory.append(min_ground_speed)
        all_trap_names.append(experiment_date+'_'+trap)
ax2.plot(np.linspace(0,2,3000), np.linspace(0,2,3000), 'k--', linewidth = 0.5)
ax2.set_xlabel('average wind speed along trajectory to trap, m/s', fontsize =12)
ax2.set_ylabel('minimum ground speed of first arrivers, m/s', fontsize =12)
ax2.tick_params(labelsize=12, top = False, right = False)
adjust_spines(ax2, spines = ['left','bottom'])
ax2.set_ylim(0, 5)
ax2.set_xlim(-2, 3)
# plt.savefig('./groundspeed_vs_windspeed_rescaled_TEST.png')
# plt.savefig('./groundspeed_vs_windspeed_rescaled.svg')
plt.show()

dictionary = {'windspeeds':all_windspeeds_along_trajectory, 'groundspeeds':all_groundspeeds_along_trajectory, 'name':all_trap_names}
with open('./machine_vision_gtraj_wtraj.json', mode = 'w') as g:
    json.dump(dictionary,g, indent = 1)
#############################################################3333
##################################################################################
def giraldo_vec_strength(alpha,axis=None):
#### vector strength = 1 - circvar
    N = len(alpha)
    R = np.sqrt(np.sum(np.sin(alpha),axis)**2 + np.sum(np.cos(alpha),axis)**2)/N
    V = 1-R
    return R

def sum_vectors_from_angles_and_mags(angle_list, magnitude_list, axis = None):
    angle_list = [x*np.pi/180. for x in angle_list]
    angle_list = [x for x in angle_list if not np.isnan(x)]
    magnitude_list = [x for x in magnitude_list if not np.isnan(x)]

    y_component = np.sum(np.sin(angle_list)*magnitude_list,axis)
    print ('y: ' +str(y_component))
    x_component = np.sum(np.cos(angle_list)*magnitude_list,axis)
    print ('x: ' +str(x_component))
    avg_vector_mag = np.sqrt(y_component**2 + x_component**2)/len(angle_list)
    avg_vector_dir = np.arctan2(y_component, x_component)
    return avg_vector_mag, avg_vector_dir

def calc_weighted_vec_strength(direction_list, speed_list):
    scaled_vector_list = np.empty([len(direction_list),2])
    direction_unit_vector_list = [ [np.cos(dir), np.sin(dir)] for dir in direction_list]
    for index, dir_vector in enumerate(direction_unit_vector_list):
        dir_vector = np.array(dir_vector)
        scaled_vector_list[index]=(speed_list[index]*dir_vector)
    scaled_vector_list = np.array(scaled_vector_list)
    vector_sum = np.sum(scaled_vector_list, axis = 0)
    scalar_sum = np.sum(speed_list)
    return np.linalg.norm(vector_sum)/float(scalar_sum)

cmap = plt.cm.Greys
norm = matplotlib.colors.Normalize(vmin=0.65 ,vmax=1)
dotcolors = plt.cm.ScalarMappable(norm, cmap)

# count = 0
# all_windspeeds_along_trajectory = []
# all_groundspeeds_along_trajectory = []
# # all_vector_averaged_windspeeds = []
# vec_avg_dict = {}
# for experiment_date in d_scored_by_eye:
#     plot_color = color_list[count]
#     count +=1
#     print (experiment_date)
#     with open(experiment_date+'/field_parameters.json') as f:
#         field_parameters = json.load(f)
#     trap_catch_dictionary = field_parameters["trap counts"]  # dictionary, keyed by trap name (e.g. trap_A)
#     second_fly_arrival_dict = d_scored_by_eye[experiment_date]['2nd_fly']
#     first_fly_arrival_dict = d_scored_by_eye[experiment_date]['1st_fly']
#
# # STEP ONE: for each trap, determine arrival time of first fraction of arrivers
# #     using distance to trap, calculate the minimal ground speed (assuming straight flight to trap, and assuming minimal plume tracking of these lucky first few flies).
#
#     moving_avg_window_size = 10
#     #for trap in second_fly_arrival_dict:
#     for trap in first_fly_arrival_dict:
#         print ()
#         print (trap)
#         print ('second fly arrival noted at '+str(first_fly_arrival_dict[trap])+' seconds post release')
#         arrival_wave_time = first_fly_arrival_dict[trap]
#         min_ground_speed = d_scored_by_eye[experiment_date]['trap_radius_meters']/arrival_wave_time
#         print ('minimal ground speed of first arrivers: %02f m/s' %(min_ground_speed))
#
#         # plt.show()
#
#         # STEP TWO: using the arrival time  calculated in the previous step, generate a summary of the anemometer data from time zero to arrival time:
#         #     - in say 10-second bins, calculate magnitude of the component of the wind parallel to fly's trajectory to trap.
#         #     (Will generate a list of magnitudes, in mph, whose variabity is related to the variability in windspeed and wind direction)
#         anemometer_analyzer = anemometer_analysis.AnemometerAnalyzer(directory = experiment_date, ax_handle = 'blah', time_shift = 0, minutes_to_vector_average_wind = 40)
#         wind_file_list = anemometer_analyzer.get_filenames(experiment_date+'/weather_data/metone_anemometer_data/', contains ='.txt', does_not_contain=['~', '.pyc'])
#         winddirection_list = []
#         windspeed_list =[]
#         timestamp_list = []
#         # print ("MAGNETIC DECLINATION: " +str(anemometer_analyzer.mag_declin_deg))
#         with open(wind_file_list[0], 'r') as wind_file:
#             for line in wind_file:
#                 winddirection_list.append(float(line.split(' ')[1])-anemometer_analyzer.mag_declin_deg) # in degrees from TRUE north
#                 timestamp_list.append(float(line.split(' ')[0])) # unix epoch
#                 windspeed_list.append(float(line.split(' ')[2])) # in meters per second
#
#         sec_since_release_list =[]
#         for timestamp in timestamp_list:
#             sec_since_release_list.append(anemometer_analyzer.get_time_since_release_from_epoch_timestamp(timestamp))
#
#         indices_to_plot=[]
#         for index, sec_since_release in enumerate(sec_since_release_list):
#             if sec_since_release - anemometer_analyzer.time_shift< 0:
#                 continue
#             if sec_since_release - anemometer_analyzer.time_shift > arrival_wave_time:
#                 print ('breaking')
#                 break
#             indices_to_plot.append(int(index))
#
#         bin_duration = 20. #always specify in seconds
#         binned_wind_dir = []
#         binned_windspeeds = []
#         # print ('length of winddirection list:' +str(len(winddirection_list)))
#         direction_sublist = winddirection_list[indices_to_plot[0]:indices_to_plot[-1]] # in degrees
#         speed_sublist = windspeed_list[indices_to_plot[0]:indices_to_plot[-1]] # in m/s
#         sec_since_release_sublist = sec_since_release_list[indices_to_plot[0]:indices_to_plot[-1]] # in seconds
#
#         start_time = sec_since_release_sublist[0]
#         bin_count = 0
#         #total_bin_number = np.ceil(len(speed_sublist)/(20*bin_duration))
#         # print ('length of speed sublist: '+str(np.ceil(len(speed_sublist))))
#         # print ('length of direction sublist: ' +str(len(direction_sublist)))
#         #######################################3
#         #total_seconds = int((arrival_wave_time/60.)*bin_duration)
#         total_seconds = int(arrival_wave_time)
#         print ('total seconds of data to analyze: ' +str(total_seconds))
#         total_bin_number = arrival_wave_time/bin_duration
#         print ('total bin number: ' +str(total_bin_number))
#         while True:
#             start_index = anemometer_analyzer.find_nearest(sec_since_release_sublist, start_time)
#             end_time = start_time+bin_duration
#             end_index = anemometer_analyzer.find_nearest(sec_since_release_sublist, end_time)
#
#             if bin_count < total_bin_number:
#                 direction_slice = direction_sublist[start_index:end_index]
#                 speed_slice = speed_sublist[start_index:end_index] #still in m/s
#             else:
#                 print ('this last bin is not a full %d seconds in length ' %(bin_duration))
#                 binned_wind_dir.append(circmean(direction_sublist[start_index:-1], high = 360, low =0))
#                 binned_windspeeds.append(np.mean(speed_sublist[start_index:-1]))
#                 break
#             binned_wind_dir.append(circmean(direction_slice, high = 360, low =0))
#             binned_windspeeds.append(np.mean(speed_slice)) # still in m/s
#             start_time = end_time
#             bin_count += 1
#
#         #########
#         #in here, I'll calculate the vector-averaged wind from t = 0 to the time of first arrival at each trap, independent of trajectoryself.
#         #this will be used to generate a distribution of vector-averaged windspeeds for the plumeless models I'm running, to make these models' pdfs more comparable to the set of field data...
#
#         vector_average_mag, vector_average_dir = sum_vectors_from_angles_and_mags(angle_list = binned_wind_dir, magnitude_list = binned_windspeeds)
#         print ('vector averaged windspeed: ' + str(vector_average_mag) + ' m/s')
#         # all_vector_averaged_windspeeds.append(vector_average_mag)
#         vec_avg_dict[experiment_date + '_'+trap] = vector_average_mag
#         #########
#
#         # now calculating the angle of the path from the release site to the trap in question
#
#         trap_number = d_scored_by_eye[experiment_date]['trap_num']
#         trap_lettering_list = ['A','B','C','D','E','F','G','H','I','J']
#         trap_index = trap_lettering_list.index(trap.split('_')[1])
#         angle_to_trap = (360./trap_number)*trap_index
#         print ('angle to trap: ' + str(angle_to_trap))
#
#         list_of_headwind_magnitudes = [] # in meters/second
#         for index, wind_angle in enumerate(binned_wind_dir[:-1]):
#             wind_speed = binned_windspeeds[index]
#             if np.isnan(angle_to_trap-wind_angle):
#                 # print (str(wind_angle))
#                 continue
#             headwind_mag =  -1 * np.cos((angle_to_trap - wind_angle)*np.pi/180)*wind_speed #defined so that, if the wind is pushing a fly along its trajectory, we'll have a positive value
#             list_of_headwind_magnitudes.append(headwind_mag)
#
#         headwind_mag_mean = np.mean(list_of_headwind_magnitudes)
#         headwind_mag_std = np.std(list_of_headwind_magnitudes)
#         #wind_vector_strength = giraldo_vec_strength(np.array([x*np.pi/180. for x in direction_sublist]),axis=None)
#         weighted_vector_strength = calc_weighted_vec_strength(np.array([x*np.pi/180. for x in direction_sublist]), np.array(speed_sublist))
#         print ('wind weighted vector strength: ' +str(weighted_vector_strength))
# # step three: plot the result from step one on y axis, and that from step two on x-axis (maybe with horiz bars to indicate variance in wind)
# #     find some way to retain the identity of each scatter point - color might be distracting
#         dotcolor = dotcolors.to_rgba(weighted_vector_strength)
#         if experiment_date == '2019_05_08':
#             if trap == 'trap_C':
#                 dotcolor = 'red'
#         if experiment_date == '2019_06_11':
#             if trap == 'trap_G':
#                 dotcolor = 'green'
#         ax2.scatter(headwind_mag_mean,min_ground_speed, label = experiment_date, color = dotcolor, edgecolor = 'black', s = 40)
#
#         all_windspeeds_along_trajectory.append(headwind_mag_mean)
#         all_groundspeeds_along_trajectory.append(min_ground_speed)
#
#         #ax2.errorbar(headwind_mag_mean, min_ground_speed, xerr=headwind_mag_std, fmt='o', color = plot_color,label = experiment_date)
# #ax2.legend(loc =2, prop={'size': 10},scatterpoints = 1)
#
# ax2.plot(np.linspace(0,2,3000), np.linspace(0,2,3000), 'k--', linewidth = 1)
# ax2.set_xlabel('average wind speed along trajectory to trap, m/s', fontsize =14)
# ax2.set_ylabel('minimum ground speed of first arrivers, m/s', fontsize =14)
# ax2.tick_params(labelsize=14, top = False, right = False, length=5, width=1)
# # ax2.spines['right'].set_visible(False)
# # ax2.spines['top'].set_visible(False)
# adjust_spines(ax_handle=ax2, spines = ['left','bottom'])
# ax2.set_ylim(0, 5)
# ax2.set_xlim(-2, 3)
# # ax2.set_xlim(-2,2)
# # ax2.axis('equal')
# plt.plot()
#
# #plt.savefig('./groundspeed_vs_windspeed_2nd_fly_arrival_dates_labeled.png')
# # plt.savefig('./groundspeed_vs_windspeed__1st_fly_arrival_wind_vec_strength.svg', transparent = True)
# # plt.savefig('./groundspeed_vs_windspeed__1st_fly_arrival_wind_vec_strength.png', transparent = True)
# plt.show()
#
# dictionary = {'windspeeds':all_windspeeds_along_trajectory, 'groundspeeds':all_groundspeeds_along_trajectory}
# # with open('./first_fly_arrival_groundspeeds_vs_windspeeds.json', mode = 'w') as f:
# #     json.dump(dictionary,f, indent = 1)
# with open('./first_fly_arrival_vector_avg_windspeeds.json', mode = 'w') as g:
#     json.dump(vec_avg_dict,g, indent = 1)
#
# fig = plt.figure(figsize=(5,5))
# ax = plt.subplot(111)
# vector_avg_winds = []
# for key in vec_avg_dict:
#     vector_avg_winds.append(vec_avg_dict[key])
# ax.hist(vector_avg_winds, bins = 10)
# plt.show()
##############################################################################################
################################################################################################
# for experiment_date in d:
#     print (experiment_date)
#     with open(experiment_date+'/field_parameters.json') as f:
#         field_parameters = json.load(f)
#     trap_catch_dictionary = field_parameters["trap counts"]  # dictionary, keyed by trap name (e.g. trap_A)
#     with open(experiment_date+'/all_traps_final_analysis_output.json') as g:
#         trap_arrival_dictionary = json.load(g)
#
# # STEP ONE: for each trap, determine arrival time of first fraction of arrivers
# #     using distance to trap, calculate the minimal ground speed (assuming straight flight to trap, and assuming minimal plume tracking of these lucky first few flies).
#
#     moving_avg_window_size = 10
#     for trap in d[experiment_date]['traps_with_decent_machine_vision']:
#         print ()
#         print (trap)
#         seconds_since_release = trap_arrival_dictionary[trap]['seconds since release:']
#         flies_on_trap_over_time = trap_arrival_dictionary[trap]['flies on trap over time:']
#         latency_list = calculate_latency_to_arrival_wave(flies_on_trap = flies_on_trap_over_time, sec_since_release=seconds_since_release, baseline_cutoff_seconds = 500, sigma_cutoff= 2.5, window_size = moving_avg_window_size)
#
#         print ('wave calculated at '+str(latency_list[0])+' seconds post release')
#
#         fig = plt.figure(figsize=(20,5))
#         ax = plt.subplot(1,1,1)
#         ax.set_title(experiment_date +' '+ trap)
#
#         #above_thresh_indices = [x+moving_avg_window_size-1 for x in latency_list[4]]
#         #above_thresh_indices = [x for x in latency_list[4]]
#         thresholded_data = latency_list[4]
#         #ax.scatter([seconds_since_release[value] for value in above_thresh_indices],[flies_on_trap_over_time[value] for value in above_thresh_indices], color = 'yellow')
#
#         for index,datapoint in enumerate(flies_on_trap_over_time):
#             scattercolor = 'black'
#             if thresholded_data[index-moving_avg_window_size] ==1:
#                 scattercolor = 'green'
#             ax.scatter(seconds_since_release[index],datapoint, color = scattercolor)
#         #ax.scatter(seconds_since_release[x -moving_avg_window_size for x in latency_list[3]:], flies_on_trap_over_time[x -moving_avg_window_size for x in latency_list[3]:])
#         ax.plot(seconds_since_release[moving_avg_window_size:], latency_list[1], color ='red')
#         ax.set_xlabel('seconds since release')
#         ax.axvline(latency_list[0])
#         ax.axhline(latency_list[2])
#         ax.axhline(latency_list[3])
#
#         arrival_wave_time = latency_list[0]
#         min_ground_speed = d[experiment_date]['trap_radius_meters']/arrival_wave_time
#         print ('minimal ground speed of first arrivers: %02f m/s' %(min_ground_speed))
#
#         # plt.show()
#
#         # STEP TWO: using the arrival time  calculated in the previous step, generate a summary of the anemometer data from time zero to arrival time:
#         #     - in say 10-second bins, calculate magnitude of the component of the wind parallel to fly's trajectory to trap.
#         #     (Will generate a list of magnitudes, in mph, whose variabity is related to the variability in windspeed and wind direction)
#         anemometer_analyzer = anemometer_analysis.AnemometerAnalyzer(directory = experiment_date, ax_handle = 'blah', time_shift = 0)
#         wind_file_list = anemometer_analyzer.get_filenames(experiment_date+'/weather_data/metone_anemometer_data/', contains ='.txt', does_not_contain=['~', '.pyc'])
#         winddirection_list = []
#         windspeed_list =[]
#         timestamp_list = []
#         # print ("MAGNETIC DECLINATION: " +str(anemometer_analyzer.mag_declin_deg))
#         with open(wind_file_list[0], 'r') as wind_file:
#             for line in wind_file:
#                 winddirection_list.append(float(line.split(' ')[1])-anemometer_analyzer.mag_declin_deg) # in degrees from TRUE north
#                 timestamp_list.append(float(line.split(' ')[0])) # unix epoch
#                 windspeed_list.append(float(line.split(' ')[2])) # in meters per second
#
#         sec_since_release_list =[]
#         for timestamp in timestamp_list:
#             sec_since_release_list.append(anemometer_analyzer.get_time_since_release_from_epoch_timestamp(timestamp))
#
#         indices_to_plot=[]
#         for index, sec_since_release in enumerate(sec_since_release_list):
#             if sec_since_release - anemometer_analyzer.time_shift< 0:
#                 continue
#             if sec_since_release - anemometer_analyzer.time_shift > arrival_wave_time:
#                 print ('breaking')
#                 break
#             indices_to_plot.append(int(index))
#
#         bin_duration = 5. #always specify in seconds
#         binned_wind_dir = []
#         binned_windspeeds = []
#         # print ('length of winddirection list:' +str(len(winddirection_list)))
#         direction_sublist = winddirection_list[indices_to_plot[0]:indices_to_plot[-1]] # in degrees
#         speed_sublist = windspeed_list[indices_to_plot[0]:indices_to_plot[-1]] # in m/s
#         sec_since_release_sublist = sec_since_release_list[indices_to_plot[0]:indices_to_plot[-1]] # in seconds
#
#         start_time = sec_since_release_sublist[0]
#         bin_count = 0
#         #total_bin_number = np.ceil(len(speed_sublist)/(20*bin_duration))
#         # print ('length of speed sublist: '+str(np.ceil(len(speed_sublist))))
#         # print ('length of direction sublist: ' +str(len(direction_sublist)))
#         #######################################3
#         #total_seconds = int((arrival_wave_time/60.)*bin_duration)
#         total_seconds = int(arrival_wave_time)
#         print ('total seconds of data to analyze: ' +str(total_seconds))
#         total_bin_number = arrival_wave_time/bin_duration
#         print ('total bin number: ' +str(total_bin_number))
#         while True:
#             start_index = anemometer_analyzer.find_nearest(sec_since_release_sublist, start_time)
#             end_time = start_time+bin_duration
#             end_index = anemometer_analyzer.find_nearest(sec_since_release_sublist, end_time)
#
#             if bin_count < total_bin_number:
#                 direction_slice = direction_sublist[start_index:end_index]
#                 speed_slice = speed_sublist[start_index:end_index] #still in m/s
#             else:
#                 print ('this last bin is not a full %d seconds in length ' %(bin_duration))
#                 binned_wind_dir.append(circmean(direction_sublist[start_index:-1], high = 360, low =0))
#                 binned_windspeeds.append(np.mean(speed_sublist[start_index:-1]))
#                 break
#             binned_wind_dir.append(circmean(direction_slice, high = 360, low =0))
#             binned_windspeeds.append(np.mean(speed_slice)) # still in m/s
#             start_time = end_time
#             bin_count += 1
#
#         # now calculating the angle of the path from the release site to the trap in question
#
#         trap_number = d[experiment_date]['trap_num']
#         trap_lettering_list = ['A','B','C','D','E','F','G','H','I','J']
#         trap_index = trap_lettering_list.index(trap.split('_')[1])
#         angle_to_trap = (360./trap_number)*trap_index
#         print ('angle to trap: ' + str(angle_to_trap))
#
#         list_of_headwind_magnitudes = [] # in meters/second
#         for index, wind_angle in enumerate(binned_wind_dir[:-1]):
#             wind_speed = binned_windspeeds[index]
#             if np.isnan(angle_to_trap-wind_angle):
#                 # print (str(wind_angle))
#                 continue
#             headwind_mag =  -1 * np.cos((angle_to_trap - wind_angle)*np.pi/180)*wind_speed #defined so that, if the wind is pushing a fly along its trajectory, we'll have a positive value
#             list_of_headwind_magnitudes.append(headwind_mag)
#
#         headwind_mag_mean = np.mean(list_of_headwind_magnitudes)
#         headwind_mag_std = np.std(list_of_headwind_magnitudes)
# # step three: plot the result from step one on y axis, and that from step two on x-axis (maybe with horiz bars to indicate variance in wind)
# #     find some way to retain the identity of each scatter point - color might be distracting
#
#         #ax2.scatter(headwind_mag_mean,min_ground_speed, label = experiment_date +'  ' +trap, color = 'black', s = 30)
#
#         ax2.errorbar(headwind_mag_mean, min_ground_speed, xerr=headwind_mag_std, fmt='o', label = experiment_date +'  ' +trap)
#         ax2.legend(loc =2, prop={'size': 10},scatterpoints = 1)
#         ax2.plot(np.linspace(0,2,3000), np.linspace(0,2,3000), 'k--', linewidth = 0.5)
#         ax2.set_xlabel('average wind speed along trajectory to trap, m/s', fontsize =12)
#         ax2.set_ylabel('minimum ground speed of first arrivers, m/s', fontsize =12)
#         ax2.tick_params(labelsize=12, top = False, right = False)
#         ax2.spines['right'].set_visible(False)
#         ax2.spines['top'].set_visible(False)
#         #ax2.set_ylim(0,2.7)
#         # ax2.set_xlim(-1.5,2.5)
#
#         #ax2.axis('equal')
#         plt.plot()
#
# plt.savefig('./groundspeed_vs_windspeed__kalman_latency.svg')
# plt.savefig('./groundspeed_vs_windspeed__kalman_latency.png')
# plt.show()
