from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import json
from scipy.optimize import curve_fit
import trap_layout as trap_layout
import trapcam_analysis_in_progress_2019_july as t
from operator import itemgetter

#******* some figure defaults *******
font = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 7}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['xtick.labelsize']=7
plt.rcParams['ytick.labelsize']=7
plt.rcParams['xtick.major.width'] = 0.75
plt.rcParams['xtick.minor.width'] = 0.75
plt.rcParams['ytick.major.width'] = 0.75
plt.rcParams['ytick.minor.width'] = 0.75
plt.rcParams['axes.linewidth']    = 0.75
#************************************

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



d_machine_vision = {'2017_10_26':              {'trap_num':8, 'trap_radius_meters':1000., 'traps_with_decent_machine_vision': ['trap_A', 'trap_F', 'trap_G', 'trap_H']},
    '2019_04_19_first_release': {'trap_num':10,'trap_radius_meters':1000., 'traps_with_decent_machine_vision': []},
    '2019_05_08':               {'trap_num':10,'trap_radius_meters':1000., 'traps_with_decent_machine_vision': ['trap_C']},
    '2019_06_11':               {'trap_num':10,'trap_radius_meters':1000., 'traps_with_decent_machine_vision': ['trap_F', 'trap_G']},
    '2019_07_06':               {'trap_num':10,'trap_radius_meters':1000., 'traps_with_decent_machine_vision': ['trap_A', 'trap_B']}}

d_scored_by_eye ={'2017_10_26': {'trap_num':8, 'trap_radius_meters':1000., '2nd_fly':
    {'trap_A':1137,'trap_B':1194,'trap_C':1190,'trap_D':1674,'trap_E':3279,'trap_F':547,'trap_G':749,'trap_H':403},
    '1st_fly':{'trap_A':1080,'trap_B':1080,'trap_C':1140,'trap_D':1442,'trap_E':655,'trap_F':544,'trap_G':665,'trap_H':391}},
    '2019_04_19_first_release': {'trap_num':10,'trap_radius_meters':1000.,'2nd_fly':{'trap_H':488,'trap_I':482,'trap_J':735}, '1st_fly':{'trap_H':472,'trap_I':468,'trap_J':520}},
    '2019_05_08': {'trap_num':10,'trap_radius_meters':1000.,'2nd_fly':{'trap_C':644,'trap_D':963,'trap_E':1174,'trap_F':3640,'trap_J':4451}, '1st_fly':{'trap_C':395,'trap_D':844,'trap_E':1142,'trap_F':3640,'trap_J':4414}},
    '2019_06_11':               {'trap_num':10,'trap_radius_meters':1000.,'2nd_fly':{'trap_A':933,'trap_B':912,'trap_C':1275,'trap_D':868,'trap_E':898,'trap_F':738,'trap_G':1122,'trap_H':1099,'trap_I':1735,'trap_J':955},'1st_fly':
     {'trap_A':674,'trap_B':856,'trap_C':968,'trap_D':844,'trap_E':689,'trap_F':587,'trap_G':839,'trap_H':968,'trap_I':426,'trap_J':929} },
    '2019_07_06':{'trap_num':10,'trap_radius_meters':1000.,'2nd_fly':{'trap_A':662,'trap_B':464,'trap_I':1860,'trap_J':2214},
    '1st_fly':{'trap_A':573,'trap_B':382,'trap_I':1641,'trap_J':1081} }}

fig = plt.figure(figsize=(0.8,0.8))
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1)
##################################################################################
sigma = 2.0
example_dict = {'2019_05_08':['trap_C',[174,113,6]],'2019_06_11':['trap_G',[173,47,146]]}
for experiment_date in d_scored_by_eye:
    with open(experiment_date+'/field_parameters.json') as f:
        field_parameters = json.load(f)
    with open(experiment_date+'/all_traps_final_analysis_output.json') as g:
        trap_arrival_dictionary = json.load(g)
    trap_catch_dictionary = field_parameters["trap counts"]  # dictionary, keyed by trap name (e.g. trap_A)
    second_fly_arrival_dict = d_scored_by_eye[experiment_date]['2nd_fly']
    first_fly_arrival_dict = d_scored_by_eye[experiment_date]['1st_fly']

    for trap in first_fly_arrival_dict:
        print ()
        print (experiment_date + '    ' + trap)
        first_arrival_time = first_fly_arrival_dict[trap]
        print ('first fly arrival noted at '+str(first_arrival_time)+' seconds post release')


        if trap in d_machine_vision[experiment_date]['traps_with_decent_machine_vision']:
            seconds_since_release = trap_arrival_dictionary[trap]['seconds since release:']
            flies_on_trap_over_time = trap_arrival_dictionary[trap]['flies on trap over time:']
            latency_list = calculate_latency_to_arrival_wave(flies_on_trap = flies_on_trap_over_time, sec_since_release=seconds_since_release, baseline_cutoff_seconds = 200, sigma_cutoff= sigma, window_size = 10)
            print ('arrival wave calculated at '+str(latency_list[0])+' seconds post release')
            ax.scatter(latency_list[0]/60., first_arrival_time/60., s = 8, color = 'k',zorder = 20)
            ax.scatter(latency_list[0]/60., first_arrival_time/60., s = 20, color = 'white',zorder = 19)
            ax.text(latency_list[0]/60., first_arrival_time/60., experiment_date+'_'+trap)
            #ax.hlines(y=first_arrival_time/60. ,xmin = 0, xmax = latency_list[0]/60., color = 'k',linewidth = 0.3)
            if experiment_date in example_dict:
                if trap == example_dict[experiment_date][0]:
                    color = [x/255. for x in example_dict[experiment_date][1]]
                    print (color)
                    ax.scatter(latency_list[0]/60., first_arrival_time/60., s = 50, color = color,zorder = 11)



ax.set_ylim([4,22])
ax.set_xlim([4,22])
plt.xticks(np.linspace(5, 20, 4, endpoint = True))
ax.spines['bottom'].set_bounds(5,20)
plt.yticks(np.linspace(5, 20, 4, endpoint = True))
ax.spines['left'].set_bounds(5,20)
ax.plot(np.linspace(0,25,100), np.linspace(0,25,100), '--k')

ax.set_ylabel('first arr. (min)')
ax.set_title('arrival wave: time at which flies-on-trap exceeds ' +str(sigma)+' sigma')
ax.set_xlabel('wave arr. (min)')
adjust_spines(ax_handle = ax, spines = ['bottom','left'])

#plt.savefig('./first_arrivals_vs_arrival_waves_'+str(sigma)+' sigma_2019_07_06_trap_J_with_beetles_omitted.svg')
plt.show()
