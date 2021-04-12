from __future__ import print_function
import os
import sys
import scipy
import scipy.signal
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib
import random
import json
from collections import Counter
import pandas as pd
import seaborn as sns
"""

"""
#******* some figure defaults *******
font = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 10}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['xtick.labelsize']=10
plt.rcParams['ytick.labelsize']=10
plt.rcParams['xtick.major.width'] = 0.75
plt.rcParams['xtick.minor.width'] = 0.75
plt.rcParams['ytick.major.width'] = 0.75
plt.rcParams['ytick.minor.width'] = 0.75
plt.rcParams['axes.linewidth']    = 0.75
#************************************
def adjust_spines(ax_handle, spines):
    for loc, spine in ax_handle.spines.items():
        if loc in spines:
            spine.set_position(('outward', 4))  # outward by 10 points
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
def normal(x, mean, sigma):
    return (sigma*np.sqrt(2*np.pi))**-1 * np.exp(-(x-mean)**2/(2*sigma**2))

def normal_kernel_2D (x_array,y_array, mean, sigma_x, sigma_y):
    return np.array(normal(x_array, mean[0], sigma_x)*normal(y_array,mean[1],sigma_y))

def extract_all_windspeeds_and_groundspeeds_from_dictionary(dictionary):
    all_windspeeds_along_trajectories=[]
    all_groundspeeds_along_trajectories=[]
    for windspeed in dictionary:
        for wind_dir in dictionary[windspeed]:
            try:
                if np.isnan(dictionary[windspeed][wind_dir]['windspeeds along trajectories'][0]) == False:
                    all_windspeeds_along_trajectories.append(dictionary[windspeed][wind_dir]['windspeeds along trajectories'][0])
                    all_groundspeeds_along_trajectories.append(dictionary[windspeed][wind_dir]['groundspeeds along trajectories'][0])
            except:
                continue
    if np.isnan(all_windspeeds_along_trajectories).any():
        print ('w')
    if np.isnan(all_groundspeeds_along_trajectories).any():
        print ('g')
    return all_windspeeds_along_trajectories, all_groundspeeds_along_trajectories

def generate_pdf (path_to_data):
    with open(path_to_data) as f:
        dictionary = json.load(f)
    w_sublist, g_sublist =extract_all_windspeeds_and_groundspeeds_from_dictionary(dictionary)
    kernel_sigma = 0.05
    all_gaussians_summed = sum(normal_kernel_2D(X,Y, point, kernel_sigma,kernel_sigma) for point in zip(w_sublist, g_sublist))
    area_under_all_gaussians = float(np.sum(all_gaussians_summed))
    pdf_normalized = all_gaussians_summed/area_under_all_gaussians
    pdf_array = pdf_normalized.tolist()
    pdf_dict = {'pdf array' : pdf_array}
    pdf_name = path_to_data.split('.json')[0] +'_pdf.json'
    ###### saving data ######
    with open(pdf_name, mode = 'w') as f:
        json.dump(pdf_dict,f, indent = 1)
    #######################
def find_nearest_index(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

def log_likelihood_score_non_ratiometric (path_to_pdf, groundspeeds, windspeeds, groundspeeds_outliers_omitted, windspeeds_outliers_omitted):
    with open(path_to_pdf) as f:
        d =  json.load(f)
    pdf = np.asarray(d['pdf array'])
    resample_number = 10000
    list_of_log_likelihood_scores = np.zeros(resample_number)
    list_of_log_likelihood_scores_outliers_omitted = np.zeros(resample_number)
    for i in range(resample_number):
        resampled_indices = np.random.choice(len(groundspeeds),len(groundspeeds), replace = True)
        running_sum_of_log_likelihood_ratios = 0
        for idx in resampled_indices:
            x_index = find_nearest_index(x, windspeeds[idx])
            y_index = find_nearest_index(y, groundspeeds[idx])
            likelihood_value = pdf[y_index, x_index]
            running_sum_of_log_likelihood_ratios += (likelihood_value)#trying without the log, just to make this metric more intuitive for now.
        list_of_log_likelihood_scores[i] = running_sum_of_log_likelihood_ratios

        resampled_indices = np.random.choice(len(groundspeeds_outliers_omitted),len(groundspeeds_outliers_omitted), replace = True)
        running_sum_of_log_likelihood_ratios = 0
        for idx in resampled_indices:
            x_index = find_nearest_index(x, windspeeds_outliers_omitted[idx])
            y_index = find_nearest_index(y, groundspeeds_outliers_omitted[idx])
            running_sum_of_log_likelihood_ratios += pdf[y_index, x_index] #trying without the log, just to make this metric more intuitive for now.
        list_of_log_likelihood_scores_outliers_omitted[i] = running_sum_of_log_likelihood_ratios
    return list_of_log_likelihood_scores, list_of_log_likelihood_scores_outliers_omitted
def get_filenames(path, contains, does_not_contain=['~', '.pyc']):
    cmd = 'ls ' + '"' + path + '"'
    ls = os.popen(cmd).read()
    all_filelist = ls.split('\n')
    try:
        all_filelist.remove('')
    except:
        pass
    filelist = ['']*(len(all_filelist))
    filename_count = 0
    for i, filename in enumerate(all_filelist):
        if contains in filename:
            fileok = True
            for nc in does_not_contain:
                if nc in filename:
                    fileok = False
            if fileok:
                #filelist.append( os.path.join(path, filename) )
                filelist[filename_count] = str(os.path.join(path,filename))
                filename_count +=1
    filelist_trimmed = filelist[0:filename_count-1]
    return filelist_trimmed

def plot_pdf (path_to_pdf, groundspeeds, windspeeds):
    with open(path_to_pdf) as f:
        d =  json.load(f)
    pdf = np.asarray(d['pdf array'])
    fig = plt.figure(figsize= (5,5))
    ax = plt.subplot(111)
    ax.set_ylim(-0.5, 4.6)
    ax.set_xlim(-2, 3.1)
    plt.xticks([-2,-1,0,1,2,3])
    ax.spines['bottom'].set_bounds(-2, 3)
    ax.spines['left'].set_bounds(0, 4)
    plt.yticks([0,1,2,3,4])
    ax.set_xlabel('windspeed along trajectory, m $\mathregular{s^{-1}}$')
    ax.set_ylabel('groundspeed along trajectory, m $\mathregular{s^{-1}}$')
    plt.imshow(pdf, cmap='binary', interpolation='nearest', origin='lower', extent = (-2, 3.5, -0.5, 5)) #cmap was 'binary'
    for idx, windspeed in enumerate(windspeeds):
        ax.scatter(windspeed, groundspeeds[idx], color = 'black',edgecolor = 'white', zorder = 19, s = 30)
    adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
    figname = path_to_pdf.split('_pdf.json')[0]+'.svg'
    plt.savefig(figname)
    plt.close()
###########3
with open('./first_fly_arrival_groundspeeds_vs_windspeeds.json') as f:
    field_data_dictionary = json.load(f)
with open('./first_fly_arrival_groundspeeds_vs_windspeeds_TWO_OUTLIERS_OMITTED.json') as f:
    field_data_dictionary_outliers_omitted = json.load(f)

groundspeeds = field_data_dictionary['groundspeeds']
windspeeds = field_data_dictionary['windspeeds']
groundspeeds_outliers_omitted = field_data_dictionary_outliers_omitted['groundspeeds']
windspeeds_outliers_omitted = field_data_dictionary_outliers_omitted['windspeeds']
x = np.arange(-2, 3.5, 0.05) #0.05
y = np.arange(-0.5, 5, 0.05)
X, Y = np.meshgrid(x, y)

path_to_model = './model_3'#raw_input ('enter path to the model: ')
files_to_analyze = get_filenames(path= path_to_model, contains = '.json', does_not_contain = ['.pyc', '_pdf.json'])
for file in files_to_analyze:
    generate_pdf(file)

    log_likelihoods, log_likelihoods_outliers_omitted = log_likelihood_score_non_ratiometric(file.split('.json')[0]+'_pdf.json', groundspeeds, windspeeds, groundspeeds_outliers_omitted, windspeeds_outliers_omitted)

    avg_likelihood = np.mean(log_likelihoods)
    avg_likelihood_outliers_omitted = np.mean(log_likelihoods_outliers_omitted)

    plot_pdf(file.split('.json')[0]+'_pdf.json', groundspeeds, windspeeds)
#####################################################


#filename = 'density_vs_w_and_g_merged_py2.pkl'
filename = 'density_vs_w_and_g_D30_D1000_tend.pkl'

smoothed = False

# Load data from file
with open(filename, 'rb') as f:
    done = False
    data_list = []
    while not done:
        try:
            data = pickle.load(f)
        except EOFError as err:
            done = True
            break
        data_list.append(data)

gs_n = 4
gs_m = 4
fig = plt.figure(1, figsize = (9,6))
gs = gridspec.GridSpec(gs_n,  gs_m)

# Extract data
for i, data in enumerate(data_list):

    w_traj = data['w_traj']
    g_traj = data['g_traj']
    density = data['density']
    diffusion_coeff = data['diffusion_coeff']
    print('{}/{}, diffusion coeff = {}'.format(i+1, len(data_list), diffusion_coeff))

    # Smooth data using smoothing kernel
    n = 5
    kernel = scipy.ones((n,n),dtype=scipy.float64)
    kernel = kernel/kernel.sum()
    density_smooth = scipy.signal.convolve2d(density,kernel,mode='same')

    ax = plt.subplot(gs[i])
    plt.axis('on')

    m = i%gs_m  + 1
    n = int(i/gs_m) + 1
    if m != 1:
        ax.set_yticklabels([])
    if n != gs_n:
        ax.set_xticklabels([])

    if not smoothed:
        plt.pcolormesh(w_traj, g_traj, density, cmap='binary')

    else:
        plt.pcolormesh(w_traj, g_traj, density_smooth, cmap='binary')

    max_w = w_traj.max()
    min_w = w_traj.min()
    text_x= min_w + 0.02*(max_w - min_w)

    max_g = g_traj.max()
    min_g = g_traj.min()
    text_y = min_g + 0.88*(max_g - min_g)


    plt.text(text_x, text_y, 'D = {:0.0f}'.format(diffusion_coeff))
    if i == 0:
        plt.text(text_x, text_y-0.5, "something's fishy",color='red')

fig.text(0.5, 0.05, 'w_traj (m/s)', ha='center')
fig.text(0.08, 0.5, 'g_traj (m/s)', va='center', rotation='vertical')

plt.show()
