from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import random
import json
import matplotlib.gridspec as gridspec
from collections import Counter
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as gridspec
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
    return pdf_normalized

def find_nearest_index(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

def log_likelihood_score(pdf, ref_pdf, groundspeeds, windspeeds, groundspeeds_outliers_omitted, windspeeds_outliers_omitted):
    resample_number = 10000
    list_of_log_likelihood_scores = np.zeros(resample_number)
    list_of_log_likelihood_scores_outliers_omitted = np.zeros(resample_number)
    for i in range(resample_number):
        resampled_indices = np.random.choice(len(groundspeeds),len(groundspeeds), replace = True)
        running_sum_of_log_likelihood_ratios = 0
        for idx in resampled_indices:
            x_index = find_nearest_index(x, windspeeds[idx])
            y_index = find_nearest_index(y, groundspeeds[idx])
            likelihood_value = (np.log(ref_pdf[y_index, x_index]/pdf[y_index, x_index]))
            running_sum_of_log_likelihood_ratios += (likelihood_value)
        list_of_log_likelihood_scores[i] = running_sum_of_log_likelihood_ratios

        resampled_indices = np.random.choice(len(groundspeeds_outliers_omitted),len(groundspeeds_outliers_omitted), replace = True)
        running_sum_of_log_likelihood_ratios = 0
        for idx in resampled_indices:
            x_index = find_nearest_index(x, windspeeds_outliers_omitted[idx])
            y_index = find_nearest_index(y, groundspeeds_outliers_omitted[idx])
            likelihood_value = (np.log(ref_pdf[y_index, x_index]/pdf[y_index, x_index]))
            running_sum_of_log_likelihood_ratios += (likelihood_value)
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

def plot_pdf (path_to_pdf, groundspeeds, windspeeds, text_str):
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
    ax.set_title(text_str, fontsize = 8)
    plt.imshow(pdf, cmap='binary', interpolation='nearest', origin='lower', extent = (-2, 3.5, -0.5, 5)) #cmap was 'binary'
    for idx, windspeed in enumerate(windspeeds):
        ax.scatter(windspeed, groundspeeds[idx], color = 'black',edgecolor = 'white', zorder = 19, s = 30)
    adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
    figname = path_to_pdf.split('_pdf.json')[0]+'.png'
    print (figname)
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



#model and its specific iteration to be used as a reference for the Bayes factor calculation
dict = {'./model_1': './model_1/airmin_-0.200_airmax_2.000_gpref_1.250.json' }
# dict = {'./model_2': './model_2/forward_airspeed_1.5.json' }
# dict = {'./model_3': './model_3/airmin_-0.2_airmax_2.25_gpref_1.5.json' }
# dict = {'./model_4': './model_4/forward_airspeed_1.5.json' }

results_dict = {}
for key in dict:
    path_to_model = key
    print (path_to_model)
    reference_file = dict[key]
    ref_pdf = generate_pdf(reference_file)
    files_to_analyze = get_filenames(path= path_to_model, contains = '.json', does_not_contain = ['.pyc', '_pdf.json', 'dictionary'])
    # files_to_analyze = get_filenames(path= path_to_model, contains = 'airmax_2.7_gpref_0.125.json', does_not_contain = ['.pyc', '_pdf.json', 'dictionary'])
    all_files_in_directory = get_filenames(path= path_to_model, contains = '.json', does_not_contain = [])

    for num, file in enumerate(files_to_analyze):
        if file.split('.json')[0]+'_pdf.json' in all_files_in_directory:
            print ('already analyzed ')
            continue
        print (str(num)+ ' out of ' + str(len(files_to_analyze))+ ' ' + file)
        pdf = generate_pdf(file)

        log_likelihoods, log_likelihoods_outliers_omitted = log_likelihood_score(pdf, ref_pdf, groundspeeds, windspeeds, groundspeeds_outliers_omitted, windspeeds_outliers_omitted)

        B = np.mean(log_likelihoods)
        B_outliers_omitted = np.mean(log_likelihoods_outliers_omitted)

        # text_str = ('relative to '+reference_file.split('.json')[0].split('/')[2]+ ' B = %0.4f, without outliers = %0.4f') %(B, B_outliers_omitted)
        # plot_pdf(file.split('.json')[0]+'_pdf.json', groundspeeds, windspeeds, text_str)

        results_dict[file] = {'Bayes': B, 'Bayes outliers omitted': B_outliers_omitted}


    results_dict_name = key +'/results_dictionary.json'
    try: #if the dictionary already exisits
        with open(results_dict_name) as f:
            data = json.load(f)
        data.update(results_dict)
        with open(results_dict_name, mode = 'w') as f:
            json.dump(data, f, indent =1)
    except:
        with open(results_dict_name, mode = 'w') as f:
            json.dump(results_dict,f, indent = 1)
