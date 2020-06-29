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
    figname = path_to_pdf.split('_pdf.json')[0]+'.svg'
    plt.savefig(figname)
    plt.close()

#  #model and its iteration that best fits the field data (based on outliers omitted)
# dict = {'./model_1': './model_1/airmin_-0.2_airmax_1.75_gpref_1.0.json',
#         './model_2': './model_2/forward_airspeed_1.0.json',
#         './model_3': './model_3/airmin_-0.2_airmax_1.75_gpref_1.0.json'
#         './model_4': './model_4/forward_airspeed_1.0.json'}

model_list = ['./model_3']
for model_path in model_list:
    results_dict_name = model_path +'/results_dictionary.json'
    with open(results_dict_name) as f:
        d = json.load(f)

    airmin_list = []
    min_B = 1
    max_B = 0
    for file in d:
        string = file.split('/')[-1].split('.json')[0]
        airmin = float(string.split('_')[1])
        if airmin not in airmin_list:
            airmin_list.append(airmin)
        if d[file]['Bayes outliers omitted'] < min_B:
            min_B = d[file]['Bayes outliers omitted']
        if d[file]['Bayes outliers omitted'] > max_B:
            max_B = d[file]['Bayes outliers omitted']
    print (airmin_list)

    for airmin_value in airmin_list:
        data = []
        for file in d:
            B_outliers_omitted =    d[file]['Bayes outliers omitted']
            B =                     d[file]['Bayes']
            string = file.split('/')[-1].split('.json')[0]
            gpref  = float(string.split('_')[-1])
            airmax = float(string.split('_')[3])
            airmin = float(string.split('_')[1])
            if airmin != airmin_value:
                continue
            data.append(
                {'airmax': airmax,
                'gpref': gpref,
                'airmin':airmin,
                'Bayes outliers omitted': B_outliers_omitted})

        dataframe = pd.DataFrame(data)
        result = dataframe.pivot(index='airmax', columns='gpref', values='Bayes outliers omitted')
        fig = plt.figure()
        ax = plt.subplot(111)
        sns.heatmap(data=result, annot = True, fmt = '0.0f', vmin = 0, vmax = 500, cbar_kws={'label': 'log Bayes, within-model, outliers omitted', })

        adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
        figname = model_path+'sensitivity_analysis_airmin_'+str(airmin_value)

        plt.savefig(figname+'.png')
        plt.savefig(figname+'.svg')

#########
# model_list = ['./model_2','./model_4']
#
# for model_path in model_list:
#     results_dict_name = model_path +'_results_dictionary.json'
#     with open(results_dict_name) as f:
#         d = json.load(f)
#
#     forward_airspeed_list = []
#     bayes_outliers_omitted_list = []
#     for file in d:
#         B_outliers_omitted =    d[file]['Bayes outliers omitted']
#         B =                     d[file]['Bayes']
#         string = file.split('/')[-1].split('.json')[0]
#         forward_airspeed = float(string.split('_')[-1])
#
#         forward_airspeed_list.append(forward_airspeed)
#         bayes_outliers_omitted_list.append(B_outliers_omitted)
#
#     fb = zip(forward_airspeed_list, bayes_outliers_omitted_list)
#     fb.sort()
#     b_sorted = [b for f, b in fb]
#     f_sorted = [f for f, b in fb]
#
#     dataframe = pd.DataFrame({'Bayes outliers omitted': b_sorted},
#                 index = f_sorted)
#     # print (dataframe)
#     # result = dataframe.pivot(index = 'forward airspeed', values='Bayes outliers omitted')
#     # print (result)
#     dataframe.sort_index(axis=0, ascending=False)
#     fig = plt.figure()
#     ax = plt.subplot(111)
#     sns.heatmap(data=dataframe, annot = True, cbar_kws={'label': 'log Bayes, within-model, outliers omitted', })
#     adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
#     figname = model_path+'sensitivity_analysis'
#     plt.savefig(figname+'.png')
#     plt.savefig(figname+'.svg')
