from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import random
import json
import cv2
from matplotlib.path import Path
import matplotlib.patches as patches
import scipy
import statsmodels as sm
from decimal import Decimal

""" adapted for Levy simulations, with simpler output dictionaries"""
#******* some figure defaults *******
font = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 9}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['xtick.labelsize']=9
plt.rcParams['ytick.labelsize']=9
plt.rcParams['xtick.major.width'] = 0.75
plt.rcParams['xtick.minor.width'] = 0.75
plt.rcParams['ytick.major.width'] = 0.75
plt.rcParams['ytick.minor.width'] = 0.75
plt.rcParams['axes.linewidth']    = 0.75
#************************************
def calculate_levels_for_kde_plot(list_of_confidence_intervals, normalized_pdf):
    #https://stackoverflow.com/questions/35225307/set-confidence-levels-in-seaborn-kdeplot
    # Take histogram bin membership as proportional to Likelihood
    # This is true when data comes from a Markovian process # <----- need to evaluate this; ask Kellan
    def objective(limit, target):
        w = np.where(normalized_pdf > limit)
        count = normalized_pdf[w]
        return count.sum() - target
    # Find levels by summing histogram to objective
    list_of_levels = []
    for target in list_of_confidence_intervals:
        lev = scipy.optimize.bisect(objective, normalized_pdf.min(), normalized_pdf.max(), args=(target,))
        list_of_levels.append(lev)
    highest_level=normalized_pdf.max()
    list_of_levels.append(highest_level)
    return list_of_levels

def find_nearest_index(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

def normal(x, mean, sigma):
    return (sigma*np.sqrt(2*np.pi))**-1 * np.exp(-(x-mean)**2/(2*sigma**2))

def normal_kernel_2D (x_array,y_array, mean, sigma_x, sigma_y):
    return np.array(normal(x_array, mean[0], sigma_x)*normal(y_array,mean[1],sigma_y))

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

def extract_all_windspeeds_and_groundspeeds_from_dictionary(dictionary):
    all_windspeeds_along_trajectories=[]
    all_groundspeeds_along_trajectories=[]
    for windspeed in dictionary:
        all_windspeeds_along_trajectories.extend(dictionary[windspeed]['w_traj'])
        all_groundspeeds_along_trajectories.extend(dictionary[windspeed]['g_traj'])
    if np.isnan(all_windspeeds_along_trajectories).any():
        print ('w')
    if np.isnan(all_groundspeeds_along_trajectories).any():
        print ('g')
    return all_windspeeds_along_trajectories, all_groundspeeds_along_trajectories

def plot_model_kde_with_conf_intervals_and_field_data(model_name, list_of_conf_int, number):
    with open(model_name+'pdf.json') as f:
        pdf_dict = json.load(f)
    for key in pdf_dict: # this is sloppy but it works
        pdf = np.asarray(pdf_dict[key])
    print ('sum of pdf: '+str(np.sum(pdf)))
    list_of_levels = calculate_levels_for_kde_plot(list_of_conf_int, pdf)
    fig = plt.figure(figsize= (5,5))
    ax = plt.subplot(111)
    ax.set_ylim(-0.5, 4.6)
    ax.set_xlim(-2, 3.1)
    plt.xticks([-2,-1,0,1,2,3])
    ax.spines['bottom'].set_bounds(-2, 3)
    ax.spines['left'].set_bounds(0, 4)
    plt.yticks([0,1,2,3,4])
    ax.contour(X, Y, pdf, levels = list_of_levels, colors ='black' , zorder = 10)
    ax.set_xlabel('windspeed along trajectory, m $\mathregular{s^{-1}}$')
    ax.set_ylabel('groundspeed along trajectory, m $\mathregular{s^{-1}}$')
    plt.imshow(pdf, cmap='binary', interpolation='nearest', origin='lower', extent = (-2, 3.5, -0.5, 5)) #cmap was 'viridis_r'
    for idx, windspeed in enumerate(windspeeds):
        ax.scatter(windspeed, groundspeeds[idx], color = 'black',edgecolor = 'white', zorder = 19, s = 30)
    adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
    figname = model_name + 'MACHINE_VISION_continuous_pdf.svg'
    plt.savefig(figname)
    plt.show()

#-------------------------------------------------------
with open('/home/kate/Documents/coyote_lake/coyote_lake_simulations/plumeless_simulations/coyote_lake_simple_plumeless_simulations/first_fly_arrival_groundspeeds_vs_windspeeds.json') as f:
    field_data_dictionary = json.load(f)
with open('/home/kate/Documents/coyote_lake/coyote_lake_simulations/plumeless_simulations/coyote_lake_simple_plumeless_simulations/first_fly_arrival_groundspeeds_vs_windspeeds_TWO_OUTLIERS_OMITTED.json') as f:
    field_data_dictionary_outliers_omitted = json.load(f)
groundspeeds_outliers_omitted = field_data_dictionary_outliers_omitted['groundspeeds']
windspeeds_outliers_omitted = field_data_dictionary_outliers_omitted['windspeeds']
groundspeeds = field_data_dictionary['groundspeeds']
windspeeds = field_data_dictionary['windspeeds']
x = np.arange(-2, 3.5, 0.05) #0.05
y = np.arange(-0.5, 5, 0.05)
X, Y = np.meshgrid(x, y)



#now for Bayes factor calculation
fig = plt.figure(figsize= (1.5,2.66))

with open('/home/kate/Documents/coyote_lake/coyote_lake_simulations/plumeless_simulations/coyote_lake_simple_plumeless_simulations/four_models_normalized_pdfs_FIELD_WINDS_flexible_airspeeds_model3_improved.json')as g:
    four_original_models = json.load(g)
sideslip_pdf = np.asarray(four_original_models['./simulated_groundspeeds_vs_windspeeds_along_trajectories_1.0ms_1.8ms_FIELD_WINDS.json'])

with open ('./PDFs_by_mu_and_kappa.json') as f:
    pdfs_by_mu_and_kappa = json.load(f)

for mu_and_kappa in pdfs_by_mu_and_kappa:
    print ('mu = '+mu_and_kappa.split('_')[0][2:])
    print ('   kappa = '+mu_and_kappa.split('_')[1][5:])
    alternate_pdf = np.asarray(pdfs_by_mu_and_kappa[mu_and_kappa])
    resample_number = 40000
    list_of_log_Bayes_factors = np.zeros(resample_number)
    list_of_log_Bayes_factors_outliers_omitted = np.zeros(resample_number)
    for i in range(resample_number):
        resampled_indices = np.random.choice(len(groundspeeds),len(groundspeeds), replace = True)
        list_of_log_likelihood_ratios =[]
        for idx in resampled_indices:
            x_index = find_nearest_index(x, windspeeds[idx])
            y_index = find_nearest_index(y, groundspeeds[idx])
            list_of_log_likelihood_ratios.append(np.log(sideslip_pdf[y_index, x_index]/alternate_pdf[y_index, x_index]))
        log_Bayes_factor = np.sum(list_of_log_likelihood_ratios)
        list_of_log_Bayes_factors[i] = log_Bayes_factor

        resampled_indices = np.random.choice(len(groundspeeds_outliers_omitted),len(groundspeeds_outliers_omitted), replace = True)
        list_of_log_likelihood_ratios_outliers_omitted =[]
        for idx in resampled_indices:
            x_index = find_nearest_index(x, windspeeds_outliers_omitted[idx])
            y_index = find_nearest_index(y, groundspeeds_outliers_omitted[idx])
            list_of_log_likelihood_ratios_outliers_omitted.append(np.log(sideslip_pdf[y_index, x_index]/alternate_pdf[y_index, x_index]))
        log_Bayes_factor_outliers_omitted = np.sum(list_of_log_likelihood_ratios_outliers_omitted)
        list_of_log_Bayes_factors_outliers_omitted[i] = log_Bayes_factor_outliers_omitted

    fig = plt.figure(figsize= (1,1.5))
    ax = plt.subplot(111)
    hist_factor = -2.
    binwidth = 30

    binlist = range(int(min(list_of_log_Bayes_factors)), int(max(list_of_log_Bayes_factors))+binwidth, binwidth)
    ax.hist(list_of_log_Bayes_factors, bins = binlist, color = 'k', weights=(10**hist_factor)*np.ones_like(list_of_log_Bayes_factors), histtype = 'step')

    binlist = range(int(min(list_of_log_Bayes_factors_outliers_omitted)), int(max(list_of_log_Bayes_factors_outliers_omitted))+binwidth, binwidth)
    ax.hist(list_of_log_Bayes_factors_outliers_omitted, bins = binlist, color = 'gray', weights=(10**hist_factor)*np.ones_like(list_of_log_Bayes_factors_outliers_omitted), histtype = 'step')
    print ('Bayes factor mean:'+str(np.mean(list_of_log_Bayes_factors)))
    print ('outliers omitted:' +str(np.mean(list_of_log_Bayes_factors_outliers_omitted)))
    if np.mean(list_of_log_Bayes_factors) > 100:
        ax.set_xlim(0, 700)
        plt.xticks(np.linspace(0, 600, 2, endpoint = True))
        ax.spines['bottom'].set_bounds(0, 600)
        ax.text(600,45, 'ln (L1/L%s)'%(mu_and_kappa))
    if np.mean(list_of_log_Bayes_factors) < 100:
        ax.set_xlim(-300, 400)
        plt.xticks(np.linspace(-300, 300, 2, endpoint = True))
        ax.spines['bottom'].set_bounds(-300, 300)
        ax.text(300,45, 'ln (L1/L%s)'%(mu_and_kappa))


    # ax.text(600,105, str(np.mean(list_of_log_Bayes_factors)) )
    # ax.text(600,75,str(np.mean(list_of_log_Bayes_factors_outliers_omitted)) )
    adjust_spines(ax_handle=ax, spines= ['bottom'])
    print ('')
    figname = './bayes_hist_'+mu_and_kappa +'.svg'
    plt.savefig(figname, transparent = True)
    # plt.show()
