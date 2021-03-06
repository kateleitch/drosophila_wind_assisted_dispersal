from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import random
import json
# import seaborn as sns
import cv2
from matplotlib.path import Path
import matplotlib.patches as patches
import scipy
import statsmodels as sm
from decimal import Decimal

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

def hex2rgb(hexval):
    r = int(hexval[1:3], 16)/255.
    g = int(hexval[3:5], 16)/255.
    b = int(hexval[5:], 16)/255.
    return [r,g,b]

def color_fade(c, n):
    c_list = [c]
    for i in range(n-1):
        c_list.append([1-(1-x)*3/4 for x in c_list[-1]])
    return c_list[::-1]

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

def plot_model_kde_with_conf_intervals_and_field_data(model_name, list_of_conf_int, number, flexible_airspeeds):
    with open(model_name) as f:
        d =  json.load(f)
    pdf = np.asarray(d)
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
    plt.imshow(pdf, cmap='viridis_r', interpolation='nearest', origin='lower', extent = (-2, 3.5, -0.5, 5)) #cmap was 'binary'
    for idx, windspeed in enumerate(windspeeds):
        ax.scatter(windspeed, groundspeeds[idx], color = 'black',edgecolor = 'white', zorder = 19, s = 30)
    adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
    if flexible_airspeeds:
        figname = 'MACHINE_VISION_model' +str(number)+'FIELD_WINDS_flexible_airspeeds_color_heatmap.svg'
    else:
        #figname = 'model' +str(number)+'FIELD_WINDS_color_heatmap.svg'
        figname = 'MACHINE_VISION_model'+str(number)+'_linear_wind_array.svg'
    plt.savefig(figname)
    plt.show()



# model_list = ['./model_1.json', './model_1_infinite_airmin.json']
model_list = ['./model_2.json', './model_3.json', 'model_3_12000_test_theta.json', 'model_4_flies_simply_hold.json','model_4_no_airspeed_flexibility.json', 'model_4_flies_kick_into_1.8ms.json']

with open('./first_fly_arrival_groundspeeds_vs_windspeeds.json') as f:
    field_data_dictionary = json.load(f)
# with open('./machine_vision_gtraj_wtraj.json') as f:
#     field_data_dictionary = json.load(f)
with open('./first_fly_arrival_groundspeeds_vs_windspeeds_TWO_OUTLIERS_OMITTED.json') as f:
    field_data_dictionary_outliers_omitted = json.load(f)


groundspeeds = field_data_dictionary['groundspeeds']
windspeeds = field_data_dictionary['windspeeds']
x = np.arange(-2, 3.5, 0.05) #0.05
y = np.arange(-0.5, 5, 0.05)
X, Y = np.meshgrid(x, y)


generate_pdfs_from_wtraj_gtraj_json = raw_input('Do you want to generate PDFs from wtraj/gtraj files? (y/n)  ')
if generate_pdfs_from_wtraj_gtraj_json == 'y':
    for model in model_list:
        with open(model) as f:
            dictionary = json.load(f)
        all_windspeeds_along_trajectories, all_groundspeeds_along_trajectories=extract_all_windspeeds_and_groundspeeds_from_dictionary(dictionary)

        f, ax= plt.subplots(figsize=(10,10))
        #ax.scatter(all_windspeeds_along_trajectories, all_groundspeeds_along_trajectories, color = [0,0,0,0.02], zorder = 25)

        kernel_sigma = 0.05
        contour_number_for_plot = 20
        greys = color_fade(hex2rgb("#000000"), contour_number_for_plot)
        w_sublist = all_windspeeds_along_trajectories
        g_sublist = all_groundspeeds_along_trajectories
        all_gaussians_summed = sum(normal_kernel_2D(X,Y, point, kernel_sigma,kernel_sigma) for point in zip(w_sublist, g_sublist))
        area_under_all_gaussians = float(np.sum(all_gaussians_summed))
        pdf_normalized = all_gaussians_summed/area_under_all_gaussians
        print (np.sum(pdf_normalized)) # should sum to 1.0
        # HERE I WANT TO SAVE THE PDF_SUM ARRAY FOR COMPARISONS AT THE END OF THE SCRIPT
        ax.contour(X, Y, pdf_normalized, contour_number_for_plot, colors ='gray' , zorder = 10)

        for idx, windspeed in enumerate(windspeeds):
            x_index = find_nearest_index(x, windspeed)
            y_index = find_nearest_index(y, groundspeeds[idx])
            ax.scatter(windspeed, groundspeeds[idx], color = 'black',zorder = 19)
            probability_string = '%.0E' %(pdf_normalized[y_index,x_index])
            ax.text(windspeed+0.03, groundspeeds[idx]+0.03, probability_string,zorder = 20)
        print (np.max(pdf_normalized))

        ax.set_ylim(-0.5, 5)
        ax.set_xlim(-2, 3.5)
        adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
        ax.set_ylabel('groundspeed along trajectory, m/s')
        ax.set_xlabel('windspeed along trajectory, m/s')
        figname = model[:-5]+'FIELD_WINDS.svg'
        plt.savefig(figname)
        plt.show()
        pdf_array = pdf_normalized.tolist()
        pdf_dict = {model : pdf_array}
        pdf_name = model.split('.json')[0] +'_pdf.json'
        ###### saving data ######
        with open(pdf_name, mode = 'w') as f:
            json.dump(pdf_dict,f, indent = 1)
        #######################

else:
    print ('OK, going straight to continuous PDF plotting, and log Bayes calculation.')

    groundspeeds = field_data_dictionary['groundspeeds']
    windspeeds = field_data_dictionary['windspeeds']

    groundspeeds_outliers_omitted = field_data_dictionary_outliers_omitted['groundspeeds']
    windspeeds_outliers_omitted = field_data_dictionary_outliers_omitted['windspeeds']

    for number, model in enumerate(model_list):
        plot_model_kde_with_conf_intervals_and_field_data(model_name=model, list_of_conf_int=[0.99], number = number +1, flexible_airspeeds = True)


    #now for Bayes factor calculation

    alternate_model_list = ['./model_2.json',
                        './model_3.json',
                        './model_4_flies_simply_hold.json']
    reference_model_name = './model_1.json'
    with open(reference_model_name) as f:
        ref_dict =  json.load(f)
    for key in ref_dict:
        sideslip_pdf = np.asarray(ref_dict[key])
    fig = plt.figure(figsize= (1.5,2))
    for number, model in enumerate(alternate_model_list):
        print ('alt. model: ' + model)
        with open(model) as f:
            model_dict =  json.load(f)
        for key in model_dict:
            alternate_pdf = np.asarray(model_dict[key])
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

        ax = plt.subplot(3,1,number+1)
        hist_factor = -2.
        binwidth = 30

        binlist = range(int(min(list_of_log_Bayes_factors)), int(max(list_of_log_Bayes_factors))+binwidth, binwidth)
        ax.hist(list_of_log_Bayes_factors, bins = binlist, color = 'k', weights=(10**hist_factor)*np.ones_like(list_of_log_Bayes_factors), histtype = 'step')

        ax.set_ylim(0, 130)
        ax.set_xlim(-750, 1300)
        plt.xticks(np.linspace(-600, 1200, 4, endpoint = True))

        adjust_spines(ax_handle=ax, spines= ['bottom'])
        ax.text(800,45, 'ln (L1/L%d)'%(number+2))
        print (np.mean(list_of_log_Bayes_factors))
        print ('')
    figname = 'log_bayes_factor_distributions.svg'
    plt.savefig(figname)
    # plt.show()
