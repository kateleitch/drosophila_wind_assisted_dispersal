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
import pickle

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

def normal_kernel_2D (x_array,y_array,mean,density,sigma_x, sigma_y):
    array = np.array(normal(x_array, mean[0], sigma_x)*normal(y_array,mean[1],sigma_y))
    scaled_array = density*array
    return (scaled_array)

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

def plot_model_kde_with_conf_intervals_and_field_data(pdf, list_of_conf_int, figname):
    # with open(model_name) as f:
    #     d =  json.load(f)
    # pdf = np.asarray(d)
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
    plt.imshow(pdf, cmap='binary', interpolation='nearest', vmax = 0.0025, origin='lower', extent = (-2, 3.5, -0.5, 5)) #cmap was 'binary'
    print (np.max(pdf))
    for idx, windspeed in enumerate(windspeeds):
        ax.scatter(windspeed, groundspeeds[idx], color = 'black',edgecolor = 'white', zorder = 19, s = 30)
    adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])

    #     #figname = 'model' +str(number)+'FIELD_WINDS_color_heatmap.svg'
    #     figname = 'MACHINE_VISION_model'+str(number)+'_linear_wind_array.svg'

    plt.savefig(figname)
    # plt.show()


#filename = 'density_vs_w_and_g_merged_py2.pkl'
#filename2 = 'new_pkl_file.pkl'
# filename = 'density_vs_w_and_g_D30_D1000_tend.pkl'
examine_outlier_omission = True
radius = int(raw_input('what radius are you analyzing? (1000 or 250)'))
if radius == 1000:
    filename = 'density_vs_w_and_g_D30_D1000_tend_v2.pkl'
    with open('./first_fly_arrival_groundspeeds_vs_windspeeds.json') as f:
        field_data_dictionary = json.load(f)
    with open('./first_fly_arrival_groundspeeds_vs_windspeeds_TWO_OUTLIERS_OMITTED.json') as f:
        field_data_dictionary_outliers_omitted = json.load(f)
if radius == 250:
    filename = 'density_vs_w_and_g_D10_D1000_small.pkl'
    with open('./first_fly_arrival_groundspeeds_vs_windspeeds_250m_2017_04_30.json') as f:
        field_data_dictionary = json.load(f)
    examine_outlier_omission = False

groundspeeds = field_data_dictionary['groundspeeds']
windspeeds = field_data_dictionary['windspeeds']
contour_number_for_plot = 25
x_grid = np.arange(-2, 3.5, 0.05)
y_grid = np.arange(-0.5, 5, 0.05)
X, Y = np.meshgrid(x_grid, y_grid)
kernel_sigma = 0.05

generate_pdfs_from_wtraj_gtraj_json = raw_input('Do you want to generate PDFs from wtraj/gtraj files? (y/n)  ')
if generate_pdfs_from_wtraj_gtraj_json == 'y':
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

    #Append data from additional file
    # data_list = []
    # with open (filename2, 'rb') as g:
    #     done = False
    #     while not done:
    #         try:
    #             data2 = pickle.load(g)
    #         except EOFError as err:
    #             done = True
    #             break
    #         data_list.append(data2)

    # Extract data
    for i, data in enumerate(data_list):
        diffusion_coeff = data['diffusion_coeff']
        if diffusion_coeff <0:
            continue
        w_traj = data['w_traj']
        g_traj = data['g_traj']
        y = g_traj[:,0] #length 200
        x = w_traj[0] # length 200
        density = data['density']
        f, ax= plt.subplots(figsize=(10,10))
        print('{}/{}, diffusion coeff = {}'.format(i+1, len(data_list), diffusion_coeff))
        x_indices = np.arange(0,len(x))
        y_indices = np.arange(0,len(y)) #array length 200

        all_gaussians_summed = np.zeros((110,110))
        for x_idx in x_indices:
            for y_idx in y_indices:
                x_point = x[x_idx]
                y_point = y[y_idx]
                all_gaussians_summed += normal_kernel_2D(X,Y, (x_point, y_point), (density[y_idx][x_idx]), kernel_sigma,kernel_sigma)  # shape (110, 110)
        #all_gaussians_summed = density
        area_under_all_gaussians = float(np.sum(all_gaussians_summed))
        # print (area_under_all_gaussians)
        pdf_normalized = all_gaussians_summed/area_under_all_gaussians # shape (110, 110)
        print ('sum of pdf: ' + str(np.sum(pdf_normalized))) # should sum to 1.0

        ax.contour(X, Y, pdf_normalized, contour_number_for_plot, colors ='gray' , zorder = 10)
        plt.pcolormesh(w_traj, g_traj, density, cmap='binary')
        #ax.contour(w_traj, g_traj, pdf_normalized, contour_number_for_plot, colors ='gray' , zorder = 1)

        for idx, windspeed in enumerate(windspeeds):
            #x_index = find_nearest_index(x, windspeed)
            #y_index = find_nearest_index(y, groundspeeds[idx])
            x_index = find_nearest_index(x_grid, windspeed)
            y_index = find_nearest_index(y_grid, groundspeeds[idx])
            ax.scatter(windspeed, groundspeeds[idx], color = 'black',zorder = 19)
            probability_string = '%.0E' %(pdf_normalized[y_index,x_index])
            ax.text(windspeed+0.03, groundspeeds[idx]+0.03, probability_string,zorder = 20)

        ax.set_ylim(-0.5, 5)
        ax.set_xlim(-2, 3.5)
        adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
        ax.set_ylabel('groundspeed along trajectory, m/s')
        ax.set_xlabel('windspeed along trajectory, m/s')
        # plt.show()
        pdf_array = pdf_normalized.tolist()
        pdf_dict = {'diff_coeff_'+str(diffusion_coeff) : pdf_array}
        if radius == 250:
            figname = '250_meters_diff_coeff_'+str(diffusion_coeff)+'_FIELD_WINDS.svg'
            pdf_name = '250_meters_diff_coeff_'+str(diffusion_coeff)+'_pdf.json'
        else:
            pdf_name = 'diff_coeff_'+str(diffusion_coeff)+'_pdf.json'
            figname = 'diff_coeff_'+str(diffusion_coeff)+'_FIELD_WINDS.svg'
        plt.savefig(figname)
        ###### saving data ######
        with open(pdf_name, mode = 'w') as f:
            json.dump(pdf_dict,f, indent = 1)
        #######################

else:
    print ('OK, going straight to continuous PDF plotting, and log Bayes calculation.')

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
    # Extract data
    for i, data in enumerate(data_list):
        diffusion_coeff = data['diffusion_coeff']
        if diffusion_coeff < 0:
            continue
        w_traj = data['w_traj']
        g_traj = data['g_traj']
        y = g_traj[:,0] #length 200
        x = w_traj[0] # length 200


    groundspeeds = field_data_dictionary['groundspeeds']
    windspeeds = field_data_dictionary['windspeeds']

    if examine_outlier_omission:
        groundspeeds_outliers_omitted = field_data_dictionary_outliers_omitted['groundspeeds']
        windspeeds_outliers_omitted = field_data_dictionary_outliers_omitted['windspeeds']



    #now for Bayes factor calculation
    diffusion_coeff_list = [10.0, 20.0, 30.0, 40.0 , 50.0 , 60.0 , 70.0, 80.0, 90.0,200.0, 300.0 , 400.0, 500.0, 600.0]#, 700.0, 800.0, 900.0, 1000.0]
    if radius == 250:
        alternate_model_list = ['250_meters_diff_coeff_' +str(diff)+'_pdf.json' for diff in diffusion_coeff_list]
        reference_model_name = '250_meters_diff_coeff_100.0_pdf.json'
    if radius == 1000:
        alternate_model_list = ['diff_coeff_' +str(diff)+'_pdf.json' for diff in diffusion_coeff_list]
        reference_model_name = 'diff_coeff_100.0_pdf.json'
    with open(reference_model_name) as f:
        ref_dict =  json.load(f)
    for key in ref_dict:
        ref_pdf = np.asarray(ref_dict[key])

    list_of_avg_log_Bayes_factors =[]
    for number, model in enumerate(alternate_model_list):

        print ('alt. model: ' + model)
        with open(model) as f:
            model_dict =  json.load(f)
        for key in model_dict:
            alternate_pdf = np.asarray(model_dict[key])
        plot_model_kde_with_conf_intervals_and_field_data(pdf = alternate_pdf, list_of_conf_int=[0.99], figname =model.split('.')[0]+'.svg')
        resample_number = 10000
        list_of_log_Bayes_factors = np.zeros(resample_number)
        list_of_log_Bayes_factors_outliers_omitted = np.zeros(resample_number)
        for i in range(resample_number):
            resampled_indices = np.random.choice(len(groundspeeds),len(groundspeeds), replace = True)
            list_of_log_likelihood_ratios =[]
            for idx in resampled_indices:
                x_index = find_nearest_index(x_grid, windspeeds[idx])
                y_index = find_nearest_index(y_grid, groundspeeds[idx])
                list_of_log_likelihood_ratios.append(np.log(ref_pdf[y_index, x_index]/alternate_pdf[y_index, x_index]))
            log_Bayes_factor = np.sum(list_of_log_likelihood_ratios)
            list_of_log_Bayes_factors[i] = log_Bayes_factor

            if examine_outlier_omission:
                resampled_indices = np.random.choice(len(groundspeeds_outliers_omitted),len(groundspeeds_outliers_omitted), replace = True)
                list_of_log_likelihood_ratios_outliers_omitted =[]
                for idx in resampled_indices:
                    x_index = find_nearest_index(x_grid, windspeeds_outliers_omitted[idx])
                    y_index = find_nearest_index(y_grid, groundspeeds_outliers_omitted[idx])
                    list_of_log_likelihood_ratios_outliers_omitted.append(np.log(ref_pdf[y_index, x_index]/alternate_pdf[y_index, x_index]))
                log_Bayes_factor_outliers_omitted = np.sum(list_of_log_likelihood_ratios_outliers_omitted)
                list_of_log_Bayes_factors_outliers_omitted[i] = log_Bayes_factor_outliers_omitted

        # ax = plt.subplot(4,4,number+1)
        # hist_factor = -2.
        # binwidth = 30
        #
        # binlist = range(int(np.nanmin(list_of_log_Bayes_factors)), int(np.nanmax(list_of_log_Bayes_factors))+binwidth, binwidth)
        # ax.hist(list_of_log_Bayes_factors, bins = binlist, color = 'k', weights=(10**hist_factor)*np.ones_like(list_of_log_Bayes_factors), histtype = 'step')

        # if examine_outlier_omission:
        #     binlist = range(int(np.nanmin(list_of_log_Bayes_factors_outliers_omitted)), int(np.nanmax(list_of_log_Bayes_factors_outliers_omitted))+binwidth, binwidth)
        #     ax.hist(list_of_log_Bayes_factors_outliers_omitted, bins = binlist, color = 'gray', weights=(10**hist_factor)*np.ones_like(list_of_log_Bayes_factors_outliers_omitted), histtype = 'step')

        # ax.set_ylim(0, 40)
        # ax.set_xlim(-750, 1300)
        # plt.xticks(np.linspace(-600, 1200, 4, endpoint = True))
        # adjust_spines(ax_handle=ax, spines= ['bottom'])
        # ax.text(800,30, 'ln (L1/L%d)'%(number+2))
        print (np.mean(list_of_log_Bayes_factors))
        list_of_avg_log_Bayes_factors.append(np.mean(list_of_log_Bayes_factors))
        if examine_outlier_omission:
            print (np.mean(list_of_log_Bayes_factors_outliers_omitted))
        print ('')
    if radius == 1000:
        figname = '1000_meter_mean_bayes_factors_wrt_D.svg'
    if radius == 250:
        figname = '250_meter_mean_bayes_factors_wrt_D.svg'
    fig = plt.figure(figsize= (0.9,0.9))
    ax = plt.subplot(111)
    ax.plot (diffusion_coeff_list, list_of_avg_log_Bayes_factors, color = 'k')
    ax.scatter(diffusion_coeff_list+[100.0], list_of_avg_log_Bayes_factors+[0], color = 'k', s = 10)
    #ax.scatter([100.0], [0], color = 'k') #because a diffusivity of 100 was used as the reference model
    ax.set_ylabel ('mean log Bayes factor')
    ax.set_xlabel ('diffusion coefficient')
    adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
    plt.savefig(figname)
    # plt.show()
