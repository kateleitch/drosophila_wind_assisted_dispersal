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
def normal_kernel_1D (x_array, mean, sigma_x):
    return np.array(normal(x_array, mean, sigma_x))
def generate_windspeed_list_from_field_wind_distribution(field_wind_list, output_number, kernel_sigma):
    xmin = 0
    xmax = 2.8
    x = np.linspace(xmin, xmax, 2000, endpoint = True)
    all_gaussians_summed = sum(normal_kernel_1D(x, speed, kernel_sigma) for speed in field_wind_list)
    area_under_all_gaussians = float(np.sum(all_gaussians_summed))
    pdf_normalized = all_gaussians_summed/area_under_all_gaussians
    windspeed_list = np.random.choice(x, output_number, p=pdf_normalized)
    return windspeed_list
def scalar_projection(a,b): #projects a onto b, yields magnitude in direction of b
    return np.dot(a,b) / np.linalg.norm(b)
def vector_projection(a,b): #projects a onto b, yields scaled vector in direction of b
    return b * np.dot(a,b) / np.dot(b,b)
def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)
def scatter_with_stacked_symbols(ax_handle, angles, bin_num, hist_start,hist_spacing):
    bins = np.linspace(0, 2*np.pi, bin_num+1, endpoint = True)
    digitized = np.digitize(angles, bins)
    x = Counter(digitized)
    for bin_index, angle in enumerate(bins):
        count = x[bin_index]
        bin_spacing = 2*np.pi/bin_num
        bin_center = angle - (bin_spacing/2)
        if count >0:
            ax_handle.scatter([bin_center]*count, np.linspace(hist_start, hist_start+(hist_spacing*count),  count, endpoint = True), marker = '.', s = 15, color = 'k')
def rotate_vector(vector, ccw_radians_to_rotate):
    c, s = np.cos(ccw_radians_to_rotate), np.sin(ccw_radians_to_rotate)
    R = np.array(((c, -s), (s, c)))
    vector = np.array(vector)
    vector_vertical = vector.reshape(2,1)
    rotated_vector_vertical = np.dot(R,vector_vertical)
    rotated_vector = rotated_vector_vertical.reshape(1,2)[0]
    return rotated_vector
def play_with_vector_rotation():
    x = float(raw_input('starting vector, x: '))
    y = float(raw_input('starting vector, y: '))
    rad_to_rotate = np.pi/180.*float(raw_input('degrees to rotate: '))
    starting_vector = [x,y]
    rotated_vector = rotate_vector(starting_vector, rad_to_rotate)
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot([0, starting_vector[0]],[0, starting_vector[1]], color = 'black')
    ax.plot([0, rotated_vector[0]],[0, rotated_vector[1]], color ='green')
    plt.axis('equal')
    plt.show()

def plot_model_kde_with_conf_intervals_and_current_frame_values(ax, pdf, current_wtraj, current_gtraj):
    ax.set_ylim(-0.5, 4.6)
    ax.set_xlim(-2, 3.1)
    plt.xticks([0])
    ax.spines['bottom'].set_bounds(-2, 3)
    ax.spines['left'].set_bounds(0, 4)
    plt.yticks([1])
    ax.set_xlabel('wtraj')
    ax.set_ylabel('gtraj')
    ax.vlines(x= 0, ymin = 0, ymax = 4, color = 'gray', zorder = 10, linewidth = 0.25)
    ax.hlines(y = 1,xmin=-1.5, xmax=3, color = 'gray', zorder = 10, linewidth = 0.25)
    ax.scatter(current_wtraj, current_gtraj, s = 5, color = 'orange', zorder = 20)
    ax.imshow(pdf, cmap='binary', interpolation='nearest', origin='lower', extent = (-2, 3.5, -0.5, 5), zorder = 9)

#############################################333333



list_of_model_data = ['../model_1.json', '../model_2.json', '../model_3.json','../model_4_flies_simply_hold.json']

# list_of_model_data = ['../model_1_infinite_airmin.json']
list_of_model_data = ['../model_3_intended_trajectory_noted.json']
factor = 0.25
for model_num, model in enumerate(list_of_model_data):
    print (model)
    with open(model) as f:
        d =  json.load(f)

    with open(model.split('.json')[0]+'_pdf.json') as g:
        pdf_dictionary = json.load(g)
    for key in pdf_dictionary:
        pdf = np.asarray(pdf_dictionary[key])

    for windspeed in d:
        trajectory_list = []
        wind_list =[]
        traj_vector_list = []
        heading_unit_vector_list = []
        w = []
        g = []
        for wind_dir in d[windspeed]:
            try:
                traj = d[windspeed][wind_dir]['trajectories'][0]
                wind = d[windspeed][wind_dir]['wind vectors'][0]
                traj_vector = d[windspeed][wind_dir]['trajectory_vectors'][0]
                wind_vector = d[windspeed][wind_dir]["wind vectors"][0]
                heading_unit_vector = d[windspeed][wind_dir]["heading_unit_vectors"][0]
                wtraj =  d[windspeed][wind_dir]["windspeeds along trajectories"][0]
                gtraj = d[windspeed][wind_dir]["groundspeeds along trajectories"][0]
                intended_traj_angle = d[windspeed][wind_dir]["intended trajectory"][0]
            except:
                wind_list.append(None)
                trajectory_list.append(None)
                traj_vector_list.append(None)
                heading_unit_vector_list.append(None)
                continue
            wind_angle = np.arctan2(wind[1],wind[0])

            angle = (traj - wind_angle +2*np.pi)%(2*np.pi)  #WILL PLOT WITH WIND FIXED EAST
            xlabel = 'initial fly azimuthal rule'

            # angle = (traj +2*np.pi)%(2*np.pi)              #WILL PLOT WITH FLY AZIMUTHAL RULE FIXED EAST
            # xlabel = 'wind blowing to'

            trajectory_vector_rotated_with_respect_to_wind_blowing_east = rotate_vector(traj_vector, -1*wind_angle)
            heading_unit_vector_rotated_with_respect_to_wind_blowing_east = rotate_vector(heading_unit_vector, -1*wind_angle)
            wind_list.append(wind_angle)
            trajectory_list.append(angle)
            traj_vector_list.append(trajectory_vector_rotated_with_respect_to_wind_blowing_east)
            heading_unit_vector_list.append(heading_unit_vector_rotated_with_respect_to_wind_blowing_east)
            w.append(wtraj)
            g.append(gtraj)
        fig = plt.figure(figsize=(14, 10))
        gs = gridspec.GridSpec(nrows=20, ncols=28)
        ax2 = fig.add_subplot(gs[0:20, 8:28], projection= 'polar')
        ax_pdf = fig.add_subplot(gs[12:20, 0:8])
        ax_lakebed = fig.add_subplot(gs[0:20, 8:28])
        ax_notes = fig.add_subplot(gs[0:10, 0:8])

        #ax2.scatter(0,0, color = 'red')
        scatter_with_stacked_symbols(ax_handle=ax2, angles=trajectory_list, bin_num=135, hist_start=0.85,hist_spacing=0.012)
        ax2.set_ylim([0,1])
        #adjust_spines(ax_handle=ax2, spines=[])
        #plt.gca().set_position([0, 0, 1, 1])
        ax2.axis("off")
        ax2.margins(0)


        adjust_spines(ax_handle=ax_pdf, spines=['bottom', 'left'])
        plot_model_kde_with_conf_intervals_and_current_frame_values(ax= ax_pdf, pdf = pdf, current_wtraj = w, current_gtraj = g)

        ax_lakebed.patch.set_visible(False)
        ax_lakebed.axis('equal')
        #ax_lakebed.set_xlim(-8.5,8.5)
        ax_lakebed.set_xlim(-4.25,4.25)
        ax_lakebed.set_ylim(-4.25,4.25)
        #plt.gca().set_position([0, 0, 1, 1])
        adjust_spines(ax_handle=ax_lakebed, spines=[])
        ax_lakebed.scatter(0,0, color = 'black')

        for i, traj_vector in enumerate(traj_vector_list):
            if traj_vector == None:
                continue
            heading_unit_vector = heading_unit_vector_list[i]
            ax_lakebed.scatter(traj_vector[0], traj_vector[1], marker = '.', s = 10,color ='black')
            ax_lakebed.plot([traj_vector[0],traj_vector[0]+factor*heading_unit_vector[0]],[traj_vector[1],traj_vector[1]+factor*heading_unit_vector[1]], linewidth = 0.5, color = 'orange')

            intended_traj_x = factor* np.cos(intended_trajectory_angle)
            intended_traj_y = factor* np.sin(intended_trajectory_angle)
            ax_lakebed.plot([traj_vector[0],traj_vector[0]+intended_traj_x],[traj_vector[1],traj_vector[1]+intended_traj_y], linewidth = 0.5, color = 'black')
            if i%7 == 0:
                ax_lakebed.plot([0,traj_vector[0]],[0,traj_vector[1]], linewidth = 0.5)

        # ax_lakebed.text( -4.2, 4.0, 'wind, %0.2f m s-1 from west'%(float(windspeed)), fontsize = 10)
        # ax_lakebed.text( -4.2, 4.2, model,  fontsize = 10)
        adjust_spines(ax_handle=ax_notes, spines=[])
        ax_notes.patch.set_visible(False)
        ax_notes.set_xlim(-2,2)
        ax_notes.set_ylim(-2,2)
        ax_notes.text( -2, 1.5, 'wind, %0.2f m s-1 from west'%(float(windspeed)), fontsize = 12)
        ax_notes.text( -2, 1.9, model,  fontsize = 12)

        script_dir = os.path.dirname(__file__)
        model_folder_name = model.split('/')[-1].split('.json')[0]+'/'
        results_dir = os.path.join(script_dir, model_folder_name)
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        filename = '/wind_%7d'%(int(float(windspeed)*10000))+'.png'
        plt.savefig(results_dir + filename, bbox_inches = 'tight')


#    ffmpeg -framerate 25 -start_number 000183076 -i ./wind_%*.jpg -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
#ffmpeg -framerate 25 -i ./wind_%*.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
