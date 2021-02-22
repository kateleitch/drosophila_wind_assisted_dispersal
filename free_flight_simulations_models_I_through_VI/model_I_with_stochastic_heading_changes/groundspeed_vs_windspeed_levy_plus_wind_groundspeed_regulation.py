from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import random
import json
import matplotlib.gridspec as gridspec
from collections import Counter
from scipy.stats import vonmises
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
def generate_run_length_distribution(min_length, max_length, mu = 2):
    """ mu = 2 taken from Reynolds and Frye (experimental data, D. melanogaster)"""
    lengths = np.linspace(min_length, max_length, 100)
    probs = [l**(-1*mu) for l in lengths]
    run_length_pdf = probs/np.sum(probs) #area under curve sums to 1.0
    return lengths, run_length_pdf
def draw_from_run_length_distribution(lengths, run_length_pdf, number_of_draws):
    run_length_list = np.random.choice(lengths, number_of_draws, p=run_length_pdf)
    return run_length_list
# def draw_saccade_angles (min_angle, max_angle, number_of_draws, turn_bias = False):
#     pos_angles = np.linspace(min_angle, max_angle, 20)
#     neg_angles = np.linspace(-1*min_angle, -1*max_angle, 20)
#     available_angles = np.append(pos_angles , neg_angles)
#     saccade_angle_list = np.random.choice(available_angles, number_of_draws)
#     return saccade_angle_list
def plot_simulated_trajectory (ax_handle, run_length_list, saccade_angle_list, color):
    x_start = 0
    y_start = 0
    ang = 0
    for idx, saccade_angle in enumerate(saccade_angle_list):
        ang += saccade_angle
        len = run_length_list[idx]
        ax_handle.plot([x_start,x_start+ np.cos(ang)*len], [y_start, y_start + np.sin(ang)*len], color = color)
        x_start += np.cos(ang)*len
        y_start += np.sin(ang)*len
def calculate_trajectory_vector (wind_speed, wind_dir,
                            fly_heading, preferred_groundspeed_along_body_axis,
                            forward_airspeed_limit, reverse_airspeed_limit):
        wind_vector = np.array([np.cos(wind_dir)*wind_speed, np.sin(wind_dir)*wind_speed]) #length in m/s
        fly_heading_unit_vector = np.array([np.cos(fly_heading), np.sin(fly_heading)]) # length  not meaningful.
        parallel_wind_vector = vector_projection(wind_vector, fly_heading_unit_vector) #length in m/s, negative values allowed
        perpendicular_wind_vector = wind_vector - parallel_wind_vector #length of this vector is in units of m/s

        attempted_groundspeed_vector_along_body_axis = fly_heading_unit_vector * preferred_groundspeed_along_body_axis
        attempted_airspeed_vector_along_body_axis = attempted_groundspeed_vector_along_body_axis - parallel_wind_vector
        a = scalar_projection(attempted_airspeed_vector_along_body_axis,fly_heading_unit_vector)
        if reverse_airspeed_limit <= a <= forward_airspeed_limit:
            groundspeed_vector_along_body_axis = attempted_groundspeed_vector_along_body_axis
        elif a > forward_airspeed_limit: # fly cannot thrust enough to achieve preferred groundspeed along body axis
            actual_airspeed_vector_along_body_axis = forward_airspeed_limit*fly_heading_unit_vector
            groundspeed_vector_along_body_axis = actual_airspeed_vector_along_body_axis + parallel_wind_vector
            if scalar_projection(groundspeed_vector_along_body_axis, fly_heading_unit_vector) <0: #fly would be pushed backwards along her body axis
                return (None, None)
        elif a < reverse_airspeed_limit: # fly cannot brake enough to achieve preferred groundspeed along body axis
            # print ('fly cannot brake enough')
            actual_airspeed_vector_along_body_axis = reverse_airspeed_limit*fly_heading_unit_vector
            groundspeed_vector_along_body_axis = actual_airspeed_vector_along_body_axis + parallel_wind_vector
        trajectory_vector = groundspeed_vector_along_body_axis + perpendicular_wind_vector
        return trajectory_vector

def run_simulation(mu, kappa, wind_direction, wind_speed, default_airspeed, s_angles_zero_mean, s_pdf, min_run_length = 2.0, max_run_length = 50.0, fly_number = 20):
    wind_velocity_vector = np.array([np.cos(wind_direction)*wind_speed, np.sin(wind_direction)*wind_speed]) #m/s
    # color_list = [[102,194,165],[252,141,98],[141,160,203],[231,138,195],[166,216,84]]
    # fig = plt.figure(figsize=(4,8))
    # gs = gridspec.GridSpec(nrows=16, ncols=8)#(nrows=16, ncols=8)
    # ax = fig.add_subplot(gs[2:7, 1:8])
    # ax2 = fig.add_subplot(gs[9:17, 0:8])
    # ax3 = fig.add_subplot(gs[9:17, 0:8], projection = 'polar')
    # miniax = fig.add_subplot(gs[0:2, 6:8])
    total_flight_duration_list = []
    ending_angular_position_list = []
    lengths, run_length_pdf = generate_run_length_distribution(min_length = min_run_length, max_length = max_run_length, mu = mu)
    trap_radius = 1000.
    number_to_plot = 4
    groundspeeds_along_trajectory = []
    windspeeds_along_trajectory   = []
    for fly in range(fly_number):
        discard_fly = False
        position = np.array([0.0,0.0])
        heading_ang = np.random.choice(np.linspace(-1*np.pi,np.pi,360,endpoint = False), 1)[0] #INITIAL HEADING, HELD BY MENOTAXIS
        s_angles = [heading_ang + s for s in s_angles_zero_mean] # <---- 2021 CHANGE HERE; HEADING ANGLES BECOME CENTERED ON INIT. HEADING
        run_length_list =[] # in time, seconds
        # run_distance_list = [] # in meters
        saccade_angle_list = []
        position_list = [np.array([0.0,0.0])]
        windless_x_pos = 0
        windless_y_pos = 0
        windless_x_list = [0]
        windless_y_list = [0]
        #while np.sqrt(position[0]**2 + position[1]**2) < trap_radius: # while fly is still within 1-km radius

        continue_count = 0
        while np.linalg.norm(position) < trap_radius:
            heading_ang = np.random.choice(s_angles, 1, p = s_pdf)[0] #*(np.pi/180.)
            saccade_angle_list.append(heading_ang)
            #heading_ang += ang  # <------ HERE'S WHAT I'M CHANGING; INSTEAD OF SACCADES, HEADING FOR EACH RUN IS DRAWN FROM A CELESTIAL-CENTERED distribution
            trajectory_vector = calculate_trajectory_vector (wind_speed=wind_speed, wind_dir=wind_direction,
                                        fly_heading = heading_ang, preferred_groundspeed_along_body_axis = 1.0,
                                        forward_airspeed_limit=1.8, reverse_airspeed_limit=-0.2) #trajectory vector in m/s
            if trajectory_vector[0] == None:
                #print ('          None trajectory vector at '+str(wind_speed)+ ' m/s')
                continue_count +=1
                if continue_count < 100:
                    continue #causes fly to re-draw a new saccade angle from the distribution.
                else:
                    discard_fly = True
                    break

            len = draw_from_run_length_distribution(lengths = lengths, run_length_pdf = run_length_pdf, number_of_draws = 1)[0] # RUN LENGTH IN SECONDS
            run_length_list.append(len)
            this_segment_in_meters = len*trajectory_vector#runlength (sec) * traj (meters/sec)
            position = position + this_segment_in_meters
            position_list.append(position)
            if fly < number_to_plot:
                windless_x_pos += np.cos(heading_ang)*len*default_airspeed
                windless_y_pos += np.sin(heading_ang)*len*default_airspeed
                windless_x_list.append(windless_x_pos )
                windless_y_list.append(windless_y_pos )
        if discard_fly == True:
            print ('                        fly discarded')
            continue
        #here, trimming final segment so it does not extend beyond the trap radius.
        last_x_velocity = (position_list[-1][0]-position_list[-2][0])/run_length_list[-1]
        last_y_velocity = (position_list[-1][1]-position_list[-2][1])/run_length_list[-1]

        overshoot = np.sqrt(position[0]**2 + position[1]**2) - trap_radius
        overall_traj_angle = np.arctan2(position[1],position[0])
        final_traj_angle = np.arctan2(trajectory_vector[1], trajectory_vector[0])
        theta = overall_traj_angle - final_traj_angle
        overshoot_along_final_heading = np.abs(overshoot/(np.cos(theta))) #this does act as the hypotenuse, always larger
        d = np.sqrt((position_list[-1][0]-position_list[-2][0])**2 +(position_list[-1][1]-position_list[-2][1])**2)
        overshoot_corr_factor = (d-overshoot_along_final_heading)/d
        run_length_list[-1] = (overshoot_corr_factor*run_length_list[-1])
        # if fly < number_to_plot:
        #     untrimmed_last_position = np.array(position_list[-1])
        #     ax2.scatter(untrimmed_last_position[0], untrimmed_last_position[1], s = 30, color = 'red',zorder =25)
        position_list[-1][0] =  position_list[-2][0] + run_length_list[-1]*last_x_velocity
        position_list[-1][1] =  position_list[-2][1] + run_length_list[-1]*last_y_velocity
        ######################################################
        total_flight_dur = np.sum(run_length_list)
        ending_ang_pos = np.arctan2(position_list[-1][1],position_list[-1][0])
        groundspeeds_along_trajectory.append(trap_radius/total_flight_dur) #in meters per second
        total_flight_duration_list.append(total_flight_dur)
        ending_angular_position_list.append(ending_ang_pos)
        w_traj = scalar_projection(wind_velocity_vector, np.array([position_list[-1][0],position_list[-1][1]]))
        windspeeds_along_trajectory.append(w_traj)# in meters per second
    #     if fly <number_to_plot: # just plot a few flies
    #         ax2.plot([p[0] for p in position_list],[p[1] for p in position_list], color = 'black', zorder = 20)
    #         ax2.scatter(position_list[-1][0], position_list[-1][1], color = [a/255. for a in color_list[fly]], s = 40, zorder = 100)
    #         ax.scatter(np.sum(run_length_list)/60., 0, color = [a/255. for a in color_list[fly]], s = 40)
    #         ax2.plot(windless_x_list, windless_y_list, color = 'gray', zorder = 10)
    # histogram = ax.hist([p/60. for p in total_flight_duration_list], bins = np.linspace(0,250,50), histtype = 'step', color = 'black')
    # ax.set_xlabel('total flight duration (min)')
    # ax.set_ylabel('count')
    # ax.set_xlim(0, 250)
    # ax.set_ylim(-2, fly_number)
    # ax.spines['bottom'].set_bounds(0,250)
    # adjust_spines(ax_handle=ax, spines= ['bottom', 'left'])
    #
    # miniax.plot(np.log(lengths), np.log(run_length_pdf), color = 'black')
    # miniax.set_xlabel('log(l)', size = 8)
    # miniax.set_ylabel('log(p(l))', size = 8)
    # miniax.set_ylim(-8, 0)
    # miniax.set_xlim(0, 7.0)
    # plt.xticks([0, 7], size = 8)
    # plt.yticks([-8, 0], size = 8)
    # miniax.spines['bottom'].set_bounds(0,7)
    # miniax.spines['left'].set_bounds(-8, 0)
    # min = int(min_run_length)
    # max = int(max_run_length)
    # miniax.text(0, 1.5,'mu = '+str(mu)+'\nmin = %d, max = %d'%(min, max)+'_muijres_angles_'+'\nwind %0.1f m s-1 from N'%(wind_speed), size = 8)
    # adjust_spines(ax_handle=miniax, spines= ['bottom', 'left'])
    #
    # ending_angular_position_list = np.array([(x +2*np.pi)%(2*np.pi) for x in ending_angular_position_list])
    # scatter_with_stacked_symbols(ax_handle=ax3, angles=ending_angular_position_list, bin_num=2000, hist_start=2100,hist_spacing=50)
    # adjust_spines(ax_handle = ax3, spines = [])
    # ax3.patch.set_alpha(0)
    # ax3.set_ylim(0,2800)
    #
    # adjust_spines(ax_handle = ax2, spines = [])
    # an = np.linspace(0, 2 * np.pi, 100)
    # ax2.plot(trap_radius * np.cos(an), trap_radius * np.sin(an), color = 'gray')
    # ax2.axis('equal')
    # ax2.set_ylim(-1800, 1800)
    # ax2.set_xlim(-1800,1800)
    #
    # figname = './mu' +str(mu)+'_kappa'+str(kappa)+'_windspeed_'+str(wind_speed)+'_min_run_'+str(min_run_length)+'_gsreg.png'
    # plt.savefig(figname)

    return np.array(total_flight_duration_list),groundspeeds_along_trajectory, windspeeds_along_trajectory, ending_angular_position_list

def scatter_with_stacked_symbols(ax_handle, angles, bin_num, hist_start,hist_spacing):
    bins = np.linspace(0, 2*np.pi, bin_num+1, endpoint = True)
    digitized = np.digitize(angles, bins)
    x = Counter(digitized)
    for bin_index, angle in enumerate(bins):
        count = x[bin_index]
        bin_spacing = 2*np.pi/bin_num
        bin_center = angle - (bin_spacing/2)
        if count >0:
            ax_handle.scatter([bin_center]*count, np.linspace(hist_start, hist_start+(hist_spacing*count),  count, endpoint = True), marker = '+', s = 10, color = 'k')

#############################################
mu_list = [1, 1.5,2]#[1, 1.5, 2.0, 2.5]
kappa_list = [0.1, 1, 10, 500]

dictionary_keyed_by_mu = {}
for mu in mu_list:

    dictionary_keyed_by_kappa ={}
    for kappa in kappa_list:
        min_run_length = 1.0
        max_run_length = 1000.0 # very generous

        h_angles = np.linspace(-np.pi,np.pi, 361,endpoint = True)
        h_pdf_not_normalized = vonmises.pdf(h_angles, kappa)
        h_pdf = h_pdf_not_normalized/np.sum(h_pdf_not_normalized)

        preferred_groundspeed_along_body_axis = 1.0 #meters per second.
        forward_airspeed_limit = 1.8 #meters per second.
        reverse_airspeed_limit = -0.2 #meters per second.
        default_airspeed = 1.0 # m/s
        wind_direction = 3*np.pi/2. #wind blowing TO the south
        with open('/home/kate/Documents/coyote_lake/coyote_lake_field_data/first_fly_arrival_vector_avg_windspeeds.json') as f:
            field_wind_dict =  json.load(f)
        field_wind_list = []
        for key in field_wind_dict:
            field_wind_list.append(field_wind_dict[key])
        wind_speed_list = generate_windspeed_list_from_field_wind_distribution(field_wind_list=field_wind_list, output_number=261, kernel_sigma= 0.1)

        fly_number = 180 # this number taken from the length of the "wind direction list" in my earlier simulations.

        new_dict_entry = {}
        for wind_speed in wind_speed_list:
            print ('mu: ' +str(mu))
            print ('     kappa: ' +str(kappa))
            print ('                     '+ str(wind_speed))
            flight_duration_list, groundspeeds_along_trajectory, windspeeds_along_trajectory, ending_angular_position_list= run_simulation(mu, kappa, wind_direction, wind_speed, default_airspeed, h_angles, h_pdf, min_run_length, max_run_length, fly_number)

            new_dict_entry[wind_speed] = {'g_traj': groundspeeds_along_trajectory, 'w_traj': windspeeds_along_trajectory, 'angular_positions': list(ending_angular_position_list)}
            plt.close()
        dictionary_keyed_by_kappa[kappa] = new_dict_entry
    dictionary_keyed_by_mu[mu] = dictionary_keyed_by_kappa

with open('./gtraj_wtraj_by_mu_and_kappa.json', mode = 'w') as f:
    json.dump(dictionary_keyed_by_mu,f, indent = 1)
