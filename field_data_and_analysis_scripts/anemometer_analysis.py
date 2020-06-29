#! /usr/bin/python
from __future__ import print_function
import cv2 # opencv
import numpy as np
import os, sys
import matplotlib.pyplot as plt
import itertools
import time
import json
from pylab import *
from scipy.optimize import curve_fit
import random
import seaborn as sns
import datetime
from scipy.stats import circmean

# ONLY COMPATIBLE WITH DATA PROCESSED THROUGH TEENSY -- > https://github.com/willdickson/m1_wind_sensor; ONLY FOR WIND DATA ACQUIRED 2018 AND LATER

# winddirection transformation is based on the anemometer's "lower notch" oriented SOUTH, which is my new convention in the field.
# any time "wind direction" is mentioned, this means "direction the wind is coming FROM"
# the anemometer's output wraps when the wind is coming from the NORTH

#Will's m1_wind_sensor script saves speed data in meters per second, and direction data in degrees

class AnemometerAnalyzer:
    def __init__(self, directory, ax_handle, time_shift, minutes_to_vector_average_wind, annie_sim_data = False, turn_off_text = False, turn_off_scalebar = False, orient_circmean_from_north = False, turn_off_title = False, bin_duration = 60., declare_fixed_scale = 0, turn_off_wind = False, plot_vectors_as_function_of_time = False, savefig = True):
        self.directory = directory
        if annie_sim_data == False:
            with open(directory+'/field_parameters.json') as f:
                field_parameters = json.load(f)
            self.release_time = field_parameters["time_of_fly_release"]
            self.mag_declin_deg = field_parameters["magnetic_declination_degrees_cw_from_true_N"]
            print (self.release_time)
        else:
            self.release_time = '16:00:00' # <---- in my timezone, unix epoch is 1969-12-31 16:00:00
        # with open(directory+'/trap_layout_parameters.json') as f:
        #     self.trap_layout_parameters = json.load(f)
        self.ax_handle = ax_handle
        self.time_shift = time_shift
        self.annie_sim_data = annie_sim_data
        self.turn_off_text = turn_off_text
        self.turn_off_scalebar = turn_off_scalebar
        self.orient_circmean_from_north = orient_circmean_from_north
        self.turn_off_title = turn_off_title
        self.bin_duration = bin_duration
        self.declare_fixed_scale = declare_fixed_scale
        self.turn_off_wind = turn_off_wind
        self.plot_vectors_as_fnc_of_time = plot_vectors_as_function_of_time
        self.minutes_to_vector_average_wind = minutes_to_vector_average_wind
        self.savefig = savefig
    def get_filenames(self, path, contains, does_not_contain=['~', '.pyc']):
        cmd = 'ls ' + '"' + path + '"'
        ls = os.popen(cmd).read()
        all_filelist = ls.split('\n')
        try:
            all_filelist.remove('')
        except:
            pass
        filelist = []
        for i, filename in enumerate(all_filelist):
            if contains in filename:
                fileok = True
                for nc in does_not_contain:
                    if nc in filename:
                        fileok = False
                if fileok:
                    filelist.append( os.path.join(path, filename) )
        return filelist

    def get_time_since_release_from_epoch_timestamp (self, epoch_timestamp):
        release_time_seconds = int(self.release_time.split(':')[0])*3600 +int(self.release_time.split(':')[1])*60 + int(self.release_time.split(':')[2])
        epoch_datetime_string = str(datetime.datetime.fromtimestamp(epoch_timestamp)) # comes out formatted like 2019-03-14 21:47:24.071131

        t = epoch_datetime_string.split(' ')[1]
        epoch_seconds = int(t.split(':')[0])*3600 +   int(t.split(':')[1])*60 + float(t.split(':')[2])

        time_elapsed = epoch_seconds - release_time_seconds
        return time_elapsed

    def adjust_spines(self, ax_handle, spines):
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

    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = np.argmin((np.abs(array - value)))
        return idx

    def run(self):
        if self.annie_sim_data == False:
            wind_filename_list = self.get_filenames(self.directory+'/weather_data/metone_anemometer_data/', contains = 'm1_wind_data')
        else:
            wind_filename_list = self.get_filenames(self.directory+'/simulated_wind_data/', contains = 'simulated_wind_data')

        if len(wind_filename_list) > 1:
            for index, wind_file_name in enumerate(wind_filename_list):
                response = raw_input('Found %d files in this dir; want to analyze file: %s? (y/n) ' %(len(wind_filename_list), wind_file_name))
                if response == 'y':
                    print ('OK, analyzing ' +wind_file_name)
                    break
                if response == 'n':
                    continue
        else:
            wind_file_name = wind_filename_list[0]

        winddirection_list = []
        windspeed_list = []
        timestamp_list = []
        with open(wind_file_name, 'r') as wind_direction_file:
            for line in wind_direction_file:
                winddirection_list.append(float(line.split(' ')[1])-self.mag_declin_deg) # in degrees from TRUE north
                timestamp_list.append(float(line.split(' ')[0])) # unix epoch
                windspeed_list.append(float(line.split(' ')[2])) # in meters per second

        sec_since_release_list =[]
        for timestamp in timestamp_list:
            sec_since_release_list.append(self.get_time_since_release_from_epoch_timestamp(timestamp))

        #   THE FOLLOWING IS ENTIRELY NOT LEGIT; JUST TO ALLOW PLOTTING OF ONE DATE'S WIND DATA AGAINST ANOTHER DATE'S FLY DATA
        #if self.annie_sim_data == False:
            #sec_since_release_list = sec_since_release_list - np.mean(sec_since_release_list) # just to ensure a zero-crossing in the time-since-release list, THIS SHOULD BE DELETED WHEN DOING REAL ANALYSES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111111111!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111111111!!!!!!!!!!!1111111

        plot_how_many_minutes_pre_release = 0   # < --------- NEEDS TO BE MADE FLEXIBLE
        plot_how_many_minutes_post_release = self.minutes_to_vector_average_wind


        indices_to_plot=[]
        for index, sec_since_release in enumerate(sec_since_release_list):

            if sec_since_release - self.time_shift< -1*plot_how_many_minutes_pre_release*60.:
                continue
            if sec_since_release - self.time_shift > plot_how_many_minutes_post_release*60.:
                print ('breaking')
                break
            indices_to_plot.append(int(index))
            # print ('min since release: ' +str(sec_since_release/60.))
        bin_duration = self.bin_duration #please always specify in seconds
        binned_wind_dir = []
        binned_windspeeds = []
        print ('length of winddirection list:' +str(len(winddirection_list)))
        print (indices_to_plot[0])
        print (indices_to_plot[-1])
        direction_sublist = winddirection_list[indices_to_plot[0]:indices_to_plot[-1]] # in degrees
        speed_sublist = windspeed_list[indices_to_plot[0]:indices_to_plot[-1]] # in m/s
        sec_since_release_sublist = sec_since_release_list[indices_to_plot[0]:indices_to_plot[-1]] # in seconds

        start_time = sec_since_release_sublist[0]
        bin_count = 0
        #total_bin_number = np.ceil(len(speed_sublist)/(20*bin_duration))
        print ('length of speed sublist: '+str(np.ceil(len(speed_sublist))))
        print ('length of direction sublist: ' +str(len(direction_sublist)))
        total_seconds = (plot_how_many_minutes_pre_release+plot_how_many_minutes_post_release)*bin_duration
        print ('total seconds of data to analyze: ' +str(total_seconds))
        total_bin_number = (plot_how_many_minutes_pre_release+plot_how_many_minutes_post_release)*60.0/bin_duration
        #total_bin_number = np.ceil(len(speed_sublist)/((plot_how_many_minutes_pre_release+plot_how_many_minutes_post_release)*bin_duration)) #2019_06_11 this does not make sense to me
        print ('total bin number: ' +str(total_bin_number))
        degrees_to_rotate = 0
        if self.orient_circmean_from_north:
            print ('ORIENTING CIRCMEAN FROM NORTH: ' + str(degrees_to_rotate))
            degrees_to_rotate = circmean(direction_sublist, high = 360, low =0)
        while True:
            start_index = self.find_nearest(sec_since_release_sublist, start_time)
            end_time = start_time+bin_duration
            end_index = self.find_nearest(sec_since_release_sublist, end_time)
            if bin_count < total_bin_number:
                direction_slice = direction_sublist[start_index:end_index]
                speed_slice = speed_sublist[start_index:end_index] #still in m/s
            else:
                print ('this last bin is not a full %d seconds in length ' %(bin_duration))
                binned_wind_dir.append(circmean(direction_sublist[start_index:-1], high = 360, low =0) - degrees_to_rotate)
                binned_windspeeds.append(np.mean(speed_sublist[start_index:-1]))
                break
            binned_wind_dir.append(circmean(direction_slice, high = 360, low =0)- degrees_to_rotate)
            binned_windspeeds.append(np.mean(speed_slice)) # still in m/s
            start_time = end_time
            bin_count += 1
        ax = self.ax_handle
        number_of_pre_release_bins = plot_how_many_minutes_pre_release * (60/bin_duration)
        if self.plot_vectors_as_fnc_of_time == False:
            wind_vector_points = [[0,0]]  #in cartesian coordinates, start point of head-to-tail vectors
            if self.annie_sim_data == False:
                for i in range (len(binned_wind_dir)):
                    x = wind_vector_points[i][0]- (np.cos(np.pi/2 - binned_wind_dir[i]*np.pi/180)* binned_windspeeds[i])
                    y = wind_vector_points[i][1]- (np.sin(np.pi/2 - binned_wind_dir[i]*np.pi/180)* binned_windspeeds[i])
                    wind_vector_points.append([x,y])
                # print (wind_vector_points)
            else:
                for i in range (len(binned_wind_dir)-1):
                    x = wind_vector_points[i][0]- (np.cos(np.pi/2 - binned_wind_dir[i]*np.pi/180)* binned_windspeeds[i])
                    y = wind_vector_points[i][1]- (np.sin(np.pi/2 - binned_wind_dir[i]*np.pi/180)* binned_windspeeds[i])
                    wind_vector_points.append([x,y])
                # print (wind_vector_points)
            for i in range(len(wind_vector_points)-1):
                color = 'black'
                if i < number_of_pre_release_bins:
                    color = 'gray'
                ax.scatter(wind_vector_points[i][0], wind_vector_points[i][1], s = 0) #dummy plot; not sure why I need this, but empirically I do
                if not self.turn_off_wind:
                    ax.annotate('', xy=(wind_vector_points[i+1][0],wind_vector_points[i+1][1]),
                                xytext=(wind_vector_points[i][0],wind_vector_points[i][1]),
                                arrowprops=dict(arrowstyle="simple, head_width=1", linewidth = 2, color = color))
            wind_vector_points = np.array(wind_vector_points)
            x_range = max(wind_vector_points[:,0]) - min(wind_vector_points[:,0])
            y_range = max(wind_vector_points[:,1]) - min(wind_vector_points[:,1])
            if self.declare_fixed_scale>0:
                x_range = self.declare_fixed_scale
                y_range = self.declare_fixed_scale
            margin = 0
            set_range = max(x_range, y_range) + margin
            x_center = (max(wind_vector_points[:,0]) + min(wind_vector_points[:,0]))/2.
            y_center = (max(wind_vector_points[:,1]) + min(wind_vector_points[:,1]))/2.
            ax.set_ylim(y_center -set_range/2, y_center +set_range/2)
            ax.set_xlim(x_center- set_range/2 ,x_center +set_range/2)
            plt.axis('off')

        else:
            time_bin_list = np.linspace(-1*plot_how_many_minutes_pre_release*60, plot_how_many_minutes_post_release*60, total_bin_number+1)
            print (len(time_bin_list), len(binned_windspeeds))
            #ymax = np.nanmax(binned_windspeeds)*1.1
            ymax = 5.0
            ymin = -0.5
            xmax = np.max(time_bin_list)*1.05/60.
            #xmin = -160
            xmin= -3

            column_num =2.9 #<--- I wanted this to be a principled thing, but because of the matplotlib axis layout I needed to empirically
            #choose a value here to allow the plotted arrows to have the same length regardless of their orientation
            y_scale = 1
            x_scale = (xmax-xmin)/(ymax-ymin) *(1/column_num)
            ax.set_ylim([ymin, ymax])
            ax.set_xlim([xmin, xmax])
            ax.plot([t/60. for t in time_bin_list], binned_windspeeds, 'gray')
            for index, time in enumerate(time_bin_list):
                speed = binned_windspeeds[index] #speed, in meters per second
                dir = binned_wind_dir[index]*np.pi/180. #direction, here converted into radians
                time = time/60.
                ax.scatter(time,speed, s=0)
                ax.annotate('', xy=(time -(np.cos(np.pi/2 - dir)*x_scale), speed-np.sin(np.pi/2 - dir)*y_scale),
                            xytext=(time, speed),
                            arrowprops=dict(arrowstyle="simple, head_width=1", color = (0.5,0.5,0.5,1), alpha = 1))

        if not self.turn_off_title:
            #ax.set_title('wind dir/speed for %d min pre, %d min post release' %(plot_how_many_minutes_pre_release, plot_how_many_minutes_post_release))
            ax.set_title('avg %01.1f m/s over %d min post-release' %(np.mean(binned_windspeeds[0:-1]), plot_how_many_minutes_post_release), fontsize = 10)
        if not self.turn_off_scalebar:
            scalebar_position = [x_center - set_range/2., y_center+set_range/6.]
            scalebar_mph = 4
            scalebar_mps = scalebar_mph/2.23694
            ax.annotate(#str(scalebar_mph)+' m.p.h.',
                        '',
                        xy =scalebar_position,
                        xytext=(scalebar_position[0]+scalebar_mps, scalebar_position[1]),
                        arrowprops = dict(arrowstyle= '|-|', linewidth =3, color = 'black'))
            ax.text(scalebar_position[0]+scalebar_mph/2., scalebar_position[1] +y_range/12., str(scalebar_mph)+' mph', ha = 'center')

        if not self.turn_off_text:
            #ax.text(scalebar_position[0], scalebar_position[1] - 6, ("%.1f" % avg_speed)+' mph average speed')
            for vector_index in [-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13]:
                if np.isnan(wind_vector_points[vector_index][0]):
                    continue
                else:
                    mean_dx =  wind_vector_points[0][0] - wind_vector_points[vector_index][0]
                    print ('mean dx: '+str(mean_dx))
                    mean_dy = wind_vector_points[0][1] - wind_vector_points[vector_index][1]
                    print ('mean dy: '+str(mean_dy))
                    break
            mean_dir = np.arctan2(mean_dy, mean_dx)
            mean_dir = (2*np.pi + mean_dir)
            mean_dir = divmod(mean_dir, 2*np.pi)[1]

            vector_sum_speed = np.sqrt(mean_dy**2 + mean_dx**2)
            number_of_vectors_summed = len(wind_vector_points[0:vector_index])
            avg_speed = vector_sum_speed/number_of_vectors_summed  # last index can contain NaN
            avg_speed_mph = avg_speed * 2.23694
            print ('average speed, meters per second: ' +str(avg_speed))
            print ('mean direction FROM which wind blew calculated as: '+str(mean_dir))
            print (self.directory)
            radians_at_16_compass_points = np.linspace(0,2*np.pi,16, endpoint=False)
            idx = self.find_nearest(radians_at_16_compass_points,mean_dir)
            compass_points = ['E','ENE', 'NE','NNE',
                              'N', 'NNW', 'NW', 'WNW',
                              'W', 'WSW', 'SW', 'SSW',
                              'S', 'SSE', 'SE', 'ESE']
            compass_point = compass_points[idx]
            ax.text(scalebar_position[0]+scalebar_mph/2., scalebar_position[1] - y_range/4,
                    'avg wind ' +('%.1f'%avg_speed_mph)+ ' mph \n' + 'out of ' + compass_point, ha = 'center')
            # ax.axis('equal')

        savefig = self.savefig
        if savefig:
            plt.savefig(self.directory+'/weather_data/metone_anemometer_data/' + 'wind_direction_and_speed.png', transparent = False)
            plt.savefig(self.directory+'/weather_data/metone_anemometer_data/' + 'wind_direction_and_speed.svg')
        return degrees_to_rotate, mean_dir, avg_speed
        #plt.show() # <---- COMMENT THIS OUT WHEN INTEGRATING THIS WITH OVERARCHING SUMMARY FIGURE
