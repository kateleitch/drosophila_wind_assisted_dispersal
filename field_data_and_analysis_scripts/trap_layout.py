#!/usr/bin/python

from __future__ import print_function
import cv2
import numpy as np
import os, sys
import matplotlib.pyplot as plt
import json

####
#directory = sys.argv[1] #specify experiment date directory with respect to current working directory
#planned_or_actual = sys.argv[2]

class TrapLayoutVisualizer:
    def __init__(self, directory, planned_or_actual, ax_handle, turn_off_text, return_histogram = False):  # <---- instances of this class will specify the directory, most likely using directory = sys.argv[1]
        self.directory = directory
        self.planned_or_actual = planned_or_actual
        self.ax_handle = ax_handle
        self.turn_off_text = turn_off_text
        self.return_histogram = return_histogram

    def run(self, dot_colors_and_sizes, trap_catch_dictionary, radians_to_rotate = 0, annie_sim_data = False):
        annotate_with_text = True
        if self.turn_off_text:
            annotate_with_text = False
        experiment_date = self.directory.split('/')[-1]
        output_name = self.directory+'/trap_layout_'+ self.planned_or_actual+'.png'
        print (self.directory+'/trap_layout_parameters.json')
        with open(self.directory+'/trap_layout_parameters.json') as f:
            plotting_parameters = json.load(f)

        figure_size = plotting_parameters['figsize']
        inner_trap_radius = plotting_parameters['inner_trap_radius']
        outer_trap_radius = plotting_parameters['outer_trap_radius']
        inner_trap_number = plotting_parameters['inner_trap_number']
        outer_trap_number = plotting_parameters['outer_trap_number']
        are_traps_ordered_clockwise = plotting_parameters['are_traps_ordered_clockwise']
        inner_ring_first_trap_rad_cw_from_n = plotting_parameters['inner_ring_first_trap_deg_cw_from_n']*np.pi/180
        outer_ring_first_trap_rad_cw_from_n = plotting_parameters['outer_ring_first_trap_deg_cw_from_n']*np.pi/180

        # fig = plt.figure(figsize=(figure_size,figure_size))
        # ax = fig.add_subplot(111, projection='polar')
        if annie_sim_data == False:
            self.ax_handle.set_theta_offset(np.pi/2)
            if are_traps_ordered_clockwise:
                self.ax_handle.set_theta_direction(-1)

        trap_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        if annie_sim_data:
            trap_names = '0123456789'
        trap_name_list = [trap_names[x] for x in range(0, len(trap_names))] #this list most likely ends up longer than the total number of traps, which is fine

        theta = np.linspace(0, 2*np.pi, 360)
        r1 = np.ones(360) *outer_trap_radius
        r2 = np.ones(360) *inner_trap_radius
        if inner_trap_number != 0:
            self.ax_handle.plot(theta, r2, 'k', zorder=1, linewidth = 0.5)
        if outer_trap_number != 0:
            self.ax_handle.plot(theta, r1, 'k', zorder=1, linewidth = 0.5)

        inner_trap_angles = np.linspace(0,2*np.pi, inner_trap_number, endpoint = False) + inner_ring_first_trap_rad_cw_from_n - radians_to_rotate
        outer_trap_angles = np.linspace(0,2*np.pi, outer_trap_number, endpoint = False) + outer_ring_first_trap_rad_cw_from_n - radians_to_rotate

        self.ax_handle.scatter(inner_trap_angles,np.ones(inner_trap_number)*inner_trap_radius, s = 1, color = 'black', zorder = 3, edgecolor = 'white')
        self.ax_handle.scatter(outer_trap_angles,np.ones(outer_trap_number)*outer_trap_radius, s = 1, color = 'black', zorder = 3, edgecolor = 'white')

        if self.return_histogram:
            histo_inner_angle_list = []
            histo_outer_angle_list = []

        for i in range(inner_trap_number+outer_trap_number):
            letter = trap_name_list[i]
            if annie_sim_data:
                catches = trap_catch_dictionary[letter]
            else:
                catches = trap_catch_dictionary['trap_'+letter]
            if i < inner_trap_number:
                self.ax_handle.scatter(inner_trap_angles[i],inner_trap_radius, s = dot_colors_and_sizes[letter]['dotsize'], color = dot_colors_and_sizes[letter]['dotcolor'], zorder =2, edgecolor = 'white')
                if annotate_with_text:
                    self.ax_handle.text(inner_trap_angles[i], inner_trap_radius*1.3, letter, va = 'center', ha = 'center')
                    self.ax_handle.text(inner_trap_angles[i], inner_trap_radius*0.7, catches, va = 'center', ha = 'center')
                if self.return_histogram:
                    histo_inner_angle_list.extend([inner_trap_angles[i]]*catches)
            else:
                self.ax_handle.scatter(outer_trap_angles[i-inner_trap_number],outer_trap_radius, s = dot_colors_and_sizes[letter]['dotsize'], color = dot_colors_and_sizes[letter]['dotcolor'], zorder =2, edgecolor = 'white')
                if annotate_with_text:
                    self.ax_handle.text(outer_trap_angles[i-inner_trap_number], outer_trap_radius*1.3, letter, va = 'center', ha = 'center')
                    self.ax_handle.text(outer_trap_angles[i-inner_trap_number], outer_trap_radius*0.7, str(catches), va = 'center', ha = 'center')
                if self.return_histogram:
                    histo_outer_angle_list.extend([outer_trap_angles[i-inner_trap_number]]*catches)

        self.ax_handle.scatter(0,0, color = 'black')
        print ('outer trap radius: '+str(outer_trap_radius))
        self.ax_handle.set_ylim(0,outer_trap_radius*1.6)

        if annotate_with_text:
            self.ax_handle.text(0, outer_trap_radius*1.3, experiment_date+' ' +self.planned_or_actual +'\n %d traps at %d m, and %d traps at %d m from release site'%(inner_trap_number, inner_trap_radius, outer_trap_number, outer_trap_radius),va = 'center', ha = 'center')

        self.ax_handle.grid(False)
        self.ax_handle.axis("off")
        if self.return_histogram:
            return (histo_inner_angle_list, histo_outer_angle_list)

    def plot_separate_scatterpoints_for_scale(self, dot_colors_and_sizes, trap_catch_dictionary, radians_to_rotate = 0, annie_sim_data = False, turn_off_circles = True):
        annotate_with_text = True
        if self.turn_off_text:
            annotate_with_text = False
        experiment_date = self.directory.split('/')[-1]
        output_name = self.directory+'/trap_layout_'+ self.planned_or_actual+'.png'
        print (self.directory+'/trap_layout_parameters.json')
        with open(self.directory+'/trap_layout_parameters.json') as f:
            plotting_parameters = json.load(f)

        figure_size = plotting_parameters['figsize']
        inner_trap_radius = 250#plotting_parameters['inner_trap_radius']
        outer_trap_radius = 1000 #plotting_parameters['outer_trap_radius']
        inner_trap_number  = 1
        outer_trap_number  = 1

        are_traps_ordered_clockwise = plotting_parameters['are_traps_ordered_clockwise']
        inner_ring_first_trap_rad_cw_from_n = 180*np.pi/180
        outer_ring_first_trap_rad_cw_from_n = 0*np.pi/180

        if annie_sim_data == False:
            self.ax_handle.set_theta_offset(np.pi/2)
            if are_traps_ordered_clockwise:
                self.ax_handle.set_theta_direction(-1)

        trap_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        if annie_sim_data:
            trap_names = '0123456789'
        trap_name_list = [trap_names[x] for x in range(0, len(trap_names))] #this list most likely ends up longer than the total number of traps, which is fine

        inner_trap_angles = np.linspace(0,2*np.pi, inner_trap_number, endpoint = False) + inner_ring_first_trap_rad_cw_from_n - radians_to_rotate
        outer_trap_angles = np.linspace(0,2*np.pi, outer_trap_number, endpoint = False) + outer_ring_first_trap_rad_cw_from_n - radians_to_rotate

        self.ax_handle.scatter(inner_trap_angles,np.ones(inner_trap_number)*inner_trap_radius, s = 1, color = 'black', zorder = 3, edgecolor = 'white')
        self.ax_handle.scatter(outer_trap_angles,np.ones(outer_trap_number)*outer_trap_radius, s = 1, color = 'black', zorder = 3, edgecolor = 'white')

        if self.return_histogram:
            histo_inner_angle_list = []
            histo_outer_angle_list = []

        for i in range(inner_trap_number+outer_trap_number):
            letter = trap_name_list[i]
            if annie_sim_data:
                catches = trap_catch_dictionary[letter]
            else:
                catches = trap_catch_dictionary['trap_'+letter]
            if i < inner_trap_number:
                self.ax_handle.scatter(inner_trap_angles[i],inner_trap_radius, s = dot_colors_and_sizes[letter]['dotsize'], color = dot_colors_and_sizes[letter]['dotcolor'], zorder =2, edgecolor = 'white')
                if annotate_with_text:
                    #self.ax_handle.text(inner_trap_angles[i], inner_trap_radius*1.3, letter, va = 'center', ha = 'center')
                    self.ax_handle.text(inner_trap_angles[i] -0.25*np.pi, inner_trap_radius, str(catches)+' flies', va = 'center', ha = 'center')
                if self.return_histogram:
                    histo_inner_angle_list.extend([inner_trap_angles[i]]*catches)
            else:
                self.ax_handle.scatter(outer_trap_angles[i-inner_trap_number],outer_trap_radius, s = dot_colors_and_sizes[letter]['dotsize'], color = dot_colors_and_sizes[letter]['dotcolor'], zorder =2, edgecolor = 'white')
                if annotate_with_text:
                    #self.ax_handle.text(outer_trap_angles[i-inner_trap_number], outer_trap_radius*1.3, letter, va = 'center', ha = 'center')
                    self.ax_handle.text(outer_trap_angles[i-inner_trap_number]-0.18*np.pi, outer_trap_radius, str(catches)+' flies', va = 'center', ha = 'center')
                if self.return_histogram:
                    histo_outer_angle_list.extend([outer_trap_angles[i-inner_trap_number]]*catches)

        #self.ax_handle.scatter(0,0, color = 'black')
        print ('outer trap radius: '+str(outer_trap_radius))

        self.ax_handle.set_ylim(0,outer_trap_radius*1.6) #was 1.35

        # self.ax_handle.tight_layout()
        # if annotate_with_text:
        #     self.ax_handle.text(0, outer_trap_radius*1.3, experiment_date+' ' +self.planned_or_actual +'\n %d traps at %d m, and %d traps at %d m from release site'%(inner_trap_number, inner_trap_radius, outer_trap_number, outer_trap_radius),va = 'center', ha = 'center')

        self.ax_handle.grid(False)
        self.ax_handle.axis("off")

        if self.return_histogram:
            return (histo_inner_angle_list, histo_outer_angle_list)

    def run_planned(self ):
        annotate_with_text = True
        if self.turn_off_text:
            annotate_with_text = False
        experiment_date = self.directory.split('/')[-1]
        output_name = self.directory+'/trap_layout_'+ self.planned_or_actual+'.png'
        print (self.directory+'/trap_layout_parameters.json')
        with open(self.directory+'/trap_layout_parameters.json') as f:
            plotting_parameters = json.load(f)

        try:
            with open(self.directory+'/trap_coords.json') as f:
                lat_lon_coords = json.load(f)["release_site_1"]
        except:
            print ('fyi, looks like you have not yet generated a kml file for this experiment')

        figure_size = plotting_parameters['figsize']
        inner_trap_radius = plotting_parameters['inner_trap_radius']
        outer_trap_radius = plotting_parameters['outer_trap_radius']
        inner_trap_number = plotting_parameters['inner_trap_number']
        outer_trap_number = plotting_parameters['outer_trap_number']
        traps_w_mini_vanes = plotting_parameters["traps_with_mini_wind_vanes"]
        are_traps_ordered_clockwise = plotting_parameters['are_traps_ordered_clockwise']
        inner_ring_first_trap_rad_cw_from_n = plotting_parameters['inner_ring_first_trap_deg_cw_from_n']*np.pi/180
        outer_ring_first_trap_rad_cw_from_n = plotting_parameters['outer_ring_first_trap_deg_cw_from_n']*np.pi/180

        # fig = plt.figure(figsize=(figure_size,figure_size))
        # ax = fig.add_subplot(111, projection='polar')
        self.ax_handle.set_theta_offset(np.pi/2)
        if are_traps_ordered_clockwise:
            self.ax_handle.set_theta_direction(-1)

        trap_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        trap_name_list = [trap_names[x] for x in range(0, len(trap_names))] #this list most likely ends up longer than the total number of traps, which is fine
        lat_lon_coords
        theta = np.linspace(0, 2*np.pi, 360)
        r1 = np.ones(360) *outer_trap_radius
        r2 = np.ones(360) *inner_trap_radius
        if inner_trap_number != 0:
            self.ax_handle.plot(theta, r2, 'k', zorder=1)
        if outer_trap_number != 0:
            self.ax_handle.plot(theta, r1, 'k', zorder=1)

        inner_trap_angles = np.linspace(0,2*np.pi, inner_trap_number, endpoint = False) + inner_ring_first_trap_rad_cw_from_n
        outer_trap_angles = np.linspace(0,2*np.pi, outer_trap_number, endpoint = False) + outer_ring_first_trap_rad_cw_from_n

        self.ax_handle.scatter(inner_trap_angles,np.ones(inner_trap_number)*inner_trap_radius, s = 1, color = 'black', zorder = 3, edgecolor = 'white')
        self.ax_handle.scatter(outer_trap_angles,np.ones(outer_trap_number)*outer_trap_radius, s = 1, color = 'black', zorder = 3, edgecolor = 'white')
        for i in range(inner_trap_number+outer_trap_number):
            letter = trap_name_list[i]
            coords = lat_lon_coords[i]
            dotcolor = (0.7,0.7,0.7,1)
            if letter in traps_w_mini_vanes:
                dotcolor = 'orange'
            if i < inner_trap_number:
                self.ax_handle.scatter(inner_trap_angles[i],inner_trap_radius, s = 100, color = dotcolor, zorder =2, edgecolor = 'white')
                self.ax_handle.text(inner_trap_angles[i], inner_trap_radius*1.1, letter+'\n'+str(coords.split(',')[0]+'\n'+coords.split(',')[1]), va = 'center', ha = 'right')
            else:
                self.ax_handle.scatter(outer_trap_angles[i],outer_trap_radius, s = 100, color = dotcolor, zorder =2, edgecolor = 'white')
                self.ax_handle.text(outer_trap_angles[i-inner_trap_number], outer_trap_radius*1.1, letter+'\n'+str(coords.split(',')[0]+'\n'+coords.split(',')[1]), va = 'center', ha = 'right')
        self.ax_handle.scatter(0,0, color = 'black')
        self.ax_handle.text (np.pi, outer_trap_radius*0.1,  'release site'+'\n'+str(coords.split(',')[0]+'\n'+coords.split(',')[1]), va = 'center', ha = 'center')
        print ('outer trap radius: '+str(outer_trap_radius))
        self.ax_handle.set_ylim(0,outer_trap_radius*1.6) #was *1.3

        if annotate_with_text:
            self.ax_handle.text(0, outer_trap_radius*1.4, experiment_date+' ' +self.planned_or_actual +'\n %d traps at %d m, and %d traps at %d m from release site'%(inner_trap_number, inner_trap_radius, outer_trap_number, outer_trap_radius),va = 'center', ha = 'center') #was *1.4

            self.ax_handle.text(0, outer_trap_radius*1.3, 'orange traps outfitted with mini-vanes',va = 'center', ha = 'center') #was *1.3

        # if radians_to_rotate !=0:
        #     # if we're rotating the whole thing by some amount, I think this necessitates the addition of a compass arrow so the viewer can easily tell what's going on
        #     self.ax_handle.annotate('North', xy=(wind_vector_points[i+1][0],wind_vector_points[i+1][1]),
        #                 xytext=(wind_vector_points[i][0],wind_vector_points[i][1]),
        #                 arrowprops=dict(arrowstyle="simple, head_width=1", linewidth = 2, color = 'black'))

        self.ax_handle.grid(False)
        self.ax_handle.axis("off")
        plt.savefig(output_name, transparent = True)
        plt.show()
