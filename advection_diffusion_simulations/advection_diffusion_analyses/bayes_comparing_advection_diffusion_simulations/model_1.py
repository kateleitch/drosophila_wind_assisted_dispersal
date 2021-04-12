from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import random
import json

def normal(x, mean, sigma):
    return (sigma*np.sqrt(2*np.pi))**-1 * np.exp(-(x-mean)**2/(2*sigma**2))
def normal_kernel_1D (x_array, mean, sigma_x):
    return np.array(normal(x_array, mean, sigma_x))
def generate_windspeed_list_from_field_wind_distribution(field_wind_list, output_number, kernel_sigma):
    xmin = 0
    xmax = 2.8
    x = np.linspace(xmin, xmax, 261, endpoint = True)
    all_gaussians_summed = sum(normal_kernel_1D(x, speed, kernel_sigma) for speed in field_wind_list)
    area_under_all_gaussians = float(np.sum(all_gaussians_summed))
    pdf_normalized = all_gaussians_summed/area_under_all_gaussians
    windspeed_list = np.random.choice(x, output_number, p=pdf_normalized)
    return windspeed_list

def scalar_projection(a,b): #projects a onto b, yields magnitude in direction of b
    return np.dot(a,b) / np.linalg.norm(b)

def vector_projection(a,b): #projects a onto b, yields scaled vector in direction of b
    return b * np.dot(a,b) / np.dot(b,b)

def calculate_trajectory_vector (wind_speed, wind_dir,
                            fly_heading, preferred_groundspeed_along_body_axis,
                            forward_airspeed_limit, reverse_airspeed_limit):
        wind_vector = np.array([np.cos(wind_dir)*wind_speed, np.sin(wind_dir)*wind_speed]) #length of this vector is in units of m/s
        fly_heading_unit_vector = np.array([np.cos(fly_heading), np.sin(fly_heading)]) # length of this vector is not meaningful.
        parallel_wind_vector = vector_projection(wind_vector, fly_heading_unit_vector) #length of this vector is in units of m/s, negative values allowed
        perpendicular_wind_vector = wind_vector - parallel_wind_vector #length of this vector is in units of m/s
        attempted_groundspeed_vector_along_body_axis = fly_heading_unit_vector * preferred_groundspeed_along_body_axis
        attempted_airspeed_vector_along_body_axis = attempted_groundspeed_vector_along_body_axis - parallel_wind_vector
        a = scalar_projection(attempted_airspeed_vector_along_body_axis,fly_heading_unit_vector)
        #print (a)
        if reverse_airspeed_limit <= a <= forward_airspeed_limit:
            groundspeed_vector_along_body_axis = attempted_groundspeed_vector_along_body_axis
        elif a > forward_airspeed_limit: # fly cannot thrust enough to achieve preferred groundspeed along body axis
            actual_airspeed_vector_along_body_axis = forward_airspeed_limit*fly_heading_unit_vector
            groundspeed_vector_along_body_axis = actual_airspeed_vector_along_body_axis + parallel_wind_vector
            if scalar_projection(groundspeed_vector_along_body_axis, fly_heading_unit_vector) <0: #fly would be pushed backwards along her body axis
                print ('maintaining this heading is not possible in this wind')
                return (None, None)
        elif a < reverse_airspeed_limit: # fly cannot brake enough to achieve preferred groundspeed along body axis
            print ('fly cannot brake enough')
            actual_airspeed_vector_along_body_axis = reverse_airspeed_limit*fly_heading_unit_vector
            groundspeed_vector_along_body_axis = actual_airspeed_vector_along_body_axis + parallel_wind_vector
        trajectory_vector = groundspeed_vector_along_body_axis + perpendicular_wind_vector
        return [wind_vector, trajectory_vector, perpendicular_wind_vector, parallel_wind_vector, fly_heading_unit_vector]

def run(preferred_groundspeed_along_body_axis, forward_airspeed_limit, reverse_airspeed_limit):
    dictionary = {}
    for wind_speed in wind_speed_list:
        dictionary[str(wind_speed)] = {}
        for wind_dir in wind_direction_list:
            experiment_dictionary = {'trajectories':[],'wind vectors': [], 'perp_vectors':[], 'par_vectors':[], 'heading_unit_vectors':[],'groundspeeds along trajectories':[], 'windspeeds along trajectories':[], 'trajectory_vectors':[]}

            for fly_heading in fly_headings:
                return_list = calculate_trajectory_vector(wind_speed, wind_dir,
                                            fly_heading, preferred_groundspeed_along_body_axis,
                                            forward_airspeed_limit, reverse_airspeed_limit)
                if return_list[0]==None:
                    continue
                wind_vector = return_list[0]
                trajectory_vector = return_list[1]

                trajectory_angle = np.arctan2(trajectory_vector[1],trajectory_vector[0]) # trajectory angle in radians
                trajectory_angle = (trajectory_angle+2*np.pi)%(2*np.pi) #recast from 0 to 2pi
                windspeed_along_trajectory = scalar_projection(wind_vector, trajectory_vector)
                groundspeed_along_trajectory = np.linalg.norm(trajectory_vector)
                experiment_dictionary['trajectories'].append(trajectory_angle)
                experiment_dictionary['groundspeeds along trajectories'].append(groundspeed_along_trajectory)
                experiment_dictionary['windspeeds along trajectories'].append(windspeed_along_trajectory)
                experiment_dictionary['wind vectors'].append(list(wind_vector))
                experiment_dictionary['perp_vectors'].append(list(return_list[2]))
                experiment_dictionary['par_vectors'].append(list(return_list[3]))
                experiment_dictionary['heading_unit_vectors'].append(list(return_list[4]))
                experiment_dictionary['trajectory_vectors'].append(list(trajectory_vector))
            dictionary[str(wind_speed)][str(wind_dir)] = experiment_dictionary

    json_name = 'airmin_'+str(reverse_airspeed_limit)+'_airmax_'+str(forward_airspeed_limit)+'_gpref_'+str(preferred_groundspeed_along_body_axis)

    with open('./model_1/'+json_name+'.json', mode = 'w') as f:
        json.dump(dictionary,f, indent = 1)

with open('/home/kate/Documents/coyote_lake_field_data/first_fly_arrival_vector_avg_windspeeds.json') as f:
    field_wind_dict =  json.load(f)
field_wind_list = []
for key in field_wind_dict:
    field_wind_list.append(field_wind_dict[key])
wind_speed_list = generate_windspeed_list_from_field_wind_distribution(field_wind_list=field_wind_list, output_number=261, kernel_sigma= 0.1)

wind_direction_list = np.linspace(0,2*np.pi,180, endpoint = False)#in radians, the direction TO WHICH the wind is blowing. NOTE THAT THIS IS PI RADIANS OFF MY TYPICAL CONVENTION

fly_number = 1
fly_headings = list(np.linspace(0, np.pi*2, fly_number, endpoint = False))

script_dir = os.path.dirname(__file__)
model_folder_name = 'model_1/'
model_dir = os.path.join(script_dir, model_folder_name)
if not os.path.isdir(model_dir):
    os.makedirs(model_dir)

airmin_list = [-0.2]#[-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]
airmax_list = [2.6, 2.7, 2.8, 2.9, 3.0]#[0.5, 0.75, 1.0, 1.25,1.5, 1.75, 2.0, 2.25, 2.5]
gpref_list = [0.5, 0.75, 1.0, 1.25, 1.5]

for forward_airspeed_limit in airmax_list:
    for reverse_airspeed_limit in airmin_list:
        for preferred_groundspeed_along_body_axis in gpref_list:
            if preferred_groundspeed_along_body_axis > forward_airspeed_limit:
                continue
            else:
                run(preferred_groundspeed_along_body_axis, forward_airspeed_limit, reverse_airspeed_limit)
