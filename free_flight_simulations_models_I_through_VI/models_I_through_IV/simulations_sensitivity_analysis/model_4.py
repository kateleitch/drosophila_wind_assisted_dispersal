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

def calculate_heading_vector (wind_speed, wind_dir,
                            fly_trajectory, forward_airspeed):
        print()
        wind_vector = np.array([np.cos(wind_dir)*wind_speed, np.sin(wind_dir)*wind_speed]) #length of this vector is in units of m/s
        fly_trajectory_unit_vector = np.array([np.cos(fly_trajectory), np.sin(fly_trajectory)]) # length of this vector is not meaningful.
        #in this script, "parallel wind vector" is the wind component parallel to the TRAJECTORY (which is fixed)
        parallel_wind_vector = vector_projection(wind_vector, fly_trajectory_unit_vector) #length of this vector is in units of m/s, negative values allowed
        #in this script, "perpendicular wind vector" is the wind component perpendicular to the TRAJECTORY (which is fixed)
        perpendicular_wind_vector = wind_vector - parallel_wind_vector #length of this vector is in units of m/s
        perpendicular_airspeed_vector = -1*perpendicular_wind_vector #in a fixed trajectory model, must oppose perp component of wind.

        perp_mag = np.linalg.norm(perpendicular_airspeed_vector) #airspeed magnitude perpendicular to trajectory. always positive
        if allow_flexible_airspeeds:
            if perp_mag > forward_airspeed_limit:
                print ('perp_mag exceeds forward airspeed limit') #immediately excludes some cases, whe
                return(None, None)
            elif perp_mag > forward_airspeed:
                if np.sign(scalar_projection(wind_vector, fly_trajectory_unit_vector)) == -1: # wind opposing intended traj
                    wind_opposing = scalar_projection(wind_vector, fly_trajectory_unit_vector)
                    par_mag = np.abs(wind_opposing) #this will give us a fly with a zero groundspeed vector, because perp_mag is being set to cancel the perp component of wind.
                    forward_airspeed = np.sqrt(par_mag**2 + perp_mag**2)
                    print ('setting forward airspeed to '+str(forward_airspeed))
                    if forward_airspeed > forward_airspeed_limit:
                        print (forward_airspeed)
                        return (None,None)
                else:
                    forward_airspeed = perp_mag
                    par_mag = 0.0
            else:
                par_mag = np.sqrt(forward_airspeed**2 - perp_mag**2)
        else:
            if perp_mag > forward_airspeed:
                #print ('maintaining this trajectory is not possible in this wind')
                return(None, None)
            else:
                par_mag = np.sqrt(forward_airspeed**2 - perp_mag**2)   # KJL 06Feb fixed sign error here. #airspeed magnitude parallel to trajectory.

        parallel_airspeed_vector = par_mag*fly_trajectory_unit_vector
        total_airspeed_vector = parallel_airspeed_vector + perpendicular_airspeed_vector
        total_trajectory_vector = total_airspeed_vector + wind_vector
        if np.sign(scalar_projection(total_trajectory_vector, fly_trajectory_unit_vector)) == -1:
            print ('backwards')
            return(None,None)
        #this is just a check:
        achieved_trajectory_angle = np.arctan2(total_trajectory_vector[1], float(total_trajectory_vector[0]))
        error_ang = fly_trajectory - achieved_trajectory_angle
        if error_ang > 0:
            print ('error: ' + error_ang)

        fly_heading_unit_vector = total_airspeed_vector/np.linalg.norm(total_airspeed_vector)
        return[wind_vector, total_trajectory_vector, perpendicular_wind_vector, parallel_wind_vector, fly_heading_unit_vector]

def run(forward_airspeed):
    dictionary = {}
    for wind_speed in wind_speed_list:
        dictionary[str(wind_speed)] = {}
        for wind_dir in wind_direction_list:
            experiment_dictionary = {'trajectories':[],'wind vectors': [], 'perp_vectors':[], 'par_vectors':[], 'heading_unit_vectors':[],'groundspeeds along trajectories':[], 'windspeeds along trajectories':[], 'trajectory_vectors':[]}

            for fly_trajectory in fly_trajectories:
                return_list = calculate_heading_vector(wind_speed, wind_dir,
                                            fly_trajectory,
                                            forward_airspeed)
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
    json_name = 'forward_airspeed_'+str(forward_airspeed)

    with open('./model_4/'+json_name+'.json', mode = 'w') as f:
        json.dump(dictionary,f, indent = 1)

with open('/home/kate/Documents/coyote_lake_field_data/first_fly_arrival_vector_avg_windspeeds.json') as f:
    field_wind_dict =  json.load(f)
field_wind_list = []
for key in field_wind_dict:
    field_wind_list.append(field_wind_dict[key])
wind_speed_list = generate_windspeed_list_from_field_wind_distribution(field_wind_list=field_wind_list, output_number=261, kernel_sigma= 0.1)

wind_direction_list = np.linspace(0,2*np.pi,180, endpoint = False)#in radians, the direction TO WHICH the wind is blowing. NOTE THAT THIS IS PI RADIANS OFF MY TYPICAL CONVENTION

fly_number = 1
fly_trajectories = list(np.linspace(0, np.pi*2, fly_number, endpoint = False))

script_dir = os.path.dirname(__file__)
model_folder_name = 'model_4/'
model_dir = os.path.join(script_dir, model_folder_name)
if not os.path.isdir(model_dir):
    os.makedirs(model_dir)

allow_flexible_airspeeds = False

airmin_list = [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]
airmax_list = [0.5, 0.75, 1.0, 1.25,1.5, 1.75, 2.0, 2.25, 2.5]
gpref_list = [0.125, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]


forward_airspeed_list = [1.5, 1.75, 2.0, 2.25, 2.5]#[0.125, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
for f in forward_airspeed_list:
    run(forward_airspeed = f)
