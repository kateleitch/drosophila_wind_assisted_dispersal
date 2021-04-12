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

def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_best_theta (wind_vector, perpendicular_airspeed_vector, preferred_groundspeed_along_body_axis):
    """ wind_speed is a scalar.
    a_perp is the perpendicular component of the airspeed vector"""
    sign = scalar_projection(perpendicular_airspeed_vector, np.array([0.0,1.0]))
    if sign == 0:
        test_theta_list = np.array([0.0])
    if sign > 0: #this means that perpendicular_airspeed_vector is above the x-axis, that is, has some northerly component.
        test_theta_list = np.linspace(0, np.pi/2., 300, endpoint = False)[1:]
    if sign < 0: #this means that perpendicular_airspeed_vector is below the x-axis, that is, has some southerly component.
        test_theta_list = np.linspace(3*np.pi/2., 2*np.pi, 300, endpoint = False)[1:]

    candidate_groundspeed_along_body_axis_scalar_list = []
    candidate_theta_list = []
    candidate_airspeed_along_body_axis_vector_list = []

    intended_trajectory_vector = np.array([1.0, 0.0]) # the convention for all my simulations is to have the fixed azimuthal behavior oriented exactly east.
    for i, test_theta in enumerate(test_theta_list):
        test_theta_unit_vector = np.array([np.cos(test_theta), np.sin(test_theta)]) # has a unitary magnitude of 1.0 meter per second along body axis (confirmed that this operation yields a unit vector, 18 February)
        projecting_1ms_onto_perp_airspeed_vector = scalar_projection(test_theta_unit_vector, perpendicular_airspeed_vector) #should range from [1 to 0]; confirmed this 18 February
        required_airspeed_along_body_axis = (np.linalg.norm(perpendicular_airspeed_vector))/projecting_1ms_onto_perp_airspeed_vector #required to achieve the airspeed perpendicular to the trajectory
        required_airspeed_vector_along_body_axis = required_airspeed_along_body_axis * test_theta_unit_vector

        if required_airspeed_along_body_axis < -0.2:
            continue
        if required_airspeed_along_body_axis > 1.8:
            continue

        wind_along_body_axis = vector_projection(wind_vector, test_theta_unit_vector)
        groundspeed_along_body_axis_vector = wind_along_body_axis + required_airspeed_vector_along_body_axis
        #groundspeed_along_body_axis = np.linalg.norm(groundspeed_along_body_axis_vector) # <------ kjl 11 Feb -- THIS WAS A HUGE PROBLEM. it forces an absolute value...
        groundspeed_along_body_axis = scalar_projection(groundspeed_along_body_axis_vector, test_theta_unit_vector) # <---- kjl 11 Feb allowing negative values
        if groundspeed_along_body_axis < 0:
            #print ('fly would get backwards optic flow along her body axis')
            continue
        ###
        total_trajectory_vector = required_airspeed_vector_along_body_axis + wind_vector
        if np.sign(scalar_projection(total_trajectory_vector, intended_trajectory_vector)) < 0: #if trajectory is antiparallel with intended trajectory
            continue
        ###
        candidate_groundspeed_along_body_axis_scalar_list.append(groundspeed_along_body_axis)
        candidate_theta_list.append(test_theta)
        candidate_airspeed_along_body_axis_vector_list.append(required_airspeed_vector_along_body_axis)
    try:
        best_index = find_nearest_index(candidate_groundspeed_along_body_axis_scalar_list, preferred_groundspeed_along_body_axis) # <---- 12 Feb this is working as expected.
    except:
        return (None, None, None)

    achieved_groundspeed_along_body_axis = candidate_groundspeed_along_body_axis_scalar_list[best_index]
    if achieved_groundspeed_along_body_axis < 0:
        print ('fly flying backwards, which I think means the code was problematic somewhere above')
        return (None, None, None)
    else:
        theta = candidate_theta_list[best_index]
        achieved_airspeed_vector_along_body_axis = candidate_airspeed_along_body_axis_vector_list[best_index]
        return (theta, achieved_groundspeed_along_body_axis, achieved_airspeed_vector_along_body_axis)

def calculate_heading_vector (wind_speed, wind_dir,
                            fly_trajectory, preferred_groundspeed_along_body_axis,forward_airspeed_limit):

        wind_vector = np.array([np.cos(wind_dir)*wind_speed, np.sin(wind_dir)*wind_speed]) #length of this vector is in units of m/s
        fly_trajectory_unit_vector = np.array([np.cos(fly_trajectory), np.sin(fly_trajectory)]) # length of this vector is not meaningful.
        #in this script, "parallel wind vector" is the wind component parallel to the TRAJECTORY (which is fixed)
        parallel_wind_vector = vector_projection(wind_vector, fly_trajectory_unit_vector) #length of this vector is in units of m/s, negative values allowed
        #in this script, "perpendicular wind vector" is the wind component perpendicular to the TRAJECTORY (which is fixed)
        perpendicular_wind_vector = wind_vector - parallel_wind_vector #length of this vector is in units of m/s
        perpendicular_airspeed_vector = -1*perpendicular_wind_vector #in a fixed trajectory model, must oppose perp component of wind.
        absolute_value_of_a_perp = np.linalg.norm(perpendicular_airspeed_vector)
        if absolute_value_of_a_perp > forward_airspeed_limit:
            return(None, None)
        theta, achieved_groundspeed_scalar_along_body_axis, airspeed_vector_along_body_axis = find_best_theta (wind_vector, perpendicular_airspeed_vector, preferred_groundspeed_along_body_axis)
        if theta == None:
            return(None,None)
        fly_heading_angle = theta
        fly_heading_unit_vector = np.array([np.cos(fly_heading_angle), np.sin(fly_heading_angle)])

        total_trajectory_vector =airspeed_vector_along_body_axis + wind_vector
        if np.abs(scalar_projection(airspeed_vector_along_body_axis, perpendicular_airspeed_vector) - np.linalg.norm(perpendicular_airspeed_vector)) > 0.000001:
            print ('something is wrong')
            print (scalar_projection(airspeed_vector_along_body_axis, perpendicular_airspeed_vector) - np.linalg.norm(perpendicular_airspeed_vector))
            print (total_trajectory_vector)
            return (None, None)
        return[wind_vector, total_trajectory_vector, perpendicular_wind_vector, parallel_wind_vector, fly_heading_unit_vector, airspeed_vector_along_body_axis, achieved_groundspeed_scalar_along_body_axis]

def run(preferred_groundspeed_along_body_axis, forward_airspeed_limit, reverse_airspeed_limit):
    dictionary = {}
    for wind_speed in wind_speed_list:
        dictionary[str(wind_speed)] = {}
        for wind_dir in wind_direction_list:
            experiment_dictionary = {'trajectories':[],'wind vectors': [], 'perp_vectors':[], 'par_vectors':[], 'heading_unit_vectors':[],'groundspeeds along trajectories':[], 'windspeeds along trajectories':[], 'trajectory_vectors':[], 'airspeeds along body axes':[], 'groundspeeds along body axes':[]}

            for fly_trajectory in fly_trajectories:
                return_list = calculate_heading_vector(wind_speed, wind_dir,
                                            fly_trajectory,
                                            preferred_groundspeed_along_body_axis,forward_airspeed_limit)

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
                airspeed_scalar_along_body_axis = scalar_projection(return_list[5],return_list[4])
                experiment_dictionary['airspeeds along body axes'].append(airspeed_scalar_along_body_axis)
                experiment_dictionary['groundspeeds along body axes'].append(return_list[6])
            dictionary[str(wind_speed)][str(wind_dir)] = experiment_dictionary

    json_name = 'airmin_'+str(reverse_airspeed_limit)+'_airmax_'+str(forward_airspeed_limit)+'_gpref_'+str(preferred_groundspeed_along_body_axis)

    with open('./model_3/'+json_name+'.json', mode = 'w') as f:
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
model_folder_name = 'model_3/'
model_dir = os.path.join(script_dir, model_folder_name)
if not os.path.isdir(model_dir):
    os.makedirs(model_dir)

airmin_list = [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]
airmax_list = [1.75, 2.0, 2.25, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]#[0.5, 0.75, 1.0, 1.25,1.5, 1.75, 2.0, 2.25, 2.5]
gpref_list = [1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2]#[0.125, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]

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
