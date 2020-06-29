import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import circmean
import math
import datetime

path_to_original_data = '/home/kate/Documents/coyote_lake_field_data/2017_10_26/weather_data/metone_anemometer_data/'
direction_file_name = 'wind_direction_20171026_101348.csv'
speed_file_name = 'windspeed_20161128_211313.csv' # this file was date-stamped by a raspberry pi that had been unaware of the actual date.
#wind speed recording started at 10:11:00 (real time) flies released at 10:16:00 (real time)
compressed_speed_file_name = '2017_10_26_wind_speed_compressed.csv'

def get_epoch_timestamp_from_time_since_release(seconds_since_release):
    # release_time_seconds = int(release_time.split(':')[0])*3600 +int(release_time.split(':')[1])*60 + int(release_time.split(':')[2])
    release_time_epoch_seconds_GMT = (datetime.datetime(2017,10,26,10,16) - datetime.datetime(1970,1,1)).total_seconds() #this is the unix timestamp of the release time, in GMT
    release_time_epoch_seconds = release_time_epoch_seconds_GMT + 7*3600
    return release_time_epoch_seconds +seconds_since_release

def get_time_since_release_from_epoch_timestamp (self, epoch_timestamp):
    release_time_seconds = int(self.release_time.split(':')[0])*3600 +int(self.release_time.split(':')[1])*60 + int(self.release_time.split(':')[2])
    epoch_datetime_string = str(datetime.datetime.fromtimestamp(epoch_timestamp)) # comes out formatted like 2019-03-14 21:47:24.071131

    t = epoch_datetime_string.split(' ')[1]
    epoch_seconds = int(t.split(':')[0])*3600 +   int(t.split(':')[1])*60 + float(t.split(':')[2])

    time_elapsed = epoch_seconds - release_time_seconds
    return time_elapsed

def convert_direction_data_to_radians(maximum, winddirection_data, east_as_zero = True, north_as_zero = False):
    if east_as_zero:
        shift = 3*np.pi/2
    if north_as_zero:
        shift = np.pi
    winddirection_radians = [(shift - (2*np.pi*x/maximum))%(2*np.pi) for x in winddirection_data]
    return winddirection_radians

def convert_direction_data_to_degrees(maximum, winddirection_data, east_as_zero = True, north_as_zero = False):
    if east_as_zero:
        shift = 270
    if north_as_zero:
        shift = 180
    #winddirection_degrees = [(shift - (360*x/maximum))%(360) for x in winddirection_data]
    winddirection_degrees = [((360*x/maximum))%(360)-shift for x in winddirection_data]

    return winddirection_degrees

def find_nearest_index(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def pulse_rate_to_miles_per_hour(pulse_rate):
    mph = (pulse_rate*0.08946) + 0.6
    return mph

def calc_seconds_from_timestamp(timestamp):
    day = int(timestamp.split(' ')[0].split('-')[2])
    time = timestamp.split(' ')[1]
    split_time = time.split(':')
    hour = int(split_time[0])
    minute = int(split_time[1])
    second = float(split_time[2].split(',')[0])
    total_seconds = 86400*day + 3600*hour + 60*minute + second
    return total_seconds

def calc_binned_windspeeds(pulse_array, min_per_bin, bin_num): #assumes pulse_array is in seconds, and seconds have been zeroed relative to fly release
    pulse_seconds_from_release = np.array(pulse_array)
    binned_windspeeds = []
    for i in range (bin_num):
        start_index = find_nearest_index(pulse_seconds_from_release, i*min_per_bin*60)
        end_index = find_nearest_index(pulse_seconds_from_release, (i+1)*min_per_bin*60)
        p = pulse_seconds_from_release[start_index:end_index]
        p_diff = np.diff(p)
        windspeed = pulse_rate_to_miles_per_hour(1/np.mean(p_diff))
        binned_windspeeds.append(windspeed)
    return binned_windspeeds

def calc_binned_wind_dir (winddirection, bin_length_minutes, bins_to_analyze, release_index, direction_seconds_array):
#     release_index = 19581 # 1964 second lag between starting wind direction recording and fly release; effective samp rate 9.97 Hz
#     bin_length_minutes = 1
#     bins_to_analyze =20
    effective_samp_rate =(len(direction_seconds_array)/(max(direction_seconds_array)-min(direction_seconds_array)))
    print effective_samp_rate

    binned_wind_dir = []
    for i in range(bins_to_analyze):
        start = release_index+int(i*bin_length_minutes*60*effective_samp_rate)
        end = start+int(bin_length_minutes*60*effective_samp_rate)
        binned_wind_dir.append(circmean(winddirection[start:end], high = 360))
    return binned_wind_dir

############33--------------------------------------###############################333333---------------

winddirection = []
timestamps = []
with open(path_to_original_data+direction_file_name, 'r') as wind_direction_file:
    for line in wind_direction_file:
        direction = float(line.split(',')[1])
        winddirection.append(direction)
        timestamps.append(line.split(',')[0])
start_time = calc_seconds_from_timestamp(timestamps[0])
direction_seconds_array = [calc_seconds_from_timestamp(x)- start_time for x in timestamps]
winddirection_north_zero = convert_direction_data_to_degrees(max(winddirection),winddirection,east_as_zero = False, north_as_zero = True)



seconds_to_skip = 300 #skipping 5 minutes of windspeed data
pulse_seconds_from_release = []
with open(path_to_original_data+compressed_speed_file_name, 'r') as compressed_speed:
    for i, line in enumerate(compressed_speed):
        if i ==0:
            start_seconds = calc_seconds_from_timestamp(line.split(',')[0])
        secondstamp = calc_seconds_from_timestamp(line.split(',')[0])
        if secondstamp - start_seconds > seconds_to_skip:
            pulse_seconds_from_release.append(secondstamp - start_seconds - seconds_to_skip)
        if i >1000000000:
            break

bin_length_minutes = 0.1
bin_num = 200

binned_windspeeds = calc_binned_windspeeds(pulse_array = pulse_seconds_from_release,
                                           min_per_bin = bin_length_minutes,
                                           bin_num =bin_num) #yields a list of windspeeds in units of miles per hour
print (binned_windspeeds) #in miles per hour!

binned_wind_dir = calc_binned_wind_dir(winddirection =winddirection_north_zero,
                                       bin_length_minutes =bin_length_minutes,
                                       bins_to_analyze = bin_num,
                                       release_index = 68616, #80 seconds between direction recording and release; samp rate of 857 Hz,
                                       direction_seconds_array = direction_seconds_array)
print (binned_wind_dir)


# OK, now that I've made these arrays (wind speed and wind direction), let's use them to write-out a file in the same format as Will Dickson's metone anemometer logger
binned_epoch_timestamps = []
for i in range(bin_num):
    binned_epoch_timestamps.append(get_epoch_timestamp_from_time_since_release(i*bin_length_minutes*60))
print binned_epoch_timestamps


with open(path_to_original_data+'m1_wind_data_rewritten_in_wbd_format.txt', mode='w') as outfile:
    # for i in range(int((plot_how_many_minutes_post_release+5)*60*sampleHz)):
    #     list_to_write=[]
    #     list_to_write.append(str(unix_timestamp_release_time+ i/sampleHz))
    #     list_to_write.append(str(a_degrees_FROM_cw_from_north))
    #     list_to_write.append(str(avg_speed_mph))
    #     string_to_write = ' '.join(list_to_write)
    #     outfile.write(string_to_write + '\n')

    # 1526325596.004 103.2728 2.644859
    for i in range(bin_num):
        line_to_write = str(binned_epoch_timestamps[i])+' '+str(binned_wind_dir[i])+' ' +str(binned_windspeeds[i]/2.23694)
        outfile.write(line_to_write + '\n')
