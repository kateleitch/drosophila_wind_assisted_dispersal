from __future__ import print_function
import os
import sys
import scipy
import scipy.signal
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

filename = 'density_vs_w_and_g_merged_py2.pkl'
filename = 'density_vs_w_and_g_D30_D1000_tend.pkl'
smoothed = False

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

gs_n = 4
gs_m = 4
fig = plt.figure(1, figsize = (9,6))
gs = gridspec.GridSpec(gs_n,  gs_m)

# Extract data
for i, data in enumerate(data_list):

    w_traj = data['w_traj']
    g_traj = data['g_traj']
    density = data['density']
    diffusion_coeff = data['diffusion_coeff']
    print('{}/{}, diffusion coeff = {}'.format(i+1, len(data_list), diffusion_coeff))

    #Smooth data using smoothing kernel
    n = 5
    kernel = scipy.ones((n,n),dtype=scipy.float64)
    kernel = kernel/kernel.sum()
    density_smooth = scipy.signal.convolve2d(density,kernel,mode='same')

    ax = plt.subplot(gs[i])
    plt.axis('on')

    m = i%gs_m  + 1
    n = int(i/gs_m) + 1
    if m != 1:
        ax.set_yticklabels([])
    if n != gs_n:
        ax.set_xticklabels([])

    if not smoothed:
        plt.pcolormesh(w_traj, g_traj, density, cmap='binary')

    else:
        plt.pcolormesh(w_traj, g_traj, density_smooth, cmap='binary')

    max_w = w_traj.max()
    min_w = w_traj.min()
    text_x= min_w + 0.02*(max_w - min_w)

    max_g = g_traj.max()
    min_g = g_traj.min()
    text_y = min_g + 0.88*(max_g - min_g)


    plt.text(text_x, text_y, 'D = {:0.0f}'.format(diffusion_coeff))
    # if i == 0:
    #     plt.text(text_x, text_y-0.5, "something's fishy",color='red')

fig.text(0.5, 0.05, 'w_traj (m/s)', ha='center')
fig.text(0.08, 0.5, 'g_traj (m/s)', va='center', rotation='vertical')

plt.show()
