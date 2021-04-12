import sys
import scipy
import pickle
import matplotlib.pyplot as plt

data_file_dict = {
        'cum_density_D300_R1000_velo_list.pkl' : {'D': 300.0, 'R': 1000, 'color': 'b', 'fignum': 1},
        'cum_density_D70_R250_velo_list.pkl'   : {'D': 70.0,  'R': 250,  'color': 'r', 'fignum': 2},
        }
save_figs = False

for in_data_file, param in data_file_dict.items():

    fignum = param['fignum']
    plt.figure(fignum,figsize=(10,7))

    with open(in_data_file,'rb') as f:
        data_list = pickle.load(f)

    for i, data in enumerate(data_list):
        angle = data['angle']
        velocity =  data['velocity']
        cum_density = data['cum_density']
        print('{}/{}: velocity = {:0.2f}'.format(i+1, len(data_list), velocity))

        # Transform to degrees for plotting
        angle_deg = scipy.rad2deg(angle)
        cum_density_deg  = cum_density/scipy.rad2deg(1)

        plt.plot(scipy.rad2deg(angle_deg), cum_density_deg,'b')

    plt.grid(True)
    plt.xlabel('angle (deg)')
    plt.ylabel('(prob/deg)')
    plt.title('exit probability density, D={:0.0f}, R={:0.0f}'.format(param['D'], param['R']))
    if save_figs:
        plt.savefig('exit_prob_vs_angle.png', bbox_inches='tight')

plt.show()

        
