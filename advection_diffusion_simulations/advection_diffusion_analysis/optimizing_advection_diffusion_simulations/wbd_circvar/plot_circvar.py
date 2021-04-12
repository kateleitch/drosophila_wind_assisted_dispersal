import sys
import scipy
import pickle
import matplotlib
import matplotlib.pyplot as plt

#******* some figure defaults *******
font = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 8}
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


def circvar_from_density(angle, density):
    """
    Calculates the circular variance from probability density specified in terms
    of an array of angles and an array of probabiltiy densities.

    Arguments:
      angle    = array of angles (radians)
      density  = array of probability densities

    Return
      cirvar   = circular variance
    """
    dangle = angle[1] - angle[0]
    cirvar =  1.0 - scipy.absolute((density*scipy.exp(1.0*1j*angle)).sum()*dangle)
    return cirvar

# -------------------------------------------------------------------------------------
if __name__ == '__main__':

    data_file_dict = {
            'cum_density_D300_R1000_velo_list.pkl' : {'D': 300.0, 'R': 1000, 'color': 'k'},
            'cum_density_D70_R250_velo_list.pkl'   : {'D': 70.0,  'R': 250,  'color': 'gray'},
            }

    save_figs = False
    save_circvar = False

    plt.figure(1,figsize=(3.5,3.5))
    # plt.figure(2,figsize=(10,7))
    ax = plt.subplot(111)

    for in_data_file, param in data_file_dict.items():
        with open(in_data_file,'rb') as f:
            data_list = pickle.load(f)
        plotcolor = param['color']
        circvar_list = []
        velocity_list = []
        for i, data in enumerate(data_list):
            angle = data['angle']
            velocity =  data['velocity']
            cum_density = data['cum_density']
            print('{}/{}: velocity = {:0.2f}'.format(i+1, len(data_list), velocity))
            circvar = circvar_from_density(angle, cum_density)
            circvar_list.append(circvar)
            velocity_list.append(velocity)
            #ax.scatter(velocity_list, circvar_list, color = plotcolor)
            ax.plot (velocity_list, circvar_list, color = plotcolor)

        ax.set_ylim([0,1.2])
        ax.set_xlim([0,3.0])
        ax.set_xticks([0,1,2,3])
        ax.set_yticks([0,0.5,1])
        ax.set_ylabel('circular variance')
        ax.set_xlabel('windspeed, m s-1')
        ax.spines['bottom'].set_bounds(0,3)
        ax.spines['left'].set_bounds(0, 1)
        adjust_spines(ax, spines = ['bottom','left'])


    filename = 'wbd_advec_diffusion_circvar_wrt_windspeed.svg'
    plt.savefig(filename, bbox_inches = 'tight')
    plt.show()
