#COURTESY OF NICK OSBORNE

import numpy as np

# map upwards via FC table and normalise:
def get_upper(ground=None, T=5000, FC_table_name =  'FCFs_2004_Fantz.npy',DE = np.array([0, 4303, 8436, 12402, 16204, 19841, 23314, 26622, 29765, 32785, 35548, 38184, 40642, 42917, 44998, 46876, 48533, 49953, 51122, 52001, 52560])
):
    import os

    FC_table_name = os.path.join(os.path.dirname(__file__), FC_table_name)
    # load Frank-Condon factor table for D2
    FC_table = np.load(FC_table_name, allow_pickle=True)
    FC_sums = np.sum(FC_table, axis=1)

    #get ground state distribution
    if type(ground)==type(None):
        ground = vibr_boltzmann(T,DE=DE)

    n_ground = len(ground)

    if n_ground>20:
        n_ground=20

    upper = np.zeros((n_ground))

    # map upwards via FC table and normalise:
    for i in range(0, n_ground):
        upper[i] = sum([ground[j] * FC_table[i][j]/FC_sums[i] for j in range(0, n_ground)])

    return upper/upper[0]

def vibr_boltzmann(T,DE=np.array([0, 4303, 8436, 12402, 16204, 19841, 23314, 26622, 29765, 32785, 35548, 38184, 40642, 42917, 44998, 46876, 48533, 49953, 51122, 52001, 52560])):
    return np.exp(-DE/T)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider

    # Insert arbitrary ground state here:
    arb = False
    arb_ground = np.array(
        [1, 0.089, 0.046, 0.034, 0.026, 0.024, 0.02, 0.017, 0.014, 0.012, 0.01, 0.008, 0.007, 0.006, 0, 0, 0, 0, 0, 0])

    # update plot
    def update(val):
        vib_upper.set_ydata(get_upper(T=vib_T_slider.val)[0:6])
        vib_ground.set_ydata(vibr_boltzmann(vib_T_slider.val)[0:6])
        fig.canvas.draw_idle()

    # initialise plot
    fig, ax22 = plt.subplots()
    plt.subplots_adjust(left=0.2)
    ax22.set_xlabel('V')
    ax22.set_xlabel('Normalised pop.')
    ax22.set_title('Vibrational Distribution', y=1.0, pad=-14)
    ax22.set_ylim(0.0, 2)
    if arb:
        vib_ground, = ax22.plot(arb_ground[0:6], '--d')
        vib_upper, = ax22.plot(np.arange(0, 6), get_upper(ground=arb_ground)[0:6], '--d')
        plt.legend(['ground', 'upper'])
    else:
        vib_ground, = ax22.plot(vibr_boltzmann(5000)[0:6], '--d')
        vib_upper, = ax22.plot(np.arange(0, 6), get_upper(T=5000)[0:6], '--d')
        plt.legend(['ground', 'upper'])
        vib_T_axis = plt.axes([0.05, 0.1, 0.02, 0.7])
        vib_T_slider = Slider(vib_T_axis, 'V_Temp/K', valmin=500, valmax=20000, valinit=5000, valstep=100, orientation='vertical')
        vib_T_slider.on_changed(update)
    plt.show()