import matplotlib.pyplot as plt

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def SetFonts():
    '''
    Set the font of the plots
    '''
    plt.rcParams['font.sans-serif'] = ['Noto Sans', 'sans-serif']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Noto Sans Math'

