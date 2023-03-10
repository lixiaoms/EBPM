import matplotlib
import socket
from ipywidgets import IntProgress
from IPython.display import display
import numpy as np


def configure_plt():
    matplotlib.rc('font', **{'family': 'sans-serif',
                             'sans-serif': ['Computer Modern Roman']})
    params = {'axes.labelsize': 10,
              'font.size': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'figure.figsize': (3, 3)}
    matplotlib.rcParams.update(params)
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'


def waitbar(i, N, h=None):
    if h is None:
        h = IntProgress(min=0, max=N)
        display(h)
    h.value = i
    return h


def deterministic(mode, state=None):
    if mode is 'on':
        state = np.random.get_state()
        np.random.seed(4251)
        return state
    else:
        if state is None:
            raise NameError('Cannot exit deterministic mode without state')
        else:
            np.random.set_state(state)
            return None
