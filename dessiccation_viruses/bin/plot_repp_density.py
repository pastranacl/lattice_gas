import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib import cm
from matplotlib import colors




def main():

    # Load data
    data = np.loadtxt("./mean_minim_lattice.dat", delimiter="\t") 
    data = np.flipud(data)
    nr,nc = data.shape
    
    yvals = np.linspace(-nr/2, nr/2, nr)
    xvals = np.linspace(-nc/2, nc/2, nc)
    
    # Plot and setting parameters of the heatmap 
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    plt.rcParams['axes.linewidth'] = 2.0 #set the value globally
    plt.rcParams['font.size'] = 17.0
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['FreeSans']
    fig = plt.figure()
    ax = fig.gca()
    
    boundaries = [0.00, 0.26, 0.51, 0.76, 1.01, 2.01]
    colorvirus=['#ffffff', '#afc6e9', '#739bd9', '#3771c8',  '#333333']
    cml = matplotlib.colors.ListedColormap(colorvirus)
    norm = matplotlib.colors.BoundaryNorm(boundaries, cml.N, clip=True)
    
    c = ax.pcolor(xvals, yvals, data, cmap=cml, norm=norm)


    ax.tick_params(axis=u'both', which=u'both',length=0)


    ax.set_xlabel(r'$x/l$',  fontsize=20)
    ax.set_ylabel(r'$y/l$',  fontsize=20)


    #plt.show()
    plt.tight_layout()

    plt.gcf().subplots_adjust(bottom=0.12)
    plt.gcf().subplots_adjust(left=0.15)
    plt.savefig('virus.pdf', dpi=300)  


if __name__ == "__main__":
    main() 
