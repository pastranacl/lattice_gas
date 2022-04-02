import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib import cm
from matplotlib import colors




def colormap_afm(name='afm_color_cm'):

    cm = np.ones((3, 4))
    
    # Color empty space
    cm[0, 0] = 1
    cm[0, 1] = 1
    cm[0, 2] = 1
    
    # Water
    cm[1, 0] = 0.215686
    cm[1, 1] = 0.443137
    cm[1, 2] = 0.784313
    
    # Surfaces
    cm[2, 0] = 0.2
    cm[2, 1] = 0.2
    cm[2, 2] = 0.2
    
    new_mp =  matplotlib.colors.ListedColormap(cm)
    plt.register_cmap(name, cmap=new_mp)


    return new_mp


def main():

    # Load data
    data = np.loadtxt("./mean_minim_lattice.dat", delimiter="\t", dtype ="int") 
    data = np.flipud(data)
    nr,nc = data.shape
    
    yvals = np.linspace(-nr/2, nr/2, nr)
    xvals = np.linspace(-nc/2, nc/2, nc)
    
    colorsafm = colormap_afm(name='afm_color_cm')
    
    # Plot and setting parameters of the heatmap 
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    plt.rcParams['axes.linewidth'] = 2.0 #set the value globally
    plt.rcParams['font.size'] = 17.0
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['FreeSans']
    fig = plt.figure()
    ax = fig.gca()
    c = ax.pcolor(xvals,yvals,data, cmap=colorsafm)


    #ax.xaxis.set_ticks(np.linspace(60, 85, 6))
    #ax.yaxis.set_ticks(np.linspace(0, 0.50, 6))
    ax.tick_params(axis=u'both', which=u'both',length=0)


    ax.set_xlabel(r'$x/l$',  fontsize=20)
    ax.set_ylabel(r'$y/l$',  fontsize=20)


    #plt.show()
    plt.tight_layout()

    plt.gcf().subplots_adjust(bottom=0.12)
    plt.gcf().subplots_adjust(left=0.15)
    plt.savefig('afm.pdf', dpi=300)  


if __name__ == "__main__":
    main() 
