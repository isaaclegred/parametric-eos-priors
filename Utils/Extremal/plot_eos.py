import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import argparse
import seaborn as sns
sns.set(font="Arial",
        rc={
    'axes.axisbelow': False,
            'axes.edgecolor': 'lightgrey',
            'axes.facecolor': 'None',
            'axes.grid': False,
            'axes.labelcolor': 'dimgrey',
            'axes.spines.right': False,
            'axes.spines.top': False,
            'figure.facecolor': 'white',
            'lines.solid_capstyle': 'round',
            'patch.edgecolor': 'w',
            'patch.force_edgecolor': True,
            'text.color': 'dimgrey',
            'xtick.bottom': False,
            'xtick.color': 'dimgrey',
            'xtick.direction': 'out',
            'xtick.top': False,
            'ytick.color': 'dimgrey',
            'ytick.direction': 'out',
            'ytick.left': False,
            'ytick.right': False})




parser = argparse.ArgumentParser(description='Get needed info')                    
parser.add_argument("--infile-prefix", type=str, dest="file_name", default="eos-draw-000000")                                            
parser.add_argument("--plot-index", type=int, dest="plot_index", default=2)                                          


if __name__ == "__main__":
    args = parser.parse_args()
    file_name = args.file_name + ".csv"
    plot_index = args.plot_index
    data = np.loadtxt(file_name, skiprows=1, delimiter=",")
    #data2 = np.loadtxt("../../css.csv", skiprows=1, delimiter=",")
    data_names = ["pressure", "energy density", "baryon_density"]
    other_data_name  =data_names[plot_index]
    ps = data[:, 0]
    #ps2= data2[:,0]
    other = data[:, plot_index]
    #other2 = data2[:, plot_index]
    ax = plt.figure()
    ax = plt.loglog(other, ps, label=args.file_name)
    #plt.loglog(other2, ps2, color="r", label="phil's")
    plt.ylabel(r"$p/c^2$ (cgs)")
    plt.xlabel(other_data_name + "(cgs)" )
    plt.ylim((10**10, 10**16))
    plt.xlim((10**13, 4*10.0**15))
    plt.legend()
    plt.show()
    plt.savefig(args.file_name + ".png", bbox_inches="tight")
