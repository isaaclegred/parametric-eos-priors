import numpy as np
from matplotlib import pyplot as plt
import argparse
import seaborn as sns
sns.set(font="Palatino",
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
sns.set_context("notebook", rc={"font.size":16,
                                "axes.titlesize":20,
                                "axes.labelsize":18})

parser = argparse.ArgumentParser(description='Get needed info')                    
parser.add_argument("--infile-prefix", type=str, dest="file_name", default="macro-draw-000000")                                            
parser.add_argument("--plot-index", type=int, dest="plot_index", default=2)
parser.add_argument("--eos-dir", type=str, dest="eos_dir", default=".")
parser.add_argument("--outdir", type=str, dest="outdir", default=".")


# if somebody leaves in a .csv we can handle it
def remove_ending(file_name):
    if file_name[-4:] == ".csv":
        return file_name[:-4]
    return file_name
if __name__ == "__main__":
    args = parser.parse_args()
    local_file_name = remove_ending(args.file_name)
    file_name =args.eos_dir + "/" +  local_file_name + ".csv"
    plot_index = args.plot_index
    data = np.loadtxt(file_name, skiprows=1, delimiter=",")
#    data2 = np.loadtxt("../../macro-css.csv", skiprows=1, delimiter=",")
    data_names = ["central density(cgs)", "radius(km)", "mass (solar masses)"]
    other_data_name  =data_names[plot_index]
    
    radii = data[:, 1]
 #   radii2 = data2[:,2]
    other = data[:, plot_index]
  #  other2 = data2[:,1]
    #plt.plot(radii2, other2, label="phil's", color="r")
    plt.plot(radii, other, label=file_name)
    plt.ylabel(other_data_name)
    plt.xlabel( data_names[1])
    plt.ylim((.2,3))
    plt.xlim((5, 21))
    plt.legend()
    plt.show()
    plt.savefig(args.outdir + "/"+ args.file_name + ".png")
