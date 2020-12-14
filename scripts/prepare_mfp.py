import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pressure", default=5, help="gas pressure in Pa [def = 5]")
    parser.add_argument("-o", "--out_file", default="", help="out file [def = '']")
    parser.add_argument("-i", "--in_file",  help="input file [energy eV, P [1/(cm*Tor)]]")
    parser.add_argument("-S", "--show_fig", action="store_true", default=False, help="show plot of the input file [def = False]")
    args = parser.parse_args()
    
    data = np.loadtxt(args.in_file)
    P = float(args.pressure)
    data = np.column_stack((data[:,0], 0.01/(data[:,1]*P*0.0075)))  #0.0075 to transform Pa into Torr
    #now data has [energy J, mfp, m]
    filename = "mfp_P_%.2f_Pa.txt" % P
    if len(args.out_file)>0:
        filename = args.out_file
    np.savetxt(filename, data, fmt="%.10e", delimiter="\t", header="energy, eV\tmfp, m")
    
    if args.show_fig:
        plt.figure()
        plt.plot(data[:,0], data[:,1])
        plt.grid()
        plt.show()