import numpy as np
import matplotlib.pyplot as plt

def maxwell(energy, T):
    return np.sqrt(energy)/T**1.5*np.exp(-energy/T)

def integral(data):
    integral = 0
    for i in range(len(data[:,0])-1):
        integral+=0.5*(data[i+1][0] - data[i][0])*(data[i+1][1] + data[i][1])
    return integral


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-T", "--temperature", default=10, help="temperature in eV [def = 10]")
    parser.add_argument("-E_min", "--E_min", default=1e-4, help="minimum energy in eV [def = 1e-4]")
    parser.add_argument("-E_max", "--E_max", default=100, help="maximum energy in eV [def = 100]")
    parser.add_argument("-o", "--out_file", default="", help="out file name [def = '' - wil be derived from temperature]")
    parser.add_argument("-S", "--show_fig", action="store_true", default=False, help="show plot of the distribution [def = False]")
    args = parser.parse_args()
    
    Te = float(args.temperature)        #eV
    E = np.linspace(float(args.E_min), float(args.E_max), 10000)
    data = np.column_stack((E, maxwell(E, Te)))
    filename = "Maxwel_T_%.2feV.txt" % Te
    if len(args.out_file)>0:
        filename = args.out_file
    np.savetxt(filename, data, fmt="%.10e", delimiter='\t', header="energy, eV\tprobability, a.u.")
    
    if args.show_fig:
        plt.figure()
        plt.plot(data[:,0], data[:,1])
        plt.grid()
        plt.show()
    