
import numpy as np
import matplotlib.pyplot as plt
import sys
# parameters file
input_file=sys.argv[1]

# read values
def read_params(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

        i = 0
        for line in lines:
            if "Number of Steps:" in line:
                # steps = number of time steps
                steps = int(lines[i +1])

            if "Time Evolution:" in line:
                file1 = lines[i + 1]
                # file1_split = name of the file that contains the observables time evolution
                file1_split = file1.split('"')[1]
                file2 = lines[i + 3]
                # file2_split = name of the file that contains the MSD evolution.
                file2_split = file2.split('"')[1]
                file3 = lines[i + 7]
                # file3_split = name of the file that contains the g(r) histogram
                file3_split = file3.split('"')[1]
                break
            i = i + 1

    return steps,file1_split,file2_split,file3_split

# obtain the file names
steps, file1_split, file2_split, file3_split = read_params(input_file)

# this function reads the time, kinetic, potential, total (energies), temperature and pressure evolution from the input file
def read_values(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    time = []
    kin = []
    pot = []
    tot = []
    temp = []
    press = []

    for line in lines:
        if "#" in line:
            continue
        
        line_split = line.split()

        time.append(float(line_split[0]))
        kin.append(float(line_split[1]))
        pot.append(float(line_split[2]))
        tot.append(float(line_split[3]))
        temp.append(float(line_split[4]))
        press.append(float(line_split[5]))

    
    return time, kin, pot, tot, temp, press

# create the arrays that contain the observables evolution
time, kin, pot, tot, temp, press = read_values(file1_split)

# function that reads the evolution of MSD
def read_msd(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    time = []
    sqr_dist = []

    for line in lines:
        if "#" in line:
            continue
        
        line_split = line.split()

        time.append(float(line_split[0]))
        sqr_dist.append(float(line_split[1]))

    return time, sqr_dist

# create the arrays that contain the MSD evolution
time_dist, sqr_dist = read_msd(file2_split)

# this function read the values of r and g(r) to generate the gofr histogram
def read_gofr(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    position = []
    gofr = []

    for line in lines:
        if "#" in line:
            continue
        
        line_split = line.split()

        position.append(float(line_split[0]))
        gofr.append(float(line_split[1]))

    return position, gofr

# # create the arrays that contain the gofr histogram
position, gofr = read_gofr(file3_split)

# generate figures
plt.xlim(0,time[-1])
plt.xlabel("Time (ps)")
plt.ylabel("Energy (kJ/mol)")
plt.plot(time,kin, label="Kinetic Energy")
plt.plot(time,pot, label="Potential Energy")
plt.plot(time,tot, label="Total Energy")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=3)
plt.savefig(f"{file1_split}_Energy.png", bbox_inches="tight", dpi=300)
plt.clf()

plt.xlim(0,time[-1])
plt.plot(time,press, label="Pressure")
plt.ylabel("Pressure (Pa)")
plt.xlabel("Time (ps)")
plt.savefig(f"{file1_split}_Pressure.png", bbox_inches="tight", dpi=300)
plt.clf()

plt.xlim(0,time[-1])
plt.plot(time,temp, label="Temperature")
plt.ylabel("Temperature (K)")
plt.xlabel("Time (ps)")
plt.savefig(f"{file1_split}_Temperature.png", bbox_inches="tight", dpi=300)
plt.clf()

plt.xlim(0,time_dist[-1])
plt.plot(time_dist,sqr_dist, label="Mean Square Distance")
plt.ylabel("Mean Square Distance (${\AA}^{2}$)")
plt.xlabel("Time (ps)")
plt.savefig(f"{file2_split}_MSD.png", bbox_inches="tight", dpi=300)
plt.clf()

plt.xlim(0,position[-1])
plt.plot(position,gofr, label="Radial Distribution Function")
plt.xlabel("r (${\AA}$)")
plt.ylabel("g(r)")
plt.savefig(f"{file3_split}_Gofr.png", bbox_inches="tight", dpi=300)
plt.clf()

print("Time evolution and g(r) figures were generated.")