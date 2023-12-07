################################################################################
# Name: map_numerics_1spec_si.py
#
# Purpose: Numerics for 1 species memory-advection-diffusion model, with periodic
#          boundary conditions
#
# Usage: python map_numerics_1spec_si.py 10000 100000 200 100 0.1 5 1 temp.inp temp.out
################################################################################

import sys 

# Get parameters from command line
curr_arg = 1
max_time = int(sys.argv[curr_arg])
curr_arg += 1
time_granularity = int(sys.argv[curr_arg])
curr_arg += 1
space_granularity = int(sys.argv[curr_arg])
curr_arg += 1
out_time_res = int(sys.argv[curr_arg])
curr_arg += 1
delta = float(sys.argv[curr_arg])
curr_arg += 1
gamma = float(sys.argv[curr_arg])
curr_arg += 1
D_val = float(sys.argv[curr_arg])
curr_arg += 1
infile = open(sys.argv[curr_arg],'r')
curr_arg += 1
outfile = open(sys.argv[curr_arg],'w')

# Derived parameters
delta_t = 1/float(time_granularity)
delta_x = 1/float(space_granularity)

# Set up arrays to store distributions and set-up initial conditions
old_u_array = []
new_u_array = []
smooth_u_array = []
u_tot = 0.0

# Box goes from -1 to 1
box_width = 2.0

# Flag whether finished
finish = 0

# Get initial conditions from file
line = infile.readline()
split_line = line.rsplit()
if(len(split_line) != space_granularity):
    sys.stderr.write("Warning: input array of length %i; should be length %i\n" % (len(split_line),space_granularity))

for counter in range(space_granularity):
    old_u_array += [float(split_line[counter])]
    new_u_array += [0]
    smooth_u_array += [0]
    u_tot += old_u_array[counter]

# Normalise    
for counter in range(space_granularity):
    old_u_array[counter] *= space_granularity/(u_tot*box_width)

# Output ICs
for counter in range(space_granularity):
    if counter == space_granularity - 1:
        outfile.write('%f\n' % (old_u_array[counter]))
    else:
        outfile.write('%f\t' % (old_u_array[counter]))

# Loop through time, solving PDE using finite difference method
u_tot = 0.0
for time in range(max_time):
    # Smoothed arrays
    for space in range(space_granularity):
        minx = int(space-delta*space_granularity)
        maxx = int(space+delta*space_granularity)+1
        smooth_u_array[space] = 0
        for z_val in range(minx,maxx):
            smooth_u_array[space] += old_u_array[z_val % space_granularity]/(maxx-minx)

    for space in range(space_granularity):
        # Change u-values
        new_u_array[space] = old_u_array[space] + (delta_t/(delta_x**2))*(
                                     D_val*(old_u_array[(space+1) % space_granularity]-2*old_u_array[space]+old_u_array[(space-1) % space_granularity]) - 
                                     gamma*(((smooth_u_array[(space+2) % space_granularity] - smooth_u_array[space])*old_u_array[(space+1) % space_granularity]) - 
                                          ((smooth_u_array[space] - smooth_u_array[(space-2) % space_granularity])*old_u_array[(space-1) % space_granularity]))/4) 

    # Ensure no negative values
    for space in range(space_granularity):
        if new_u_array[space] < 0:
            sys.stderr.write("WARNING: gone below zero at position %i\n" % space)
            finish = 1
            break
    if finish == 1:
        # We have gone below zero; break
        break

    # Copy new arrays into old    
    u_tot = 0.0
    for space in range(space_granularity):
      old_u_array[space] = new_u_array[space]
      u_tot += new_u_array[space]
    if time % max(out_time_res,100) == 0:
      sys.stderr.write('time: %i\t%f\n' % (time,u_tot))
      # Put out resulting distribution
    if time % out_time_res == 0:
      for counter in range(space_granularity):
        if counter == space_granularity - 1:
          outfile.write('%f\n' % (new_u_array[counter]))
        else:
          outfile.write('%f\t' % (new_u_array[counter]))
