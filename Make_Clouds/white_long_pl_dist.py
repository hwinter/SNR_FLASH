#This code does the following:
# Computes a random x,y,z coordinate based on uniform deviates
# Computes a radius based on 


#import scipy
import numpy as np
#import matplotlib.pyplot as plt

###########################################################
#Input parameters
#Define the number of "clouds".
n_clouds=1000
#Index of the power law
pl_index=3
#Output file name
filename='clouds.txt'
separator=','

#Scale of the computational domain in pc.
xrange=np.array([0,200])
yrange=np.array([0,200])
zrange=np.array([0,200])

#Scale of the cloud size in pc.
min_size_cloud=0.001
max_size_cloud=1.0
#Cloud density in g/cc
rho=10
###########################################################
parsec2cm=1#3.09E+18
xrange=xrange*parsec2cm
yrange=yrange*parsec2cm
zrange=zrange*parsec2cm

min_size_cloud=min_size_cloud*parsec2cm
max_size_cloud=max_size_cloud*parsec2cm
###########################################################
x_scale_comp_domain=xrange[1]-xrange[0]
y_scale_comp_domain=yrange[1]-yrange[0]
z_scale_comp_domain=zrange[1]-zrange[0]
###########################################################
#Define a subroutine function that generates a random x,y,z coordinate and a
# random radius
def generate_pos_rad():
    x=np.random.uniform()*x_scale_comp_domain
    y=np.random.uniform()*y_scale_comp_domain
    z=0 #np.random.uniform()*z_scale_comp_domain
    size_counter=0
    while size_counter==0 :
        r=np.random.power(pl_index)*max_size_cloud
        size_counter=test_cloud_size(r)

    return x,y,z,r ;

###########################################################
#Define a subroutine function that detects an overlap between an input x, y, 
def test_overlap(x,y,z,r,x_pos,y_pos,z_pos,cloud_rad):
    overlap=0
    #The following will not run if x_pos is an empty array
    #if x_pos.size:
    for iii in range( x_pos.size):
        dist=np.sqrt(((x-x_pos[iii])**2)+((y-y_pos[iii])**2)+((z-z_pos[iii])**2))
        if dist <= r:
            overlap=1
    ###
    #Check to see of we overlap the xrange
    dist=np.sqrt((x-xrange[0])**2)
    if dist <= r:
        overlap=1
    dist=np.sqrt((x-xrange[1])**2)
    if dist <= r:
        overlap=1
    #Check to see of we overlap the yrange
    dist=np.sqrt((y-yrange[0])**2)
    if dist <= r:
        overlap=1
    dist=np.sqrt((y-yrange[1])**2)
    if dist <= r:
        overlap=1
    #Check to see of we overlap the zrange
    dist=np.sqrt((z-zrange[0])**2)
    #if dist <= r:
    #    overlap=1
    dist=np.sqrt((z-zrange[1])**2)
    #if dist <= r:
    #    overlap=1
    ###
    return overlap;
###########################################################
def test_cloud_size(r):
    #1 for pass, 0 for fail.
    pass_fail=1
    if r > max_size_cloud:
        pass_fail=0
    
    if r < min_size_cloud:
        pass_fail=0

    return pass_fail;
###########################################################

#Open text file for writing only. (Not updating or appending)
out_file = open(filename, "w")

#out_file.write( 'X')
#out_file.write(separator)
#out_file.write('   Y')
#out_file.write(separator)
#out_file.write('    Z')
#out_file.write(separator)
#out_file.write('         Rad')
#out_file.write('\n')



x_pos=np.zeros(0)
y_pos=np.zeros(0)
z_pos=np.zeros(0)

cloud_rad=np.zeros(0)
cloud_rho=np.zeros(0)
#Write out the number of clouds to the file
out_file.write(str(n_clouds))
out_file.write('\n')
#For every cloud do the following
for ind in range(0, n_clouds):
    while_counter=1
    while while_counter == 1 :
        #Generate a x,y,x uniform random coordinate and
        # a radius of the cloud from a powerlaw
        x,y,z,r=generate_pos_rad()
        #Check to see if the cloud overlaps
        while_counter=test_overlap(x,y,z,r,x_pos,y_pos,z_pos,cloud_rad)
        
        #If the cloud overlaps, choose new random variables and
        #  retest.  If not exit while loop
        
    # Append the new coordinates and radius to the arrays
    x_pos=np.append(x_pos, x)
    y_pos=np.append(y_pos, y)
    z_pos=np.append(z_pos, z)
    cloud_rad=np.append(cloud_rad, r)
    cloud_rho=np.append(cloud_rho, rho)
        
    #Write data out to a text file
    out_file.write(str(x_pos[ind]))
    out_file.write(separator)
    out_file.write(str(y_pos[ind]))
    out_file.write(separator)
    out_file.write(str(z_pos[ind]))
    out_file.write(separator)
    out_file.write(str(cloud_rad[ind]))
    out_file.write(separator)
    out_file.write(str(cloud_rho[ind]))
    out_file.write('\n')

#Close the file
out_file.close()
#Done
