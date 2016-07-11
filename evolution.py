from visual import *
from solid import solid
from surface import surface
import os
import sys

n_sol = int(sys.argv[1])
path = 'solids'
file_evol = 'evolution_' + str(n_sol) + '.data'
n_iter = 5000
fps = 1000
rho = 2500

def get_size(n_sol, path):
    """
    Method which gets the size of each solid.
    Input: folder path.
    Output: list of size
    """
    list_of_size = []

    if (path == ''):
        file_size = open('size_' + str(n_sol) + '.data','r')
    else:
        file_size = open(path + os.sep + 'size_' + str(n_sol) + '.data','r')
    
    list_of_lines = file_size.readlines()
    file_size.close()  
    
    for line in list_of_lines:
        line = line.strip()
        list_of_size.append([float(val) for val in line.split()])
        
    return list_of_size

def get_initial_position(n_sol, path):
    """
    Method which gets the initial position of each solid.
    Input: folder path, amount of solids.
    Output: list of initial position, list of initial axis.
    """                
    list_of_pos, list_of_axis = [],[]
    #list_of_aux = []

    if (path == ''):
        f = open('position_' + str(n_sol) + '.data','r')
    else:
        f = open(path + os.sep + 'position_' + str(n_sol) + '.data','r')
        
    for i in range(n_sol):
        axis = []
        for j in range(6):
            line = f.readline()       
            line = line.strip()
            #Index 0 is the position
            if (j == 0):
                list_of_pos.append([float(val) for val in line.split()])
            #Index 3,4 and 5 are the axis.
            elif (j >= 3):
                axis.append([float(val) for val in line.split()])
        list_of_axis.append(axis)
    f.close()
    return list_of_pos, list_of_axis

def create_solids(n_sol,list_of_size,list_of_pos,list_of_axis):
    """
    Method which creates the solids.
    Input: amount of solids, list of solids size, list of solids position, list of solids axis.
    Output: list of solids
    """
    list_of_solids = []
    
    for i in range(n_sol):
            list_of_solids.append(solid(list_of_size[i],list_of_pos[i],list_of_axis[i], i))        
    return list_of_solids

scene.title = "Solid Collision Simulation"
scene.background = (1., 1., 1.)
scene.fullscreen = True    
scene.forward = (0.,1.,0.)
scene.up = (0.,0.,1.)
scene.center = (0.,-8.,8.)
#scene.fov = pi/2.
scene.range = 20

#curve(pos = [[0.,0.,0.],[1.,0.,0.]],color = color.red, radius = 0.05)
#curve(pos = [[0.,0.,0.],[0.,1.,0.]],color = color.green, radius = 0.05)
#curve(pos = [[0.,0.,0.],[0.,0.,1.]],color = color.blue, radius = 0.05)

#Surface
surface({'point':(-20.,-20.),'delta':(40.,40.)},.05)

#Getting size of solids.
list_of_size = get_size(n_sol, path)

#Getting initial position of solids.
list_of_pos, list_of_axis = get_initial_position(n_sol, path)

#Creating the solids.
list_of_solids = create_solids(n_sol,list_of_size,list_of_pos,list_of_axis)

#Drawing each solid.
for solid in list_of_solids:
    solid.draw()

if (path == ''):
    try:
        f = open(file_evol,'r')
    except IOError:
        print 'File evolution_' + str(n_sol) + '.data not found'
        raise
else:
    try:
        f = open(path + os.sep + file_evol,'r')
    except IOError:
        print 'File ' + path + os.sep + 'evolution_' + str(n_sol) + '.data not found'
        raise

flag = True
for count in range(1,n_iter):    
    if flag:
        rate(fps)
        print 'Iteration: ',count
        #Reading for each solid.
        for i in range(n_sol):
            axis = []
            for j in range(6):
                line = f.readline()
                line = line.strip()
                if len(line) == 0:
                    flag = False
                    break
                #Reading the position line.
                if (j == 0):                                
                    pos = [float(val) for val in line.split()]					
                    #str_pos += ','.join(map(str,pos)) + '\n'
                if (j == 1):                                
                    vel = [float(val) for val in line.split()]                    
                    #str_vel += ','.join(map(str,vel)) + ','
                if (j == 2):                                
                    omg = [float(val) for val in line.split()]                    
                    #str_vel += ','.join(map(str,omg)) + ';\n'
                #Reading the current axis lines.
                elif (j >= 3):                     
                     axis_aux = [float(val) for val in line.split()]
                     axis.append(axis_aux)
            if flag:
                #Updating position and axis of the solid.
                list_of_solids[i].set_position(pos, axis)
                #Updating the edges.            
                #list_of_solids[i].set_edges()                
            else:
                break
    else:
        break
f.close()

#Getting the directory path.
pathname = os.getcwd()
#Getting the list of files in the directory path.
list_of_files = os.listdir(pathname)
#Remove the python file compiled.
for filename in list_of_files:
    if 'pyc' in filename: os.remove(filename)

