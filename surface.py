from visual import *
from numpy import arange

#Input: reference = {'point':(x,y),'delta':(delta_x,delta_y)}, delta_h
def surface(reference,delta_h):
    #Surface
    
    for i in arange(reference['point'][0],reference['point'][0] + reference['delta'][0] + 2*delta_h,delta_h):
        curve(pos = [(i,reference['point'][1],0.),(i,reference['point'][1] + reference['delta'][1] + delta_h,0.)],color=color.black)

    for j in arange(reference['point'][1],reference['point'][1] + reference['delta'][1] + 2*delta_h,delta_h):
        curve(pos = [(reference['point'][0],j,0.),(reference['point'][0] + reference['delta'][0] + delta_h,j,0.)],color=color.black)
    
    '''
    alpha = pi/4.
    gx = reference['point'][0] + .5*reference['delta'][0]
    gy = reference['point'][1]
    ground = box(pos=(gx,gy,-0.025),axis = (1.,0.,0.),up = (0.,1.,0.),length = reference['delta'][0],width = 0.05, height = reference['delta'][1], color = color.gray(0.9))

    sx = reference['point'][0]
    sy = reference['point'][1]
    hz = .5*sin(alpha)*reference['delta'][0]/2.
    hx = .5*cos(alpha)*reference['delta'][0]/2.
    print hz
    slope = box(pos=(sx,sy,-0.025),axis = (1.,0.,0.),up = (0.,1.,0.),length = .5*reference['delta'][0],width = 0.05, height = reference['delta'][1], color = color.gray(0.5))
    slope.rotate(angle = alpha,axis = (0.,1.,0.), origin = (0.,0.,0.))
    slope.pos = slope.pos + (-hx,0.,hz)
    '''
