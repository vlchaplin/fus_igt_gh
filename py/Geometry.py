# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:27:01 2016

@author: caskeylab
"""

import numpy as np
from math import cos,sin,tan

class quaternion:
    """ Object implementation of 3-space (length 4) quaternions"""
    def __init__(self):
        self.qs=1.0
        self.qx=0.0
        self.qy=0.0
        self.qz=0.0
        self.__rot=None
        return
    
    def __init__(self,q0xyz):
        self.qs=q0xyz[0]
        self.qx=q0xyz[1]
        self.qy=q0xyz[2]
        self.qz=q0xyz[3]
        self.__rot=None
        return
        
    def rot(self):
        if self.__rot is None:
            self.__rot=np.zeros([3,3])
           
        R = self.__rot
        
        R[0,0] =  self.qx*self.qx - self.qy*self.qy - self.qz*self.qz + self.qs*self.qs
        R[1,1] = -self.qx*self.qx + self.qy*self.qy - self.qz*self.qz + self.qs*self.qs
        R[2,2] = -self.qx*self.qx - self.qy*self.qy + self.qz*self.qz + self.qs*self.qs
        
        R[1,0] = 2 * ( self.qx*self.qy + self.qz*self.qs )
        R[2,0] = 2 * ( self.qx*self.qz - self.qy*self.qs )
        R[0,1] = 2 * ( self.qx*self.qy - self.qz*self.qs )
        R[2,1] = 2 * ( self.qx*self.qs + self.qy*self.qz )
        R[0,2] = 2 * ( self.qx*self.qz + self.qy*self.qs )
        R[1,2] = 2 * ( self.qy*self.qz - self.qx*self.qs )
        return R