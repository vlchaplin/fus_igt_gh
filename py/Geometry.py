# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:27:01 2016

@author: caskeylab
"""

import numpy as np
from numpy.linalg import inv
from math import cos,sin,tan,acos

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
    
    def from_mat(self,M):
        T=np.trace(M)+1
        
        if T>0:
            self.qs=np.sqrt(T)/2
            self.qx=(M[2,1] - M[1,2])/(4*self.qs)
            self.qy=(M[0,2] - M[2,0])/(4*self.qs)
            self.qz=(M[1,0] - M[0,1])/(4*self.qs)
        else:
            raise KeyError('Trace <-1')
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
        
    def rotinv(self):
        R=self.rot()
        return inv(R)
        
    def conjugate(self):
        return quaternion([self.qs, -self.qx, -self.qy, -self.qz])
        
    def magnitude(self):
        return np.sqrt( np.sum(np.array([self.qs, self.qx, self.qy, self.qz])**2) )
        
    
    def inverse(self):
        """
        Very ineffcient route. Need to code the more direct algorithm.
        """
        q=quaternion([1,0,0,0])
        Rq=self.rot()
        q.from_mat(inv(Rq))
        return q
        
    def dot(self, other):
        
        return self.qs*other.qs + self.qx*other.qx + self.qy*other.qy + self.qz*other.qz
        
    
def q_linterp(Q1, Q2, unit_time, ret_array=True):
    """
    Compute linear interpolation of components between Q1 and Q2. Linear interpolation is fast
    but introduces a false acceleration between Q1 and Q2. Only use if the orientation is slowly
    varying.
    
    unit_time is the fractional distance between. See q_slerp doc for other info.
    
    """
    if type(unit_time) == np.ndarray:    
        #return N quaternions
    
        q_sxyz = np.outer( 1. - unit_time, np.array( [Q1.qs, Q1.qx, Q1.qy, Q1.qz]) ) + np.outer(unit_time, np.array( [Q2.qs, Q2.qx, Q2.qy, Q2.qz]) )
        
        q_sxyz /= np.sqrt( np.sum(q_sxyz**2,axis=1) )
        
        if ret_array:
            return q_sxyz
        else:
            N = q_sxyz.shape[0]
            return [quaternion(q_sxyz[i]) for i in range(N) ]
    else:
        
        q_sxyz = (1. - unit_time) * np.array( [Q1.qs, Q1.qx, Q1.qy, Q1.qz]) + (unit_time)*np.array( [Q2.qs, Q2.qx, Q2.qy, Q2.qz]) 
        
        q_sxyz /= np.sqrt( np.sum(q_sxyz**2) )
        
        if ret_array:
            return q_sxyz
        else:
            return quaternion(q_sxyz)

def q_slerp(Q1, Q2, unit_time, ret_array=True):
    """
    Compute spherical linear interpolation of quaternions Q1 and Q2. SLERP is a little slower
    than simple linear interpolation but does not introduce second-derivative artifacts.
    Assumes intertial frame between the two points.
    
    Inputs:
    Q1 and Q2 - quaternion objects.
    unit_time - scalar or numpy array. Time at which the interpolated value will be computed, as a number between 0 and 1.0. 
        Dimensionless time from Q1 to Q2 (the fraction of time the total interval elapsed).
        
    Return one or more quaternion objects at the interpolated times, or a the raw quaternion components if ret_array=True
    """
    angle = acos ( Q1.dot(Q2) / (Q1.magnitude()*Q2.magnitude() ) )
    
    if type(unit_time) == np.ndarray:    
        #return N quaternions
    
        q_sxyz = np.outer( np.sin(angle*(1. - unit_time)), np.array( [Q1.qs, Q1.qx, Q1.qy, Q1.qz]) ) + np.outer( np.sin(angle*(unit_time)), np.array( [Q2.qs, Q2.qx, Q2.qy, Q2.qz]) )
        
        q_sxyz /= angle
        
        if ret_array:
            return q_sxyz
        else:
            N = q_sxyz.shape[0]
            return [quaternion(q_sxyz[i]) for i in range(N) ]
    else:
        
        q_sxyz = np.sin(angle*(1. - unit_time)) * np.array( [Q1.qs, Q1.qx, Q1.qy, Q1.qz]) + np.sin(angle*(unit_time))*np.array( [Q2.qs, Q2.qx, Q2.qy, Q2.qz]) 
        
        q_sxyz /= angle
        
        if ret_array:
            return q_sxyz
        else:
            return quaternion(q_sxyz)