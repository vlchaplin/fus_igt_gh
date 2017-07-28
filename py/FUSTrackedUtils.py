# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:33:51 2017

@author: vlchaplin@gmail.com
"""

import numpy as np
import numpy.linalg 

import re;
import csv
import math



def point_register(X,Y,w=None,n_t=None):
    """
    Adapted from Fitzpatrick matlab code:
    
    [R,t,FRE] = POINT_REGISTER(X,Y,w,n_t)
    Find R,t such that R*X + t*ones(1,N) ~ Y,
    where N = number of columns of X and Y. Optional weightings may be
    specified in the vector w. n_t is the number of points which contribute
    to the calculation of the translation vector t. If n_t is zero,
    then t is all zeros. If n_t is omitted, then all points contribute to t.
    
     Author: J. Michael Fitzpatrick
     Creation: Fall 2001 (for my CS 357 course, Image Processing)
    
     Each column of X and each column of Y represents a
     K-dimensional point, where K is the number of rows.
     X and Y must have the same number of points (columns).
     FRE = RMS fiducial registration error.
     FREcomponents = R*X + t*ones(1,N) - Y.
     Uses algorithm 8.1 and notation of pp. 469-70 of
     J. Michael Fitzpatrick, Derek L. G. Hill, and Calvin R. Maurer, Jr.,
     "Image Registration", Chapter 8 of "Handbook of Medical Imaging,
     Volume 2, Medical Image Processing and Analysis",
     Milan Sonka and J. Michael Fitzpatrick, eds., SPIE Press,
     Bellingham, Wa, 2000.
    """

    #N x 3 in most cases
    (N,K) = np.shape(X)
    (Ny,Ky) = np.shape(Y)

    if w is None:
        w = np.ones(N)
    if n_t is None:
        n_t = N

    wsq = w**2
    wsq_normed = wsq/np.sum(wsq) 


    #X*np.transpose(np.repeat([wsq_normed],3,axis=0))
    Xavg = np.mean(X,axis=0)
    Yavg = np.mean(Y,axis=0)
    X0 = X - Xavg
    Y0 = Y - Yavg

    H = np.transpose(X0).dot( np.diag(wsq_normed).dot( Y0 ) )

    (U,S,V) = numpy.linalg.svd(H)
    #evidently need this when comparing to fitzpatrick code
    V = np.transpose(V)
    U = np.transpose(U)

    D = np.diag(np.append(np.ones(K-1), numpy.linalg.det(np.dot(V,U)) ))

    R = V.dot( np.dot(D,U) )
    t = Yavg - R.dot(Xavg)


    FREcomp = np.transpose( np.tensordot( R,X,axes=[1,1]) ) + t - Y
    FRE = np.sqrt( np.sum(FREcomp*FREcomp)/np.sum(wsq)  )

    return (R,t,FRE)
    
    
def write_fcsv(filename, vectors,filemode="w",end="\n"):
    """
    Write the xyz vectors into a Slicer-compatible .fcsv file. 
    """
    text_file = open(filename, filemode)
    text_file.write("# Markups fiducial file version = 4.6"+end)
    text_file.write("# CoordinateSystem = 0"+end)
    text_file.write("# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID"+end)
    for i in range(len(vectors)):
        (x,y,z)=vectors[i]
        text_file.write( ("vtkMRMLMarkupsFiducialNode_%d,%0.4f,%0.4f,%0.4f,0,0,0,1,1,1,1,absphys-%d,,vtkMRMLScalarVolumeNode1"+end)%(i,x,y,z,i))
    text_file.close()
    
def read_fcsv(filename, nmax=None, headerlength=3):
    """
    Load Slicer-compatible .fcsv into numpy array of vectors. 
    """
    csvfile=open(filename)
    reader=csv.reader(csvfile)
    text=list(reader)
    csvfile.close()
    
    for n in range(headerlength):
        text.pop(0)
        
    if nmax is not None:
        text = text[0:nmax]
    
    xyz = np.array( list(map( lambda x: x[1:4], text )), dtype='float')

    return xyz
    
def readNDI_csv(filename, nmax=None, startcol=0, delimiter=','):
    """
    Returns (quats, vecs, errors) 
    
    Loads the csv format created by "NDI Track" software.
    Quaternions, positions, and error is a returned as numpy arrays.
    """
    csvfile=open(filename)
    reader=csv.reader(csvfile, delimiter=delimiter)
    text=list(reader)
    csvfile.close()
    
    if re.search( 'Tools', text[0][0]):
        text.pop(0)
        
    if nmax is not None:
        text = text[0:nmax]
    
    err = np.array( list(map( lambda x: x[startcol + 12], text )), dtype='float')
    txyz = np.array( list(map( lambda x: x[startcol + 9:startcol + 12], text )), dtype='float')
    q0xyz = np.array( list(map( lambda x: x[startcol + 5:startcol + 9], text )), dtype='float')
    
    return (q0xyz, txyz, err)
