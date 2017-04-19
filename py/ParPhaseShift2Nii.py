#!/usr/bin/python3
"""
Created on Wed Apr 19 11:37:17 2017

@author: caskeylab
"""

import nibabel
import os.path
import MRDataAnalysis
import argparse
import re
import numpy as np

import pyqtgraph as pg


def parse_scan_parts(imageset, num_dyns=None, num_parts=1, ordering=0):
    """
    Given an NX x NY x NS x D array of images, this routine 
    parses the ordering of the last dimension into a set of
    3D/4D image stacks, one per image part.
    
    The number of image parts (e.g., 'M','P', etc.)
    is specified by num_parts, and the number of scan dynamics by num_dyns.
    
    num_parts*num_dyns must equal the length of the last dimension, D. 
    
    If ordering==0, the last dimension will be interpreted as sorted by 
    dynamic first, then by image part:
    (dyn0, part0)
    (dyn0, part1)
    (dyn1, part0)
    (dyn1, part1)
    etc.
    
    
    If ordering!=0, the opposite ordering is assumed.
    
    Returns: a tuple of 4D image stacks of length num_parts.
    """
    
    if len(imageset.shape)==4:
        (nx,ny,ns,D) = imageset.shape   
    else:
        (nx,ny,ns) = imageset.shape
        D=0

    if num_dyns is None:
        num_dyns = int(D/num_parts)
     
    if D!=(num_dyns*num_parts):
        raise ValueError("num_dyns*num_parts must equal length of last dimension")
    
    result_set = []
    
    if num_dyns >=1:    
    
        for partnum in range(num_parts):
            if ordering:
                dynamics = partnum*num_dyns + np.arange(0,num_dyns,dtype=int)
            else:
                dynamics = np.arange(partnum,num_dyns*num_parts,num_parts,dtype=int)
            
            result_set.append(imageset[:,:,:, dynamics])
    else:
        result_set.append( imageset )
        
    return tuple(result_set)


def change_file_extension(file,newext):
    """
    Remove the trailing '.xxx' and replace with newext, and return.
    newext should include '.'
    """
    path = os.path.dirname(file)
    parts=os.path.basename(file).split(".")
    (name,ext) = ( ".".join(parts[0:-1]), parts[-1] )
    
    return os.path.join(path , name + newext)
    



## Begin main ##
parser = argparse.ArgumentParser()
parser.add_argument("files", metavar="input_par_files",help="List of PAR/REC files to write as NIFTI using nibabel", type=str, nargs='+' )
parser.add_argument("-o", metavar="output_nifti",help="Optional output name for .nii file. Only valid with a single input argument.", type=str )

parser.add_argument("-dyn", metavar="dyn", required=True, help="Specify the dynamic at which to compute phase shift wrt '-b'", type=int  )
parser.add_argument("-p", metavar="inparts", required=True, help="Specify the input image channels and order. For PAR files (with 4 parts) from gstudy this should typically be '-p M R I P'. ** For exports from 7T directly this is usually '-p R I M P'.", type=str, nargs='+'  )

parser.add_argument("-b", metavar="baseline", help="Specify the baseline dynamic (default is 0)", default=0, type=int  )


parser.add_argument("--preview", action='store_true',help="Optionally preview the output phase change."  )


args = parser.parse_args()  

#if len(args.files)==0 :
#    quit()
    
    
if args.o and len(args.files)>1:
    print("Only one input file allowed when using '-o' argument")
    quit()
    
    
parse_parts={}
use_RI=True
errflag=False
if args.p:
    for i in range(len(args.p)):
        
        if (len(args.p[i])!=1) or (not re.match('[MPRI]',args.p[i])):
            print("unrecognize image part '%s'"%args.p[i])
            errflag=True
            continue
        
        parse_parts[args.p[i]] = i
        
    if errflag:
        quit()
        
have_RI = 'R' in parse_parts and 'I' in parse_parts
have_MP = 'P' in parse_parts and 'M' in parse_parts

if not (have_RI or have_MP):
    print('Must have at least M,P or R,I image parts')
    quit()
    
        
if args.preview:
    from pyqtgraph.Qt import QtCore, QtGui
    pg.mkQApp()
        
for parfile in args.files:
    
    par_obj=nibabel.load(parfile)

    if not args.o:
        outfile = change_file_extension(parfile,".nii")
    else:
        outfile = args.o
        
    parts = parse_scan_parts(par_obj.get_data(), num_parts=len(parse_parts) )
    
    if have_RI:
        R = parts[parse_parts['R']]
        I = parts[parse_parts['I']]
        
        complex_data = R + 1j*I
    else:
        M = parts[parse_parts['M']]
        P = parts[parse_parts['P']]
        
        complex_data = M * np.exp(1j*P)
        
    
    #num_dyns = complex_data.shape[3]    

    output_array = np.angle(  complex_data[:,:,:,args.b] * np.conjugate( complex_data[:,:,:,args.dyn])  )
    newImageObject = nibabel.nifti1.Nifti1Image(output_array, par_obj.affine)
       
    if args.preview:
        pg.image(output_array.transpose(), title='Preview of output array slices')

    
    nibabel.loadsave.save(newImageObject, outfile)

if args.preview:
    QtGui.QApplication.instance().exec_()