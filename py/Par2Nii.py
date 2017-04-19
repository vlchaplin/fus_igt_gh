#!/usr/bin/python3
"""
Created on Wed Apr  5 10:43:00 2017

@author: vlchaplin@gmail.com
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

parser.add_argument("-p", metavar="inparts",help="Specify the input image parts and order. For PAR files (with 4 parts) from gstudy this should typically be '-p M R I P'. ** For exports from 7T directly this is usually '-p R I M P'. Only one of these will be saved to the .nii file ('-op'). If -p isn't used, then the entire image stack gets saved.", type=str, nargs='+'  )
parser.add_argument("-op", metavar="outpart", help="Specify which part is written to the .nii file,  e.g., '-op M'. Only a single part is allowed for Slicer compatability.", type=str )

parser.add_argument("--preview", action='store_true',help="Optionally preview the output image part. Requires -op."  )


args = parser.parse_args()  

#if len(args.files)==0 :
#    quit()
    
    
if args.o and len(args.files)>1:
    print("Only one input file allowed when using '-o' argument")
    quit()
    
if args.op and not args.p:
    print("'-op' given without '-p'" )
    quit()

if args.p and not args.op:
    print("'-p' given without '-op'" )
    quit()
    
if args.preview and not args.op:
    print("'-preview' given without '-op'" )
    quit()
    
parse_parts={}
out_part=''
errflag=False
if args.p:
    for i in range(len(args.p)):
        
        if (len(args.p[i])!=1) or (not re.match('[MPRI]',args.p[i])):
            print("unrecognize image part '%s'"%args.p[i])
            errflag=True
            continue
        
        parse_parts[args.p[i]] = i
        
    if (len(args.op)!=1) or (not re.match('[MPRI]',args.op)):        
        print("unrecognize output image part '%s'"%args.p[i])
        errflag=True
    else:
        out_part= args.op
    
    if errflag:
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
        
    if len(parse_parts)==0:
        newImageObject = nibabel.nifti1.Nifti1Image(par_obj.get_data(), par_obj.affine)
    else:
        parts = parse_scan_parts(par_obj.get_data(), num_parts=len(parse_parts) )
        
        output_array = parts[parse_parts[out_part]].squeeze()
        newImageObject = nibabel.nifti1.Nifti1Image(output_array, par_obj.affine)
           
        if args.preview:
            pg.image(output_array.transpose(), title='Preview of output array slices')

    
    nibabel.loadsave.save(newImageObject, outfile)

if args.preview:
    QtGui.QApplication.instance().exec_()