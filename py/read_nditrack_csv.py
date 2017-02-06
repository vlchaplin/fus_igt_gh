# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:38:42 2016

@author: caskeylab
"""

import pandas
import numpy as np
from scipy import io
from math import *;


import matplotlib.pyplot as plt
from matplotlib import image

import sys
import glob

sys.path.append('/home/caskeylab/Vandiver/repos/fus_igt/py')

from Geometry import quaternion as quat


#csv_capture="~/DataDrive/NDI_calibrations/H115_ShortCone_Linear_20160816/ndi_tool_tracker_stylus340_body339_absolute.csv"
#csv_capture="~/DataDrive/NDI_calibrations/H115_ShortLinear_20160623/ndi_tool_tracker_stylus340_wrt_body339.csv"
#csv_capture="~/DataDrive/NDI_calibrations/H115_ShortCone_Linear_20160816/ndi_tool_tracker_stylus340_body339_0818.csv"

csv_capture="~/DataDrive/NDI_calibrations/H115_ShortLinear_20161002/stylus340_body339_calib.csv"
csv_capture="~/DataDrive/NDI_calibrations/H115_ShortCone2_20161122_Fomblin/339_340_calib_20161122_Fomblin.csv"

csv_capture="/home/caskeylab/DataDrive/NDI_calibrations/H115_ShortLinear_20170111/stylus340_339.csv"


frame= pandas.read_csv(csv_capture)

refTxyz = np.mean( frame[['Tx','Ty','Tz']].as_matrix(), axis=0)
refQsxyz = np.mean( frame[['Q0','Qx','Qy','Qz']].as_matrix(), axis=0)


stylusTxyz = np.mean( frame[['Tx.1','Ty.1','Tz.1']].as_matrix(), axis=0 )
stylusQsxyz = np.mean( frame[['Q0.1','Qx.1','Qy.1','Qz.1']].as_matrix(), axis=0)

rq = quat(refQsxyz)
sq = quat(stylusQsxyz)

offsetToStylus_Absolutexyz = stylusTxyz - refTxyz





#csv_capture="~/DataDrive/NDI_calibrations/H115_ShortCone_Linear_20160816/junk_stylus340_wrt_body339.csv"
#
#csv_capture="~/DataDrive/NDI_calibrations/H115_ShortCone_Linear_20160816/junk_340_wrt_339.csv"
#
#
#frame= pandas.read_csv(csv_capture)
#
#refTxyz = np.mean( frame[['Tx','Ty','Tz']].as_matrix(), axis=0)
#refQsxyz = np.mean( frame[['Q0','Qx','Qy','Qz']].as_matrix(), axis=0)
#
#
#stylusTxyz = np.mean( frame[['Tx.1','Ty.1','Tz.1']].as_matrix(), axis=0 )
#stylusQsxyz = np.mean( frame[['Q0.1','Qx.1','Qy.1','Qz.1']].as_matrix(), axis=0)
#
#rq = quat(refQsxyz)
#sq = quat(stylusQsxyz)
#
#offsetToStylus_rel_xyz = stylusTxyz 


print( "Offset mm (absolute xyz): %0.2f, %0.2f, %0.2f" %(tuple(offsetToStylus_Absolutexyz)))

#print( "Offset mm (wrt 339 xyz): %0.2f, %0.2f, %0.2f" %(tuple(offsetToStylus_rel_xyz)))



y=np.array([46.69, -86.80, -1091.5])

xp=np.array([-74.379, -80.758,-1086.388])

#In [149]: y-xp
#Out[149]: array([ 121.069,   -6.042,   -5.112])

probeToTip=np.array([-19.15, 0.48, -159.86])

x=np.array([-227.26,-112.49,-1047.82])

R=sq.rot()
TipPosition = stylusTxyz + R.dot(probeToTip)

offsetVec = TipPosition - refTxyz

R339=rq.rotinv()

#create a new transfrom with this final vector, and put in under the body 339 hierarchy 
offsetInBodyFrame = R339.dot(offsetVec)

print( "Offset mm (bodyf frame): %0.2f, %0.2f, %0.2f" %(tuple(offsetInBodyFrame)))


#%%
#fiducials
fidofiles = glob.glob("/home/caskeylab/DataDrive/MRscans/NHP_20161122/conefids_*csv")
fi=1
fidoutput=[]
for f in fidofiles:
    fidFrame = pandas.read_csv(f)
    stylusTxyz = np.mean( fidFrame[['Tx.1','Ty.1','Tz.1']].as_matrix(), axis=0 )
    stylusQsxyz = np.mean( fidFrame[['Q0.1','Qx.1','Qy.1','Qz.1']].as_matrix(), axis=0)
    
    print("F%d = %f,%f,%f"%(fi,stylusTxyz[0],stylusTxyz[1],stylusTxyz[2]))
    fidoutput.append(stylusTxyz)
    fi+=1
    
#%%
outf = "/home/caskeylab/DataDrive/MRscans/NHP_20161122/PhysicalFids_20161122.fcsv"
fh = open(outf,'w')
header=["# Markups fiducial file version = 4.5","# CoordinateSystem = 0",
        "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n"]
fh.writelines(header)

for fi in range(len(fidoutput)):
    (x,y,z)=fidoutput[fi][[0,1,2]]
    label="phys%d"%(fi+1)
    fh.write("vtkMRMLMarkupsFiducialNode_%d,%f,%f,%f,0,0,0,1,1,1,0,%s,,vtkMRMLScalarVolumeNode1\n"%(fi,x,y,z,label) )
fh.close()