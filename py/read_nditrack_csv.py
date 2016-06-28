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

sys.path.append('/home/caskeylab/Vandiver/repos/fus_igt/py')

from Geometry import quaternion as quat


csv_capture="~/DataDrive/NDI_calibrations/H115_ShortLinear_20160623/ndi_tool_tracker_stylus340_wrt_body339.csv"

frame= pandas.read_csv(csv_capture)

refTxyz = np.mean( frame[['Tx','Ty','Tz']].as_matrix(), axis=0)
refQsxyz = np.mean( frame[['Q0','Qx','Qy','Qz']].as_matrix(), axis=0)


stylusTxyz = np.mean( frame[['Tx.1','Ty.1','Tz.1']].as_matrix(), axis=0 )
stylusQsxyz = np.mean( frame[['Q0.1','Qx.1','Qy.1','Qz.1']].as_matrix(), axis=0)

rq = quat(refQsxyz)
sq = quat(stylusQsxyz)