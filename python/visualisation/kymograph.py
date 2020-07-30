import numpy as np
import matplotlib.pyplot as plt
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
import mpacts.io.datasave as ds
import mpacts.io.vtpdatareader as read
from scipy import stats
import time,os

storage = h5s.H5Storage( 'simulation.h5' )

a = op.AnalysisContainer()
dr = ds.DataReader( "simulation", folder='./')

def takeframe( frame ):
   return frame % 10 == 0

def bin( y, vy, bins = [] ):
   if len( bins ) == 0:
      bins = np.linspace( min(y), max(y), 100 )
   return stats.binned_statistic( y, vy, statistic =np.nanmean, bins = bins )[0]

#-------------------------------------------------------------------------------------
ti = op.GetData( a, 'time')
final_frame = dr[-1].index
frame = op.GetData(a, 'DataFrameIndex')
time = a.Recorder(ti/3600)
time_mask = (frame >= 0)#Selecting all frames for the kymograph
op.FrameIndexPrinter( a, final_frame )
va = 60./1e6/3600#Velocity in m/s
R = 10e-6#Cell radius

x = a.GetData("cells/x")
y = x[:,2]
v = a.GetData("cells/v_average")/va
vx = v[:,0]

ybins = np.linspace( -4.5e-3,4.5e-3,250 )

vb = a.Function( bin, y, vx, bins = ybins )

vb = a.Recorder( vb )

a.loop( dr, mask_function = time_mask )

ymid = nf.bin_centers( ybins )

fig=plt.figure(figsize=(3,10))
ax=fig.add_subplot(111)
plt.pcolormesh(time(), ymid*1e3, np.transpose(vb()), cmap='coolwarm', vmin=-1, vmax=1 )
ax.set_ylabel("x (mm)", fontsize=14)
ax.set_xlabel("time (h)", fontsize=14)
plt.tight_layout()

sample_idx = os.path.abspath(os.path.curdir).split("/")[-1]

plt.savefig( f"kymograph_vx_time_{sample_idx}.png", dpi=300 )

#--------------------------------------------------------------
