import numpy as np
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
import mpacts.io.datasave as ds
import mpacts.io.vtpdatareader as read
from scipy import stats
import time,os,matplotlib
import numpy.ma as ma
import DEMutilities.postprocessing.plottools.mplsphereplotter as mps
import DEMutilities.postprocessing.analysistools.spatial.structurefactor as sf

sftor = sf.StructureFactor()

storage = h5s.H5Storage( 'simulation.h5' )
#params = storage('simulation_info/settings/params')

a = op.AnalysisContainer()
dr = ds.DataReader( "simulation", folder='./')
#dr = read.VTPDataReader( "simulation", folder='./')

def computeSF( x, weights, ks=[], qs=[] ):
   x_list = list(map(tuple, x))
   w_list = list(map(tuple,weights))
   return sftor( x_list, ks, qs, weights = w_list )

def takeframe( frame ):
   return frame % 10 == 0


def create_bins( array, resolution = 20e-6):
   return np.linspace( min(array), max(array), int( (max(array)-min(array))/resolution ))

def bin( x, y, vy, bins_x, bins_y ):
   return stats.binned_statistic_2d( x, y, vy, statistic = np.nanmean, bins=[bins_x,bins_y] )[0]

#-------------------------------------------------------------------------------------
ti = op.GetData( a, 'time')
final_frame = dr[-1].index
frame = op.GetData(a, 'DataFrameIndex')
time = a.Recorder(ti/3600)
#time_mask = a.Function( takeframe, frame)
time_mask = frame==final_frame
#time_mask = (frame >= 0)
#time_mask = frame==50
op.FrameIndexPrinter( a, final_frame )
va = 60./1e6/3600
R = 10e-6

pos = a.GetData("cells/x")
x = pos[:,0]
y = pos[:,2]

v = a.GetData("cells/v_average")/va
vx = v[:,0]

bins_x = a.Function( create_bins, x )
bins_y = a.Function( create_bins, y)

vb = a.Function( bin, x, y, vx, bins_x, bins_y )

a.loop( dr, mask_function = time_mask )

import matplotlib.pyplot as plt

xmid = nf.bin_centers( bins_x() )
ymid = nf.bin_centers( bins_y() )

fig=plt.figure(figsize=(1.6,10))
ax=fig.add_subplot(111)

z = np.transpose(vb())
zm = ma.masked_invalid(z)

cmap = matplotlib.cm.coolwarm
cmap.set_bad('white',1.)

plt.pcolormesh(xmid*1e3, ymid*1e3, zm, cmap=cmap, vmin=-1, vmax=1, shading='gouraud' )
#ax.set_xlabel("x (mm)", fontsize=14)
#ax.set_ylabel("y (mm)", fontsize=14)
plt.xlim( bins_x().min(), bins_x().max() )
plt.ylim( bins_y().min(), bins_y().max())
plt.axis('equal')
plt.axis('off')
#plt.tight_layout()

sample_idx = os.path.abspath(os.path.curdir).split("/")[-1]

fname = f"snapshot_{sample_idx}.png"
plt.savefig( fname , dpi=300 )
import os
os.system(f'convert -trim {fname} {fname}')


#fig, (ax1, ax2) = plt.subplots(1, 2)

#color.cycle_cmap(len( s() ), cmap='jet', ax=ax1 )
#color.cycle_cmap(len( b() ), cmap='jet', ax=ax2 )
#for si in s():
   ##L = 1/qs
   #ax1.plot( qs, si, marker='o' )
#ax1.set_xlabel("$q$ ($1/2R$)", fontsize=15)
##ax1.set_xlabel("$L$ ($2R$)", fontsize=15)
#ax1.set_ylabel("$S(q)$", fontsize=15)

#edges = np.linspace( np.min(y()), np.max(y()), 100)
#midpoints = edges[1:]-0.5*np.diff(edges)

#for bi in b():
   #ax2.plot( midpoints , bi )
   #ax2.set_xlabel("$y$ ($2R$)",fontsize=15)
   #ax2.set_ylabel("$\\langle v (y)\\rangle / v_a$", fontsize=15)


#plt.tight_layout()
#plt.show()


#storage = h5s.H5Storage("simulation.h5", openmode='a')
#results = storage.data_section( "results", overwrite=True )
#sf = results.data_section("structure_factor", overwrite=True)

#t = sf.add_data("time",time()/3600)
#q=sf.add_data("q", qs )
#s=sf.add_data("S", s(), axes=[q,t])

#storage.close()




#--------------------------------------------------------------


