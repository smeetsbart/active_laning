import numpy as np
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
import mpacts.io.datasave as ds
import mpacts.io.vtpdatareader as read
from scipy import stats
import time,os
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

def bin( y, vy, bins = [] ):
   if len( bins ) == 0:
      bins = np.linspace( min(y), max(y), 100 )
   return stats.binned_statistic( y, vy, statistic =np.nanmean, bins = bins )[0]

#-------------------------------------------------------------------------------------
ti = op.GetData( a, 'time')
final_frame = dr[-1].index
frame = op.GetData(a, 'DataFrameIndex')
time = a.Recorder(ti/3600)
#time_mask = a.Function( takeframe, frame)
#time_mask = frame==final_frame
time_mask = (frame >= 0)
#time_mask = frame==50
op.FrameIndexPrinter( a, final_frame )
va = 60./1e6/3600
R = 10e-6

x = a.GetData("cells/x")
y = x[:,2]
v = a.GetData("cells/v_average")/va
vx = v[:,0]

ybins = np.linspace( -4.5e-3,4.5e-3,250 )

vb = a.Function( bin, y, vx, bins = ybins )

vb = a.Recorder( vb )

a.loop( dr, mask_function = time_mask )

import matplotlib.pyplot as plt

ymid = nf.bin_centers( ybins )

fig=plt.figure(figsize=(3,10))
ax=fig.add_subplot(111)
plt.pcolormesh(time(), ymid*1e3, np.transpose(vb()), cmap='coolwarm', vmin=-1, vmax=1 )
ax.set_ylabel("x (mm)", fontsize=14)
ax.set_xlabel("time (h)", fontsize=14)
plt.tight_layout()

sample_idx = os.path.abspath(os.path.curdir).split("/")[-1]

plt.savefig( f"kymograph_vx_time_{sample_idx}.png", dpi=300 )


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


