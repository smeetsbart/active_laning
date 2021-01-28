import numpy as np
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
import mpacts.io.datasave as ds
import mpacts.io.vtpdatareader as read
from mpltools import color
from scipy import stats
import time
import DEMutilities.postprocessing.plottools.mplsphereplotter as mps
import DEMutilities.postprocessing.analysistools.spatial.structurefactor as sf

sftor = sf.StructureFactor()

storage = h5s.H5Storage( 'simulation.h5' )

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
   return stats.binned_statistic( y, vy, statistic ='mean', bins = bins )[0]

#-------------------------------------------------------------------------------------
ti = op.GetData( a, 'time')
final_frame = dr[-1].index
frame = op.GetData(a, 'DataFrameIndex')
time = a.Recorder(ti)
time_mask = a.Function( takeframe, frame)
op.FrameIndexPrinter( a, final_frame )
va = 60./1e6/3600
R = 16e-6

x = a.GetData("cells/x") / ( 2*R )
y = x[:,2]
v = a.GetData("cells/v_average")/va
vx = v[:,0]

diameters_min = 2.0
diameters_max = 30.0
Nq = 100
qs = np.linspace(1./diameters_max,1./diameters_min, Nq)
#qs = 1./np.log10(np.logspace(diameters_min,diameters_max, Nq))

s = a.Recorder( a.Function( computeSF, x, v, ks=[(0,0,1)], qs=qs.tolist() ) )
b = a.Recorder( a.Function( bin, y, vx ))

a.loop( dr, mask_function = time_mask )

storage = h5s.H5Storage("simulation.h5", openmode='a')
results = storage.data_section( "results", overwrite=True )
sf = results.data_section("structure_factor", overwrite=True)

t = sf.add_data("time",time()/3600)
q=sf.add_data("q", qs )
s=sf.add_data("S", s(), axes=[q,t])

storage.close()




#--------------------------------------------------------------


