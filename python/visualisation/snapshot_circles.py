import numpy as np
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
#import mpacts.io.datasave as ds
import mpacts.io.vtpdatareader as read
from scipy import stats
import time
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.plottools.mplsphereplotter as mps



storage = h5s.H5Storage( 'simulation.h5' )
#params = storage('simulation_info/settings/params')

a = op.AnalysisContainer()
#dr = ds.DataReader( "simulation", folder='./')
dr = read.VTPDataReader( "simulation", folder='./')
Lx = 900
Ly = 9000
class Plotter():
   def __init__(self):
      self.p = mps.SpherePlotter( figsize = (2.5,int(2*(Ly/Lx)))
                                , perpendicular_to='y-'
                                , cmap_name = 'RdBu'
                                , scale = 1e6
                                , bounds = [(-Lx/2,-20,-Ly/2), (Lx/2,20,Ly/2)])

   def __call__( self, i, x, r, c ):
      self.p.plot( x, r, c, vmin=-1.0, vmax=1.0 )
      self.fname = f"frame-{i:04d}.png"
      mps.plt.axis('off')
      mps.plt.tight_layout()
      self.p.savefig( self.fname )


#-------------------------------------------------------------------------------------
plot = Plotter()
ti = op.GetData( a, 'time')
final_frame = dr[-1].index
frame = op.GetData(a, 'DataFrameIndex')
time_mask = (frame == final_frame)
op.FrameIndexPrinter( a, final_frame )

va = 60./1e6/3600
x = a.GetData("cells/x")
v = a.GetData("cells/v_average")[:,0]/va
r = a.GetData("cells/r")

a.Function( plot, frame, x, r, v )

a.loop( dr, mask_function=time_mask )


#--------------------------------------------------------------


