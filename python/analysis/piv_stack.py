import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
import DEMutilities.postprocessing.plottools.mplsphereplotter as mps
import matplotlib.pyplot as plt
import mpacts.io.datasave as ds
from mpacts.core.new_units import registry as u
import pathlib
import time, os
import numpy as np

#-----------------------------------------------------------------------------------------------------------------------

#We save the frames in a new directory 'stack', which we create if it does not exist yet
pathlib.Path('./stack').mkdir(parents=True, exist_ok=True)

#-----------------------------------------------------------------------------------------------------------------------

#Utility function to make sure I always get the SI value of a certain parameter
def fetch_param( name, settings ):
   param = settings.f['params'][name]#This is the actual h5 entry of our parameter
   unit = param.attrs['unit'].decode('utf-8')#Actually stored as 'byte' objects
   value = float(np.array(param['data']))#np array conversion is somewhat of a hack
   return (value * u(unit)).value_SI()#Conversion to quantity and next to SI value

#-----------------------------------------------------------------------------------------------------------------------

#Ensure that the color cycler always returns the 'black' value which is rgb (0,0,0)
def color_cycler():
   return (0,0,0)

#-----------------------------------------------------------------------------------------------------------------------

#The __call__ from this class will actually perform the plotting. This function is called by the analysis operators:
class Plotter():
   def __init__(self):
      self.p = mps.SpherePlotter( figsize = (Npixel_x/dpi,Npixel_y/dpi)
                                , perpendicular_to='y-'
                                , cmap_name = 'RdBu'
                                , scale = 1/um
                                , bounds = [(-Lx/2,-20,-Ly/2), (Lx/2,20,Ly/2)])

   def __call__( self, i, x, r, c ):
      self.p.reset()
      self.p.colorcycler = color_cycler
      border_color = (1.0,0.0,0.0)#Set this to something like 'r' to ensure that borders are proper size
      self.p.axes.plot( [-Lx/2,-Lx/2],[-Ly/2, Ly/2],ls='-', color=border_color, lw=1 )
      self.p.axes.plot( [ Lx/2, Lx/2],[-Ly/2, Ly/2],ls='-', color=border_color, lw=1 )
      self.p.axes.plot( [-Lx/2, Lx/2],[-Ly/2,-Ly/2],ls='-', color=border_color, lw=1 )
      self.p.axes.plot( [-Lx/2, Lx/2],[ Ly/2, Ly/2],ls='-', color=border_color, lw=1 )
      self.p.axes.scatter( x[:,0]/um,x[:,2]/um, r, color='k' )
      #NOTE: The commented line below would generate more detailed spheres with precise radius. It is much slower though.
      #self.p.plot( x, r, None, vmin=-1.0, vmax=1.0)
      plt.xlim(-Lx/2,Lx/2)
      plt.ylim(-Ly/2,Ly/2)
      self.basename = f"stack/frame-{i:04d}"
      self.pngname = f"{self.basename}.png"
      self.tiffname = f"{self.basename}.tiff"
      self.p.axes.set_aspect('auto')
      plt.subplots_adjust(0,0,1,1)
      mps.plt.axis('off')
      plt.savefig(self.pngname, bbox_inches='tight', dpi=dpi)
      os.system(f"convert -colorspace gray -depth 8 -colors 256 -trim {self.pngname} {self.pngname}")
      png_list.append(self.pngname)

#-----------------------------------------------------------------------------------------------------------------------

png_list = []
pixel_size = 2.0#Size of one pixels in um
um = 1e-6#meter in a micrometer
h = 3600.#seconds in an hour
dpi = 300#Dpi for saving. Only used for conversion purposes.
#Radius of nucleus in number of pixels
nucleus_size = 1.5#Making this value bigger increases the size of the saved points. irrespective of the resolution
marker_size = nucleus_size / ( dpi / 72. )

settings = h5s.H5Storage( 'simulation.h5' )("simulation_info/settings")
Lx = fetch_param( "width" , settings )/um#Width in micro meter
Ly = fetch_param( "height", settings )/um#Height in micro meter
va = fetch_param( "v_cell", settings )#Velocity in m/s

#-----------------------------------------------------------------------------------------------------------------------

Npixel_x = Lx / pixel_size
Npixel_y = Ly / pixel_size
print( f"- Required resolution of frames: {Npixel_x:1.0f}x{Npixel_y:1.0f}" )

#-----------------------------------------------------------------------------------------------------------------------

#An analysiscontainer is an object that keeps track of all 'instructions' that should be performed on every
#simulation frame.
a = op.AnalysisContainer()
dr = ds.DataReader( "simulation", folder='./')
plot = Plotter()
ti = op.GetData( a, 'time')
final_frame = dr[-1].index
frame = op.GetData(a, 'DataFrameIndex')
time_mask = frame >= 0
op.FrameIndexPrinter( a, final_frame )
x = a.GetData("cells/x")
#If we just plot black circles, the 'coloring' with 'v' is not actually used. Following line has no consequence
v = a.GetData("cells/v_average")[:,0]/va

r = marker_size

a.Function( plot, frame, x, r, v )

#This line performs the actual loop over the detected mpactsd frames.
a.loop( dr, mask_function=time_mask )
#-----------------------------------------------------------------------------------------------------------------------
