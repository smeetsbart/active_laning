import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
import mpacts.io.datasave as ds
from scipy import stats
import time,os,pickle
from scipy.io import savemat
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import matplotlib

simname = os.getcwd().split('/')[-1]
storage = h5s.H5Storage( 'simulation.h5' )

mm = 1e-3#one milimeter in meter
um = 1e-6#One micrometer in meter
hr = 3600.#One hour in seconds

vref = 60*um/hr#Reference 'free' velocity.
height  = 10*mm#Height of the patch
width   = 2*mm#Width of the patch
bin_pixel = 52*um#Pixel size of the final image
aspect_ratio = height/width
save_snapshots = False#If true, save a .png for every frame with a picture of x velocity

with open("params.pickle",'rb') as fob:
   params = pickle.load(fob)

#-------------------------------------------------------------------------------------
a = op.AnalysisContainer()
dr = ds.DataReader( "simulation", folder='./')

final_frame = dr[-1].index
a.FrameIndexPrinter( final_frame )
index = a.GetData("DataFrameIndex")
ti = a.GetData( 'time')
time = a.Recorder(ti)
mask = index <= final_frame
#mask = index <= 10
x = a.GetData("cells/x")
x = a.Recorder( x )
a.loop( dr, mask_function=mask )

x = x()
t = time()

vd = []
xd = []
yd = []
td = []
for i in range( 1, len(x)-1):
   ti = t[i]
   td.append(ti)
   xi = x[i][:,0]
   yi = x[i][:,2]
   dx = x[i+1] - x[i-1]
   dxmag = nf.norm( dx )
   #These are the ones that jumped the periodic. For now, we disregard these as 'invalids'
   #We could try to correct for the periodic jump, but have to be very careful not to introduce bias!!
   mask = dxmag < 0.25*width
   dt = t[i+1] - t[i-1]
   vi = dx/dt
   vd.append(vi[mask])
   xd.append(xi[mask])
   yd.append(yi[mask])

#Perfectly symmetric around zero
xbin = np.r_[-np.arange(0,width/2, bin_pixel)[::-1],np.arange(0,width/2, bin_pixel)[1:]]
ybin = np.r_[-np.arange(0,height/2,bin_pixel)[::-1],np.arange(0,height/2,bin_pixel)[1:]]

vxs = []
vys = []
for i in range( len(vd) ):
   vx,_,_,Nbix = stats.binned_statistic_2d( xd[i],yd[i],vd[i][:,0],statistic='mean', bins=[xbin,ybin] )
   vy,_,_,Nbiy = stats.binned_statistic_2d( xd[i],yd[i],vd[i][:,2],statistic='mean', bins=[xbin,ybin] )
   vxs.append(vx)
   vys.append(vy)

vxs = np.array(vxs)
vys = np.array(vys)

mdic = {'Vx' : vxs, "Vy" : vys, "t" : np.array(td), "x_edges" : xbin, "y_edges" : ybin, "pixel_size" : bin_pixel}
#Also add the parameters of this specific simulation in the dictionary to have no confusion
mdic.update( {el.split('/')[-1]:params[el] for el in params} )
savemat( f"binned_v_{simname}.mat", mdic )

if save_snapshots:
   fig=plt.figure(figsize=(2.0,aspect_ratio*2.0))
   ax=fig.add_subplot(111)

   for i in range(len(vxs)):
      #Normalized velocities
      vn = vxs[i] / vref
      z = np.transpose(vn)
      zm = ma.masked_invalid(z)

      cmap = matplotlib.cm.coolwarm
      cmap.set_bad('black',1.)
      xmid = xbin[:-1]+0.5*np.diff(xbin)
      ymid = ybin[:-1]+0.5*np.diff(ybin)

      plt.pcolormesh(xmid/mm, ymid/mm, zm, cmap=cmap, vmin=-1.5, vmax=1.5, shading='gouraud' )
      plt.xlim( xbin.min()/mm, xbin.max()/mm )
      plt.ylim( ybin.min()/mm, ybin.max()/mm )
      plt.axis('off')
      plt.tight_layout()
      filename = f'frame-{(i+1):04d}.png'
      print(f"Saving snapshot {filename}")
      plt.savefig(filename, dpi=300)
      ax.clear()
