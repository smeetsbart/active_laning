import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
from mpacts.core.new_units import registry as u
import mpacts.io.datasave as ds
from scipy import stats
import time,os,pickle
from scipy.io import savemat
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import matplotlib

def si_entry( entry ):
   value = float(entry.get_data())
   unit  = entry.get_attr("unit")
   quantity = value*u(unit)
   si_value = quantity.value_SI()
   return si_value

##Computes the 'absolute displacement' starting from x0
def absolute_disp( x, dpx, dpy ):
   xa = np.zeros( np.shape(x) )
   xa[0] = x[0]
   for i in range(1,len(x)):
      d = x[i] - x[i-1]
      #Correct periodics in x-direction
      d[:,0] = np.where( d[:,0] >  0.5*dpx, d[:,0]-dpx, d[:,0] )
      d[:,0] = np.where( d[:,0] < -0.5*dpx, d[:,0]+dpx, d[:,0] )
      #Correct periodics in the y-direction:
      d[:,2] = np.where( d[:,2] >  0.5*dpy, d[:,2]-dpy, d[:,2] )
      d[:,2] = np.where( d[:,2] < -0.5*dpy, d[:,2]+dpy, d[:,2] )
      xa[i] = xa[i-1] + d
   return xa

simname = os.getcwd().split('/')[-1]
storage = h5s.H5Storage( 'simulation.h5','a' )
ppath = 'simulation_info/settings/params'
dr = ds.DataReader( "simulation", folder='./')

#Get thse values from the original simulation:
vref   = si_entry(storage(f'{ppath}/v0'))
height = si_entry(storage(f'{ppath}/height'))
width  = si_entry(storage(f'{ppath}/width') )
Rcell  = si_entry(storage(f'{ppath}/R_cell'))

#Periodic y translation is a function of initial positions (hex grid). Do not change this formula:
x0 = np.array(dr[0]('cells/x'))
length_y = np.max(x0[:,2])-np.min(x0[:,2])+Rcell

#Periodic displacement in the x-direction:
dpx = width-3**0.5*Rcell
#Periodic displacement in the y-direction:
dpy = length_y - 0.5*Rcell

image_info = storage.data_section( 'image_info', overwrite=True )

#Settings for the image generation. We will save these in the database file as well
bin_pixel = (52*u("um")).value_SI()#Pixel size of the final image
Rcell = (10*u('um')).value_SI()
px_2 = bin_pixel / 2.
aspect_ratio = height/width
save_snapshots = False#If true, save a .png for every frame with a picture of x velocity

image_info.add_data("bin_pixel", bin_pixel, unit='m', description='Bin shifting window size (2x actual pixel size)')
image_info.add_data("real_pixel", px_2, unit='m', description='Real pixel size of the resulting image')
image_info.add_data("aspect_ratio", aspect_ratio, unit='', description='Aspect ratio of resulting image')

with open("params.pickle",'rb') as fob:
   params = pickle.load(fob)

#-------------------------------------------------------------------------------------
a = op.AnalysisContainer()


final_frame = dr[-1].index
a.FrameIndexPrinter( final_frame )
index = a.GetData("DataFrameIndex")
ti = a.GetData( 'time')
time = a.Recorder(ti)
mask = index <= final_frame
#mask = index <= 5
x = a.GetData("cells/x")
x = a.Recorder( x )

v = a.GetData("cells/v")
v = a.Recorder( v )
a.loop( dr, mask_function=mask )

x = x()
xabs = absolute_disp( x, dpx, dpy )
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
   #Use absolute positions to not have problems with periodic jumps:
   dx = xabs[i+1] - xabs[i-1]
   dxmag = nf.norm( dx )
   #These are the ones that jumped the periodic. For now, we disregard these as 'invalids'
   #We could try to correct for the periodic jump, but have to be very careful not to introduce bias!!
   mask = dxmag < 0.25*width
   if not np.all(mask):
      raise ValueError("Something went wrong with the periodic correction. Jumps should not be that high")

   dt = t[i+1] - t[i-1]
   vi = dx/dt
   vd.append(vi[mask])
   xd.append(xi[mask])
   yd.append(yi[mask])

#Perfectly symmetric around zero
px4 = px_2/2.
xbin0 = np.r_[-np.arange(0,width/2-Rcell , bin_pixel)[::-1],np.arange(0,width/2-Rcell, bin_pixel)[1:]]
ybin0 = np.r_[-np.arange(0,height/2-Rcell, bin_pixel)[::-1],np.arange(0,height/2-Rcell,bin_pixel)[1:]]
xbin = np.interp( np.arange(2*len(xbin0)-1),2*np.arange(len(xbin0)),xbin0)[:-1]+px4
ybin = np.interp( np.arange(2*len(ybin0)-1),2*np.arange(len(ybin0)),ybin0)[:-1]+px4

vxs = []
vys = []
for i in range( len(vd) ):
   vaxis = []
   for axi in [0,2]:
      vx1,_,_,Nbix = stats.binned_statistic_2d( xd[i],yd[i],vd[i][:,axi],statistic='mean', bins=[xbin0-px4,ybin0-px4] )
      vx2,_,_,Nbix = stats.binned_statistic_2d( xd[i],yd[i],vd[i][:,axi],statistic='mean', bins=[xbin0-px4,ybin0+px4] )
      vx3,_,_,Nbix = stats.binned_statistic_2d( xd[i],yd[i],vd[i][:,axi],statistic='mean', bins=[xbin0+px4,ybin0-px4] )
      vx4,_,_,Nbix = stats.binned_statistic_2d( xd[i],yd[i],vd[i][:,axi],statistic='mean', bins=[xbin0+px4,ybin0+px4] )

      vx1 = np.repeat(np.repeat(vx1, 2, axis=0), 2, axis=1)[1: ,1: ]
      vx2 = np.repeat(np.repeat(vx2, 2, axis=0), 2, axis=1)[1: ,:-1]
      vx3 = np.repeat(np.repeat(vx3, 2, axis=0), 2, axis=1)[:-1, 1:]
      vx4 = np.repeat(np.repeat(vx4, 2, axis=0), 2, axis=1)[:-1,:-1]

      #We use nanmean, to not propagate nans in the entire field when one of them was empty
      vx = np.nanmean( [ vx1,vx2,vx3,vx4 ],axis=0 )
      vaxis.append(vx)

   vxs.append(vaxis[0])
   vys.append(vaxis[1])

vxs = np.array(vxs)
vys = np.array(vys)

mdic = {'Vx' : vxs, "Vy" : vys, "t" : np.array(td), "x_edges" : xbin, "y_edges" : ybin, "pixel_size" : bin_pixel}
#Also add the parameters of this specific simulation in the dictionary to have no confusion
mdic.update( {el.split('/')[-1]:params[el] for el in params} )

res_x = len(xbin)-1
res_y = len(ybin)-1
matname = f"binned_v_{simname}.mat"
image_info.add_data("mat_name", matname, description='Name of the .mat file that contains all binned info')
image_info.add_data("resolution", f"{res_x}x{res_y}", description="Resolution of the binned images")

savemat(matname, mdic )

if save_snapshots:
   fs = 2.
   fig=plt.figure(figsize=(fs,aspect_ratio*fs))
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
      plt.savefig(filename, dpi=int( np.shape(vxs)[0]/fs ))
      os.system(f"convert -trim {filename} {filename}")
      ax.clear()

storage.close()
