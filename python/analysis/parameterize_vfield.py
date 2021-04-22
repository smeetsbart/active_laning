import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.postprocessing.representation.cmap_tools as cmpt
from mpacts.core.new_units import registry as u
import DEMutilities.numeric_functions as nf
import mpacts.io.datasave as ds
from scipy import stats
import time,os,pickle
from scipy.io import savemat, loadmat
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import matplotlib

def corx( Vx, dj ):
   xval = np.shape( Vx )[1]
   R = np.zeros( (np.shape(Vx)[0], np.shape(Vx)[2]) )
   for j in range(xval-dj):
      product = Vx[:,j,:]*Vx[:,j+dj,:]
      R += product / np.abs( product )
   R /= xval-dj
   return np.mean( R, axis=1 )
   #inner_prod = Vx[:,:,:-dj]*Vx[:,:,dj:]
   #sum = np.sum(np.sum( inner_prod / np.abs(inner_prod),axis=1),axis=1)
   #return sum / ( np.shape(Vx)[1]*(np.shape(Vx)[2]-dj) )

def si_entry( entry ):
   value = float(entry.get_data())
   unit  = entry.get_attr("unit")
   quantity = value*u(unit)
   si_value = quantity.value_SI()
   return si_value

#-----------------------------------------------------------------------------------------------------------------------

storage = h5s.H5Storage( 'simulation.h5', 'a' )

#The name of the mat file is stored directly in the original database file, so we cannot make a mistake here:
matname = storage("image_info/mat_name").value()

#Get thse values from the original simulation:
ppath = 'simulation_info/settings/params'
vref   = si_entry(storage(f'{ppath}/v0'))
height = si_entry(storage(f'{ppath}/height'))
width  = si_entry(storage(f'{ppath}/width') )
Rcell  = si_entry(storage(f'{ppath}/R_cell'))
bin_pixel = si_entry(storage("image_info/bin_pixel"))
px_2 = si_entry(storage("image_info/real_pixel"))

line_width = 5#For computing the vx vy profile in the middle of the field. This number of pixels
N_dbins = 20

storage("image_info").add_data( "mid_pixels", line_width )
storage("image_info").add_data( "autocorrelation_pixels", N_dbins)

#-----------------------------------------------------------------------------------------------------------------------
data = loadmat(matname)

um = 1e-6#One micrometer in meter
hr = 3600.#One hour in seconds

Vx = data['Vx']/um*hr
Vy = data['Vy']/um*hr
Vmag = np.sqrt( Vx**2 + Vy**2 )#Magnitude of the velocity, in micro-meter / hr
Vmag_mean = np.mean( np.mean(Vmag, 1), 1)#Averaging over x and y, but not over time.

d = [ px_2*i for i in range(1,N_dbins) ]
d = np.array(d)

Cdx = np.array([ corx(Vx,dx+1) for dx in range(len(d)) ])
Cdx = np.transpose(Cdx)#Time major as the other ones

#We compute the nearest value for the correlation at 500 here:
n_pix_500um = int( 500*um / px_2 )
Cd500 = corx(Vx, n_pix_500um + 1 )

x_midpoints = data['x_edges'][0,:][:-1]+0.5*np.diff(data['x_edges'][0,:])
y_midpoints = data['y_edges'][0,:][:-1]+0.5*np.diff(data['y_edges'][0,:])
Vx_line = np.mean(Vx[:,np.argsort(np.abs(x_midpoints))[:line_width],:],1)
Vy_line = np.mean(Vy[:,np.argsort(np.abs(x_midpoints))[:line_width],:],1)

sum_vy = np.nansum(np.nansum(np.abs(Vy),axis=2),axis=1)
sum_vx = np.nansum(np.nansum(np.abs(Vx),axis=2),axis=1)

th = data['t'][0,:]/hr
ht = 1 - sum_vy / sum_vx

zero_crossings = np.array([((el[:-1] * el[1:]) < 0).sum() for el in Vx_line])
height = data['y_edges'].max() - data['y_edges'].min()
average_lane_width = height/(zero_crossings+1)/um

#for i,Cdi in enumerate(Cdx):
   #plt.plot( d, Cdi,alpha=i/len(Cdx),color='C0' )
#plt.show()

results = storage.data_section( "results/grid_stats", overwrite=True)

t = results.add_data( "time", th
                    , description = 'Time'
                    , label='$t$'
                    , unit = 'hour')

vt = results.add_data("vmag", Vmag_mean
                     , axes = [t]
                     , description = 'Magnitude of the velocity'
                     , unit='um/hr'
                     , label='$v_\\mathrm{mag}$')

vm = results.add_data("vmag_mean", np.mean(Vmag_mean)*um/hr / vref )
d = results.add_data( "d", d/um, label='$dx$', unit='um' )

results.add_data( "ht", ht
                , description = 'Global anisotropy ratio'
                , label='$h_t$'
                , axes= [t]
                , unit = '')

results.add_data( "lane_width", average_lane_width
                , label='$d_l$'
                , unit='um'
                , axes=[t])

results.add_data( "C_dx", Cdx
                , label='$C(dx)$'
                , description = "Velocity autocorrelation"
                , axes=[t,d])

results.add_data( "C500", Cd500
                , label="C(500)"
                , description="Velocity autocorrelation at 500 um"
                , axes = [t])

storage.close()
#-----------------------------------------------------------------------------------------------------------------------
