from scipy.io import loadmat
import DEMutilities.postprocessing.h5storage as h5s
from mpacts.core.new_units import registry as u
import matplotlib.pyplot as plt
import glob, os
import numpy as np

def si_entry( entry ):
   value = float(entry.get_data())
   unit  = entry.get_attr("unit")
   quantity = value*u(unit)
   si_value = quantity.value_SI()
   return si_value

def load_data():
   matfile = glob.glob("*.mat")[0]
   D = loadmat(matfile)
   return D


def save_snapshot( data, i ):
   D = data
   simname = os.getcwd().split('/')[-1]
   storage = h5s.H5Storage( 'simulation.h5','a' )
   ppath = 'simulation_info/settings/params'
   vref   = si_entry(storage(f'{ppath}/v0'))
   storage.close()

   xe = D['x_edges'][0,:]
   ye = D['y_edges'][0,:]

   Vx = D['Vx'][i]/vref

   Y,X = np.meshgrid(ye,xe)

   plt.pcolormesh(X,Y,Vx, cmap='coolwarm', vmax = 1, vmin=-1)
   plt.xlim( X.min(), X.max() )
   plt.ylim( Y.min(), Y.max())
   plt.axis('equal')
   plt.axis('off')

   figname = f"{simname}-{i:04d}.png"
   print(f" - Created snapshot {figname}")
   plt.savefig(figname)
   os.system(f"convert -trim {figname} {figname}")

def save_movie( data ):
   Nframes = len(data['Vx'])
   for i in range(Nframes):
      save_snapshot(data, i)

if __name__ == "__main__":
   data = load_data()
   #save_movie(data)
   save_snapshot( data, -1 )
