#!/usr/bin/env python-mpacts
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import matplotlib.pyplot as plt
import glob,matplotlib
import numpy as np
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from scipy import interpolate
import matplotlib.colors as colors

def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def compute_Lx( d, C, t):
   xdata = d
   ydata = np.mean( C[t>t.max()-50,:],axis=0 )
   G = -np.gradient(ydata,xdata)[1:]
   return 1./G.max()

def heatmap( name, xname, yname
           , label = None
           , xlabel = None
           , ylabel = None
           , levels=None
           , func = lambda x:x
           , **kwargs ):

   if label is None:
      label = name
   if xlabel is None:
      xlabel = xname
   if ylabel is None:
      ylabel = yname

   path = f"results/{name}"

   flist = glob.glob("data/grouped_sample_*/grouped_simulation.h5")
   flist.sort()

   reader = h5s.HdfReader(flist)
   a = op.AnalysisContainer()
   dirs   = a.Recorder( a.GetData( 'simulation_info/working_directory'))

   C = a.GetData( path )
   d = a.GetData( path+"/axis_1")
   tf = a.GetData( "results/tf")
   Lx = a.Function( compute_Lx, d, C, tf )
   v = a.Recorder(Lx)

   #a.FrameIndexPrinter( len(flist))
   a.loop(reader)

   writer = h5s.H5Storage('data/pstudy.h5')

   results = writer.data_section("results/analysis", overwrite=True)

   results.copy_all_data_recursive( '/doe/values/params'#Source
                                  , '/results/params'#Destination (must not exist, will be automatically created)
                                  , samples = dirs() #Filter to 're-determine the order'
                                  )

   axis0 = results(f'/results/params/{xname}').scale( 1. ).set_attrs( unit = '', axis_label='$k$')
   axis1 = results(f"/results/params/{yname}").scale( 1. ).set_attrs( unit='', axis_label= '$f$')

   cell_radius = 10#Radius of cell in um
   func = lambda x : x/cell_radius

   v = results.add_data( name, func(v()), axes=[axis0, axis1] )

   x = v("axis_0").get_data()
   y = v("axis_1").get_data()
   v = v.get_data()

   fig = plt.figure(figsize=(0.8*4.5,0.8*4.))

   zoomf = 2
   #This is just used for the approximating transition line
   ipol = interpolate.interp2d(x,y,v, bounds_error=False, fill_value=None,kind='cubic')
   xi = np.linspace(x[0]-0.5*(x[1]-x[0]),x[-1]+0.5*(x[1]-x[0]),len(x)*zoomf)
   yi = np.linspace(y[0]-0.5*(y[1]-y[0]),y[-1]+0.5*(y[1]-y[0]),len(y)*zoomf)
   vzoom = ipol(xi,yi)

   f = lambda x,y: v[int(y),int(x) ]
   g = np.vectorize(f)

   xx = np.linspace(0,v.shape[1], v.shape[1]*100)
   yy = np.linspace(0,v.shape[0], v.shape[0]*100)
   X, Y= np.meshgrid(xx[:-1],yy[:-1])
   Z = g(X[:-1],Y[:-1])

   xvalues = np.linspace( y[0]-0.5*np.diff(y)[0], y[-1]+0.5*np.diff(y)[0], len(y)+1)
   yvalues = np.linspace( x[0]-0.5*np.diff(x)[0], x[-1]+0.5*np.diff(x)[0], len(x)+1 )

# Alternatively, you can manually set the levels
# and the norm:
   # lev_exp = np.arange(np.floor(np.log10(z.min())-1),
   #                    np.ceil(np.log10(z.max())+1))
   # levs = np.power(10, lev_exp)
   # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())

   plt.pcolormesh( xvalues, yvalues, v
                 , norm=colors.LogNorm(vmin=v.min(), vmax=v.max())
                 , cmap='BrBG',alpha=0.65,**kwargs)

   xmag = np.linspace( xvalues.min(),xvalues.max(), Z.shape[1] )
   ymag = np.linspace( yvalues.min(),yvalues.max(), Z.shape[0] )
   if levels is not None:
      plt.contour( xmag, ymag, Z, levels, colors='k', linewidths=[2], linestyles=['--'], alpha=0.5)
      #plt.contour( yi, xi, vzoom, levels=levels, colors=[(0.3,0.2,0.2)], linestyles=['--'], linewidths=2)
   #for i in range(len(x)):
      #for j in range(len(y)):
         #s = f"{v[i,j]:1.2f}"
         #ci = 'k' if v[i,j] < 0.6 and v[i,j] > 0.11 else 'w'
         #plt.annotate(  s, (y[j]+0.05,x[i]-0.30), ha='center', va='center', color=ci, fontsize=9
                     #, alpha=max(0.35,v[i,j]))
   cbar = plt.colorbar()
   cbar.set_label('x correlation length (cell radii)',fontsize=12, rotation=270, labelpad=20)
   roman_labels = False
   if roman_labels:
      #for symbol, location in zip(['I',"II",'III',"IV"],[(1,3.0),(1,10.5),(2.5,6.0),(4,9.0)]):
      for symbol, location in zip(['I',"II",'III'],[(1,3.0),(1,10.5),(4,9.0)]):
         if symbol != 'II':
            plt.annotate( symbol, location, (location[0]+0.05,location[1]+0.05)
                        , ha='left',va='bottom', color='k', fontsize=15, fontweight="bold")
         else:
            plt.annotate( symbol, location, (location[0]+0.05,location[1]-0.05)
                        , ha='left',va='top', color='k', fontsize=15, fontweight="bold")

         #if symbol != 'III':
            #plt.annotate( symbol, location, (location[0],location[1]+0.05)
                        #,ha='right',va='bottom', color='k', fontsize=15, fontweight="bold")
         #else:
            #plt.annotate( symbol, location, (location[0]+0.05,location[1]-0.15)
                        #,ha='left',va='top', color='k', fontsize=15, fontweight="bold")
         plt.plot( [location[0]],[location[1]],ls='None',marker='o', color='k', ms=6)

   max_ks = 11.25
   max_hs = 4.75

   #hs = np.linspace(0.75,4.75,10000)
   #power = 1.5

   ##fp = (c_1-c_inf)/(hs**power) + c_inf
   #fp = 6.49/hs**power + 2.86

   transition = np.loadtxt('phase_transition_fit.txt')
   hs_fit = transition[:,0]
   ks_fit = transition[:,1]

   fitting = ( ks_fit < max_ks ) * ( hs_fit < max_hs)
   hs_fit = hs_fit[fitting]
   ks_fit = ks_fit[fitting]

   plt.plot(hs_fit, ks_fit, '-k', lw=2)

   hs = np.linspace(0.8,4.75,20)
   plt.plot(hs, [6 for _ in hs],'--k')
   plt.plot([3],[6], marker='X', mec='k', mfc='w', ms=8)

   #plt.annotate( '$\sim {h_s}^{-1/2}$'
               #, xy = ( 1.7,1.6*4. )
               #, xytext = (1.65,1.6*4.85)
               #, xycoords='data'
               #, textcoords='data'
               #, ha='left'
               #, arrowprops=dict(arrowstyle="->", connectionstyle="arc3")
               #, fontsize=15 )

   plt.xticks([1,2,3,4])
   plt.yticks([0,3,6,9])

   plt.xlabel(ylabel, fontsize=14)
   plt.ylabel(xlabel, fontsize=14)
   #plt.title( label , fontsize=15)
   #plt.axis('equal')
   #plt.ylim(x[0]-0.5*(x[1]-x[0]),x[-1]+0.5*(x[1]-x[0]))
   #plt.xlim(y[0]-0.5*(y[1]-y[0]),y[-1]+0.5*(y[1]-y[0]))
   plt.tight_layout()
   plt.savefig(f"fig4c_phasediagram-simulation.png",dpi=300)
   plt.savefig(f"fig4c_phasediagram-simulation.pdf")
   return x,y,v


if __name__ == "__main__":

   x,y,v = heatmap( "C_mean", 'f_n', 'h'
              , label='$\phi$'
              , xlabel='velocity-polarity coupling'
              , ylabel="friction anisotropy"
              , levels = None
              , func = lambda x : x )
