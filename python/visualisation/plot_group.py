#!/usr/bin/env python-mpacts
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import matplotlib.pyplot as plt
import glob
import numpy as np

def heatmap( name, xname, yname
           , label = None
           , xlabel = None
           , ylabel = None ):

   if label is None:
      label = name
   if xlabel is None:
      xlabel = xname
   if ylabel is None:
      ylabel = yname

   path = f"results/{name}"

   flist = glob.glob("grouped_sample_*/grouped_simulation.h5")
   flist.sort()

   reader = h5s.HdfReader(flist)
   a = op.AnalysisContainer()
   dirs   = a.Recorder( a.GetData( 'simulation_info/working_directory'))

   v = a.Recorder(a.GetData( path ))

   #a.FrameIndexPrinter( len(flist))
   a.loop(reader)

   writer = h5s.H5Storage('pstudy.h5')

   results = writer.data_section("results/analysis", overwrite=True)

   results.copy_all_data_recursive( '/doe/values/params'#Source
                                  , '/results/params'#Destination (must not exist, will be automatically created)
                                  , samples = dirs() #Filter to 're-determine the order'
                                  )

   axis0 = results(f'/results/params/{xname}').scale( 1. ).set_attrs( unit = '', axis_label='$k$')
   axis1 = results(f"/results/params/{yname}").scale( 1. ).set_attrs( unit='', axis_label= '$f$')

   v = results.add_data( name, v(), axes=[axis0, axis1] )

   x = v("axis_0").get_data()
   y = v("axis_1").get_data()
   v = v.get_data()

   fig = plt.figure(figsize=(5,5))

   plt.pcolormesh( np.linspace( x[0]-0.5*np.diff(x)[0], x[-1]+0.5*np.diff(x)[0], len(x)+1 )
                 , np.linspace( y[0]-0.5*np.diff(y)[0], y[-1]+0.5*np.diff(y)[0], len(y)+1 )
                 , np.transpose(v)
                 , edgecolor = 'k'
                 , cmap='coolwarm')

   for i in range(len(x)):
      for j in range(len(y)):
         s = f"{v[i,j]:1.2f}"
         plt.annotate(  s, (x[i],y[j]), ha='center', va='center', color='k', fontsize=10)

   plt.xticks(x)
   plt.yticks(y)

   plt.xlabel(xlabel, fontsize=15)
   plt.ylabel(ylabel, fontsize=15)
   plt.title( label , fontsize=15)
   plt.savefig(f"map_{name}_{xname}_{yname}.png",dpi=300)
   writer.close()


if __name__ == "__main__":

#g.add_statistic( "v_av", "results/particle_stats/vmag_mean", [], statistics = stats )
#g.add_statistic( "v0_f", "results/particle_stats/v0_fit_rel", [], statistics = stats )
#g.add_statistic( "Dr_f", "results/particle_stats/D_fit_rel", [], statistics = stats )
#g.add_statistic( "Deff_f", "results/particle_stats/Deff_fit_rel", [], statistics = stats )

   heatmap( "v_av_mean", 'k_n', 'h'
          , label='$\\langle v_p \\rangle/v_0$'
          , xlabel='$k_c$'
          , ylabel="$h$" )

   heatmap( "vg_av_mean", 'k_n', 'h'
          , label='$\\langle v_g \\rangle/v_0$'
          , xlabel='$k_c$'
          , ylabel="$h$" )

   heatmap( "v0_f_mean", 'k_n', 'h'
          , label='$v_{0,f}/v_0$'
          , xlabel='$k_c$'
          , ylabel="$h$" )

   heatmap( "Dr_f_mean", 'k_n', 'h'
          , label='$D_{\\varphi,f}/D_{\\varphi}^*$'
          , xlabel='$k_c$'
          , ylabel="$h$" )

   #heatmap( "Deff_f", 'k_n', 'h'
          #, label='$D_{\\mathrm{eff},f}$'
          #, xlabel='$k_c$'
          #, ylabel="$h$" )
