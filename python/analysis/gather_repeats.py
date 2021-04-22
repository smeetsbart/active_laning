#!/usr/bin/env python-mpacts
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
from DEMutilities.postprocessing.operators.gather_doe_repeats import RepeatStatistics
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize
from scipy import stats
import glob, os
import numpy as np

#-----------------------------------------------------------------------------------------------------------------------

def plot( time, m, s ):
   plt.plot( t,m,lw=2, color='C0' )
   plt.fill_between( t,m-s,m+s,color='C0',alpha=0.2 )
   #plt.xscale('log')
   plt.show()

#Example statistic: The interquartile range:
class IQR:
   def __init__(self):
      self.p = [72,25]
   def __call__(self, data, axis):
      q75, q25 = np.percentile(data, self.p, axis=axis )
      return q75 - q25

#-----------------------------------------------------------------------------------------------------------------------

stats={'mean' : np.mean, "std" : np.std, "median" : np.median, "iqr" : IQR()}

g = RepeatStatistics( "simulation*.h5", 'grouped_simulation.h5' )
tf  = g.add_axis( 'tf',  'results/grid_stats/time' )
tp  = g.add_axis( 'tp',  'results/particle_stats/time')
tl  = g.add_axis( 'tl',  'results/particle_stats/lag_time')
vb  = g.add_axis( "vb",  "results/particle_stats/v_midpoints")
vbx = g.add_axis( "vbx", "results/particle_stats/vx_midpoints")

g.add_statistic( "v",    "results/grid_stats/vmag", [tf], statistics = stats  )
g.add_statistic( "ht",   "results/grid_stats/ht", [tf], statistics = stats  )
g.add_statistic( "C500", "results/grid_stats/C500", [tf], statistics = stats  )
g.add_statistic( "lw",   "results/grid_stats/lane_width", [tf], statistics = stats  )

g.add_statistic( "vp",   "results/particle_stats/vmag", [tp], statistics = stats  )
g.add_statistic( "msd",  "results/particle_stats/msd", [tl], statistics = stats  )
g.add_statistic( "pv",   "results/particle_stats/p_v", [vb], statistics = stats  )
g.add_statistic( "pvx",  "results/particle_stats/p_vx", [vbx], statistics = stats  )
g.add_statistic( "pvy",  "results/particle_stats/p_vy", [vbx], statistics = stats  )
g.add_statistic( "v_av", "results/particle_stats/vmag_mean", [], statistics = stats )
g.add_statistic( "v0_f", "results/particle_stats/v0_fit_rel", [], statistics = stats )
g.add_statistic( "Dr_f", "results/particle_stats/D_fit_rel", [], statistics = stats )
g.add_statistic( "Deff_f", "results/particle_stats/Deff_fit_rel", [], statistics = stats )

g.add_statistic( "vg_av", "results/grid_stats/vmag_mean", [], statistics = stats)

g.execute()
g.close()
#-----------------------------------------------------------------------------------------------------------------------
