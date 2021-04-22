#!/usr/bin/env python-mpacts
import DEMutilities.multisim_layout as ms
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.doe_samples as doe
import numpy as  np

#Simple example of a design of experiments script where we generate a full factorial design
#and make the according simulation layout.

f = h5s.H5Storage( "pstudy.h5", openmode='wa' )
section = f.data_section( "doe", overwrite=False)

sampler = doe.FullFactorialDesign( section )

sampler.add_parameter("params/k_n"  , doe.ValueList( [ 1.,2.,3.,4.,5.,6.,7.] ) )
sampler.add_parameter("params/h"    , doe.ValueList( [ 1.0, 2.0, 3.0, 4.0 ] ) )

#Three random repeats for EACH case
sampler.add_full_repeat( "params/rng_seed", doe.ValueList(range( 3 )), unit='-', axis_label='sample' )
sample_section = sampler.generate()

ms.generate_layout( sample_section )
