#! /usr/bin/env python-mpacts

"""
   Script to simulate cells behaving like an active polar liquid, driven by contact inhibition of locomotion
   and guided by a global director field that reorients the polarization of the cells.

   This script is based on the same framework described in:

   Smeets et al. "Emergent structures and dynamics of cell colonies
   by contact inhibition of locomotion" (2016) PNAS 113(51), 14621-14626.

   To run this script, you need to have a recent (2020) version of the Mpacts software and available in Python.
   This script only works in Python 3.
"""

#-----------------------------------------------------------------------------------------------------------------------
#If a user module exists, load this. Otherwise, load the generic specified option from mpacts itself.
def load_user_module( umod_name, alternative, cmdname ) -> None:
   import importlib
   try:
      mod = importlib.import_module( f"mpacts_usermodules.{umod_name}" )
      cmd = mod.__dict__[cmdname]
      print(f" - Loading user module: {cmd.__name__}")
   except ImportError as e:
      mod = importlib.import_module(f"{alternative}")
      cmd = mod.__dict__[cmdname]
      print(f" - Loading generic module: {cmd.__name__}")
   globals()[cmdname] = cmd

#-----------------------------------------------------------------------------------------------------------------------

#The user model will contain the open source code of all (physically relevant) parts of the simulation.
#These can be separately compiled by any reviewer / collaborator to verify the implementation.
load_user_module( "active_force"#User module
                , "mpacts.commands.force.body"#Alternative
                , "PersistentRandomForceCommand" )
load_user_module( "contractilerepulsive"#User module
                , "mpacts.contact.models.collision.contractilerepulsive"#Alternative
                , "ContractileRepulsiveMatrix")

from mpacts.commands.onarrays.linearcombination import LinearCombinationCommand
import mpacts.commands.monitors.progress as pp
from mpacts.commands.onarrays.setvalue import SetValueCommand
from mpacts.commands.time_evolution.integration import ForwardEuler_Generic
import mpacts.commands.misc.stresscalculation as stress_calculation
from mpacts.contact.detectors.multigrid import MultiGridContactDetector
from mpacts.contact.models.collision.linearforce.linearforce_matrix import LinearForceMatrix
from mpacts.commands.onarrays.average import ExponentialMovingAverageCommand

from mpacts.contact.matrix.conjugategradient import DefaultConjugateGradientSolver
from mpacts.contact.matrix.cmtypes import ConjugateGradientSolver_old, ConjugateGradientSolver
from mpacts.commands.time_evolution.integration import ForwardEuler_UncoupledOverdamped
from mpacts.core.valueproperties import Variable, VariableFunction
from mpacts.contact.models.collision.jkr.jkr_matrix import JKRMatrix
from mpacts.io.vtk_writer import VTKSerialWriter
import mpacts.geometrygenerators.pointgenerators as pointgen
import mpacts.geometrygenerators.quadgeometries as quadgen
from mpacts.io.datasave import DataSaveCommand
from mpacts.core.arrays import create_array
from mpacts.core.simulation import Simulation
import DEMutilities.enumerate_properties as ep
import DEMutilities.postprocessing.analysistools.clip as clip
import mpacts.tools.random_seed as randseed
from mpacts.tools.load_parameters import load_parameters
import DEMutilities.postprocessing.h5storage as h5s
from mpacts.core.new_units import registry as u
import mpacts.particles as prt
import os, sys, glob
import mpacts.core.command as cmd
import generate_track as gt
import numpy as np

#-----------------------------------------------------------------------------------------------------------------------------
#Utility function to compute correct hexagonal packing:
#Aspect ratio is defined as width / height
#N is the number of cells, when the size would be a square of width*width size!
def Nxz_from_N( N, aspect_ratio ) -> tuple:
   Nz = round((1 + (1+(8*N)/(3**0.5))**0.5)/(4.*aspect_ratio))
   Nx = round(3**0.5 * Nz * aspect_ratio)
   return (Nx,Nz)

#Just a short hand to get access to SI values of Variables and VariableFunctions
def si(variable) -> float:
    return variable.get_value().value_SI()

#Solver for equation of motion, reasonably chosen based on gamma.
#For final simulations, probably best to keep one fixed solver!
def automatic_solver(params) -> str:
    gamma = si(gamma_n_n)
    if gamma == 0:#No solver needed
        return "Explicit"
    elif gamma < 1.0:#More efficient for few iteration steps since no copying of data
        return "Internal"
    else:#Highly optimized, so more efficient for many iteration steps
        return "Eigen"

#-----------------------------------------------------------------------------------------------------------------------------
#First all the properties that might change in a parameter study:
mysim = Simulation( "simulation" )#We will set the actual timestep later!
params=mysim("params")

storage = h5s.H5Storage(f"{mysim.name()}.h5",openmode='wa')

#Geometry of the track:

R_inner    = Variable("R_inner", params, value=120*u('um'))
R_outer    = Variable("R_outer", params, value=250*u("um"))
L_track    = Variable("L_track", params, value=500*u('um'))

#NOTE that for the racetrack, we obviously will assume no anisotropy, so we keep the value of 'h' to 1.
h = Variable("h", params, value = 1.0 )
k_perpf    = VariableFunction("k_perp_factor", params, function='1.0')#set this to 1/$h$ for anisotropic reinforcement
g_perpf    = VariableFunction("gamma_perp_factor", params, function='$h$' )#set this to $h$ for anisotropic friction

H          = Variable("H"      , params, value=(1,0,0,0,0,0)*u('dimensionless'))
I          = Variable("I"      , params, value=(1,0,1,0,0,1)*u('dimensionless'))
Y          = Variable("Y"      , params, value=(0,0,1,0,0,0)*u("dimensionless"))

densprop  = Variable( "density"   , params, value = 1.2 * u("dimensionless") )#Cell density
gamma_n_n = Variable( 'gamma_n_n' , params, value = 0.0 * u("dimensionless") )#Viscosity
f_a_n     = Variable( 'f_a_n'     , params, value = 2.5 * u("dimensionless") )#Dimensionless rate of cell-substrate adhesion
f_c_n     = Variable( 'f_c_n'     , params, value = 0.0 * u("dimensionless") )#Dimensionless rate of cell-cell adhesion
Dr_n      = Variable( 'Dr_n'      , params, value = 0.45 * u("dimensionless") )#Rotational diffusivity rate
f_n       = Variable( "f_n"       , params, value = 4.0 * u('dimensionless'))#Acceleration 'rate' on self-reinforcement of velocity


R_cell     = Variable( "R_cell"   , params, value = 10 * u("um"))#Spread out cell size. Sets the units of length
R_std      = Variable( "R_std"    , params, value = 0.5 * u('um') )#Only to prevent crystal-like effects
v_cells    = Variable( "v_cell"   , params, value = 60* u("um/hour"))#Cell speed. Sets the units of time
xi         = Variable( "xi"       , params, value = 6 * u("Pa*hour/um"))#Just sets the units of force. Taken from Duclos et al SI.
dt_n       = Variable( "dt_n"     , params, value = 0.02 *u('dimensionless') )
N_out      = Variable( "N_out"    , params, value = 500 )#Number of output frames
rng_seed   = Variable( 'rng_seed' , params, value = randseed.produce_random_seed() )
rhv        = Variable( "rhv"      , params, value = 1e8*u("dimensionless"))
cg_solver  = Variable( "CG_solver", params, value = "Explicit")#Options: Internal/Eigen/Explicit
cg_tol     = Variable( 'CG_tol'   , params, value = 1e-4 )
cd_ue      = Variable( "cd_ue"    , params, value = 10 )#Do contact detection every cd_ue steps

cell_area  = VariableFunction("area"   , params, function='math.pi*$R_cell$**2')
gamma_s    = VariableFunction("gamma_s", params, function="$xi$*$area$")
F_cell     = VariableFunction("F_cell" , params, function='$gamma_s$*$v_cell$')
Wm         = VariableFunction("Wm"     , params, function="2*$F_cell$*$R_cell$")
tscale     = VariableFunction("tscale" , params, function="2*$R_cell$/$v_cell$")
endtime    = VariableFunction("endtime"  , params,function='250*$tscale$' )

Dr         = VariableFunction("Dr"     , params, function="$Dr_n$/$tscale$")
f          = VariableFunction("f"      , params, function="$f_n$/$tscale$")
f_a        = VariableFunction("f_a"    , params, function="$f_a_n$/$tscale$")
f_c        = VariableFunction("f_c"    , params, function="$f_c_n$/$tscale$")

#These are the timescales:
tau_Dr     = VariableFunction("tau_Dr" , params, function="1/$Dr$")
tau        = VariableFunction("tau"    , params, function='1/$f$')

k          = VariableFunction("k"      , params, function='$gamma_s$/$tau$')
b          = VariableFunction("b"      , params, function='$k$/( $v_cell$**2*$gamma_s$**3 )')
D          = VariableFunction("D"      , params, function='$F_cell$**2/$tau_Dr$')
Wa         = VariableFunction("Ws"     , params, function="$gamma_s$*$R_cell$**2*$f_a$")
Wc         = VariableFunction("Wc"     , params, function="$gamma_s$*$R_cell$**2*$f_c$")
R0         = VariableFunction("R0"     , params, function="$R_cell$*(2*$Ws$+$Wc$)/($Ws$+$Wc$)")
gamma_c    = VariableFunction("gamma_c", params, function="$gamma_n_n$*$gamma_s$")
dt         = VariableFunction("dt"     , params, function='$dt_n$*$tscale$')
cd_kd      = VariableFunction("cd_kd"  , params, function="3*$v_cell$*$dt$*$cd_ue$")
out_int    = VariableFunction("out_int", params, function="$endtime$/float($N_out$)")
alpha_s    = VariableFunction("alpha_s", params, function="2./(2*$out_int$/$dt$+1)")

width     = VariableFunction( "width"     , params, function="2*$R_outer$+$L_track$")#width of the strip (x direction)
height    = VariableFunction( "height"    , params, function="2*$R_outer$")#height of the strip (y direction)

K          = VariableFunction( "K"    , params, function='$k$*$H$+$k_perp_factor$*$k$*($I$-$H$)')
Gamma      = VariableFunction( "Gamma", params, function='$gamma_s$*$H$+$gamma_perp_factor$*$gamma_s$*($I$-$H$)+$rhv$*$gamma_s$*$Y$')

if os.path.exists('params.pickle'):
   load_parameters( mysim )

randseed.set_random_seed( rng_seed.get_value() )
mysim.set_timestep(dt.get_value())

storage.add_simulation_info( mysim )

#-----------------------------------------------------------------------------------------------------------------------------

surface_substrate = si(height)*si(width)
surface_cells = surface_substrate * si(densprop)

surface_width = si(width)**2
surface_cells_width = surface_width * si(densprop)
R = si(R_cell)
N_width = int( surface_cells_width / ( np.pi*R**2 ) )

Nx, Nz = Nxz_from_N( N_width, 1. )

#Make the cells particle container cells:
cells = prt.ParticleContainer('cells',prt.Sphere0, mysim)
box = prt.ParticleContainer("box", prt.RigidBody.compose( (prt.RigidQuad, "RigidQuads")
                                                        , (prt.RigidCapsule, "capsules")), mysim)

box.VTKWriterCmd( executeInterval=out_int, select_all=True )#Only execute once, since the box does not change
box('RigidQuads').QuadNormalsCmd( x=box('controlPoints')['x'], executeOnce=True)

energy = create_array("Scalar",'virial_energy',cells)
v_av = create_array("Vector",'v_average',cells)

#Array with combined body forces:
Fb = create_array( "Vector", "Fb", cells )
#Array with combined contact forces: Computed as F - Fb + Ff
Fc = create_array( "Vector", "Fc", cells )
#Array with active forces. Part of body forces. Without other external forces, Fa = Fb
Fa = create_array( "Vector", "Fa", cells )
#Array with friction forces. Will be added to the contact force array
Ff = create_array( "Vector", "Ff", cells )
#Total force F = Fb + Fc

x_idx = create_array("Index", "x_idx", cells)
y_idx = create_array("Index", "y_idx", cells)

L = si(width) - si(R0)
distance = max(L / float(Nx-1), L / float(3**0.5 * (Nz-1)) )

xs = pointgen.hexagonal_lattice_in_rectangle( distance, (Nx,Nz), crop_xtop = True, crop_ztop = True )
xs_np = np.array(xs)

pcyl1 = ( -si(L_track)/2, 0, 0 )
pcyl2 = (  si(L_track)/2, 0, 0 )

print(f"Starting with: {len(xs_np)}")
#Clip outside of track
xs_np = xs_np[clip.in_capsule( xs_np, pcyl1, pcyl2, si(R_outer)-0.5*si(R_cell))]
print(f"Removing out: {len(xs_np)}")
#Clip inside of track
xs_np = xs_np[clip.invert( clip.in_capsule( xs_np, pcyl1, pcyl2, si(R_inner)+0.5*si(R_cell) ) )]
print(f"Removing in: {len(xs_np)}")

xs = list(map(tuple,xs_np))

rs = [ np.random.normal(si(R_cell), si(R_std) ) for xi in xs ]

length_x = np.max(xs_np[:,0])-np.min(xs_np[:,0])+si(R0)
length_y = np.max(xs_np[:,2])-np.min(xs_np[:,2])+si(R0)


#-----------------------------------------------------------------------------------------------------------------------------

x_idx = np.argsort(np.argsort( np.array(xs)[:,0]))
y_idx = np.argsort(np.argsort( np.array(xs)[:,2]))

cells.resize( len(xs) )
cells["x"].set_array(xs)
cells["r"].set_array(rs)
cells["v_average"].set_array((0,0,0))
cells['Fa'].set_array((0,0,0))
cells['Fb'].set_array((0,0,0))
cells['Fc'].set_array((0,0,0))
cells['Ff'].set_array((0,0,0))
cells["virial_energy"].set_array(0.)
cells["x_idx"].set_array([int(el) for el in x_idx])
cells["y_idx"].set_array([int(el) for el in y_idx])

cells.SetContactMatrixDiagonalCmd(visc = Gamma)

#We calculate a (virial) stress-tensor in the contact forces.
#Its eigen values give the principal stress values
stress_calculation.AddStressCalculation( mysim, cells, evals = True, evecs = False)

actual_domain = length_x * length_y
actual_cell_area = np.pi * R**2 * len(xs)
print(f" - corrected density: {actual_cell_area / actual_domain:1.3f}")
print(f"Nx: {Nx}")
print(f"Nz: {Nz}")
print(f"N cells: {cells.size()}")

df = PersistentRandomForceCommand("PersistentRandomForceCommand", mysim, pc=cells
                                 , K = K, b=b, D_force= D
                                 , F_active = cells['Fa'] )

#-----------------------------------------------------------------------------------------------------------------------------
# Cell-cell contact model

cm = ContractileRepulsiveMatrix( pc1=cells, pc2=cells
                               , Wa = Wa, Wc = Wc
                               , gamma_normal = gamma_c
                               , gamma_tangential = gamma_c
                               , allow_dewet = False
                               , implicitness = 1. )

cmw = LinearForceMatrix( pc1=box("RigidQuads"), pc2=cells
                      , k = 1.0
                      , gamma_normal = gamma_c
                      , gamma_tangential = gamma_c
                      , implicitness = 1. )


cd = MultiGridContactDetector("CD", mysim, cmodel=cm, update_every=cd_ue, keep_distance=cd_kd)
cdw = MultiGridContactDetector("CD_wall", mysim, cmodel=cmw, update_every=cd_ue, keep_distance=cd_kd)

#-----------------------------------------------------------------------------------------------------------------------------

dtrack = 2*si(R_cell)
htrack = 2*si(R_cell)

pts = np.c_[gt.discrete_track( si(L_track), si(R_outer), d=dtrack, join=False )]
points, quads = quadgen.quads_from_curve( pts, htrack, invert=True )
points = list(points)

boxo = box.add_particle()
boxo.controlPoints.add_and_set_particles( xBF = points )
boxo.RigidQuads.add_and_set_particles( vertexIndices = quads )
#boxi.capsules.add_and_set_particles( vertexIndices = vilsegments, r = si(R_cell))
boxo.x = (0,0,0)
boxo.q = (1,0,0,0)

pts = np.c_[gt.discrete_track( si(L_track), si(R_inner), d=dtrack, join=False )]
points, quads = quadgen.quads_from_curve( pts, htrack, invert=False )
points = list(points)

boxi = box.add_particle()
boxi.controlPoints.add_and_set_particles( xBF = points )
boxi.RigidQuads.add_and_set_particles( vertexIndices = quads )
#boxi.capsules.add_and_set_particles( vertexIndices = vilsegments, r = si(R_cell))
boxi.x = (0,0,0)
boxi.q = (1,0,0,0)

#-----------------------------------------------------------------------------------------------------------------------------

if cg_solver.get_value() == "Explicit":
    #Faster way if we do not use conjugate gradient: Especially the copy method with Eigen is very expensive
    ForwardEuler_UncoupledOverdamped("TimeIntegration", mysim, pc=cells, gamma_array=cells["ContactMatrixDiagonal"])
else:
    solvers = { "Eigen"    : ConjugateGradientSolver
              , "Internal" : ConjugateGradientSolver_old }
    ConjugateGradient = DefaultConjugateGradientSolver( mysim, tolerance=cg_tol, reset_x=False, F_fric = "Ff"
                                                      , CG_solver = solvers[cg_solver.get_value()])

    #time integration
    ForwardEuler_Generic( "integrate_v_cells", mysim, pc = cells, x=cells["x"], dx=cells["v"])

#This prevents that fast velocity oscillations cause strong flicker when looking at v in sparse phases (dilute
ExponentialMovingAverageCommand("AverageVelocityCommand", mysim("loop_cmds/post_contact_cmds")
                               , array = cells["v"]
                               , result = cells["v_average"]
                               , alpha=alpha_s )

#Add a command that prints some information to the screen at a certain time interval, so we know what the simulation is doing:
printer_list = [ pp.DefaultProgressPrinter(mysim, sim_time_unit='h') ]

printer = pp.PrinterChain(printer_list)

pp.ProgressIndicator("PrintProgress", mysim, printer, print_interval=5.)

#-----------------------------------------------------------------------------------------------------------------------------

DataSaveCommand("DataSaveCommand", mysim, executeInterval = out_int )
wrt_c = cells.VTKWriterCmd(executeInterval = out_int, select_all = True)

mysim.run_until( endtime.get_value() )
storage.simulation_completed()
print("\n *** Simulation finished succesfully! ***")
