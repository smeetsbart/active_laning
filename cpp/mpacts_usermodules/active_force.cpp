#include <boost/python.hpp>
#include <ETility/BaseObject.h>
#include <ETility/PYDEMETER_CLASS.h>
#include <ETility/Command.h>
#include <ETility/Random.h>
#include <ETility/TreeRelationship.h>
#include <ArLite/ArrayManager.h>
#include <ArLite/LoopFunctorCommand.h>
#include <ArLite/FindAndSetArrayPointer.h>
#include <DEMeter/Sim/Simulation.h>
#include <PrimitiveTypes/R3.h>

namespace User
{

//----------------------------------------------------------------------------------------------------------------------

    //This the loop functor. This one will actually do the work. LoopFunctorCommand calls this functor for very particle.
    class PersistentRandomForceLoopFunctor : public Ar::ParallelCompatible
    {
    public:
        ET::VectorVariateGenerator< R3::Vector_t, boost::normal_distribution<R3::Scalar_t> > vvg_;
        ET_MAKE_MEMBER( R3::Scalar_t,  D,   public);
        ET_MAKE_MEMBER( R3::SymmetricMatrix_t,  K,    public);
        ET_MAKE_MEMBER( R3::Scalar_t,  b,     public);
        ET_MAKE_MEMBER( R3::Scalar_t,  gamma,     public);
        ET_MAKE_MEMBER( R3::Scalar_t,  dt,      public );
        ET_MAKE_MEMBER( R3::Vector_t,  dim, public);

        PersistentRandomForceLoopFunctor() : vvg_( R3::Vector_t(0.,0.,0.), R3::Vector_t( 1.,1.,1.) ) {}

        R3::Vector_t generate_frand(){
            R3::Vector_t F_rand = R3::Vector_t(0.,0.,0.);
            vvg_(F_rand);
            return F_rand.cwiseProduct(dim_);
        }

        void operator() ( size_t index
                        , R3::Vector_t &F
                        , const R3::Vector_t &v
                        , R3::Vector_t &Fa )
        {
            R3::Scalar_t sb = gamma_ > 0 ? b_*gamma_*gamma_ *  v.squaredNorm()
                                         : b_               * Fa.squaredNorm();
            Fa += dt_ * ( K_ * v - sb * Fa ) + sqrt( 2. * dt_ * D_ ) * generate_frand();
            F += Fa;
        }
    };

//----------------------------------------------------------------------------------------------------------------------

 //This is the LoopFunctorCommand that interfaces the module with the Python API. Provides input (parameters)
 //and accesses the necessary arrays. It puts the command by default under 'body forces' in the simulation tree.
 //=============================================================================
    class PersistentRandomForceCommand : public Ar::LoopFunctorCommand
 //=============================================================================
    {
    public:
        PersistentRandomForceCommand( std::string name, boost::intrusive_ptr<ET::BaseObject> parent = NULL )
            : Ar::LoopFunctorCommand( Default_Name( "PersistentRandomForceCmd", name )
                                    , parent, "loop_cmds/body_force_cmds" )
            , sim_(this)
        {
            addProperty( ftor_.D(), "D_force"
            , "Random force diffusivity"
            , ET::required, ET::Units::newton*ET::Units::newton / ET::Units::second ) ;
            addProperty( ftor_.K(), "K"
            , "Stiffness coefficient tensor (N/m) that drives the persistence in the direction of the velocity"
            , ET::required, ET::Units::newton / ET::Units::meter ) ;
            addProperty( Fa_, "F_active"
            , "Array in which the state of the active force will be stored. "
              "Do not forget to ensure that this array starts with reasonable value (e.g. 0,0,0)"
            , ET::required );
            addProperty( ftor_.b(), "b"
            , "Restoring coefficient of quadratic term."
            , ET::required/*, ET::Units::second / ET::Units::meter*/ );
            addProperty( ftor_.gamma(), "gamma"
            , "If given, this will use b*gamma^2*v^2 instead of b*Fa^2 for the restorative term. Both are identical for free cells"
            , ET::optional, 0. );
            addProperty( ftor_.dim(), "dimensions"
            , "Vector with 1 for dimensions where the random process should be simulated. E.g. 2d xy = (1,1,0)"
            , ET::optional, R3::Vector_t(1.,1.,1.) );
        }

        virtual void afterInit()
        {
            Ar::FindAndSetArrayPointer( pc(), "F", F_, ET_THROW_INFO(this));
            Ar::FindAndSetArrayPointer( pc(), "v", v_, ET_THROW_INFO(this));
        }

        virtual void execute()
        {
            ftor_.dt() = sim_->getTimestep();
            Loop( ftor_, F_, v_, Fa_);
        }

     private:
        ET::TreeRelationship<DEMeter::Simulation> sim_;
        PersistentRandomForceLoopFunctor ftor_;
        boost::intrusive_ptr<Ar::Array<R3::Vector_t> > F_;
        boost::intrusive_ptr<Ar::Array<R3::Vector_t> > v_;
        boost::intrusive_ptr<Ar::Array<R3::Vector_t> > Fa_;
    };
}

//----------------------------------------------------------------------------------------------------------------------

//This generates the python module using Boost Python: The name to load from python will be "PersistentRandomForceCommand"

using namespace boost::python;

BOOST_PYTHON_MODULE(active_force) //IMPORTANT the filename and the module name here should be identical!
{
    PYDEMETER_CLASS( User::PersistentRandomForceCommand, "PersistentRandomForceCommand"
            , (ET::Command)(ET::BaseObject)
            , "Persistent active force with velocity-alignment, restoring term and random fluctuations.");
}

//----------------------------------------------------------------------------------------------------------------------
