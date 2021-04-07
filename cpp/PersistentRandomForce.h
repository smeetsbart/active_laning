//------------------------------------------------------------------------------
// DEMeter10 Library
// Copyright (c) 2009 Engelbert CM Tijskens
//------------------------------------------------------------------------------

#ifndef DEMeter_PersistentRandomForce_h_INCLUDED
#define DEMeter_PersistentRandomForce_h_INCLUDED

#include <PrimitiveTypes/R3.h>
#include <ETility/Property.h>
#include <ArLite/LoopFunctorCommand.h>
#include <ArLite/FindAndSetArrayPointer.h>
#include <ETility/Random.h>
#include <DEMeter/Sim/Simulation.h>
#include <ETility/TreeRelationship.h>

namespace Ar = ArLite;

namespace DEMeter
{//=============================================================================
    class PersistentRandomForceLoopFunctor : public Ar::ParallelCompatible
    {
    public:
        ET::VectorVariateGenerator< R3::Vector_t, boost::normal_distribution<R3::Scalar_t> > vvg_;
        ET_MAKE_MEMBER( R3::Scalar_t,  D,   public);
        ET_MAKE_MEMBER( R3::SymmetricMatrix_t,  K,    public);
        ET_MAKE_MEMBER( R3::Scalar_t,  b,     public);
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
            Fa += dt_ * ( K_ * v - b_*Fa.squaredNorm()*Fa ) + sqrt( 2. * dt_ * D_ ) * generate_frand();
            F += Fa;
        }
    };
 //=============================================================================
 // Command to compute a random force with a random direction and a given standard deviation
    class PersistentRandomForceCommand : public Ar::LoopFunctorCommand
 //=============================================================================
    {
    public:
        PersistentRandomForceCommand( std::string name, boost::intrusive_ptr<ET::BaseObject> parent = NULL )
            : Ar::LoopFunctorCommand( Default_Name( "PersistentRandomForceCmd", name ), parent, "loop_cmds/body_force_cmds")
            , sim_(this)
        {
            addProperty( ftor_.D(), "D_force"
            , "Random force diffusivity"
            , ET::required, ET::Units::newton*ET::Units::newton / ET::Units::second ) ;
            addProperty( ftor_.K(), "K"
            , "Stiffness coefficient tensor (N/m) that drives the persistence in the direction of the velocity"
            , ET::required, ET::Units::newton / ET::Units::meter ) ;
//             addProperty( ftor_.mean(), "mean_force"
//             , "Mean force. The active force is a deviation from this force", ET::optional, R3::Vector_t(0.,0.,0.)) ;
            addProperty( Fa_, "F_active"
            , "Array in which the state of the active force will be stored. "
              "Do not forget to ensure that this array starts with reasonable value (e.g. 0,0,0)"
            , ET::required );
            addProperty( ftor_.b(), "b"
            , "Restoring coefficient of quadratic term."
            , ET::required/*, ET::Units::second / ET::Units::meter*/ );
            addProperty( ftor_.dim(), "dimensions"
            , "Vector with 1 for dimensions where the random process should be simulated. E.g. 2d xy = (1,1,0)"
            , ET::optional, R3::Vector_t(1.,1.,1.) );
//             addProperty( ftor_.H(), "H"
//             , "Director"
//             , ET::required/*, ET::Units::second / ET::Units::meter*/ );

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
        ET::TreeRelationship<Simulation> sim_;
        PersistentRandomForceLoopFunctor ftor_;
//         boost::intrusive_ptr<Ar::Array<R3::Vector_t> > Fc_;
        boost::intrusive_ptr<Ar::Array<R3::Vector_t> > F_;
        boost::intrusive_ptr<Ar::Array<R3::Vector_t> > v_;
        boost::intrusive_ptr<Ar::Array<R3::Vector_t> > Fa_;
    };
//=============================================================================
}// namespace DEMeter
#endif // ifndef DEMeter_PersistentRandomForce_h_INCLUDED
