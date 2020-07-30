#include <boost/python.hpp>
#include <ETility/BaseObject.h>
#include <ETility/PYDEMETER_CLASS.h>
#include <ETility/Command.h>
#include <ArLite/ArrayManager.h>
#include <ArLite/LoopFunctorCommand.h>
#include <ArLite/FindAndSetArrayPointer.h>
#include <DEMeter/ContactModels/ContactModelBase.h> //Only required for contactmodels.
#include <PrimitiveTypes/R3.h>
#include <DEMeter/Sim/Simulation.h>
#include <ETility/TreeRelationship.h>

namespace User
{
    //------------------------------------------------------------------------------------------------------------------
      void add_gradient( R3::Scalar_t k_normal
                       , R3::Scalar_t k_tangential
                       , const R3::Vector_t &normalUnitVector
                       , R3::SymmetricMatrix_t &tensor )
      {
         R3::SymmetricMatrix_t g = normalUnitVector.dyadic_self();
         R3::SymmetricMatrix_t I = R3::TypeTraits::oneElement(R3::SymmetricMatrix_t());
         tensor += k_tangential * (I-g) + k_normal * g;
      }


     struct ContactMatrixStateBase
     {
        ContactMatrixStateBase() : element_( R3::SymmetricMatrix_t::Zero() )
        {}
        R3::SymmetricMatrix_t element_;
     };

    class ContractileRepulsiveMatrix : public DEMeter::CM::ContactModelBase
    {
        public:
            ContractileRepulsiveMatrix( std::string const & name, boost::intrusive_ptr<ET::BaseObject> parent )
                : DEMeter::CM::ContactModelBase( name, parent)
                , sim_(this)
            {
                //This tells the contactdetector how it should provide memory for the contactstate (ie memory reserved for every contact candidate it finds)
                this->CS() = boost::intrusive_ptr<DEMeter::CM::ContactStateDescriptorBase>(
                        new DEMeter::CM::ContactStateDescriptor< ContactMatrixStateBase, DEMeter::CM::ContactStateSetter>("", NULL));
                addProperty( Wa_, "Wa", "Cell-substrate adhesion energy", ET::required );
                addProperty( Wc_, "Wc", "Contractile energy", ET::required );
                addProperty( gamma_normal_, "gamma_normal", "Normal friction coefficient (kg/s)", ET::required );
                addProperty( gamma_tangential_, "gamma_tangential", "Tangential friction coefficient (kg/s)", ET::required );
                addProperty( allow_dewet_, "allow_dewet", "Must be false, maintained for compatibility", ET::optional, false );
                addProperty( lambda_, "implicitness", "Choose the semi-implicit integration method. "
                                                      "(0 for explicit Euler, 1 for implicit Euler, 0.5 for Cranck-Nicholson)."
                           , ET::optional, 0.0 );


            }
            virtual R3::Scalar_t doContactModel( unsigned int p1, unsigned int p2, void* contactstate, bool newcontact )
            {
                //This guaranteed to be safe as long as the ContactstateDescriptor above is passed the same class we cast to here.
                ContactMatrixStateBase* cs = static_cast<ContactMatrixStateBase*>(contactstate);
                //Note: the constructor of ExampleContactstate is NOT called automatically. This can be done like this:
                if(newcontact)
                {
                    cs= new(cs) ContactMatrixStateBase();//use of placement new is absolutely required, since the memory is already located.
                }
                //check if the two "things" here assumed to be spheres collide
                R3::Vector_t x12 = pos_pc2_[p2]- pos_pc1_[p1];
                R3::Scalar_t d12 = x12.norm();

                R3::Scalar_t overlap = radius_pc1_[p1] + radius_pc2_[p2] - d12;
                R3::Vector_t nuv = x12 / d12;

                //Reset the contact force to zero:
                Fc_ = R3::Vector_t(0.,0.,0.);

                if( overlap > 0)
                {
                    k_ = 0.;
                    R3::Scalar_t rmean = 0.5 * (radius_pc1_[p1] + radius_pc2_[p2]);

                    //Note that using this definition: d12 = distance_12 - radius
                    d12 = rmean - overlap;

                    R3::Scalar_t inv_rmean= 1./rmean;

                    R3::Scalar_t f = 2*inv_rmean*(Wa_ - d12 * inv_rmean * (Wa_ + Wc_ ));

                    Fc_ += f * nuv;
                    k_ = 2*(Wa_ + Wc_)*inv_rmean*inv_rmean;

                    add_gradient( gamma_normal_, gamma_tangential_, nuv, gradF_v_ );
                    add_gradient( k_, k_, nuv, gradF_x_);

                    F1_pc1_[p1] = F1_pc1_[p1] - Fc_;
                    F2_pc2_[p2] = F2_pc2_[p2] + Fc_;

                    contactMatrixElement_ = gradF_v_ + lambda_ * dt_ * gradF_x_;
                    cs->element_ = - contactMatrixElement_;
                }
                else{
                    cs->element_ *= 0.;
                }
                return overlap;
            }

            virtual void doAtStartOfContactLoop()
            {
                DEMeter::CM::ContactModelBase::doAtStartOfContactLoop();
                dt_ = sim_->getTimestep();
            }

            virtual void call_init()
            {
                Ar::FindAndSetArrayPointer( pc1(), "r", radius_pc1_, ET_THROW_INFO(this));
                Ar::FindAndSetArrayPointer( pc2(), "r", radius_pc2_, ET_THROW_INFO(this));
                Ar::FindAndSetArrayPointer( pc1(), "x", pos_pc1_, ET_THROW_INFO(this));
                Ar::FindAndSetArrayPointer( pc2(), "x", pos_pc2_, ET_THROW_INFO(this));
                Ar::FindAndSetArrayPointer( pc1(), "F", F1_pc1_, ET_THROW_INFO(this));
                Ar::FindAndSetArrayPointer( pc2(), "F", F2_pc2_, ET_THROW_INFO(this));
            }

            virtual void afterInit( )
            {
                samePC() = false;
                if( pc1() == pc2())
                {
                    samePC() = true;
                }
                //and we need to register for changes in the layout/size of the pc's
                pc1()->add_observer( this);
                pc2()->add_observer( this);
                //start the chain of inits
                call_init();
            }

            virtual bool is_parallel_compatible() { return false; } //Note: if your implementation is safe to be executed in parallel feel free to change this.
            virtual std::unique_ptr<ContactModelBase> clone()
            {
                return std::unique_ptr<ContactModelBase>(new  ContractileRepulsiveMatrix( *this )); //returns a copy of this object.
            }
        public:
            ET::TreeRelationship<DEMeter::Simulation> sim_;
            R3::Point_t* pos_pc1_; //ALL positions of pc1
            R3::Point_t* pos_pc2_;
            R3::Scalar_t* radius_pc1_;//ALL radii of pc2
            R3::Scalar_t* radius_pc2_;
            R3::Vector_t* F1_pc1_;
            R3::Vector_t* F2_pc2_;
            R3::Vector_t Fc_;

            R3::Scalar_t lambda_;
            R3::Scalar_t dt_;
            R3::Scalar_t k_;
            R3::Scalar_t Wa_;
            R3::Scalar_t Wc_;
            bool allow_dewet_;
            R3::Scalar_t gamma_normal_, gamma_tangential_;
            R3::SymmetricMatrix_t gradF_x_;//Gradient of force wrt position, i.e. stiffness tensor
            R3::SymmetricMatrix_t gradF_v_;//Gradient of force wrt velocity, i.e. friction tensor (gamma)
            R3::SymmetricMatrix_t contactMatrixElement_;
    };

}

using namespace boost::python;

BOOST_PYTHON_MODULE(contractilerepulsive) //IMPORTANT the filename and the module name here should be identical!
{
    PYDEMETER_CLASS( User::ContractileRepulsiveMatrix, "ContractileRepulsiveMatrix"
            , (DEMeter::CM::ContactModelBase)(ET::BaseObject)
            , "Linear contractility up to a user-specified, overlap, after which "
              "(hard) repulsion follows a power-law slope.");
}
