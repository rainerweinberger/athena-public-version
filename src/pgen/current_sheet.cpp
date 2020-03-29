//========================================================================================
// Current sheet problem generator for Athena++ astrophysical MHD code
//
//========================================================================================
//! \file current_sheet.cpp
//  \brief Problem generator for magnetic current sheet (Gardiner&Stone 2005)
//

// C headers

// C++ headers
#include <algorithm>  // min, max
#include <cmath>      // log
#include <cstring>    // strcmp()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../defs.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/utils.hpp"

namespace {
Real vflow;
int iprob;
Real PassiveDyeEntropy(MeshBlock *pmb, int iout);
} // namespace

Real threshold;
int RefinementCondition(MeshBlock *pmb);

//----------------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin) {
  iprob = 1;  // option for new variants

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the current sheet

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::int64_t iseed = -1 - gid;
  Real gm1 = peos->GetGamma() - 1.0;

  //--- iprob=1. sin velocity perturbation with specified amplitude
  // density = 1.0, pressure 0.1, B = b0 or -b0
  if (iprob == 1) {
    // Read problem parameters
    Real amp = pin->GetReal("problem","amp");
      std::cout << "setting hydro quantities" << std::endl;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          
          phydro->u(IDN,k,j,i) = 1.0;
          phydro->u(IM1,k,j,i) = amp*sin( 3.1415 * (pcoord->x2v(j) + 1.0) );
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
              0.1/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                             SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
            }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
        std::cout << "setting B-field" << std::endl;
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = -b0;
            if(std::abs(pcoord->x1v(i)) < 0.5)
              pfield->b.x2f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  return;
}


// refinement condition: none

int RefinementCondition(MeshBlock *pmb) {
  return 0;
}

namespace {
Real PassiveDyeEntropy(MeshBlock *pmb, int iout) {
  Real total_entropy = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &r = pmb->pscalars->r;
  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
      for (int i=is; i<=ie; i++) {
        // no loop over NSCALARS; hardcode assumption that NSCALARS=1
        Real specific_entropy = -r(0,k,j,i)*std::log(r(0,k,j,i));
        total_entropy += volume(i)*w(IDN,k,j,i)*specific_entropy;  // Lecoanet (2016) eq 5
      }
    }
  }
  return total_entropy;
}
} // namespace
