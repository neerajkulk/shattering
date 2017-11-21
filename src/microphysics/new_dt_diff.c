#include "../copyright.h"
/*============================================================================*/
/*! \file new_dt_diff.c
 *  \brief Computes stability constraint on timestep for all diffusive
 *   processes currently implemented in code.
 *
 *  These include:
 *     - Ohmic dissipation, Hall effect, ambipolar diffusion
 *     - Navier-Stokes and Braginskii viscosity
 *     - isotropic and anisotropic thermal conduction
 *  The function returns maximum inverse of dt for all Domains at all Levels.
 *  With MPI, this value is calculated only for the Grids being updated on
 *  this processor.  The calling function new_dt() is responsible for finding
 *  the maximum over all processors.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - new_dt_diff()  - computes maximum inverse of dt */
/*============================================================================*/

#include <stdio.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/* private function */
Real max_kappa(DomainS *pD);
Real max_nu(DomainS *pD);

/*----------------------------------------------------------------------------*/
/*! \fn Real diff_dt(MeshS *pM)
 *  \brief Computes diffusion timestep */
Real new_dt_diff(MeshS *pM)
{
  Real max_dti_diff=(TINY_NUMBER);
  Real dxmin,qa;
  int nl,nd;
#ifdef RESISTIVITY
  int i,j,k;
  GridS *pG;

/* Calculate the magnetic diffusivity array */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {

        pG=pM->Domain[nl][nd].Grid;

        get_eta(pG);
      }
    }
  }
#endif
#ifdef THERMAL_CONDUCTION
  double kappa, kappa_max = TINY_NUMBER;
#endif
#ifdef VISCOSITY
  double nu, nu_max = TINY_NUMBER;
#endif

/* Calculate minimum dx.  Always given by Grid on highest level of refinement */

  dxmin = pM->dx[0]/pow(2,((pM->NLevels)-1));
  if (pM->Nx[1] > 1) dxmin = MIN( dxmin, (pM->dx[1]/pow(2,((pM->NLevels)-1))) );
  if (pM->Nx[2] > 1) dxmin = MIN( dxmin, (pM->dx[2]/pow(2,((pM->NLevels)-1))) );

  qa = (dxmin*dxmin)/4.0;
  if (pM->Nx[1] > 1) qa = (dxmin*dxmin)/8.0;
  if (pM->Nx[2] > 1) qa = (dxmin*dxmin)/16.0;

#ifdef THERMAL_CONDUCTION
  /* Calculate the maximum kappa */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        /* Note: don't put function calls in MAX statements! */
        kappa     = max_kappa(&(pM->Domain[nl][nd]));
        kappa_max = MAX(kappa_max, kappa);
      }
    }
  }

  max_dti_diff = MAX(max_dti_diff,(kappa_max/qa));
#endif

#ifdef VISCOSITY
  /* Calculate the maximum nu */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        /* Note: don't put function calls in MAX statements! */
        nu     = max_nu(&(pM->Domain[nl][nd]));
        nu_max = MAX(nu_max, nu);
      }
    }
  }

  max_dti_diff = MAX(max_dti_diff,(nu_max/qa));
#endif

#ifdef RESISTIVITY
/* Since resistivities can vary from cell to cell, must loop over all cells */
  for (nl=pM->NLevels-1; nl>=0; nl--){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pG = pM->Domain[nl][nd].Grid;

        dxmin = pG->dx1;
        if (pG->Nx[1] > 1) dxmin = MIN( dxmin, (pG->dx2) );
        if (pG->Nx[2] > 1) dxmin = MIN( dxmin, (pG->dx3) );

        qa = (dxmin*dxmin)/4.0;
        if (pG->Nx[1] > 1) qa = (dxmin*dxmin)/8.0;
        if (pG->Nx[2] > 1) qa = (dxmin*dxmin)/16.0;

        for (k=pG->ks; k<=pG->ke; k++) {
        for (j=pG->js; j<=pG->je; j++) {
        for (i=pG->is; i<=pG->ie; i++) {

          max_dti_diff = MAX( max_dti_diff, ((pG->eta_Ohm[k][j][i] +
                              pG->eta_AD[k][j][i])/qa) );

        }}}
        if (Q_Hall > 0.0) {
          for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
          for (i=pG->is; i<=pG->ie; i++) {

            max_dti_diff = MAX( max_dti_diff, fabs(pG->eta_Hall[k][j][i])/qa);

          }}}
        }
      }
    }
  }
#endif

  return max_dti_diff;
}

#ifdef THERMAL_CONDUCTION
Real max_kappa(DomainS *pD)
{
  GridS *pG = (pD->Grid);

  int i, il, iu, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;

  Real d, T;
  Real x1, x2, x3;

  Real kappa, kappa_max = (TINY_NUMBER); /* NOT zero */

  il = is - 1;  iu = ie + 1;
  if (pG->Nx[1] > 1){
    jl = js - 1;  ju = je + 1;
  } else {
    jl = js;  ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 1;  ku = ke + 1;
  } else {
    kl = ks;  ku = ke;
  }

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        cc_pos(pG, i, j, k, &x1, &x2, &x3);

        d = pG->U[k][j][i].d;

        T = pG->U[k][j][i].E -
            (0.5/d) * (SQR(pG->U[k][j][i].M1)
                       +SQR(pG->U[k][j][i].M2)
                       +SQR(pG->U[k][j][i].M3));

#ifdef MHD
        T -= (0.5)*(SQR(pG->U[k][j][i].B1c)
                    +SQR(pG->U[k][j][i].B2c)
                    +SQR(pG->U[k][j][i].B3c));
#endif
        T *= (Gamma_1/d);

        kappa = 0.0;
        if (KappaFun_i != NULL)
          kappa += (*KappaFun_i)(d, T, x1, x2, x3);
        if (KappaFun_a != NULL)
          kappa += (*KappaFun_a)(d, T, x1, x2, x3);

        kappa /= d;             /* convert to diffusion coefficient
                                   TODO: factor of (3/2)? */
        kappa_max = MAX(kappa_max, 2.0 * kappa);

      }
    }
  }

  return kappa_max;
}
#endif  /* THERMAL_CONDUCTION */


#ifdef VISCOSITY
Real max_nu(DomainS *pD)
{
  GridS *pG = (pD->Grid);

  int i, il, iu, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;

  Real d, T;
  Real x1, x2, x3;

  Real nu, nu_max = (TINY_NUMBER); /* NOT zero */

  il = is - 1;  iu = ie + 1;
  if (pG->Nx[1] > 1){
    jl = js - 1;  ju = je + 1;
  } else {
    jl = js;  ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 1;  ku = ke + 1;
  } else {
    kl = ks;  ku = ke;
  }

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        cc_pos(pG, i, j, k, &x1, &x2, &x3);

        d = pG->U[k][j][i].d;

        T = pG->U[k][j][i].E -
            (0.5/d) * (SQR(pG->U[k][j][i].M1)
                       +SQR(pG->U[k][j][i].M2)
                       +SQR(pG->U[k][j][i].M3));

#ifdef MHD
        T -= (0.5)*(SQR(pG->U[k][j][i].B1c)
                    +SQR(pG->U[k][j][i].B2c)
                    +SQR(pG->U[k][j][i].B3c));
#endif
        T *= (Gamma_1/d);

        nu = 0.0;
        if (NuFun_i != NULL)
          nu += (*NuFun_i)(d, T, x1, x2, x3);
        if (NuFun_a != NULL)
          nu += (*NuFun_a)(d, T, x1, x2, x3);

        nu_max = MAX(nu_max, 2.0 * nu);

      }
    }
  }

  return nu_max;
}
#endif  /* VISCOSITY */
