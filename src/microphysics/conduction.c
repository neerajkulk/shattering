#include "../copyright.h"
/*============================================================================*/
/*! \file conduction.c
 *  \brief Adds explicit thermal conduction term to the energy equation,
 *      dE/dt = Div(Q)
 *
 *   where
 *    - Q = kappa_iso Grad(T) + kappa_aniso([b Dot Grad(T)]b) = heat flux
 *    - T = (P/d)*(mbar/k_B) = temperature
 *    - b = magnetic field unit vector
 *
 *   Here
 *    - kappa_iso   is the   isotropic coefficient of thermal diffusion
 *    - kappa_aniso is the anisotropic coefficient of thermal diffusion
 *
 * Note the kappa's are DIFFUSIVITIES, not CONDUCTIVITIES.  Also note this
 * version uses "dimensionless units" in that the factor (mbar/k_B) is not
 * included in calculating the temperature (instead, T=P/d is used).  For cgs
 * units, kappa must be entered in units of [cm^2/s], and the heat fluxes would
 * need to be multiplied by (k_B/mbar).
 *
 * The heat flux Q is calculated by calls to HeatFlux_* functions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - conduction() - updates energy equation with thermal conduction
 * - conduction_init() - allocates memory needed
 * - conduction_destruct() - frees memory used */
/*============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef THERMAL_CONDUCTION

#ifdef BAROTROPIC
#error : Thermal conduction requires an adiabatic EOS
#endif

/* Arrays for the temperature and heat fluxes */
static Real ***Temp=NULL;
static Real ***Press=NULL;      /* for finding shocks */
static Real3Vect ***Q=NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   HeatFlux_iso   - computes   isotropic heat flux
 *   HeatFlux_aniso - computes anisotropic heat flux
 *============================================================================*/

void HeatFlux_iso(DomainS *pD);
void HeatFlux_aniso(DomainS *pD);

static Real limiter2(const Real A, const Real B);
static Real limiter4(const Real A, const Real B, const Real C, const Real D);
static Real vanleer (const Real A, const Real B);
static Real minmod  (const Real A, const Real B);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void conduction(DomainS *pD)
 *  \brief Explicit thermal conduction
 */
void conduction(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;
#ifdef STS
  Real my_dt = STS_dt;
#else
  Real my_dt = pG->dt / ncycle;
#endif
  Real dtodx1=my_dt/pG->dx1, dtodx2=0.0, dtodx3=0.0;

  if (KappaFun_a == NULL && KappaFun_i == NULL)
    return;

  if (pG->Nx[1] > 1){
    jl = js - 1;
    ju = je + 1;
    dtodx2 = my_dt/pG->dx2;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 1;
    ku = ke + 1;
    dtodx3 = my_dt/pG->dx3;
  } else {
    kl = ks;
    ku = ke;
  }

/* Zero heat heat flux array; compute temperature at cell centers.  Temperature
 * includes a factor [k_B/mbar].  For cgs units, the heat flux would have to
 * be multiplied by this factor.
 */

  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
  for (i=is-1; i<=ie+1; i++) {

    Q[k][j][i].x1 = 0.0;
    Q[k][j][i].x2 = 0.0;
    Q[k][j][i].x3 = 0.0;

    Temp[k][j][i] = pG->U[k][j][i].E - (0.5/pG->U[k][j][i].d)*
      (SQR(pG->U[k][j][i].M1) +SQR(pG->U[k][j][i].M2) +SQR(pG->U[k][j][i].M3));
#ifdef MHD
    Temp[k][j][i] -= (0.5)*(SQR(pG->U[k][j][i].B1c) +
      SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
    Temp[k][j][i] *= (Gamma_1/pG->U[k][j][i].d);

    /* also compute the pressure so we can detect shocks.  (at a
       contact discontinuity, temp can be discontinuous, but pressure
       isn't.  we want to turn off conduction at shocks, but not
       contact discontinuities.) */
    Press[k][j][i] = Temp[k][j][i] * pG->U[k][j][i].d;

  }}}

/* Compute isotropic and anisotropic heat fluxes.  Heat fluxes and temperature
 * are global variables in this file. */

  if (KappaFun_i != NULL)   HeatFlux_iso(pD);
  if (KappaFun_a != NULL) HeatFlux_aniso(pD);


/* Update energy using x1-fluxes */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[k][j][i].E += dtodx1*(Q[k][j][i+1].x1 - Q[k][j][i].x1);
    }
  }}

/* Update energy using x2-fluxes */

  if (pG->Nx[1] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx2*(Q[k][j+1][i].x2 - Q[k][j][i].x2);
      }
    }}
  }

/* Update energy using x3-fluxes */

  if (pG->Nx[2] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx3*(Q[k+1][j][i].x3 - Q[k][j][i].x3);
      }
    }}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void HeatFlux_iso(DomainS *pD)
 *  \brief Calculate heat fluxes with isotropic conduction
 */


void HeatFlux_iso(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie, iu;
  int j, js = pG->js, je = pG->je, ju;
  int k, ks = pG->ks, ke = pG->ke, ku;

  Real x1, x2, x3, kappa;
  Real dTdx, dTdy, dTdz;
  Real dPdx, dPdy, dPdz;

  const Real shock_thresh = 0.15;

  iu = ie+1;
  ju = (pG->Nx[1] > 1) ? je+1 : je;
  ku = (pG->Nx[2] > 1) ? ke+1 : ke;

  for (k=ks; k<=ku; k++) {
    for (j=js; j<=ju; j++) {
      for (i=is; i<=iu; i++) {
        /* Calculate the temperature gradient */
        dTdx = dTdy = dTdz = 0.0;
                           dTdx = (Temp[k][j][i]-Temp[k  ][j  ][i-1])/pG->dx1;
        if (pG->Nx[1] > 1) dTdy = (Temp[k][j][i]-Temp[k  ][j-1][i  ])/pG->dx2;
        if (pG->Nx[2] > 1) dTdz = (Temp[k][j][i]-Temp[k-1][j  ][i  ])/pG->dx3;

        /* Calculate kappa, but zero it at shocks */
        cc_pos(pG, i, j, k, &x1, &x2, &x3);
        kappa = (*KappaFun_i)(pG->U[k][j][i].d, Temp[k][j][i], x1, x2, x3);

        /* temp should be continuous... */
        if ((fabs(dTdx) > shock_thresh * Temp[k][j][i]/pG->dx1) ||
            (fabs(dTdy) > shock_thresh * Temp[k][j][i]/pG->dx2) ||
            (fabs(dTdz) > shock_thresh * Temp[k][j][i]/pG->dx3)) {

          /* Calculate the pressure gradient */
          dPdx = dPdy = dPdz = 0.0;
          dPdx = (Press[k][j][i]-Press[k  ][j  ][i-1])/pG->dx1;
          if (pG->Nx[1] > 1) dPdy = (Press[k][j][i]-Press[k  ][j-1][i  ])/pG->dx2;
          if (pG->Nx[2] > 1) dPdz = (Press[k][j][i]-Press[k-1][j  ][i  ])/pG->dx3;

          /* check that pressure is continuous */
          if ((fabs(dPdx) > shock_thresh * Press[k][j][i]/pG->dx1) ||
              (fabs(dPdy) > shock_thresh * Press[k][j][i]/pG->dx2) ||
              (fabs(dPdz) > shock_thresh * Press[k][j][i]/pG->dx3)) {

            /* The potential near the center of the cluster is
               unresolved, so can look like a shock.  Adding a cut P
               < 100 separates out the center of the cluster pretty
               cleanly over the whole redshift range of interest. */
            if (Press[k][j][i] < 100.0)
              kappa = 0.0;
          }
        }

                           Q[k][j][i].x1 += kappa*dTdx;
        if (pG->Nx[1] > 1) Q[k][j][i].x2 += kappa*dTdy;
        if (pG->Nx[2] > 1) Q[k][j][i].x3 += kappa*dTdz;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void HeatFlux_aniso(DomainS *pD)
 *  \brief Calculate heat fluxes with anisotropic conduction
 */

void HeatFlux_aniso(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real Bx,By,Bz,B02,dTc,dTl,dTr,lim_slope,dTdx,dTdy,dTdz,bDotGradT;

  Real d, x1, x2, x3, kappau, kappad, kappa;

  /* added by mkmcc to detect the virial shock */
  Real dPdx, dPdy, dPdz;
  const Real shock_thresh = 0.15;
  int at_shock;


#ifdef MHD
  if (pD->Nx[1] == 1) return;  /* problem must be at least 2D */

/* Compute heat fluxes in 1-direction  --------------------------------------*/

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {

/* Use unprojected gradients to determine whether we're at a shock. */
      /* Calculate the temperature gradient */
      at_shock = 0;

      dTdx = dTdy = dTdz = 0.0;
                         dTdx = (Temp[k][j][i]-Temp[k  ][j  ][i-1])/pG->dx1;
                         dTdy = (Temp[k][j][i]-Temp[k  ][j-1][i  ])/pG->dx2;
      if (pG->Nx[2] > 1) dTdz = (Temp[k][j][i]-Temp[k-1][j  ][i  ])/pG->dx3;

      /* temp should be continuous... */
      if ((fabs(dTdx) > shock_thresh * Temp[k][j][i]/pG->dx1) ||
          (fabs(dTdy) > shock_thresh * Temp[k][j][i]/pG->dx2) ||
          (fabs(dTdz) > shock_thresh * Temp[k][j][i]/pG->dx3)) {

        /* Calculate the pressure gradient */
        dPdx = dPdy = dPdz = 0.0;
                           dPdx = (Press[k][j][i]-Press[k  ][j  ][i-1])/pG->dx1;
                           dPdy = (Press[k][j][i]-Press[k  ][j-1][i  ])/pG->dx2;
        if (pG->Nx[2] > 1) dPdz = (Press[k][j][i]-Press[k-1][j  ][i  ])/pG->dx3;

        /* check that pressure is continuous */
        if ((fabs(dPdx) > shock_thresh * Press[k][j][i]/pG->dx1) ||
            (fabs(dPdy) > shock_thresh * Press[k][j][i]/pG->dx2) ||
            (fabs(dPdz) > shock_thresh * Press[k][j][i]/pG->dx3)) {

          /* The potential near the center of the cluster is
             unresolved, so can look like a shock.  Adding a cut P <
             100 separates out the center of the cluster pretty
             cleanly over the whole redshift range of interest. */
          if (Press[k][j][i] < 100.0)
           at_shock = 1;
        }
      }


/* Monotonized temperature difference dT/dy */
      dTr = 0.5*((Temp[k][j+1][i-1] + Temp[k][j+1][i]) -
                 (Temp[k][j  ][i-1] + Temp[k][j  ][i]));
      dTl = 0.5*((Temp[k][j  ][i-1] + Temp[k][j  ][i]) -
                 (Temp[k][j-1][i-1] + Temp[k][j-1][i]));
      dTc = dTr + dTl;

      dTdy = 0.0;
      if (dTl*dTr > 0.0) {
        lim_slope = MIN(fabs(dTl),fabs(dTr));
        dTdy = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized temperature difference dT/dz, 3D problem ONLY */
      if (pD->Nx[2] > 1) {
        dTr = 0.5*((Temp[k+1][j][i-1] + Temp[k+1][j][i]) -
                   (Temp[k  ][j][i-1] + Temp[k  ][j][i]));
        dTl = 0.5*((Temp[k  ][j][i-1] + Temp[k  ][j][i]) -
                   (Temp[k-1][j][i-1] + Temp[k-1][j][i]));
        dTc = dTr + dTl;

        dTdz = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdz = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx3;
        }
      }


      /* Interpolate the conductivity to the cell interface using
         harmonic averaging */
      d = pG->U[k][j][i].d;

      cc_pos(pG, i-1, j, k, &x1, &x2, &x3);
      kappad = (*KappaFun_a)(pG->U[k][j][i-1].d, Temp[k][j][i-1], x1, x2, x3);
      cc_pos(pG, i, j, k, &x1, &x2, &x3);
      kappau = (*KappaFun_a)(d, Temp[k][j][i], x1, x2, x3);

      kappa = 2.0 * (kappad * kappau) / (kappad + kappau + TINY_NUMBER);
      if (at_shock) kappa = 0.0;

/* Add flux at x1-interface, 2D PROBLEM */

      if (pD->Nx[2] == 1) {
        By = 0.5*(pG->U[k][j][i-1].B2c + pG->U[k][j][i].B2c);
        B02 = SQR(pG->B1i[k][j][i]) + SQR(By);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = pG->B1i[k][j][i]*(Temp[k][j][i]-Temp[k][j][i-1])/pG->dx1
           + By*dTdy;
        Q[k][j][i].x1 += kappa*(pG->B1i[k][j][i]*bDotGradT)/B02;

/* Add flux at x1-interface, 3D PROBLEM */

      } else {
        By = 0.5*(pG->U[k][j][i-1].B2c + pG->U[k][j][i].B2c);
        Bz = 0.5*(pG->U[k][j][i-1].B3c + pG->U[k][j][i].B3c);
        B02 = SQR(pG->B1i[k][j][i]) + SQR(By) + SQR(Bz);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = pG->B1i[k][j][i]*(Temp[k][j][i]-Temp[k][j][i-1])/pG->dx1
           + By*dTdy + Bz*dTdz;
        Q[k][j][i].x1 += kappa*(pG->B1i[k][j][i]*bDotGradT)/B02;
      }
    }
  }}

/* Compute heat fluxes in 2-direction  --------------------------------------*/

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {

/* Use unprojected gradients to determine whether we're at a shock. */
      /* Calculate the temperature gradient */
      at_shock = 0;

      dTdx = dTdy = dTdz = 0.0;
                         dTdx = (Temp[k][j][i]-Temp[k  ][j  ][i-1])/pG->dx1;
                         dTdy = (Temp[k][j][i]-Temp[k  ][j-1][i  ])/pG->dx2;
      if (pG->Nx[2] > 1) dTdz = (Temp[k][j][i]-Temp[k-1][j  ][i  ])/pG->dx3;

      /* temp should be continuous... */
      if ((fabs(dTdx) > shock_thresh * Temp[k][j][i]/pG->dx1) ||
          (fabs(dTdy) > shock_thresh * Temp[k][j][i]/pG->dx2) ||
          (fabs(dTdz) > shock_thresh * Temp[k][j][i]/pG->dx3)) {

        /* Calculate the pressure gradient */
        dPdx = dPdy = dPdz = 0.0;
                           dPdx = (Press[k][j][i]-Press[k  ][j  ][i-1])/pG->dx1;
                           dPdy = (Press[k][j][i]-Press[k  ][j-1][i  ])/pG->dx2;
        if (pG->Nx[2] > 1) dPdz = (Press[k][j][i]-Press[k-1][j  ][i  ])/pG->dx3;

        /* check that pressure is continuous */
        if ((fabs(dPdx) > shock_thresh * Press[k][j][i]/pG->dx1) ||
            (fabs(dPdy) > shock_thresh * Press[k][j][i]/pG->dx2) ||
            (fabs(dPdz) > shock_thresh * Press[k][j][i]/pG->dx3)) {

          /* The potential near the center of the cluster is
             unresolved, so can look like a shock.  Adding a cut P <
             100 separates out the center of the cluster pretty
             cleanly over the whole redshift range of interest. */
          if (Press[k][j][i] < 100.0)
            at_shock = 1;
        }
      }


/* Monotonized temperature difference dT/dx */
      dTr = 0.5*((Temp[k][j-1][i+1] + Temp[k][j][i+1]) -
                 (Temp[k][j-1][i  ] + Temp[k][j][i  ]));
      dTl = 0.5*((Temp[k][j-1][i  ] + Temp[k][j][i  ]) -
                 (Temp[k][j-1][i-1] + Temp[k][j][i-1]));
      dTc = dTr + dTl;

      dTdx = 0.0;
      if (dTl*dTr > 0.0) {
        lim_slope = MIN(fabs(dTl),fabs(dTr));
        dTdx = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized temperature difference dT/dz, 3D problem ONLY */
      if (pD->Nx[2] > 1) {
        dTr = 0.5*((Temp[k+1][j-1][i] + Temp[k+1][j][i]) -
                   (Temp[k  ][j-1][i] + Temp[k  ][j][i]));
        dTl = 0.5*((Temp[k  ][j-1][i] + Temp[k  ][j][i]) -
                   (Temp[k-1][j-1][i] + Temp[k-1][j][i]));
        dTc = dTr + dTl;

        dTdz = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdz = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx3;
        }
      }

      /* Interpolate the conductivity to the cell interface using
         harmonic averaging */
      d = pG->U[k][j][i].d;

      cc_pos(pG, i, j-1, k, &x1, &x2, &x3);
      kappad = (*KappaFun_a)(pG->U[k][j-1][i].d, Temp[k][j-1][i], x1, x2, x3);
      cc_pos(pG, i, j, k, &x1, &x2, &x3);
      kappau = (*KappaFun_a)(d, Temp[k][j][i], x1, x2, x3);

      kappa = 2.0 * (kappad * kappau) / (kappad + kappau + TINY_NUMBER);
      if (at_shock) kappa = 0.0;

/* Add flux at x2-interface, 2D PROBLEM */

      if (pD->Nx[2] == 1) {
        Bx = 0.5*(pG->U[k][j-1][i].B1c + pG->U[k][j][i].B1c);
        B02 = SQR(Bx) + SQR(pG->B2i[k][j][i]);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

        bDotGradT = pG->B2i[k][j][i]*(Temp[k][j][i]-Temp[k][j-1][i])/pG->dx2
           + Bx*dTdx;
        Q[k][j][i].x2 += kappa*(pG->B2i[k][j][i]*bDotGradT)/B02;

/* Add flux at x2-interface, 3D PROBLEM */

      } else {
        Bx = 0.5*(pG->U[k][j-1][i].B1c + pG->U[k][j][i].B1c);
        Bz = 0.5*(pG->U[k][j-1][i].B3c + pG->U[k][j][i].B3c);
        B02 = SQR(Bx) + SQR(pG->B2i[k][j][i]) + SQR(Bz);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = pG->B2i[k][j][i]*(Temp[k][j][i]-Temp[k][j-1][i])/pG->dx2
           + Bx*dTdx + Bz*dTdz;
        Q[k][j][i].x2 += kappa*(pG->B2i[k][j][i]*bDotGradT)/B02;
      }
    }
  }}

/* Compute heat fluxes in 3-direction, 3D problem ONLY  ---------------------*/

  if (pD->Nx[2] > 1) {
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

/* Use unprojected gradients to determine whether we're at a shock. */
        /* Calculate the temperature gradient */
        at_shock = 0;

        dTdx = dTdy = dTdz = 0.0;
                           dTdx = (Temp[k][j][i]-Temp[k  ][j  ][i-1])/pG->dx1;
                           dTdy = (Temp[k][j][i]-Temp[k  ][j-1][i  ])/pG->dx2;
        if (pG->Nx[2] > 1) dTdz = (Temp[k][j][i]-Temp[k-1][j  ][i  ])/pG->dx3;

        /* temp should be continuous... */
        if ((fabs(dTdx) > shock_thresh * Temp[k][j][i]/pG->dx1) ||
            (fabs(dTdy) > shock_thresh * Temp[k][j][i]/pG->dx2) ||
            (fabs(dTdz) > shock_thresh * Temp[k][j][i]/pG->dx3)) {

          /* Calculate the pressure gradient */
          dPdx = dPdy = dPdz = 0.0;
                             dPdx = (Press[k][j][i]-Press[k  ][j  ][i-1])/pG->dx1;
                             dPdy = (Press[k][j][i]-Press[k  ][j-1][i  ])/pG->dx2;
          if (pG->Nx[2] > 1) dPdz = (Press[k][j][i]-Press[k-1][j  ][i  ])/pG->dx3;

          /* check that pressure is continuous */
          if ((fabs(dPdx) > shock_thresh * Press[k][j][i]/pG->dx1) ||
              (fabs(dPdy) > shock_thresh * Press[k][j][i]/pG->dx2) ||
              (fabs(dPdz) > shock_thresh * Press[k][j][i]/pG->dx3)) {

            /* The potential near the center of the cluster is
               unresolved, so can look like a shock.  Adding a cut P <
               100 separates out the center of the cluster pretty
               cleanly over the whole redshift range of interest. */
            if (Press[k][j][i] < 100.0)
             at_shock = 1;
          }
        }


/* Monotonized temperature difference dT/dx */
        dTr = 0.5*((Temp[k-1][j][i+1] + Temp[k][j][i+1]) -
                   (Temp[k-1][j][i  ] + Temp[k][j][i  ]));
        dTl = 0.5*((Temp[k-1][j][i  ] + Temp[k][j][i  ]) -
                   (Temp[k-1][j][i-1] + Temp[k][j][i-1]));
        dTc = dTr + dTl;

        dTdx = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdx = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx1;
        }

/* Monotonized temperature difference dT/dy */
        dTr = 0.5*((Temp[k-1][j+1][i] + Temp[k][j+1][i]) -
                   (Temp[k-1][j  ][i] + Temp[k][j  ][i]));
        dTl = 0.5*((Temp[k-1][j  ][i] + Temp[k][j  ][i]) -
                   (Temp[k-1][j-1][i] + Temp[k][j-1][i]));
        dTc = dTr + dTl;

        dTdy = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdy = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx2;
        }

        /* Interpolate the conductivity to the cell interface using
           harmonic averaging */
         d = pG->U[k][j][i].d;

        cc_pos(pG, i, j, k-1, &x1, &x2, &x3);
        kappad = (*KappaFun_a)(pG->U[k-1][j][i].d, Temp[k-1][j][i], x1, x2, x3);
        cc_pos(pG, i, j, k, &x1, &x2, &x3);
        kappau = (*KappaFun_a)(d, Temp[k][j][i], x1, x2, x3);

        kappa = 2.0 * (kappad * kappau) / (kappad + kappau + TINY_NUMBER);
        if (at_shock) kappa = 0.0;

/* Add flux at x3-interface, 3D PROBLEM */

        Bx = 0.5*(pG->U[k-1][j][i].B1c + pG->U[k][j][i].B1c);
        By = 0.5*(pG->U[k-1][j][i].B2c + pG->U[k][j][i].B2c);
        B02 = SQR(Bx) + SQR(By) + SQR(pG->B3i[k][j][i]);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = pG->B3i[k][j][i]*(Temp[k][j][i]-Temp[k-1][j][i])/pG->dx3
           + Bx*dTdx + By*dTdy;
        Q[k][j][i].x3 += kappa*(pG->B3i[k][j][i]*bDotGradT)/B02;
      }
    }}
  }
#endif /* MHD */

  return;
}


/*----------------------------------------------------------------------------*/
/* limiter2 and limiter4: call slope limiters to preserve monotonicity
 */

static Real limiter2(const Real A, const Real B)
{
  /* van Leer slope limiter */
  return vanleer(A,B);

  /* monotonized central (MC) limiter */
  /* return minmod(2.0*minmod(A,B),0.5*(A+B)); */
}

static Real limiter4(const Real A, const Real B, const Real C, const Real D)
{
  return limiter2(limiter2(A,B),limiter2(C,D));
}

/*----------------------------------------------------------------------------*/
/* vanleer: van Leer slope limiter
 */

static Real vanleer(const Real A, const Real B)
{
  if (A*B > 0) {
    return 2.0*A*B/(A+B);
  } else {
    return 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/* minmod: minmod slope limiter
 */

static Real minmod(const Real A, const Real B)
{
  if (A*B > 0) {
    if (A > 0) {
      return MIN(A,B);
    } else {
      return MAX(A,B);
    }
  } else {
    return 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void conduction_init(MeshS *pM)
 *  \brief Allocate temporary arrays
 */

void conduction_init(MeshS *pM)
{
  int nl,nd,size1=1,size2=1,size3=1,Nx1,Nx2,Nx3;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  Nx1 = size1 + 2*nghost;

  if (pM->Nx[1] > 1){
    Nx2 = size2 + 2*nghost;
  } else {
    Nx2 = size2;
  }

  if (pM->Nx[2] > 1){
    Nx3 = size3 + 2*nghost;
  } else {
    Nx3 = size3;
  }
  if ((Temp = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((Press = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((Q = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
    goto on_error;
  return;

  on_error:
  conduction_destruct();
  ath_error("[conduct_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void conduction_destruct(void)
 *  \brief Free temporary arrays
 */

void conduction_destruct(void)
{
  if (Temp != NULL) free_3d_array(Temp);
  if (Press != NULL) free_3d_array(Press);
  if (Q != NULL) free_3d_array(Q);
  return;
}
#endif /* THERMAL_CONDUCTION */
