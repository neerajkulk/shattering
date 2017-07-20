#include "copyright.h"
/*==============================================================================
 * FILE: streaming2d.c
 *
 * PURPOSE: Problem generator for streaming instability test in non-stratified
 *   disks. This code works in 2D, with single particle species ONLY. Isothermal
 *   eos is assumed, and the value etavk/iso_sound is fixed. In the case of
 *   growth rate test, since one requires an integer number of wavelengths
 *   across the grid, the isothermal sound speed is re-evaluated.
 *
 * Perturbation modes:
 *    ipert = 0: no perturbation, test for the NSH equilibrium, need FARGO
 *    ipert = 1: linA of YJ07 (cold start), need FARGO
 *    ipert = 2: linB of YJ07 (cold start), need FARGO
 *    ipert = 3: random perturbation (warm start), do not need FARGO
 *    ipert = 11-13: three new modes for small particles used in Bai & Stone
 *                   in their appendix (need FARGO).
 *  Must be configured using --enable-shearing-box and --with-eos=isothermal.
 *  FARGO is need to establish the NSH equilibrium (ipert=0,1,2).
 *
 * Reference:
 *   Youdin & Johansen, 2007, ApJ, 662, 613
 *   Johansen & Youdin, 2007, ApJ, 662, 627
 *   Bai & Stone,       2010, ApJS, 190, 297
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

#ifndef SHEARING_BOX
#error : The streaming2d problem requires shearing-box to be enabled.
#endif /* SHEARING_BOX */

#ifndef PARTICLES
#error : The streaming2d problem requires particles to be enabled.
#endif /* PARTICLES */

#ifndef ISOTHERMAL
#error : The streaming2d problem requires isothermal equation of state.
#endif /* ISOTHERMAL */

/*------------------------ filewide global variables -------------------------*/
/* NSH equilibrium parameters */
Real rho0, mratio, etavk, uxNSH, uyNSH, wxNSH, wyNSH;
/* particle number related variables */
int Npar,Npar2,downsamp,Nx;
/* eigen vector for a streaming instability mode */
Real Reux,Imux,Reuy,Imuy,Reuz,Imuz,Rewx,Imwx,Rewy,Imwy,Rewz,Imwz,Rerho,Imrho,omg,s;
/* perturbation variables */
Real amp, kx, kz;
int ipert, nwave;
/* output filename */
char name[50];

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()            - random number generator
 * OutputModeAmplitude() - output the perturbation amplitude
 * ShearingBoxPot()  - shearing box tidal gravitational potential
 * pert_???()        - perturbation wave form for linear growth rate test
 * property_mybin()  - particle property selection function 
 * property_???()    - particle property selection function
 *============================================================================*/
double ran2(long int *idum);
void OutputModeAmplitude(MeshS *pM, OutputS *pOut);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real pert_even(Real fR, Real fI, Real x, Real z, Real t);
static Real pert_odd(Real fR, Real fI, Real x, Real z, Real t);
static int property_mybin(const GrainS *gr, const GrainAux *grsub);
extern Real expr_V1par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V2par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V3par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V2(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,ip,jp,interp;
  long p;
  Real x1,x2,x3,t,x1l,x1u,x2l,x2u,x1p,x2p,x3p,paramp,factor2,reduct;
  Real x1min,x1max,x2min,x2max,Lx,Lz;
  Real rhog,cs,u1,u2,u3,w1,w2,w3,Kxn,Kzn,denorm1;
  long int iseed = -1; /* Initialize on the first call to ran2 */

  if (pDomain->Nx[1] == 1) {
    ath_error("[streaming2d]: streaming2D only works for 2D problem.\n");
  }
  if (pDomain->Nx[2] > 1) {
    ath_error("[streaming2d]: streaming2D only works for 2D problem.\n");
  }

/* Initialize boxsize */
  x1min = pDomain->RootMinX[0];
  x1max = pDomain->RootMaxX[0];
  Lx = x1max - x1min;

  x2min = pDomain->RootMinX[1];
  x2max = pDomain->RootMaxX[1];
  Lz = x2max - x2min;

  Nx = pDomain->Nx[0]; /* used for particle selection function */

/* Read initial conditions */
  rho0 = 1.0;
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  amp = par_getd_def("problem","amp",0.0);
  ipert = par_geti_def("problem","ipert", 1);
  etavk = par_getd_def("problem","etavk",0.05);/* in unit of iso_sound (N.B.) */

  /* particle number */
  if (npartypes != 1)
    ath_error("[streaming2d]: This test only allows ONE particle species!\n");

  Npar  = (int)(sqrt(par_geti("particle","parnumcell")));
  Npar2 = SQR(Npar);

  pGrid->nparticle         = Npar2*pGrid->Nx[0]*pGrid->Nx[1];
  grproperty[0].num = pGrid->nparticle;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  /* set the down sampling of the particle list output
   * (by default, output 1 particle per cell)
   */
  downsamp = par_geti_def("problem","downsamp",Npar2);

  /* particle stopping time */
  tstop0[0] = par_getd("problem","tstop"); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[streaming2d]: This test works only for fixed stopping time!\n");

  /* assign particle effective mass */
#ifdef FEEDBACK
  mratio = par_getd_def("problem","mratio",0.0); /* total mass fraction */
  if (mratio < 0.0)
    ath_error("[streaming2D]: mratio must be positive!\n");
#else
  mratio = 0.0;
  if ((ipert == 1) || (ipert == 2))
    ath_error("[streaming2d]: linear growth test requires FEEDBACK to be\
 turned on!\n");
#endif
  amp = amp*mratio; /* N.B.! */

  if (ipert == 1) {
    Kxn   =  30.0;
    Kzn   =  30.0;
    Reux  = -0.1691398*amp;
    Imux  =  0.0361553*amp;
    Reuy  =  0.1336704*amp;
    Imuy  =  0.0591695*amp;
    Reuz  =  0.1691389*amp;
    Imuz  = -0.0361555*amp;
    Rerho =  0.0000224*amp;
    Imrho =  0.0000212*amp;
    Rewx  = -0.1398623*amp;
    Imwx  =  0.0372951*amp;
    Rewy  =  0.1305628*amp;
    Imwy  =  0.0640574*amp;
    Rewz  =  0.1639549*amp;
    Imwz  = -0.0233277*amp;
    omg   = -0.3480127*Omega_0;
    s     =  0.4190204*Omega_0;
    mratio=  3.0;
    tstop0[0]=  0.1/Omega_0;
    etavk = 0.05;
  }

  if (ipert == 2) {
    Kxn   =  6.0;
    Kzn   =  6.0;
    Reux  = -0.0174121*amp; 
    Imux  = -0.2770347*amp; 
    Reuy  =  0.2767976*amp;
    Imuy  = -0.0187568*amp; 
    Reuz  =  0.0174130*amp;
    Imuz  =  0.2770423*amp;
    Rerho = -0.0000067*amp;
    Imrho = -0.0000691*amp;
    Rewx  =  0.0462916*amp;
    Imwx  = -0.2743072*amp;
    Rewy  =  0.2739304*amp;
    Imwy  =  0.0039293*amp;
    Rewz  =  0.0083263*amp;
    Imwz  =  0.2768866*amp;
    omg   =  0.4998786*Omega_0;
    s     =  0.0154764*Omega_0;
    mratio=  0.2;
    tstop0[0]=  0.1/Omega_0;
    etavk = 0.05;
  }

 if (ipert == 11) {
    Kxn   =  1500.0;
    Kzn   =  1500.0;
    Reux  = -0.15987506226181683*amp;
    Imux  =  0.007966909776692505*amp;
    Reuy  =  0.11644231768989997*amp;
    Imuy  =  0.012237732175520223*amp;
    Reuz  =  0.1598750616685892*amp;
    Imuz  = -0.007966912082409646*amp;
    Rerho =  8.684871615767337e-8*amp;
    Imrho =  5.350037142664432e-7*amp;
    Rewx  = -0.15671736667638894*amp;
    Imwx  =  0.002883666233000987*amp;
    Rewy  =  0.11597822199096137*amp;
    Imwy  =  0.01611452420454858*amp;
    Rewz  =  0.15900951328796148*amp;
    Imwz  = -0.002484953589870477*amp;
    omg   =  0.1049236206315842*Omega_0;
    s     =  0.5980689646967706*Omega_0;
    mratio=  2.0;
    tstop0[0]=  0.01/Omega_0;
    etavk = 0.05;
  }

  if (ipert == 12) {
    Kxn   =  750.0;
    Kzn   =  750.0;
    Reux  = -0.05572772398711219*amp;
    Imux  = -0.026949633413065097*amp;
    Reuy  =  0.24569460629315637*amp;
    Imuy  = -0.01567907955658157*amp;
    Reuz  =  0.05572772420698503*amp;
    Imuz  =  0.026949633268483*amp;
    Rerho = -4.355385656601321e-8*amp;
    Imrho =  2.8002601621265505e-8*amp;
    Rewx  = -0.0497205491056745*amp;
    Imwx  = -0.028901299502860814*amp;
    Rewy  =  0.24608486243618063*amp;
    Imwy  = -0.006455493453903131*amp;
    Rewz  =  0.054638391601795934*amp;
    Imwz  =  0.028953557818336886*amp;
    omg   =  -0.0615243802498589*Omega_0;
    s     =  0.03919373660665014*Omega_0;
    mratio=  1.0;
    tstop0[0]=  0.01/Omega_0;
    etavk = 0.05;
  }

  if (ipert == 13) {
    Kxn   =  2000.0;
    Kzn   =  2000.0;
    Reux  = -0.1719650139257441*amp;
    Imux  =  0.07407124939002366*amp;
    Reuy  =  0.1918893205596143*amp;
    Imuy  =  0.07865187494765362*amp;
    Reuz  =  0.17196501382406704*amp;
    Imuz  = -0.07407124937574784*amp;
    Rerho =  2.9546312378863307e-7*amp;
    Imrho =  1.1413846778145968e-7*amp;
    Rewx  = -0.17158399112103911*amp;
    Imwx  =  0.07407376465545336*amp;
    Rewy  =  0.19185420849094573*amp;
    Imwy  =  0.07873714071990433*amp;
    Rewz  =  0.17196745751715414*amp;
    Imwz  = -0.07391604601712277*amp;
    omg   =  0.32248839717349675*Omega_0;
    s     =  0.3154372766709299*Omega_0;
    mratio=  2.0;
    tstop0[0]=  0.001/Omega_0;
    etavk = 0.05;
  }

#ifdef FEEDBACK
  grproperty[0].m = rho0*mratio/Npar2;
#endif

  /* Adjust code units */ 
  if ((ipert == 1) || (ipert == 2) || (ipert > 10))
  {
    if (Lx != Lz)
      ath_error("[streaming2d]: Linear test requires Lx1=Lx2!\n");
    nwave = par_geti_def("problem","nwave",1);
    kx    = 2.0*PI/Lx*nwave;
    kz    = 2.0*PI/Lz*nwave;

    /* Reset isothermal sound speed */
    cs = Kxn/kx/etavk*Omega_0;
    if (cs != Iso_csound)
      ath_pout(0, "[streaming2d]: Iso_csound=%f is modified to cs=%f.\n",
                                             Iso_csound, cs);
    Iso_csound = cs;
    Iso_csound2 = SQR(Iso_csound);

    interp = par_geti("particle","interp");
    if (interp == 3) {/* QP */
      paramp = amp*kx*pGrid->dx1/sin(kx*pGrid->dx1);
      factor2 = 0.5*tan(kx*pGrid->dx1)/(kx*pGrid->dx1);
      reduct = 1.0;
    }
    else if (interp == 2) {/* TSC */
        paramp = 4.0*amp*kx*pGrid->dx1/sin(kx*pGrid->dx1)/
                                 (3.0+cos(kz*pGrid->dx2));
        factor2 = 0.5*tan(kx*pGrid->dx1)/(kx*pGrid->dx1);
        reduct = 1.0/(1.0-0.25*SQR(kx*pGrid->dx1)); reduct=1.0;
    }
    else {
        paramp = amp;
        factor2 = 0.5;
        reduct = 1.0;
    }
  }
  etavk = etavk * Iso_csound; /* switch to code unit (N.B.!) */

  /* calculate NSH equilibrium velocity */
  denorm1 = 1.0/(SQR(1.0+mratio)+SQR(tstop0[0]*Omega_0));

  wxNSH = -2.0*tstop0[0]*Omega_0*denorm1*etavk;
  wyNSH = -(1.0+mratio)*denorm1*etavk;

  uxNSH = -mratio*wxNSH;
  uyNSH = -mratio*wyNSH;

  wyNSH += etavk;

  ath_pout(0,"etavk=%f, Iso_csound=%f\n",etavk,Iso_csound);

/* Now set initial conditions for the gas */
  t = 0.0;

  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    cc_pos(pGrid,i,j,pGrid->ks,&x1,&x2,&x3);

    if ((ipert == 1) || (ipert == 2) || (ipert > 10)) {
      rhog = rho0 * (1.0+reduct*pert_even(Rerho,Imrho,x1,x2,t));
      u1 = etavk * reduct * pert_even(Reux,Imux,x1,x2,t);
      u2 = etavk * reduct * pert_odd (Reuz,Imuz,x1,x2,t);
      u3 = etavk * reduct * pert_even(Reuy,Imuy,x1,x2,t);
    } else {
      rhog = rho0;  u1 = u2 = u3 = 0.0;
    }

    pGrid->U[pGrid->ks][j][i].d = rhog;

    pGrid->U[pGrid->ks][j][i].M1 = rhog * (uxNSH+u1);
    pGrid->U[pGrid->ks][j][i].M3 = rhog * (uyNSH+u3);

    pGrid->U[pGrid->ks][j][i].M2 = rhog * u2;
#ifndef FARGO
    pGrid->U[pGrid->ks][j][i].M3 -= qshear*rhog*Omega_0*x1;
#endif

  }}

/* Now set initial conditions for the particles */
  p = 0;
  x3p = 0.5*(pGrid->MinX[2]+pGrid->MaxX[2]);

  Lx = pGrid->Nx[0]*pGrid->dx1;
  Lz = pGrid->Nx[1]*pGrid->dx2;

  x1min = pGrid->MinX[0];
  x2min = pGrid->MinX[1]; 

  for (j=pGrid->js; j<=pGrid->je; j++)
  {
    x2l = pGrid->MinX[1] + (j-pGrid->js)*pGrid->dx2;
    x2u = pGrid->MinX[1] + (j-pGrid->js+1.0)*pGrid->dx2;

    for (i=pGrid->is; i<=pGrid->ie; i++)
    {
      x1l = pGrid->MinX[0] + (i-pGrid->is)*pGrid->dx1;
      x1u = pGrid->MinX[0] + (i-pGrid->is+1.0)*pGrid->dx1;

        for (ip=0;ip<Npar;ip++)
        {

          if (ipert != 3)
            x1p = x1l+pGrid->dx1/Npar*(ip+0.5);

          for (jp=0;jp<Npar;jp++)
          {
            if (ipert == 3)
            { /* ramdom particle position in the grid */
              x1p = x1min + Lx*ran2(&iseed);
              x2p = x2min + Lz*ran2(&iseed);
            }
            else
              x2p = x2l+pGrid->dx2/Npar*(jp+0.5);

            pGrid->particle[p].property = 0;
            pGrid->particle[p].x1 = x1p;
            pGrid->particle[p].x2 = x2p;
            pGrid->particle[p].x3 = x3p;

            if ((ipert == 1) || (ipert == 2) || (ipert > 10)) {
              pGrid->particle[p].x1 += paramp*cos(kz*x2p)*(-sin(kx*x1p)
                                      +factor2*paramp*sin(2.0*kx*x1p))/kx;
//              pGrid->particle[p].x1 += amp*cos(kz*x2p)*(-sin(kx*x1p)
//                                        +0.5*amp*sin(2.0*kx*x1p))/kx;
              w1 = etavk * pert_even(Rewx,Imwx,pGrid->particle[p].x1,x2p,t);
              w2 = etavk * pert_odd (Rewz,Imwz,pGrid->particle[p].x1,x2p,t);
              w3 = etavk * pert_even(Rewy,Imwy,pGrid->particle[p].x1,x2p,t);
            } else {
              w1 = w2 = w3 = 0.0;
            }

            pGrid->particle[p].v1 = wxNSH+w1;
            pGrid->particle[p].v3 = wyNSH+w3;

            pGrid->particle[p].v2 = w2;
#ifndef FARGO
            pGrid->particle[p].v3 -= qshear*Omega_0*x1p;
#endif
            pGrid->particle[p].pos = 1; /* grid particle */
            pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
            pGrid->particle[p].init_id = myID_Comm_world;
#endif
            p += 1;
          }
        }
    }
  }

/* enroll gravitational potential function, shearing sheet BC functions */
  ShearingBoxPot = UnstratifiedDisk;

  ShBoxCoord = xz;

  if (myID_Comm_world == 0) {
    /* flush output file */
    sprintf(name, "%s_%d_%d.dat", "Streaming2d",Nx,ipert);
    FILE *fid = fopen(name,"w");
    fclose(fid);
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_%d_%d.dat", "Streaming2d",Nx,ipert);
#endif
  }

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  fwrite(&rho0, sizeof(Real),1,fp);  fwrite(&mratio, sizeof(Real),1,fp);
  fwrite(&etavk, sizeof(Real),1,fp);
  fwrite(&uxNSH, sizeof(Real),1,fp); fwrite(&uyNSH, sizeof(Real),1,fp);
  fwrite(&wxNSH, sizeof(Real),1,fp); fwrite(&wyNSH, sizeof(Real),1,fp);
  if ((ipert==1) || (ipert==2) || (ipert > 10)) {
    fwrite(&Reux, sizeof(Real),1,fp); fwrite(&Imux, sizeof(Real),1,fp);
    fwrite(&Reuy, sizeof(Real),1,fp); fwrite(&Imuy, sizeof(Real),1,fp);
    fwrite(&Reuz, sizeof(Real),1,fp); fwrite(&Imuz, sizeof(Real),1,fp);
    fwrite(&Rewx, sizeof(Real),1,fp); fwrite(&Imwx, sizeof(Real),1,fp);
    fwrite(&Rewy, sizeof(Real),1,fp); fwrite(&Imwy, sizeof(Real),1,fp);
    fwrite(&Rewz, sizeof(Real),1,fp); fwrite(&Imwz, sizeof(Real),1,fp);
    fwrite(&Rerho, sizeof(Real),1,fp);fwrite(&Imrho, sizeof(Real),1,fp);
    fwrite(&omg, sizeof(Real),1,fp);  fwrite(&s, sizeof(Real),1,fp);
    fwrite(&Iso_csound, sizeof(Real),1,fp);
    fwrite(&kx, sizeof(Real),1,fp);   fwrite(&kz, sizeof(Real),1,fp);
  }
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]);
  GridS *pG = pD->Grid;

  ShearingBoxPot = UnstratifiedDisk;
  ShBoxCoord = xz;

  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  amp = par_getd_def("problem","amp",0.0);
  ipert = par_geti_def("problem","ipert", 1);

  Nx = pM->Nx[0];
  Npar  = (int)(sqrt(par_geti("particle","parnumcell")));
  Npar2 = SQR(Npar);
  downsamp = par_geti_def("problem","downsamp",Npar2);

  fread(&rho0, sizeof(Real),1,fp);  fread(&mratio, sizeof(Real),1,fp);
  fread(&etavk, sizeof(Real),1,fp);
  fread(&uxNSH, sizeof(Real),1,fp); fread(&uyNSH, sizeof(Real),1,fp);
  fread(&wxNSH, sizeof(Real),1,fp); fread(&wyNSH, sizeof(Real),1,fp);
  if ((ipert==1) || (ipert==2) || (ipert > 10)) {
    fread(&Reux, sizeof(Real),1,fp); fread(&Imux, sizeof(Real),1,fp);
    fread(&Reuy, sizeof(Real),1,fp); fread(&Imuy, sizeof(Real),1,fp);
    fread(&Reuz, sizeof(Real),1,fp); fread(&Imuz, sizeof(Real),1,fp);
    fread(&Rewx, sizeof(Real),1,fp); fread(&Imwx, sizeof(Real),1,fp);
    fread(&Rewy, sizeof(Real),1,fp); fread(&Imwy, sizeof(Real),1,fp);
    fread(&Rewz, sizeof(Real),1,fp); fread(&Imwz, sizeof(Real),1,fp);
    fread(&Rerho, sizeof(Real),1,fp);fread(&Imrho, sizeof(Real),1,fp);
    fread(&omg, sizeof(Real),1,fp);  fread(&s, sizeof(Real),1,fp);
    fread(&Iso_csound, sizeof(Real),1,fp);  Iso_csound2 = SQR(Iso_csound);
    fread(&kx, sizeof(Real),1,fp);   fread(&kz, sizeof(Real),1,fp);
  }

  if (myID_Comm_world == 0)
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_%d_%d.dat", "Streaming2d",Nx,ipert);
#else
    sprintf(name, "%s_%d_%d.dat", "Streaming2d",Nx,ipert);
#endif

  return;
}

/* difd */
static Real expr_rhodif(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->U[k][j][i].d - rho0;
}

/* dVx */
static Real expr_dVx(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->U[k][j][i].M1/pG->U[k][j][i].d - uxNSH;
}

/* dVy */
static Real expr_dVy(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M3/pG->U[k][j][i].d - uyNSH;
#else
  return (pG->U[k][j][i].M3/pG->U[k][j][i].d -uyNSH + qshear*Omega_0*x1);
#endif
}

/* difdpar */
static Real expr_rhopardif(const GridS *pG,
                           const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->Coup[k][j][i].grid_d - rho0*mratio;
}

/* dVxpar */
static Real expr_dVxpar(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return expr_V1par(pG,i,j,k)-wxNSH;
}

/* dVypar */
static Real expr_dVypar(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return expr_V3par(pG,i,j,k)-wyNSH;
#else
  return expr_V3par(pG,i,j,k)-wyNSH+qshear*Omega_0*x1;
#endif
}

ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"difd")==0) return expr_rhodif;
  if(strcmp(expr,"dVx")==0) return expr_dVx;
  if(strcmp(expr,"dVy")==0) return expr_dVy;
  if(strcmp(expr,"difdpar")==0) return expr_rhopardif;
  if(strcmp(expr,"dVxpar")==0) return expr_dVxpar;
  if(strcmp(expr,"dVypar")==0) return expr_dVypar;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  if(strcmp(name,"amp")==0) return OutputModeAmplitude;
  return NULL;
}

#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  if (strcmp(name,"downsamp")==0) return property_mybin;
  return NULL;
}

void gasvshift(const Real x1, const Real x2, const Real x3,
                                    Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2,
                                       const Real x3, Real w1, Real w2, Real w3)
{
  ft->x1 -= 2.0*etavk*Omega_0;
  return;
}
#endif

void Userwork_in_loop(MeshS *pM)
{
  return;
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*--------------------------------------------------------------------------- */
/* Output function */
void OutputModeAmplitude(MeshS *pM, OutputS *pOut)
{
  DomainS *pDomain = (DomainS*)&(pM->Domain[0][0]);
  GridS *pGrid = pDomain->Grid;
  FILE *fid;
  int i,j,ks=pGrid->ks;
  long p;
  GrainS *gr;
  Real dm,dparm,uxm,uym,uzm,wxm,wym,wzm;

  particle_to_grid(pDomain, property_all);

  dm=0.0; dparm=0.0; uxm=0.0; uym=0.0; uzm=0.0; wxm=0.0; wym=0.0; wzm=0.0;

  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    dm  = MAX(dm,fabs(expr_rhodif(pGrid,i,j,ks)));
    dparm=MAX(dparm,fabs(expr_rhopardif(pGrid,i,j,ks)));
    uxm = MAX(uxm,fabs(expr_dVx(pGrid,i,j,ks)));
    wxm = MAX(wxm,fabs(expr_dVxpar(pGrid,i,j,ks)));
    uym = MAX(uym,fabs(expr_dVy(pGrid,i,j,ks)));
    wym = MAX(wym,fabs(expr_dVypar(pGrid,i,j,ks)));
    uzm = MAX(uzm,fabs(expr_V2(pGrid,i,j,ks)));
    wzm = MAX(wzm,fabs(expr_V2par(pGrid,i,j,ks)));
  }}

#ifdef MPI_PARALLEL
  Real sendbuf[8], recvbuf[8];
  int err;
  sendbuf[0]=dm;        sendbuf[1]=dparm;
  sendbuf[2]=uxm;       sendbuf[3]=wxm;
  sendbuf[4]=uym;       sendbuf[5]=wym;
  sendbuf[6]=uzm;       sendbuf[7]=wzm;

  err = MPI_Reduce(sendbuf,recvbuf,8,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(err) ath_error("[streaming2d]: MPI_Reduce returned error code %d\n",err);

  if (myID_Comm_world == 0) {
    dm=recvbuf[0];      dparm=recvbuf[1];
    uxm=recvbuf[2];     wxm=recvbuf[3];
    uym=recvbuf[4];     wym=recvbuf[5];
    uzm=recvbuf[6];     wzm=recvbuf[7];
  }
#endif

  if (myID_Comm_world == 0) {
    fid = fopen(name,"a+");
    fprintf(fid,"%e   %e    %e    %e    %e    %e    %e    %e    %e\n",
               pGrid->time,dm,dparm,uxm,wxm,uym,wym,uzm,wzm);

    fclose(fid);
  }

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/* ShearingBoxPot */
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*SQR(Omega_0*x1);
#endif
  return phi;
}

/* even perturbation mode */
static Real pert_even(Real fR, Real fI, Real x, Real z, Real t)
{
  return (fR*cos(kx*x-omg*t)-fI*sin(kx*x-omg*t))*cos(kz*z)*exp(s*t);
}

/* odd perturbation mode */
static Real pert_odd(Real fR, Real fI, Real x, Real z, Real t)
{
  return -(fR*sin(kx*x-omg*t)+fI*cos(kx*x-omg*t))*sin(kz*z)*exp(s*t);
}

/* user defined particle selection function (1: true; 0: false) */
static int property_mybin(const GrainS *gr, const GrainAux *grsub)
{
  long a,b,c,d,e,ds,sp;

  sp = MAX(downsamp/Npar2,1);         /* spacing in cells */
  ds = Npar2*sp;               /* actual dowmsampling */

  a = gr->my_id/ds;
  b = gr->my_id - a*ds;

  c = gr->my_id/(Npar2*Nx);    /* column number */
  d = c/sp;
  e = c-sp*d;

  if ((e == 0) && (b == 0) && (gr->pos == 1))
    return 1;
  else
    return 0;
}

/*------------------------------------------------------------------------------
 * ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1.
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX
