#include "copyright.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#include "prob/math_functions.h"

#ifdef MPI_PARALLEL
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */

#define FOLLOW_CLOUD
#define REPORT_NANS
#define ENERGY_COOLING

static void bc_ix1(GridS *pGrid);
static void bc_ox1(GridS *pGrid);
static void bc_ix2(GridS *pGrid);
static void bc_ox2(GridS *pGrid);
static void bc_ix3(GridS *pGrid);
static void bc_ox3(GridS *pGrid);

#ifdef REPORT_NANS
static void report_nans(DomainS *pDomain);
#endif  /* REPORT_NANS */

#ifdef ENERGY_COOLING
static Real cooling_func(const Real rho, const Real P, const Real dt);
#endif  /* ENERGY_COOLING */

#ifdef FOLLOW_CLOUD
static Real cloud_mass_weighted_velocity(MeshS *pM);
static void boost_frame(DomainS *pDomain, Real dv);
static Real x_shift;

static Real hst_xshift(GridS *pG, int i, int j, int k);
static Real hst_vflow(GridS *pG, int i, int j, int k);
#endif

static Real drat, vflow, vflow0, betain, betaout,Bz;

#ifdef VISCOSITY
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3);

static Real nu;
#endif  /* VISCOSITY */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  int il,iu,jl,ju,kl,ku;
  Real x1,x2,x3, r;
  Real rho, vx, vy, dp, beta, v_turb, v_rot, r_cloud;

  int iseed;

#if (NSCALARS > 0)
  Real dye;
#endif
  x_shift = 0.;
  drat   = par_getd("problem", "drat");
  vflow  = par_getd("problem", "vflow");
  vflow0 = vflow;

  betain = par_getd("problem", "betain");
  betaout = par_getd("problem", "betaout");
  Bz = par_getd_def("problem","Bz",0.0);
  v_turb = par_getd_def("problem", "v_turb", 0.0);
  v_rot  = par_getd_def("problem", "v_rot",  0.0);

  r_cloud = par_getd_def("problem", "r_cloud", 0.25);

#ifdef VISCOSITY
  nu   = par_getd("problem","nu");
  NuFun_i = NULL;
  NuFun_a = nu_fun;
#endif

  iseed = -10;
#ifdef MPI_PARALLEL
  iseed -= myID_Comm_world;
#endif
  srand(iseed);

#ifdef FOLLOW_CLOUD
#if (NSCALARS == 0)
  ath_error("[problem]: requires NSCALARS > 0.\n");
#endif
#endif

#ifdef ENERGY_COOLING
  CoolingFunc = cooling_func;
#endif

#ifdef FOLLOW_CLOUD
  dump_history_enroll_alt(hst_xshift, "x_shift");
  dump_history_enroll_alt(hst_vflow,  "v_flow");
#endif

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        r = sqrt(x1*x1+x2*x2+x3*x3);

        rho  = 1.0;
        dp   = 0.0;
        vx   = vflow;

#if (NSCALARS > 0)
        dye = 0.0;
#endif

        if (r < r_cloud) {
          vx   = -(x2/r) * (r/r_cloud) * v_rot;
          vy   =  (x1/r) * (r/r_cloud) * v_rot;
          dp   = 0.0;
          rho *= drat;

#if (NSCALARS > 0)
          dye = 1.0;
#endif
        }

        /* write values to the grid */
        pGrid->U[k][j][i].d = rho;
        pGrid->U[k][j][i].M1 = rho * vx;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

        if (r < r_cloud) {
          pGrid->U[k][j][i].M1 += rho * RandomNormal(0.0, v_turb);
          pGrid->U[k][j][i].M2 += rho * RandomNormal(0.0, v_turb);
          if (pGrid->Nx[2] > 1)
            pGrid->U[k][j][i].M3 += rho * RandomNormal(0.0, v_turb);
        }

#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E = 1.0 + 0.5 * rho * SQR(vx);
        pGrid->U[k][j][i].E += dp / Gamma_1;
#endif  /* ISOTHERMAL */

#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = dye;
#endif
      }
    }
  }


  /* interface magnetic field */
  x1 -= 0.5 * pGrid->dx1;
  x2 -= 0.5 * pGrid->dx2;
  x3 -= 0.5 * pGrid->dx3;

  il = is;    iu = ie+1;
  jl = js;    ju = je+1;
  kl = ks;    ku = (pGrid->Nx[2] > 1) ? ke+1 : ke;

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef MHD
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        r = sqrt(x1*x1+x2*x2);

        beta = (r < r_cloud) ? betain : betaout;

        pGrid->B1i[k][j][i] = -x2/r * sqrt(2.0 * Gamma_1 / beta);
        pGrid->B2i[k][j][i] =  x1/r * sqrt(2.0 * Gamma_1 / beta);
        pGrid->B3i[k][j][i] = Bz;
#endif  /* MHD */
      }
    }
  }
  x1 += 0.5 * pGrid->dx1;
  x2 += 0.5 * pGrid->dx2;
  x3 += 0.5 * pGrid->dx3;


  /* cell-centered magnetic field */
  /*   derive this from interface field to be internally consistent
       with athena*/
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef MHD
        pGrid->U[k][j][i].B1c =
          0.5 * (pGrid->B1i[k][j][i] + pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = pGrid->U[k][j][i].B3c = 0.0;

        if (pGrid->Nx[1] > 1)
          pGrid->U[k][j][i].B2c =
            0.5 * (pGrid->B2i[k][j][i] + pGrid->B2i[k][j+1][i]);
        if (pGrid->Nx[2] > 1){
          pGrid->U[k][j][i].B3c =
            0.5 * (pGrid->B3i[k][j][i] + pGrid->B3i[k+1][j][i]);
	}
	else{
	   pGrid->U[k][j][i].B3c = Bz;
	}
#ifndef ISOTHERMAL
        /* add magnetic energy to the total energy */
        pGrid->U[k][j][i].E +=
          0.5 * (SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B2c)+
                 SQR(pGrid->U[k][j][i].B3c));
#endif  /* ISOTHERMAL */
#endif  /* MHD */
      }
    }
  }


  if (pDomain->Disp[0] == 0)
    bvals_mhd_fun(pDomain, left_x1,  bc_ix1);
  if (pDomain->Disp[1] == 0)
    bvals_mhd_fun(pDomain, left_x2,  bc_ix2);
  if (pDomain->Disp[2] == 0)
    bvals_mhd_fun(pDomain, left_x3,  bc_ix3);

  if (pDomain->MaxX[0] == pDomain->RootMaxX[0])
    bvals_mhd_fun(pDomain, right_x1, bc_ox1);
  if (pDomain->MaxX[1] == pDomain->RootMaxX[1])
    bvals_mhd_fun(pDomain, right_x2, bc_ox2);
  if (pDomain->MaxX[2] == pDomain->RootMaxX[2])
    bvals_mhd_fun(pDomain, right_x3, bc_ox3);
  /* seed a perturbation */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        pGrid->U[k][j][i].d *= (1.0 + RandomNormal(0.0, 0.05));
      }
    }
  }

  return;
}




/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
#ifdef FOLLOW_CLOUD
  fwrite(&x_shift, sizeof(Real), 1, fp);
  fwrite(&vflow,   sizeof(Real), 1, fp);
#endif

  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;

  for (nl=0; nl<(pM->NLevels); nl++) {
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++) {
      if (pM->Domain[nl][nd].Disp[0] == 0)
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x1,  bc_ix1);
      if (pM->Domain[nl][nd].Disp[1] == 0)
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x2,  bc_ix2);
      if (pM->Domain[nl][nd].Disp[2] == 0)
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x3,  bc_ix3);

      if (pM->Domain[nl][nd].MaxX[0] == pM->Domain[nl][nd].RootMaxX[0])
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x1, bc_ox1);
      if (pM->Domain[nl][nd].MaxX[1] == pM->Domain[nl][nd].RootMaxX[1])
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x2, bc_ox2);
      if (pM->Domain[nl][nd].MaxX[2] == pM->Domain[nl][nd].RootMaxX[2])
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x3, bc_ox3);
    }
  }

  drat  = par_getd("problem", "drat");
  vflow = par_getd("problem", "vflow");
  vflow0 = vflow;

  betain  = par_getd("problem", "betain");
  betaout = par_getd("problem", "betaout");
  Bz      = par_getd_def("problem", "Bz", 0.0);
#ifdef VISCOSITY
  nu   = par_getd("problem","nu");
  NuFun_i = NULL;
  NuFun_a = nu_fun;
#endif

#ifdef ENERGY_COOLING
  CoolingFunc = cooling_func;
#endif

#ifdef FOLLOW_CLOUD
  dump_history_enroll_alt(hst_xshift, "x_shift");
  dump_history_enroll_alt(hst_vflow,  "v_flow");
#endif


  /* DANGER: make sure the order here matches the order in write_restart() */
#ifdef FOLLOW_CLOUD
  fread(&x_shift, sizeof(Real), 1, fp);
  fread(&vflow,   sizeof(Real), 1, fp);
#endif

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  int nl, nd;

#ifdef REPORT_NANS
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {

        report_nans(&(pM->Domain[nl][nd]));
}
}
}
#endif


#ifdef FOLLOW_CLOUD
  Real dvx;
  dvx = cloud_mass_weighted_velocity(pM);
  ath_pout(0,"[Change dv]: %e\n", dvx);
  vflow -= dvx;
#endif

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef FOLLOW_CLOUD
        boost_frame(&(pM->Domain[nl][nd]), dvx);
#endif
       }
    }
  }

  return;
}

void Userwork_after_loop(MeshS *pM)
{
#ifdef FOLLOW_CLOUD
  ath_pout(0, "[follow_cloud]: shifted the domain by a total amount %e\n",
           x_shift);
#endif

  return;
}



/*==============================================================================
 * PHYSICS FUNCTIONS:
 * boost_frame()         - boost simulation frame by a velocity increment
 * report_nans()         - apply a ceiling and floor to the temperature
 * cloud_velocity()      - find the mass-weighted velocity of the cloud.
 * nu_fun()              - kinematic viscosity (i.e., cm^2/s)
 * cooling_func()        - cooling function for the energy equation
 *----------------------------------------------------------------------------*/
#ifdef FOLLOW_CLOUD
static void boost_frame(DomainS *pDomain, Real dvx)
{
  int i, j, k;
  int is,ie,js,je,ks,ke;

  Real d;

  GridS *pGrid = pDomain->Grid;
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        d = pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].M1 -= dvx * d;

#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E += 0.5 * d * SQR(dvx);
        pGrid->U[k][j][i].E -= dvx * pGrid->U[k][j][i].M1;
#endif  /* ISOTHERMAL */
      }
    }
  }

  x_shift -= dvx * pDomain->Grid->dt;

  return;
}
#endif  /* FOLLOW_CLOUD */


#ifdef REPORT_NANS
static void report_nans(DomainS *pDomain)
{
#ifndef ISOTHERMAL
  int i, j, k;
  int is,ie,js,je,ks,ke;

  Real KE, rho, press, temp,x1,x2,x3;
  int ntemp=0, nrho=0, nbad;
#ifdef MHD
  Real ME, beta;
  int nbeta=0;
#endif  /* MHD */

  Real scal[3];
#ifdef MPI_PARALLEL
  Real my_scal[3];
  int ierr;
#endif

  GridS *pGrid = pDomain->Grid;
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        rho = pGrid->U[k][j][i].d;

        KE = (SQR(pGrid->U[k][j][i].M1) +
              SQR(pGrid->U[k][j][i].M2) +
              SQR(pGrid->U[k][j][i].M3)) /
          (2.0 * rho);

        press = pGrid->U[k][j][i].E - KE;
#ifdef MHD
        ME = (SQR(pGrid->U[k][j][i].B1c) +
              SQR(pGrid->U[k][j][i].B2c) +
              SQR(pGrid->U[k][j][i].B3c)) * 0.5;
        press -= ME;
#endif  /* MHD */

        press *= Gamma_1;
        temp = press / rho;
#ifdef MHD
        beta = press / ME;
#endif  /* MHD */

        /* nans should also fail these tests */
        if ((temp < 1.0e-2 / drat)||(temp != temp)) {
          ntemp++;
	  if(ntemp<4){
	    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	    printf("bad temp:  %e Pos: %e %e %e\n",temp,x1,x2,x3);
	  }
          temp = 1.0e-2 / drat;
        }

        if ((rho < 1.0e-4)||(rho != rho)) {
          nrho++;
	  if(nrho<4){
	    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	    printf("bad rho:  %e Pos: %e %e %e\n",rho,x1,x2,x3);
	  }
	rho = 1.0e-4;
        }

#ifdef MHD
        if ((beta < 1.0e-2)||(beta != beta))
          nbeta++;
#endif  /* MHD */

        pGrid->U[k][j][i].E = temp*rho/Gamma_1 + KE;
#ifdef MHD
        pGrid->U[k][j][i].E += ME;
#endif  /* MHD */
      }
    }
  }

  /* synchronize over grids */
  scal[0] = ntemp;
  scal[1] = nrho;
#ifdef MHD
  scal[2] = nbeta;
#endif  /* MHD */

#ifdef MPI_PARALLEL
  my_scal[0] = scal[0];
  my_scal[1] = scal[1];
#ifdef MHD
  my_scal[2] = scal[2];
#endif  /* MHD */

  ierr = MPI_Allreduce(&my_scal, &scal, 3, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[report_nans]: MPI_Allreduce returned error %d\n", ierr);
#endif  /* MPI_PARALLEL */

  ntemp = scal[0];
  nrho  = scal[1];
#ifdef MHD
  nbeta = scal[2];
#endif  /* MHD */


  /* sum up the # of bad cells and report */
  nbad = MAX(ntemp, nrho);
#ifdef MHD
  nbad = MAX(nbad, nbeta);
#endif  /* MHD */

  if (nbad > 0) {
 #ifdef MHD
    ath_pout(0, "[report_nans]: found %d bad cells: %d temp, %d rho, %d beta.\n",
             nbad, ntemp, nrho, nbeta);
#else
    ath_pout(0, "[report_nans]: found %d bad cells: %d temp, %d rho.\n",
          nbad, ntemp, nrho);
#endif  /* MHD */
  }

#endif  /* ISOTHERMAL */

  return;
}
#endif  /* REPORT_NANS */


/* Cooling function to be called by the integrator */
/*   this actually goes into the energy equation, and is solved at 2nd
     order.  */
#ifdef ENERGY_COOLING
static Real cooling_func(const Real rho, const Real P, const Real dt)
{
  Real temp = P/rho;
  Real tcloud = Gamma_1 / drat;
  Real Edot, Edotmax;

  /* cool to tcloud, but no further */
  Edot = (temp < tcloud) ? 0.0 : 1000.0 * exp(-4.0 * temp / tcloud);

  /* limit cooling rate to keep pressure from going negative... */
  Edotmax = 0.3 * P/Gamma_1/dt;

  return MIN(Edot, Edotmax);
}
#endif  /* ENERGY_COOLING */


/* Return the center-of-mass velocity of the cloud */
/*   the scalar s obeys the same equation as density, so we weight
 *   the average as s*v. */
/*   for now, compute the average using the level 1 domain (i.e, the
 *   first refined level) */
#ifdef FOLLOW_CLOUD
static Real cloud_mass_weighted_velocity(MeshS *pM)
{
  GridS *pG;
  int i, j, k, is, ie, js, je, ks, ke;
  int nl, nd;

  Real s, d, scal[2];
#ifdef MPI_PARALLEL
  Real my_scal[2];
  int ierr;
#endif

  /* do the integral over level-1 domains, if they exist */
  nl = (pM->NLevels > 1) ? pM->NLevels - 2  : 0;

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if ((pM->Domain[nl][nd].Level == nl)
        && (pM->Domain[nl][nd].Grid != NULL)) {

      pG = pM->Domain[nl][nd].Grid;
      is = pG->is;  ie = pG->ie;
      js = pG->js;  je = pG->je;
      ks = pG->ks;  ke = pG->ke;

      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            d = pG->U[k][j][i].d;
            s = pG->U[k][j][i].s[0];

            scal[0] += s * pG->U[k][j][i].M1 / d;
            scal[1] += s;
          }
        }
      }

    }
  }

#ifdef MPI_PARALLEL
  my_scal[0] = scal[0];
  my_scal[1] = scal[1];

  ierr = MPI_Allreduce(&my_scal, &scal, 2, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[cloud_velocity]: MPI_Allreduce returned error %d\n", ierr);
#endif

  return scal[0] / scal[1];
}
#endif  /* FOLLOW_CLOUD */


#ifdef VISCOSITY
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3)
{
  return nu;
}
#endif  /* VISCOSITY */



/*==============================================================================
 * HISTORY OUTPUTS:
 *
 *----------------------------------------------------------------------------*/

#ifdef FOLLOW_CLOUD
static Real hst_xshift(GridS *pG, int i, int j, int k)
{
  return x_shift;
}
#endif

#ifdef FOLLOW_CLOUD
static Real hst_vflow(GridS *pG, int i, int j, int k)
{
  return vflow;
}
#endif




/*==============================================================================
 * BOUNDARY CONDITIONS:
 *
 *----------------------------------------------------------------------------*/

static void bc_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
  Real x1, x2, x3, r;
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
        pGrid->U[k][j][is-i].d  = 1.0;
        pGrid->U[k][j][is-i].M1 = 1.0 * vflow;
        pGrid->U[k][j][is-i].M2 = 0.0;
        pGrid->U[k][j][is-i].M3 = 0.0;
        pGrid->U[k][j][is-i].E  = 1.0 + 0.5*SQR(vflow) + Gamma_1 / betaout + 0.5*SQR(Bz);
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        r = sqrt(x1*x1+x2*x2);

        pGrid->B1i[k][j][i] = -x2/r * sqrt(2.0 * Gamma_1 / betaout);
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        r = sqrt(x1*x1+x2*x2);

        pGrid->B2i[k][j][i] =  x1/r * sqrt(2.0 * Gamma_1 / betaout);
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = Bz;
      }
    }
  }
#endif /* MHD */

  return;
}

static void bc_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;
}

static void bc_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i] = pGrid->U[k][js][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js][i];
      }
    }
  }
#endif /* MHD */

  return;
}


static void bc_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i] = pGrid->U[k][je][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je][i];
      }
    }
  }

/* j=je+1 is not a boundary condition for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je][i];
      }
    }
  }
#endif /* MHD */

  return;
}



/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x3 boundary (bc_ix3=2) */

static void bc_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ks][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks][j][i];
      }
    }
  }

/* B3i is not set at k=ks-nghost */
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2) */

static void bc_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ke][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke][j][i];
      }
    }
  }

/* k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}
