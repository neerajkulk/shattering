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
/* #define ENERGY_COOLING */

static void bc_ix1(GridS *pGrid);
static void bc_ox1(GridS *pGrid);
static void bc_ix2(GridS *pGrid);
static void bc_ox2(GridS *pGrid);

static void initial_field(DomainS *pDomain, Real r_cloud);

static void check_div_b(GridS *pGrid);

#ifdef REPORT_NANS
static void report_nans(MeshS *pM, DomainS *pDomain);
static OutputS nan_dump;
static int nan_dump_count;
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

static Real drat, vflow, vflow0, betain, betaout;

#ifdef VISCOSITY
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3);

static Real nu;
#endif  /* VISCOSITY */

static Real pro(Real r, Real rcloud)
{
  return (r/rcloud - log(cosh(r/rcloud))) / log(2);
}

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  int il,iu,jl,ju,kl,ku;
  Real x1,x2,x3, r;
  Real rho, vx, vy, dp, beta, v_turb, v_rot, r_cloud;
  Real fact;

  int iseed;

#if (NSCALARS > 0)
  Real dye;
#endif

  drat   = par_getd("problem", "drat");
  vflow  = par_getd("problem", "vflow");
  vflow0 = vflow;
  x_shift = 0.0;

  betain = par_getd("problem", "betain");
  betaout = par_getd("problem", "betaout");

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

#ifdef REPORT_NANS
  nan_dump_count = 0;
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

          rho  *= drat;

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

  initial_field(pDomain, r_cloud);

#ifdef MHD
  check_div_b(pGrid);
#endif MHD

  if (pDomain->Disp[0] == 0)
    bvals_mhd_fun(pDomain, left_x1,  bc_ix1);
  /* if (pDomain->Disp[1] == 0) */
  /*   bvals_mhd_fun(pDomain, left_x2,  bc_ix2); */

  /* if (pDomain->MaxX[0] == pDomain->RootMaxX[0]) */
  /*   bvals_mhd_fun(pDomain, right_x1, bc_ox1); */
  /* if (pDomain->MaxX[1] == pDomain->RootMaxX[1]) */
  /*   bvals_mhd_fun(pDomain, right_x2, bc_ox2); */

  /* seed a perturbation */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        fact = -1.0;
        while (fact < 0.1)
          fact = (1.0 + RandomNormal(0.0, 0.05));
        pGrid->U[k][j][i].d *= fact;
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
      /* if (pM->Domain[nl][nd].Disp[1] == 0) */
      /*   bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x2,  bc_ix2); */

      /* if (pM->Domain[nl][nd].MaxX[0] == pM->Domain[nl][nd].RootMaxX[0]) */
      /*   bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x1, bc_ox1); */
      /* if (pM->Domain[nl][nd].MaxX[1] == pM->Domain[nl][nd].RootMaxX[1]) */
      /*   bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x2, bc_ox2); */
    }
  }

  drat  = par_getd("problem", "drat");
  vflow = par_getd("problem", "vflow");
  vflow0 = vflow;

  betain  = par_getd("problem", "betain");
  betaout = par_getd("problem", "betaout");

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

#ifdef MHD
static void check_div_b(GridS *pGrid)
{
  int i,j,k;
  int is,js,ie,je,ks,ke;
  int n;
  Real divb[5];
#ifdef MPI_PARALLEL
  int ierr;
  Real my_divb[5];
#endif
  int three_d = (pGrid->Nx[2] > 1) ? 1 : 0;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;


  divb[0] = 0.0;
  n = 0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        if(three_d)
          divb[0] += fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
                          + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2
                          + (pGrid->B3i[k+1][j][i]-pGrid->B3i[k][j][i])/pGrid->dx3);
        else
          divb[0] +=
            fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
                 + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
        n++;
      }
    }
  }
  divb[0] /= n;

  /* check divB along ix1: */
  divb[1] = 0.0;
  n=0;
  k = ks;
  i = is;
  for (j=js; j<=je; j++) {
    divb[1] +=
      fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
           + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
    n++;
  }
  divb[1] /= n;

  /* check divB along ox1: */
  divb[2] = 0.0;
  n=0;
  k = ks;
  i = ie;
  for (j=js; j<=je; j++) {
    divb[2] +=
      fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
           + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
    n++;
  }
  divb[2] /= n;


  /* check divB along ix2: */
  divb[3] = 0.0;
  n=0;
  k = ks;
  j = js;
  for (i=is; i<=ie; i++) {
    divb[3] +=
      fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
           + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
    n++;
  }
  divb[3] /= n;


  /* check divB along ox2: */
  divb[4] = 0.0;
  n=0;
  k = ks;
  j = je;
  for (i=is; i<=ie; i++) {
    divb[4] +=
      fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
           + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
    n++;
  }
  divb[4] /= n;


#ifdef MPI_PARALLEL
  for (i=0; i<=4; i++)
    my_divb[i] = divb[i];
  ierr = MPI_Allreduce(&my_divb, &divb, 5, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[check_div_b]: MPI_Allreduce returned error %d\n", ierr);
#endif

  if (three_d) {
    ath_pout(0, "divb = %e\n", divb[0]);
  } else {
    ath_pout(0, "divb = %e\t%e\t%e\t%e\t%e\n",
             divb[0], divb[1], divb[2], divb[3], divb[4]);
  }

  return;
}
#endif  /* MHD */

void Userwork_in_loop(MeshS *pM)
{
  int nl, nd;
#ifdef FOLLOW_CLOUD
  Real dvx;
#endif

  /* report nans first, so we can fix them before they propagate into
     the following functions. */
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef REPORT_NANS
        report_nans(pM, &(pM->Domain[nl][nd]));
#endif
      }
    }
  }

#ifdef FOLLOW_CLOUD
  dvx = cloud_mass_weighted_velocity(pM);
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

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef MHD
        check_div_b(pM->Domain[nl][nd].Grid);
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
static void report_nans(MeshS *pM, DomainS *pDomain)
{
#ifndef ISOTHERMAL
  int i, j, k;
  int is,ie,js,je,ks,ke;

  Real KE, rho, press, temp;
  int nanpress=0, nanrho=0, nanv=0, nnan;   /* nan count */
  int npress=0,   nrho=0,   nv=0,   nfloor; /* floor count */
#ifdef MHD
  Real ME, beta;
  int nanmag=0;
  int nmag=0;
#endif  /* MHD */

  Real scal[8];
#ifdef MPI_PARALLEL
  Real my_scal[8];
  int ierr;
#endif

  Real tfloor    = 1.0e-2 / drat;
  Real tceil     = 100.0;
  Real rhofloor  = 1.0e-2;
  Real betafloor = 1.0e-2;

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

        if (press != press) {
          nanpress++;
          temp = tfloor;
        } else if (temp < tfloor) {
          npress++;
          temp = tfloor;
        } else if (temp > tceil) {
          npress++;
          temp = tceil;
        }

        if (rho != rho) {
          nanrho++;
          rho = rhofloor;
        } else if (rho < rhofloor) {
          nrho++;
          rho = rhofloor;
        }

        if (pGrid->U[k][j][i].M1 != pGrid->U[k][j][i].M1) {
          nanv++;
          pGrid->U[k][j][i].M1 = 0.0;
        }
        if (pGrid->U[k][j][i].M2 != pGrid->U[k][j][i].M2) {
          nanv++;
          pGrid->U[k][j][i].M2 = 0.0;
        }
        if (pGrid->U[k][j][i].M3 != pGrid->U[k][j][i].M3) {
          nanv++;
          pGrid->U[k][j][i].M3 = 0.0;
        }

#ifdef MHD
        if (ME != ME) {
          nanmag++;
          /* TODO: apply a fix to B? */
        } else if (beta < betafloor) {
          rho  = MAX(rho, sqrt(betafloor * ME));
          temp = betafloor * ME / rho;
        }
#endif  /* MHD */

        /* write values back to the grid */
        /* TODO: what about B??? */
        pGrid->U[k][j][i].d  = rho;
        pGrid->U[k][j][i].E = temp*rho/Gamma_1 + KE;
#ifdef MHD
        pGrid->U[k][j][i].E += ME;
#endif  /* MHD */
      }
    }
  }

  /* synchronize over grids */
#ifdef MPI_PARALLEL
  my_scal[0] = nanpress;
  my_scal[1] = nanrho;
  my_scal[2] = nanv;
#ifdef MHD
  my_scal[3] = nanmag;
#endif  /* MHD */
  my_scal[4] = npress;
  my_scal[5] = nrho;
  my_scal[6] = nv;
#ifdef MHD
  my_scal[7] = nmag;
#endif  /* MHD */

  ierr = MPI_Allreduce(&my_scal, &scal, 8, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[report_nans]: MPI_Allreduce returned error %d\n", ierr);

  nanpress = scal[0];
  nanrho   = scal[1];
  nanv     = scal[2];
#ifdef MHD
  nanmag   = scal[3];
#endif  /* MHD */
  npress   = scal[4];
  nrho     = scal[5];
  nv       = scal[6];
#ifdef MHD
  nmag     = scal[7];
#endif  /* MHD */
#endif  /* MPI_PARALLEL */


  /* sum up the # of bad cells and report */
  nnan = MAX(nanpress, nanrho);
  nnan = MAX(nnan,     nanv);
#ifdef MHD
  nnan = MAX(nnan,     nanmag);
#endif  /* MHD */

  /* sum up the # of floored cells and report */
  nfloor = MAX(npress, nrho);
  nfloor = MAX(nfloor, nv);
#ifdef MHD
  nfloor = MAX(nfloor, nmag);
#endif  /* MHD */

  if (nfloor > 0) {
#ifdef MHD
    ath_pout(0, "[report_nans]: floored %d cells: %d P, %d d, %d v, %d beta.\n",
      nfloor, npress, nrho, nv, nmag);
#else
    ath_pout(0, "[report_nans]: floored %d cells: %d P, %d d, %d v.\n",
      nfloor, npress, nrho, nv);
#endif  /* MHD */
  }


  if (nnan > 0) {
#ifdef MHD
    ath_pout(0, "[report_nans]: found %d nan cells: %d P, %d d, %d v, %d B.\n",
             nnan, nanpress, nanrho, nanv, nanmag);
#else
    ath_pout(0, "[report_nans]: found %d nan cells: %d P, %d d, %d v.\n",
             nnan, nanpress, nanrho, nanv);
#endif  /* MHD */

    nan_dump.n      = 100;
    nan_dump.dt     = HUGE_NUMBER;
    nan_dump.t      = pM->time;
    nan_dump.num    = 1000 + nan_dump_count;
    nan_dump.out    = "cons";
    nan_dump.nlevel = -1;       /* dump all levels */
    nan_dump_count++;

    dump_vtk(pM, &nan_dump);

    if (nan_dump_count > 10)
      ath_error("[report_nans]: too many nan'd timesteps.\n");
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
  Edotmax = 0.01 * P/Gamma_1/dt;

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

  Real s, d, scal[2], tmp;
#ifdef MPI_PARALLEL
  Real my_scal[2];
  int ierr;
#endif

  /* do the integral over level-1 domains, if they exist */
  nl = (pM->NLevels > 1) ? 1 : 0;

  scal[0] = scal[1] = 0.0;
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {

      pG = pM->Domain[nl][nd].Grid;
      is = pG->is;  ie = pG->ie;
      js = pG->js;  je = pG->je;
      ks = pG->ks;  ke = pG->ke;

      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            d = pG->U[k][j][i].d;
            s = pG->U[k][j][i].s[0];

            tmp = s * pG->U[k][j][i].M1 / d;
            if (tmp == tmp) {
              scal[0] += tmp;
              scal[1] += s;
            }

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
 * INITIAL CONDITION:
 *
 *----------------------------------------------------------------------------*/

static void initial_field(DomainS *pDomain, Real r_cloud)
{
  GridS *pGrid = pDomain->Grid;
  int i, j, k, i2, j2, k2;

  int is, ie, js, je, ks, ke;
  int il, iu, jl, ju, kl, ku;

  Real x1, x2, x3, x1f, x2f, x3f, r;
  int nx1, nx2, nx3;

  /* vector potential */
  Real ***az, ***ay, ***ax;

  /* Read in a random vector potential from a file */
  /* -- hard-code the size in for now */
  Real *xvals, *yvals, *zvals;
  Real ***axin, ***ayin, ***azin;
  const int  insize  = 256;
  Real inrange = 3.0 * r_cloud;
  int num;

  Real grad_a[3], dr[3];

  FILE *infile;
  char *fname = par_gets_def("problem", "vecpot_file", "vecpot.dat");



  /* populate indata[] from file */
  xvals = (Real*) calloc_1d_array(insize, sizeof(Real));
  yvals = (Real*) calloc_1d_array(insize, sizeof(Real));
  zvals = (Real*) calloc_1d_array(insize, sizeof(Real));

  axin = (Real***) calloc_3d_array(insize, insize, insize, sizeof(Real));
  ayin = (Real***) calloc_3d_array(insize, insize, insize, sizeof(Real));
  azin = (Real***) calloc_3d_array(insize, insize, insize, sizeof(Real));

  for (i=0; i<insize; i++){
    xvals[i] = yvals[i] = zvals[i] = (2.0*inrange)*i/(insize-1) - inrange;
  }

  infile = fopen(fname, "r");
  if (infile == NULL)
    ath_error("[initial_field]: could not open input file: %s\n", fname);

  for(k=0; k<insize; k++) {
    for (j=0; j<insize; j++) {
      num = fread(axin[k][j], sizeof(double), insize, infile);
      if (num != insize)
        ath_error("[initial_field]: error reading ax from input file: %s\n", fname);
    }
  }

  for(k=0; k<insize; k++) {
    for (j=0; j<insize; j++) {
      num = fread(ayin[k][j], sizeof(double), insize, infile);
      if (num != insize)
        ath_error("[initial_field]: error reading ay from input file: %s\n", fname);
    }
  }

  for(k=0; k<insize; k++) {
    for (j=0; j<insize; j++) {
      num = fread(azin[k][j], sizeof(double), insize, infile);
      if (num != insize)
        ath_error("[initial_field]: error reading az from input file: %s\n", fname);
    }
  }

  fclose(infile);



  /* sort out coordinates */
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  il = is;    iu = ie+1;
  jl = js;    ju = je+1;
  kl = ks;    ku = (pGrid->Nx[2] > 1) ? ke+1 : ke;

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  ax = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real));
  ay = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real));
  az = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real));


  /*  */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef MHD
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        x1f = x1 - 0.5*pGrid->dx1;
        x2f = x2 - 0.5*pGrid->dx2;
        x3f = (pGrid->Nx[2] > 1) ? x3 - 0.5*pGrid->dx3 : 0.0;

        r = sqrt(x1f*x1f + x2f*x2f + x3f*x3f);

        i2 = locate(xvals, x1f, insize);
        j2 = locate(yvals, x2f, insize);
        k2 = locate(zvals, x3f, insize);

        if (r > r_cloud) {
          az[k][j][i] = ay[k][j][i] = ax[k][j][i] = 0.0;
        } else {
          dr[0] = x1f-xvals[i2];
          dr[1] = x2f-yvals[j2];
          dr[2] = x3f-zvals[k2];

          /* interpolate ax from input*/
          grad_a[0] = (axin[k2  ][j2  ][i2+1]-axin[k2][j2][i2])/(xvals[i2+1]-xvals[i2]);
          grad_a[1] = (axin[k2  ][j2+1][i2  ]-axin[k2][j2][i2])/(yvals[j2+1]-yvals[j2]);
          grad_a[2] = (axin[k2+1][j2  ][i2  ]-axin[k2][j2][i2])/(zvals[k2+1]-zvals[k2]);

          ax[k][j][i] = axin[k2][j2][i2]
            + dr[0] * grad_a[0]
            + dr[1] * grad_a[1]
            + dr[2] * grad_a[2];

          /* interpolate ay from input*/
          grad_a[0] = (ayin[k2  ][j2  ][i2+1]-ayin[k2][j2][i2])/(xvals[i2+1]-xvals[i2]);
          grad_a[1] = (ayin[k2  ][j2+1][i2  ]-ayin[k2][j2][i2])/(yvals[j2+1]-yvals[j2]);
          grad_a[2] = (ayin[k2+1][j2  ][i2  ]-ayin[k2][j2][i2])/(zvals[k2+1]-zvals[k2]);

          ay[k][j][i] = ayin[k2][j2][i2]
            + dr[0] * grad_a[0]
            + dr[1] * grad_a[1]
            + dr[2] * grad_a[2];

          /* interpolate az from input*/
          grad_a[0] = (azin[k2  ][j2  ][i2+1]-azin[k2][j2][i2])/(xvals[i2+1]-xvals[i2]);
          grad_a[1] = (azin[k2  ][j2+1][i2  ]-azin[k2][j2][i2])/(yvals[j2+1]-yvals[j2]);
          grad_a[2] = (azin[k2+1][j2  ][i2  ]-azin[k2][j2][i2])/(zvals[k2+1]-zvals[k2]);

          az[k][j][i] = azin[k2][j2][i2]
            + dr[0] * grad_a[0]
            + dr[1] * grad_a[1]
            + dr[2] * grad_a[2];
        }

#endif  /* MHD */
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        pGrid->B1i[k][j][i] = (az[k][j+1][i] - az[k][j][i])/pGrid->dx2 -
          (ay[k+1][j][i] - ay[k][j][i])/pGrid->dx3;
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][j][i] = (ax[k+1][j][i] - ax[k][j][i])/pGrid->dx3 -
          (az[k][j][i+1] - az[k][j][i])/pGrid->dx1;
      }
    }
  }

  ku = (ke > ks) ? ke+1 : ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[k][j][i] = (ay[k][j][i+1] - ay[k][j][i])/pGrid->dx1 -
          (ax[k][j+1][i] - ax[k][j][i])/pGrid->dx2;
      }
    }
  }
#endif

#ifdef MHD
  /* Normalize the RMS magnitude of B */
  Real rmsB = 0.0;
  int cnt = 0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        x1f = x1 - 0.5*pGrid->dx1;
        x2f = x2 - 0.5*pGrid->dx2;
        x3f = (pGrid->Nx[2] > 1) ? x3 - 0.5*pGrid->dx3 : 0.0;

        r = sqrt(x1f*x1f + x2f*x2f + x3f*x3f);

        if (r <= r_cloud) {
          rmsB += SQR(pGrid->B1i[k][j][i])
                  + SQR(pGrid->B2i[k][j][i])
                  + SQR(pGrid->B3i[k][j][i]);
          cnt += 1;
        }
      }
    }
  }
  printf("[beta orig]: %e\n", sqrt(rmsB/cnt));
  rmsB = 9.147912e+02; //for 32 5.752040e+02; // For 16 cells per cloud radius!  3.430934e+02;//sqrt(rmsB / cnt);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        pGrid->B1i[k][j][i] *= sqrt(2.0 * Gamma_1 / betain) / rmsB;
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][j][i] *= sqrt(2.0 * Gamma_1 / betain) / rmsB;
      }
    }
  }

  ku = (ke > ks) ? ke+1 : ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[k][j][i] *= sqrt(2.0 * Gamma_1 / betain) / rmsB;
      }
    }
  }
#endif

  /* cell-centered magnetic field */
  /*   derive this from interface field to be internally consistent
       with athena */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef MHD
        pGrid->U[k][j][i].B1c =
          0.5 * (pGrid->B1i[k][j][i] + pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = pGrid->B2i[k][j][i];
        pGrid->U[k][j][i].B3c = pGrid->B3i[k][j][i];

        if (pGrid->Nx[1] > 1)
          pGrid->U[k][j][i].B2c =
            0.5 * (pGrid->B2i[k][j][i] + pGrid->B2i[k][j+1][i]);
        if (pGrid->Nx[2] > 1)
          pGrid->U[k][j][i].B3c =
            0.5 * (pGrid->B3i[k][j][i] + pGrid->B3i[k+1][j][i]);

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


  free_3d_array((void***) ax);
  free_3d_array((void***) ay);
  free_3d_array((void***) az);

  free_1d_array((void*) xvals);
  free_1d_array((void*) yvals);
  free_1d_array((void*) zvals);

  free_3d_array((void***) axin);
  free_3d_array((void***) ayin);
  free_3d_array((void***) azin);

  return;
}



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
        pGrid->U[k][j][is-i].E  = 1.0 + 0.5*SQR(vflow);
#ifdef MHD
        pGrid->U[k][j][is-i].B1c = 0.0;
        pGrid->U[k][j][is-i].B2c = sqrt(2.0 * Gamma_1 / betaout);;
        pGrid->U[k][j][is-i].B3c = 0.0;
        if(i == 1)
          pGrid->U[k][j][is-i].B1c = 0.5*pGrid->B1i[k][j][is];

        pGrid->U[k][j][is-i].E  += 0.5*(SQR(pGrid->U[k][j][is-i].B1c)
                                        +SQR(pGrid->U[k][j][is-i].B2c)
                                        +SQR(pGrid->U[k][j][is-i].B3c));
#endif
      }
    }
  }


#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=1; i<nghost; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        x1 -= 0.5 * pGrid->dx1;
        r = sqrt(x1*x1+x2*x2);

        pGrid->B1i[k][j][i] = 0.0;
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=0; i<nghost; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        x2 -= 0.5 * pGrid->dx2;
        r = sqrt(x1*x1+x2*x2);

        pGrid->B2i[k][j][i] = sqrt(2.0 * Gamma_1 / betaout);
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=0; i<nghost; i++) {
        pGrid->B3i[k][j][i] = 0.0;
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
