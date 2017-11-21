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


  Real ***az, ***ay, ***ax;
  int nx1, nx2, nx3;
  Real x1f, x2f, x3f;

  /* interface magnetic field */
  il = is;    iu = ie+1;
  jl = js;    ju = je+1;
  kl = ks;    ku = (pGrid->Nx[2] > 1) ? ke+1 : ke;

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  if ((ay = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[problem]: Error allocating memory for vector pot\n");
  }
  if ((az = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[problem]: Error allocating memory for vector pot\n");
  }
  if ((ax = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[problem]: Error allocating memory for vector pot\n");
  }

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef MHD
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        x1f = x1 - 0.5*pGrid->dx1;
        x2f = x2 - 0.5*pGrid->dx2;
        x3f = x3 - 0.5*pGrid->dx3;

        r = sqrt(x1f*x1f + x2f*x2f);
        beta = (r < r_cloud) ? betain : betaout;

        ax[k][j][i] = 0.0;
        ay[k][j][i] = 0.0;

        az[k][j][i] = -x1f * pow(pro(r, r_cloud), 10) * sqrt(2.0 * Gamma_1 / betaout)
          + pro(r, 0.25) * sqrt(2.0 * Gamma_1 / betain);
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
  check_div_b(pGrid);
#endif MHD

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


  if (pDomain->Disp[0] == 0)
    bvals_mhd_fun(pDomain, left_x1,  bc_ix1);
  if (pDomain->Disp[1] == 0)
    bvals_mhd_fun(pDomain, left_x2,  bc_ix2);

  if (pDomain->MaxX[0] == pDomain->RootMaxX[0])
    bvals_mhd_fun(pDomain, right_x1, bc_ox1);
  if (pDomain->MaxX[1] == pDomain->RootMaxX[1])
    bvals_mhd_fun(pDomain, right_x2, bc_ox2);

  /* seed a perturbation */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	//        pGrid->U[k][j][i].d *= (1.0 + RandomNormal(0.0, 0.05));
      }
    }
  }

  free_3d_array((void***)az);
  free_3d_array((void***)ay);
  free_3d_array((void***)ax);

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

      if (pM->Domain[nl][nd].MaxX[0] == pM->Domain[nl][nd].RootMaxX[0])
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x1, bc_ox1);
      if (pM->Domain[nl][nd].MaxX[1] == pM->Domain[nl][nd].RootMaxX[1])
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x2, bc_ox2);
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

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;


  divb[0] = 0.0;
  n = 0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
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

  /* ath_error("total   divb = %e\nbc_ix1: divb = %e\nbc_ox1: divb = %e\nbc_ix2: divb = %e\nbc_ox2: divb = %e\n", */
  /*          divb[0], divb[1], divb[2], divb[3], divb[4]); */

  ath_pout(0, "divb = %e\t%e\t%e\t%e\t%e\n",
           divb[0], divb[1], divb[2], divb[3], divb[4]);

  return;
}
#endif  /* MHD */

void Userwork_in_loop(MeshS *pM)
{
  int nl, nd;


  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        /* report nans first, so we can fix them before they propagate
           into the following functions. */
#ifdef REPORT_NANS
        report_nans(pM, &(pM->Domain[nl][nd]));
#endif

      }
    }
  }

#ifdef FOLLOW_CLOUD
  Real dvx;
  dvx = cloud_mass_weighted_velocity(pM);
  vflow -= dvx;
#endif

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef FOLLOW_CLOUD
        boost_frame(&(pM->Domain[nl][nd]), dvx);
#endif
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
  int nanpress=0, nanrho=0, nanv=0, nnan;
#ifdef MHD
  Real ME, beta;
  int nanmag=0;
#endif  /* MHD */

  Real scal[4];
#ifdef MPI_PARALLEL
  Real my_scal[4];
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

        if ((pGrid->U[k][j][i].M1 != pGrid->U[k][j][i].M1)
            || (pGrid->U[k][j][i].M2 != pGrid->U[k][j][i].M2)
            || (pGrid->U[k][j][i].M3 != pGrid->U[k][j][i].M3))
          nanv++;

        if (press != press)
          nanpress++;
        if (rho != rho)
          nanrho++;
#ifdef MHD
        if (ME != ME)
          nanmag++;
#endif
      }
    }
  }

  /* synchronize over grids */
#ifdef MPI_PARALLEL
  my_scal[0] = scal[0] = nanpress;
  my_scal[1] = scal[1] = nanrho;
  my_scal[2] = scal[2] = nanv;
#ifdef MHD
  my_scal[3] = scal[3] = nanmag;
#endif  /* MHD */

  ierr = MPI_Allreduce(&my_scal, &scal, 4, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[report_nans]: MPI_Allreduce returned error %d\n", ierr);
#endif  /* MPI_PARALLEL */

  nanpress = scal[0];
  nanrho   = scal[1];
  nanv     = scal[2];
#ifdef MHD
  nanmag   = scal[3];
#endif  /* MHD */



  /* sum up the # of bad cells and report */
  nnan = MAX(nanpress, nanrho);
  nnan = MAX(nnan,     nanv);
#ifdef MHD
  nnan = MAX(nnan,     nanmag);
#endif  /* MHD */

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

        if ((temp != temp) || (temp < (1.0e-2 / drat)))
          temp = 1.0e-2 / drat;
        else if (temp > 100.0)
          temp = 100.0;

        if ((rho != rho) || (rho < 1.0e-2))
          rho = 1.0e-2;

        if (pGrid->U[k][j][i].M1 != pGrid->U[k][j][i].M1)
          pGrid->U[k][j][i].M1 = 0.0;
        if (pGrid->U[k][j][i].M2 != pGrid->U[k][j][i].M2)
          pGrid->U[k][j][i].M2 = 0.0;
        if (pGrid->U[k][j][i].M3 != pGrid->U[k][j][i].M3)
          pGrid->U[k][j][i].M3 = 0.0;


#ifdef MHD
        if ((beta != beta) || (beta < 1.0e-2)) {
          rho = MAX(rho, 0.1);
          temp = 1.0e-2 * ME / rho;
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
