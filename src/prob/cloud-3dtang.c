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

static void initial_field(DomainS *pDomain, Real r_cloud);

static void check_div_b(GridS *pGrid);


/* add a force-free term to B */
void add_term(Real3Vect ***A, GridS *pG,
              Real theta, Real phi,
              Real alpha, Real beta, Real amp);

Real randomreal2(Real min, Real max);
static Real RandomNormal2(Real mu, Real sigma);

/* unit vectors */
Real3Vect get_e1(Real theta, Real phi);
Real3Vect get_e2(Real theta, Real phi);
Real3Vect get_e3(Real theta, Real phi);

/* dye-weighted hst quantities */
#if (NSCALARS > 0)
static Real hst_cE(const GridS *pG, const int i, const int j, const int k);
static Real hst_crho(const GridS *pG, const int i, const int j, const int k);

static Real hst_cMx(const GridS *pG, const int i, const int j, const int k);
static Real hst_cMy(const GridS *pG, const int i, const int j, const int k);
static Real hst_cMz(const GridS *pG, const int i, const int j, const int k);
#ifdef MHD
static Real hst_cBx(const GridS *pG, const int i, const int j, const int k);
static Real hst_cBy(const GridS *pG, const int i, const int j, const int k);
static Real hst_cBz(const GridS *pG, const int i, const int j, const int k);
#endif  /* MHD */
#endif  /* NSCALARS */



#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix);
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

static Real drat, vflow, vflow0, betain, betaout,dr,bz;

static Real tfloor, tceil, rhofloor, betafloor; /* Used in nancheck*/

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
  Real rho, vx, vy, dp,  v_turb, v_rot, r_cloud;
  Real fact;

  int iseed;
  Real3Vect ***A;    
  int nterms,nx1,nx2,nx3,ierr;
  Real theta, phi, alpha, beta, amp;
  Real scal[5];
  Real my_scal[5];

  Real x1min, x1max, x2min, x2max, x3min, x3max, tmp;


  Real jmag, bmag, JdB, JcBforce, norm, Brms, Bs, Bmax, ncells;
  Real bscale;
  int tangled;
#if (NSCALARS > 0)
  Real dye;
#endif

  drat   = par_getd("problem", "drat");
  vflow  = par_getd("problem", "vflow");
  vflow0 = vflow;
#ifdef FOLLOW_CLOUD
  x_shift = 0.0;
#endif
  betain = par_getd("problem", "betain");
  betaout = par_getd("problem", "betaout");

  v_turb = par_getd_def("problem", "v_turb", 0.0);
  v_rot  = par_getd_def("problem", "v_rot",  0.0);

  r_cloud = par_getd_def("problem", "r_cloud", 0.25);
  dr = par_getd_def("problem", "dr", 0.0);
  bz =  par_getd_def("problem", "Bz", 0.0);
  tangled = par_getd_def("problem","tangled",1);
  
  /*Real tfloor    = 1.0e-2 / drat;
  Real tceil     = 100.0;
  Real rhofloor  = 1.0e-2;
  Real betafloor = 3.0e-3;
  */
  tfloor = par_getd_def("problem", "tfloor", 1.e-2/drat);
  tceil = par_getd_def("problem", "tceil", 100.);
  rhofloor = par_getd_def("problem", "rhofloor", 1.e-2);
  betafloor = par_getd_def("problem", "betafloor", 3.e-3);
  d_MIN = par_getd_def("problem", "d_MIN", 1.e-4);


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

#if (NSCALARS > 0)
  dump_history_enroll(hst_cE,   "<c * E>");
  dump_history_enroll(hst_crho, "<c * rho>");

dump_history_enroll(hst_cMx, "<c * Mx>");
  dump_history_enroll(hst_cMy, "<c * My>");
  dump_history_enroll(hst_cMz, "<c * Mz>");
#ifdef MHD
  dump_history_enroll(hst_cBx, "<c * Bx>");
  dump_history_enroll(hst_cBy, "<c * By>");
  dump_history_enroll(hst_cBz, "<c * Bz>");
#endif /* MHD */
#endif

#ifdef REPORT_NANS
  nan_dump_count = 0;
#endif

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

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
          vx   = -(x2/r_cloud) * v_rot;
          vy   =  (x1/r_cloud) * v_rot;
          dp   = 0.0;

          rho  *= drat;

#if (NSCALARS > 0)
          dye = drat ; //1.0;
#endif
        }
	if (dr > 0.0){
	  rho = (1.0 + drat*0.5*(1.0+tanh((r_cloud-r)/(dr*r_cloud))));
	  vx = vflow*0.5*(1.0+tanh((-(1.0+dr)*r_cloud+r)/(dr*r_cloud)))/rho;
	  vx   += -(x2/r_cloud) * v_rot;
          vy   =  (x1/r_cloud) * v_rot;
	  if(r < r_cloud){
	    dye = rho;
	  }
	}
        /* write values to the grid */
        pGrid->U[k][j][i].d = rho;
        pGrid->U[k][j][i].M1 = rho * vx;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

	/* if (r < r_cloud && v_turb > 0.0) {
          pGrid->U[k][j][i].M1 += rho * RandomNormal(0.0, v_turb);
          pGrid->U[k][j][i].M2 += rho * RandomNormal(0.0, v_turb);
          if (pGrid->Nx[2] > 1)
            pGrid->U[k][j][i].M3 += rho * RandomNormal(0.0, v_turb);
        }
	*/

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
  if(tangled){
    // initial_field(pDomain, r_cloud);
    
    A = (Real3Vect***) calloc_3d_array(nx3, nx2, nx1, sizeof(Real3Vect));
    
    
    for (k=0; k<nx3; k++) {
      for (j=0; j<nx2; j++) {
	for (i=0; i<nx1; i++) {
	  A[k][j][i].x3 = A[k][j][i].x2 = A[k][j][i].x1 = 0.0;
	}
      }
    }
    
    
    
    nterms = par_getd_def("problem", "nterms",10);
    alpha = par_getd_def("problem","alpha",50.0);
    for (i=0; i<nterms; i++) {
      if (myID_Comm_world == 0) {
	theta = randomreal2(0.0, 2.0*PI);
	phi   = randomreal2(0.0, 2.0*PI);
	beta  = randomreal2(0.0, 2.0*PI);
	
	amp   = RandomNormal2(1.0, 0.25);
	
      } else {
	theta = phi = beta = amp = alpha = 0.0;
      }
      
#ifdef MPI_PARALLEL
      my_scal[0] = theta;
      my_scal[1] = phi;
      my_scal[2] = beta;
      my_scal[3] = amp;
      my_scal[4] = alpha;
      
      
      ierr = MPI_Allreduce(&my_scal, &scal, 5, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
      if (ierr)
	ath_error("[problem]: MPI_Allreduce returned error %d\n", ierr);
      
      theta = scal[0];
      phi   = scal[1];
      beta  = scal[2];
      amp   = scal[3];
      alpha = scal[4];
#endif
      
      add_term(A, pGrid, theta, phi, alpha, beta, amp);
    }
    
    for (k=0; k<nx3; k++) {
      for (j=0; j<nx2; j++) {
	for (i=0; i<nx1; i++) {
	  A[k][j][i].x1 *= sqrt(2.0*Gamma_1/nterms/betain);
	  A[k][j][i].x2 *= sqrt(2.0*Gamma_1/nterms/betain);
	  A[k][j][i].x3 *= sqrt(2.0*Gamma_1/nterms/betain);
	  
	  cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
	  r = sqrt(x1*x1 + x2*x2 + x3*x3);
	  if (r > 0.25)
	    A[k][j][i].x1 = A[k][j][i].x2 = A[k][j][i].x3 = 0.0;
	}
      }
    }
    

  } /*end of if(tangled) */


  /* take curl here */
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        pGrid->B1i[k][j][i] = 0.0;
	if(tangled){
	  pGrid->B1i[k][j][i] = (A[k][j+1][i].x3 - A[k][j][i].x3)/pGrid->dx2 -
	    (A[k+1][j][i].x2 - A[k][j][i].x2)/pGrid->dx3;
	}
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
	pGrid->B2i[k][j][i] = 0.0;
	if(tangled){
	  cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
	  bscale = 1.01/2.+(1.01/2.-.01)*tanh((-r_cloud*3.-x1)*3./r_cloud);
	  pGrid->B2i[k][j][i] = (A[k+1][j][i].x1 - A[k][j][i].x1)/pGrid->dx3 -
	    (A[k][j][i+1].x3 - A[k][j][i].x3)/pGrid->dx1 + bscale*sqrt(2.0 * Gamma_1 / betaout);
	}
	else{
	  pGrid->B2i[k][j][i] =  sqrt(2.0 * Gamma_1 / betaout);
	}
      }
    }
  }

   for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	pGrid->B3i[k][j][i] = 0.0;
	if(tangled){
	  cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
	  bscale = 1.01/2.+(1.01/2.-.01)*tanh((-r_cloud*3.-x1)*3./r_cloud);
	  pGrid->B3i[k][j][i] = (A[k][j][i+1].x2 - A[k][j][i].x2)/pGrid->dx1 -
	    (A[k][j+1][i].x1 - A[k][j][i].x1)/pGrid->dx2 + bscale*bz;
	}
	else{
	  pGrid->B3i[k][j][i] = bz;
	}
      }
    }
  }
   
 
#endif
   if(tangled){
     free_3d_array((void***) A);
   }
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

#ifdef MHD


  Brms = Bmax = 0.0;
  ncells = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
       

        Bs = (SQR(pGrid->U[k][j][i].B1c) +
              SQR(pGrid->U[k][j][i].B2c) +
              SQR(pGrid->U[k][j][i].B3c));

        if(Bs > 1.e-6){
	  Brms += Bs;
	  ncells += 1.0;
	}
        Bmax = MAX(Bmax, Bs);
      }
    }
  }

#ifdef MPI_PARALLEL
  my_scal[0] = Brms;
  my_scal[1] = ncells;

  ierr = MPI_Allreduce(&my_scal, &scal, 2, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[problem]: MPI_Allreduce returned error %d\n", ierr);

  Brms   = scal[0];
  ncells = scal[1];

  my_scal[0] = Bmax;

  ierr = MPI_Allreduce(&my_scal, &scal, 1, MPI_RL, MPI_MAX, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[problem]: MPI_Allreduce returned error %d\n", ierr);

  Bmax   = scal[0];
#endif

  Brms = sqrt(Brms/ncells);
  Bmax = sqrt(Bmax);

  ath_pout(0, "Brms = %f\tBmax = %f\n", Brms, Bmax);
  /* ath_error("Brms = %f\tBmax = %f\n", Brms, Bmax); */




#endif //MHD


#ifdef MHD
  check_div_b(pGrid);
#endif //MHD

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
	fact = -1.0;
	while (abs(fact) > 0.03)
	  fact = (RandomNormal(0.0, 0.01));
	if(abs(fact) < .03)
	  pGrid->U[k][j][i].d *= (1.0+fact);
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
  bz =  par_getd_def("problem", "Bz", 0.0);
  tfloor = par_getd_def("problem", "tfloor", 1.e-2/drat);
  tceil = par_getd_def("problem", "tceil", 100.);
  rhofloor = par_getd_def("problem", "rhofloor", 1.e-2);
  betafloor = par_getd_def("problem", "betafloor", 3.e-3);
  d_MIN = par_getd_def("problem", "d_MIN", 1.e-4);
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
  Real divb[5],divcell,x,y,z;
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
	divcell = (pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
	  + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2;
        if (three_d)
          divcell+= (pGrid->B3i[k+1][j][i]-pGrid->B3i[k][j][i])/pGrid->dx3;
	divb[0] += fabs(divcell);
#if (NSCALARS > 1)
	pGrid->U[k][j][i].s[1] = divcell;
#endif
	
	/*if(fabs(divcell) > 1.e-13){
	  cc_pos(pGrid,i,j,k,&x,&y,&z);

	  printf("bdv %e %d %d %d  %e  %e    %e %e %e\n", pGrid->time, i,j,k,  divcell, 1./pGrid->dx1, x,y,z);
	}
	*/
	n++;
      }
    }
  }
  //divb[0] /= n;

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

void Userwork_before_loop(MeshS *pM)
{
  int nl, nd, ntot;

  /* report nans first, so we can fix them before they propagate into
     the following functions. */
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef REPORT_NANS
        ntot = report_nans(pM, &(pM->Domain[nl][nd]),1);
        //if(ntot > 0)
	//report_nans(pM, &(pM->Domain[nl][nd]),1);
#endif
      }
    }
  }

  return;
}

void Userwork_in_loop(MeshS *pM)
{
  int nl, nd, ntot;
#ifdef FOLLOW_CLOUD
  Real dvx, newdvx, expt;
#endif


#ifdef FOLLOW_CLOUD
  dvx = cloud_mass_weighted_velocity(pM);
  if(abs(dvx) > .01 || dvx < 0.0){
    ath_pout(0,"[bad dvx:] %0.20e setting to 0.\n",dvx);
    dvx = 0.0;
  }
 

  if(dvx > 0.0){
    expt = floor(log10(dvx));
    newdvx = dvx / pow(10, expt);
    newdvx = floor(newdvx * 1.0e4) / 1.0e4;
    newdvx = newdvx * pow(10.0, expt);
    dvx = newdvx;
  }
  ath_pout(0,"[dvx:]  %0.20e\n",dvx);
  
  if(vflow - dvx < 0.0){
    dvx = vflow;
    vflow = 0.0;
    
  }
  else{
    vflow -= dvx;
  }

 
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

  if (pM->dt < 1.e-6){
    data_output(pM,1);
    ath_error("dt too small\n");
    
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
        
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E += 0.5 * d * SQR(dvx);
        pGrid->U[k][j][i].E -= dvx * pGrid->U[k][j][i].M1;
#endif  /* ISOTHERMAL */
	pGrid->U[k][j][i].M1 -= dvx * d;

      }
    }
  }

  //  x_shift -= dvx * pDomain->Grid->dt;
  x_shift -= vflow * pDomain->Grid->dt;

  return;
}
#endif  /* FOLLOW_CLOUD */


#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix)
{
#ifndef ISOTHERMAL
  int i, j, k;
  int is,ie,js,je,ks,ke;
  Real x1, x2, x3;
  int V=0; //verbose off
  int NO = 2;
  Real KE, rho, press, temp;
  int nanpress=0, nanrho=0, nanv=0, nnan;   /* nan count */
  int npress=0,   nrho=0,   nv=0,   nfloor; /* floor count */
#ifdef MHD
  Real ME;
  int nanmag=0;
  int nmag=0;
#endif  /* MHD */
  Real beta;
  Real scal[8];
#ifdef MPI_PARALLEL
  Real my_scal[8];
  int ierr;
#endif

  /*Real tfloor    = 1.0e-2 / drat;
  Real tceil     = 100.0;
  Real rhofloor  = 1.0e-2;
  Real betafloor = 3.0e-3;
  */
  GridS *pGrid = pDomain->Grid;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        rho = pGrid->U[k][j][i].d;
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
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
#else
	beta = abs(200.*betafloor);
#endif  /* MHD */

        if (press != press) {
          nanpress++;
	  if(V && nanpress < NO) printf("bad press %e R %e  %e %e %e  %d %d %d %e %e %e %e\n",press, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta); 
          if(fix) 
	    temp = tfloor;
        } else if (temp < tfloor) {
          npress++;
	  if(V && npress < NO) printf("bad tempF %e R %e  %e %e %e  %d %d %d %e %e %e %e\n",temp, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta); 
          if(fix) 
	    temp = tfloor;
        } else if (temp > tceil && beta > 10.0 * betafloor) {
          npress++;
	  if(V && npress < NO) printf("bad tempC %e R %e  %e %e %e  %d %d %d %e %e %e %e\n",temp, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta); 
          if(fix) 
	    temp = tceil;
        }

        if (rho != rho) {
          nanrho++;
	  if(V&&nanrho < NO) printf("bad rho %e R %e  %e %e %e  %d %d %d\n",rho, 1.0/pGrid->dx1, x1, x2, x3, i, j, k); 
          if(fix) 
	    rho = rhofloor;
        } else if (rho < rhofloor) {
          nrho++;
	  if(V&& nrho < NO) printf("bad rho %e R %e  %e %e %e  %d %d %d\n",rho, 1.0/pGrid->dx1, x1, x2, x3, i, j, k); 
          if(fix) 
	    rho = rhofloor;
        }

        if (pGrid->U[k][j][i].M1 != pGrid->U[k][j][i].M1) {
          nanv++;
          if(fix) 
	    pGrid->U[k][j][i].M1 = 0.0;
        }
        if (pGrid->U[k][j][i].M2 != pGrid->U[k][j][i].M2) {
          nanv++;
          if(fix) 
	    pGrid->U[k][j][i].M2 = 0.0;
        }
        if (pGrid->U[k][j][i].M3 != pGrid->U[k][j][i].M3) {
          nanv++;
          if(fix) 
	    pGrid->U[k][j][i].M3 = 0.0;
        }

#ifdef MHD
        if (ME != ME) {
          nanmag++;
          /* TODO: apply a fix to B? */
        } else if (beta < betafloor && betafloor > 0.0) {
	  nmag++;
	  if(V && nmag < NO) printf("bad mag %e R %e  %e %e %e  %d %d %d %e %e %e %e %e\n",beta, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta, ME); 
	  if(fix){
	    rho  = MAX(rho, sqrt(betafloor * ME));
	    temp = betafloor * ME / rho;
	    temp = MAX(temp,tfloor);
	    if(V && nmag < NO) printf("bad magf %e R %e  %e %e %e  %d %d %d %e %e %e %e %e\n",beta, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta,ME); 
	  }
	}
#endif  /* MHD */

        /* write values back to the grid */
        /* TODO: what about B??? */
	if(fix) {
	  pGrid->U[k][j][i].d  = rho;
	  	  KE = (SQR(pGrid->U[k][j][i].M1) +
		SQR(pGrid->U[k][j][i].M2) +
		SQR(pGrid->U[k][j][i].M3)) /
	    (2.0 * rho);
	  
	  pGrid->U[k][j][i].E = temp*rho/Gamma_1 + KE;
#ifdef MHD
	  pGrid->U[k][j][i].E += ME;
#endif  /* MHD */
	}
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
  nnan = nanpress+nanrho+nanv;

#ifdef MHD
  nnan += nanmag;
#endif  /* MHD */

  /* sum up the # of floored cells and report */
  nfloor = npress+nrho+nv;

#ifdef MHD
  nfloor += nmag;
#endif  /* MHD */

  // if (fix == 0) {
#ifdef MHD
    ath_pout(0, "[report_nans]: floored %d cells: %d P, %d d, %d v, %d beta.\n",
      nfloor, npress, nrho, nv, nmag);
#else
    ath_pout(0, "[report_nans]: floored %d cells: %d P, %d d, %d v.\n",
      nfloor, npress, nrho, nv);
#endif  /* MHD */
    // }


  //  if ((nnan > 0 || nfloor -nmag > 30) && fix == 0) {
    if (nnan > 0 ){// && fix == 0) {
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
    nan_dump.out    = "prim";
    nan_dump.nlevel = -1;       /* dump all levels */


    dump_vtk(pM, &nan_dump);
    if(nnan) nan_dump_count++;  
    if (nan_dump_count > 10)
      ath_error("[report_nans]: too many nan'd timesteps.\n");
  }


#endif  /* ISOTHERMAL */

  return nfloor+nnan;
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

  Real s, d, scal[2], tmp, x1, x2, x3;
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
	    cc_pos(pG,i,j,k,&x1,&x2,&x3);
            tmp = s * pG->U[k][j][i].M1 / d;
            if (tmp == tmp && x1 < 0.0) {

              scal[0] += tmp*d;
              scal[1] += s*d;
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
  Real newnu;
  return nu;
  /* newnu = nu*pow(T*1.5,2.5); */
  /* if(newnu != newnu || newnu < 0.0) */
  /*   newnu = TINY_NUMBER; */
  /* if(newnu > 3.*nu) */
  /*   newnu = 3.*nu; */

  /* return newnu; */
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

#if (NSCALARS > 0)
	pGrid->U[k][j][i].s[0] = 0.0;
#endif

        pGrid->U[k][j][is-i].d  = 1.0;
        pGrid->U[k][j][is-i].M1 = 1.0 * vflow;
        pGrid->U[k][j][is-i].M2 = 0.0;
        pGrid->U[k][j][is-i].M3 = 0.0;
        pGrid->U[k][j][is-i].E  = 1.0 + 0.5*SQR(vflow);
#ifdef MHD
        pGrid->U[k][j][is-i].B1c = 0.0;
        pGrid->U[k][j][is-i].B2c = sqrt(2.0 * Gamma_1 / betaout);;
        pGrid->U[k][j][is-i].B3c = bz;
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
    for (j=js; j<=je; j++) {
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
    for (j=js; j<=ju; j++) {
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
    for (j=js; j<=je; j++) {
      for (i=0; i<nghost; i++) {
        pGrid->B3i[k][j][i] = bz;
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
	if(pGrid->U[k][j][ie+i].M1 < 0.0){
	  pGrid->U[k][j][ie+i].E -= 0.5*SQR(pGrid->U[k][j][ie+i].M1)/pGrid->U[k][j][ie+i].d;
	  pGrid->U[k][j][ie+i].M1 = 0.0;
	}
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
	if( pGrid->U[k][js-j][i].M2 > 0.0){
	  pGrid->U[k][js-j][i].E -=  0.5*SQR(pGrid->U[k][js-j][i].M2)/ pGrid->U[k][js-j][i].d;
	  pGrid->U[k][js-j][i].M2 = 0.0; 
	}
	  
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
      if(pGrid->U[k][je+j][i].M2 < 0.0){
	pGrid->U[k][je+j][i].E -=  0.5*SQR(pGrid->U[k][je+j][i].M2)/ pGrid->U[k][je+j][i].d;
	pGrid->U[k][je+j][i].M2 = 0.0; 
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
	if(pGrid->U[ks-k][j][i].M3 > 0.0){
	  pGrid->U[ks-k][j][i].E -= 0.5*SQR(pGrid->U[ks-k][j][i].M3)/pGrid->U[ks-k][j][i].d;
	  pGrid->U[ks-k][j][i].M3 = 0.0;
	}
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
	if(pGrid->U[ke+k][j][i].M3 < 0.0){
	  pGrid->U[ke+k][j][i].E -= 0.5*SQR(pGrid->U[ke+k][j][i].M3)/pGrid->U[ke+k][j][i].d;
	  pGrid->U[ke+k][j][i].M3 = 0.0;
	}
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


void add_term(Real3Vect ***A, GridS *pG,
              Real theta, Real phi,
              Real alpha, Real beta, Real amp)
{
  int i, j ,k;
  Real phase;
  Real3Vect e1, e2, e3, kv, r;

  int is, ie, ks, ke, js, je;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  e1 = get_e1(theta, phi);
  e2 = get_e2(theta, phi);
  e3 = get_e3(theta, phi);

  /* e3 x e1 = e2 */
  /* think of e3 as x-axis, e1 as y-axis, e2 as z-axis */
  /* then use this: */
  /* A =   cos(a x)/a z^ +   sin(a x)/a y^ */
  /* B =   cos(a x)   z^ +   sin(a x)   y^ */
  /* j = a*cos(a x)   z^ + a*sin(a x)   y^ */

  /* think of k as x^, e1 as y^, e2 as z^ */
  /* k cross y = z, so this is consistent  */

  /* e3.x1 = 1.0;   e3.x2 = 0.0;   e3.x3 = 0.0; */
  /* e2.x1 = 0.0;   e2.x2 = 0.0;   e2.x3 = 1.0; */
  /* e1.x1 = 0.0;   e1.x2 = 1.0;   e1.x3 = 0.0; */
  /* alpha = 2.0 * PI; */

  kv.x1 = alpha * e3.x1;
  kv.x2 = alpha * e3.x2;
  kv.x3 = alpha * e3.x3;
  if(ie+nghost >= (ie-is)+1+2*nghost){
    ath_error("overrunning A array!\n");
  }
    

  for (k=0; k<=ke+nghost; k++) {
    for (j=0; j<=je+nghost; j++) {
      for (i=0; i<=ie+nghost; i++) {
        cc_pos(pG, i, j, k, &r.x1, &r.x2, &r.x3);
        r.x1 -= 0.5 * pG->dx1;
        r.x2 -= 0.5 * pG->dx2;
        r.x3 -= 0.5 * pG->dx3;
        phase = r.x1 * kv.x1 + r.x2 * kv.x2 + r.x3 * kv.x3 + beta;

        A[k][j][i].x1 += amp * (e2.x1 * cos(phase) + e1.x1 * sin(phase))/alpha;
        A[k][j][i].x2 += amp * (e2.x2 * cos(phase) + e1.x2 * sin(phase))/alpha;
        A[k][j][i].x3 += amp * (e2.x3 * cos(phase) + e1.x3 * sin(phase))/alpha;
      }
    }
  }

  return;
}




Real3Vect get_e1(Real theta, Real phi)
{
  Real3Vect ret;
  ret.x1 = -1.0 * sin(theta);
  ret.x2 = cos(theta) * cos(phi);
  ret.x3 = cos(theta) * sin(phi);

  return ret;
}

Real3Vect get_e2(__attribute__((unused))Real theta, Real phi)
{
  Real3Vect ret;
  ret.x1 = 0.0;
  ret.x2 = -1.0 * sin(phi);
  ret.x3 = cos(phi);

  return ret;
}

Real3Vect get_e3(Real theta, Real phi)
{
  Real3Vect ret;
  ret.x1 = cos(theta);
  ret.x2 = sin(theta) * cos(phi);
  ret.x3 = sin(theta) * sin(phi);

  return ret;
}

Real randomreal2(Real min, Real max)
{
  Real eta = ((Real)rand()/(Real)RAND_MAX);
  return min + eta * (max-min);
}

static Real RandomNormal2(Real mu, Real sigma)
/* Implements the box-muller routine.  Gives a mean of mu, and a
   standard deviation sigma.  */
{
  Real x1, x2, w, y1;
  static Real y2;
  static int use_last = 0;

  if (use_last){ /* use value from previous call */
    y1 = y2;
    use_last = 0;
  }
  else {
    do {
      x1 = randomreal2(-1.0, 1.0);
      x2 = randomreal2(-1.0, 1.0);
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }

  return (mu + y1 * sigma);
}

/* dye-weighted hst quantities */
#if (NSCALARS > 0)
static Real hst_cE(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].E * pG->U[k][j][i].d);
}

static Real hst_crho(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].s[0]);
}

static Real hst_cMx(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M1  /  pG->U[k][j][i].d);
}

static Real hst_cMy(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M2  /  pG->U[k][j][i].d);
}

static Real hst_cMz(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M3 /  pG->U[k][j][i].d);
}

#ifdef MHD
static Real hst_cBx(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].B1c);
}

static Real hst_cBy(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].B2c);
}

static Real hst_cBz(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].B3c);
}
#endif  /* MHD */
#endif  /* NSCALARS */
