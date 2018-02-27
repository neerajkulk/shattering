#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#ifdef MPI_PARALLEL
#include "mpi.h"
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else /* DOUBLE_PREC */
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */

// Global variables used

static Real mach;
static Real alpha;              // exponent for cooling function
static Real lambda_0;           // normalization of the cooling curve
static Real drat;
static int cooling_flag;
static int heating_flag;
static Real f, gm;
static Real tkh,tcool,v;
static Real tkhtcool;
static Real tfloor;

// helper functions to reduce code
//
static inline Real window(Real x, Real width, Real a);
static inline Real get_tfloor(Real cstcool);
static inline Real get_cstcool(const Real time);
static void integrate_cooling(GridS *pG);

// history outputs



static Real hst_rho_hot(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_cold(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_cool(const GridS *pG, const int i, const int j, const int k);


static Real hst_rho_v_hot(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_v_cold(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_v_cool(const GridS *pG, const int i, const int j, const int k);

static Real hst_rhosq(const GridS *pG, const int i, const int j, const int k);



static double ran2(long int *idum);
/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real x1,x2,x3;
  long int iseed = -1;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  
  mach         = par_getd("problem", "mach");
  drat         = par_getd("problem", "drat");
  tkhtcool     = par_getd("problem", "tkh/tcool");
  alpha        = par_getd("cooling", "alpha");
  tfloor       = par_getd("cooling", "floor");

  gm           = par_getd("problem", "gamma");

  cooling_flag = par_geti("cooling", "cooling");
  heating_flag = par_geti("cooling", "heating");


  //define velocity and the kh timescale

  
  //gm and f are used in cooling routines
  v = mach * sqrt(gm);
  tkh = (1.0+drat)/(sqrt(drat)*v); // rho_hot = 1.0 and rho_cold = drat
  

  //tcool = T^(2-alpha) * gm / ((gm-1)*(2-alpha)*lambda0) :linear theory result from mike's TI paper
  // Use tkh/tcool (set in athinput) to solve for lambda0...

  lambda_0 = tkhtcool/tkh;
  lambda_0 *= pow(drat,alpha-2.0);
  lambda_0 *= gm/(gm-1.0);
  lambda_0 *= 1.0/(2.0 - alpha);
  

  
  /*write initial conditions*/

  Real amp   = par_getd("problem", "amp");
  Real a     = par_getd("problem", "a");
  Real width = par_getd("problem", "width");
  Real P = 1.0;
  Real r;
  srand(iseed);
  
#ifdef MPI_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srand(iseed-rank);
#endif
    /* set hot cold kh beam with random noise in velocity */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	r = sqrt(SQR(x2)+SQR(x3));
	
	//set up hot cold kh beam
	pGrid->U[k][j][i].d = 1.0 + (drat - 1.0) * window(r, width, a);
	pGrid->U[k][j][i].M1 = (v * (1.0 - window(r, width, a)))*pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

	
	// white noise perturbations to velocity
	pGrid->U[k][j][i].M1 += amp*(((Real)rand()/RAND_MAX) - 0.5)*pGrid->U[k][j][i].d * v;
        pGrid->U[k][j][i].M2 += amp*(((Real)rand()/RAND_MAX) - 0.5)*pGrid->U[k][j][i].d * v;
      

	
	pGrid->U[k][j][i].E = P/Gamma_1 + (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) + SQR(pGrid->U[k][j][i].M3))/(2.0*pGrid->U[k][j][i].d);
	
#ifdef MHD
	pGrid->U[k][j][i].E += (SQR(pGrid->U[k][j][i].B1c)
				+ SQR(pGrid->U[k][j][i].B2c)
				+ SQR(pGrid->U[k][j][i].B3c))*0.5;
#endif
      }
    }
  }

#ifdef VISCOSITY
  extern Real nu_iso, nu_aniso;
  nu_aniso = 0.0;
  Real dx = 1.0/(par_getd("domain1","Nx1"));
  nu_iso =  par_getd_def("problem","nu_iso",0.0); // in units of dx*cs
  nu_iso *= dx * sqrt(gm);
#endif /* VISCOSITY */


  
/* enroll new history variables, only once  */
  dump_history_enroll(hst_rho_hot, "rho_hot");
  dump_history_enroll(hst_rho_v_hot, "rho_v_hot");
  dump_history_enroll(hst_rho_cold, "rho_cold");
  dump_history_enroll(hst_rho_v_cold, "rho_v_cold");
  dump_history_enroll(hst_rho_cool, "rho_cool");
  dump_history_enroll(hst_rho_v_cool, "rho_v_cool");
  dump_history_enroll(hst_rhosq, "rho^2");

  
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
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

#if (NSCALARS > 0)
/*! \fn static Real color(const GridS *pG, const int i, const int j,const int k)
 *  \brief Returns first passively advected scalar s[0] */
static Real color(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].s[0]/pG->U[k][j][i].d;
}
#endif

ConsFun_t get_usr_expr(const char *expr)
{
#if (NSCALARS > 0)
  if(strcmp(expr,"color")==0) return color;
#endif
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  int nl, nd, ntot;
  GridS *pGrid;
  
  int gnx1 = pM->Nx[0];
  int gnx2 = pM->Nx[1];
  int gnx3 = pM->Nx[2];
  
  /* loop over the mesh to find the grid on this processor */
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pGrid = pM->Domain[nl][nd].Grid;

#ifndef BAROTROPIC
        integrate_cooling(pGrid);
#endif  /* BAROTROPIC */
      }
    }
  }

  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

static inline Real get_cstcool(const Real time)
{

  return 1.0;

}


static inline Real get_tfloor(Real cstcool)
{

  return 1.0;
}

static Real window(Real x, Real width, Real a)
{
  return 0.5 * (tanh(a * (width - x)) +
                tanh(a * (width + x)));
}


static void integrate_cooling(GridS *pG)
{
  PrimS W;
  ConsS U;
  int i, j, k;
  int is, ie, js, je, ks, ke;
  Real temp;

  Real deltaE[2];
#ifdef MPI_PARALLEL
  Real deltaE_global[2];
  int ierr;
#endif  /* MPI_PARALLEL */

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  
  deltaE[0] = deltaE[1] = 0.0;
  
  if(cooling_flag){
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {

          /* first, get the temperature */
          W = Cons_to_Prim(&(pG->U[k][j][i]));
          temp = W.P/W.d;

          /* cooling law */
          temp -= (gm-1) * W.d * lambda_0 * pow(temp, alpha) * pG->dt;

          /* apply a temperature floor (nans tolerated) */
          if (isnan(temp) || temp < tfloor)
            temp = tfloor;

	  /* apply a temperature ceiling which is 1.5 times initial hot gas temperature */
	  else if (temp > 1.5)
	    temp = 1.5;

          W.P = W.d * temp;
          U = Prim_to_Cons(&W);

          deltaE[0] += pG->U[k][j][i].E - U.E;
	  deltaE[1] += 1.0;
          pG->U[k][j][i].E = U.E;
        }
      }
    }
  }


  if (heating_flag){
#ifdef MPI_PARALLEL
    ierr = MPI_Allreduce(&deltaE, &deltaE_global, 2, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
    if (ierr)
      ath_error("[integrate_cooling]: MPI_Allreduce returned error %d\n", ierr);

    deltaE[0] = deltaE_global[0];
    deltaE[1] = deltaE_global[1];
#endif  /* MPI_PARALLEL */

    deltaE[0] /= deltaE[1];

    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pG->U[k][j][i].E += deltaE[0];
        }
      }
    }
  }

  return;

}



/* ========================================================================== */
/* history outputs
 */


Real hst_rho_hot(const GridS *pG, const int i, const int j, const int k)
{
  Real tcut = 1.1;
  Real ret = 0.0;

  PrimS W;
  ConsS U;

  W = Cons_to_Prim(&(pG->U[k][j][i]));
  Real temp = W.P/W.d;

  if (temp >=tcut) {
    ret = W.d;
  }

  return ret;
}


/* TODO: DECIDE ON A PROPER VALUE FOR TCUT */



Real hst_rho_v_hot(const GridS *pG, const int i, const int j, const int k)
{
  Real tcut = 1.1;
  Real ret = 0.0;

  PrimS W;
  ConsS U;

  W = Cons_to_Prim(&(pG->U[k][j][i]));
  Real temp = W.P/W.d;

  if (temp >=tcut) {
    ret = pG->U[k][j][i].M1;
  }

  return ret;
}


/*tracks gas just above the floor temperature*/

Real hst_rho_cold(const GridS *pG, const int i, const int j, const int k)
{
  Real tcut = 2.0*tfloor;
  Real ret = 0.0;
  PrimS W;
  ConsS U;

  W = Cons_to_Prim(&(pG->U[k][j][i]));
  Real temp = W.P/W.d;

  if (temp <=tcut) {
    ret = W.d;
  }

  return ret;
}




Real hst_rho_v_cold(const GridS *pG, const int i, const int j, const int k)
{
  Real tcut = 2.0*tfloor;
  Real ret = 0.0;
  PrimS W;
  ConsS U;

  W = Cons_to_Prim(&(pG->U[k][j][i]));
  Real temp = W.P/W.d;

  if (temp <=tcut) {
    ret = pG->U[k][j][i].M1;
  }

  return ret;
}


/*tracks initial cool gas as it cools further (this gas is not necessarily at the floor temperature)*/
Real hst_rho_cool(const GridS *pG, const int i, const int j, const int k)
{
  Real tcut = (2.0/drat);
  Real ret = 0.0;
  PrimS W;
  ConsS U;

  W = Cons_to_Prim(&(pG->U[k][j][i]));
  Real temp = W.P/W.d;

  if (temp <=tcut) {
    ret = W.d;
  }

  return ret;
}




Real hst_rho_v_cool(const GridS *pG, const int i, const int j, const int k)
{
  Real tcut = (2.0/drat);
  Real ret = 0.0;
  PrimS W;
  ConsS U;

  W = Cons_to_Prim(&(pG->U[k][j][i]));
  Real temp = W.P/W.d;

  if (temp <=tcut) {
    ret = pG->U[k][j][i].M1;
  }

  return ret;
}



Real hst_rhosq(const GridS *pG, const int i, const int j, const int k)
{
  return SQR(pG->U[k][j][i].d);
}


