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

// Global variables uses

static Real mach;               // Mach number
static Real alpha;              // exponent for cooling function
static int num_steps;           // number of steps for cooling function
static Real lambda_0;           // normalization of the cooling curve
static Real n_kh;               // t_sim / t_kh
static Real end_time;
static Real f, gm;                 // gamma = cp/cv


static int cooling_flag;
static int heating_flag;


// parameters for simulation length
//
static const Real n_cool = 1.0/3.0; // t_sim / t_{cool,hot}


// helper functions to reduce code
//
static inline Real window(Real x, Real width, Real a);
static inline Real get_tfloor(Real cstcool);
static inline Real get_cstcool(const Real time);
static void integrate_cooling(GridS *pG);
// history outputs
//
//Real hst_cstcool(MeshBlock *pmb, int iout);
//Real hst_tfloor(MeshBlock *pmb, int iout);


// mass of hot and cold gas...
// ... (scalars are useless for this problem)
//Real hst_rho_hot(MeshBlock *pmb, int iout);
//Real hst_rho_cold(MeshBlock *pmb, int iout);

// mass*velocity for hot and cold gas...
// ... ratio with previous entries to get mass-weighted velocity
//Real hst_rho_v_hot(MeshBlock *pmb, int iout);
//Real hst_rho_v_cold(MeshBlock *pmb, int iout);



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
  static int frst=1;  /* flag so new history variables enrolled only once */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  
  mach         = par_getd("problem", "mach");
  alpha        = par_getd("cooling", "alpha");
  end_time     = par_getd("time", "tlim");
  n_kh         = par_getd("problem", "n_kh");
  gm           = par_getd("problem", "gamma");
  
  num_steps    = par_geti("cooling", "steps");

  cooling_flag = par_geti("cooling", "cooling");
  heating_flag = par_geti("cooling", "heating");


  //gm and f are used in cooling routines
  f = pow(gm, 1.5)/(gm-1.0)/(2.0-alpha);
  f = 1.0/f;

  //calculate lambda_0

  lambda_0=1/f;
  lambda_0 *= pow(get_cstcool(0.0), 1.0 / (4.0-2.0*alpha));
  lambda_0 *= pow(mach * n_cool / n_kh, (5.0-2.0*alpha)/(4.0-2.0*alpha));
  
  /*check that runtime is consistent with mach number*/

  Real tkh = pow(lambda_0 * f * get_cstcool(0.0), -1.0/(5.0-2.0*alpha));
  tkh = tkh / (mach * sqrt(gm));

  Real err = fabs(end_time - n_kh * tkh)/end_time;

  if (err >= 1.0e-4) {
    ath_error("#### FATAL ERROR in problem file \n Simulation tlim should equal %f\n", n_kh * tkh);
  }

  /*check to see if variables are defined correctly*/  

  ath_pout(0, "gamma =  %f\n", gm);
  ath_pout(0, "f =  %f\n", f);
  ath_pout(0, "lambda_0 = %f\n", lambda_0);
  ath_pout(0, "init_cstcool = %f\n", get_cstcool(0.0));
  

  /*write initial conditions*/

  Real amp   =  par_getd("problem", "amp");
  Real a     = par_getd("problem", "a");
  Real width = par_getd("problem", "width");

  Real vflow = mach * sqrt(gm);

  Real drat = 1.0 / get_tfloor(get_cstcool(0.0));
  Real P = 1.0;
 
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


	//set up hot cold kh beam
	pGrid->U[k][j][i].d = 1.0 + (drat-1) * window(x2, width, a);
	pGrid->U[k][j][i].M1 = (vflow * (1.0 - window(x2, width, a)))*pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

	
	// white noise perturbations to velocity
	pGrid->U[k][j][i].M1 += amp*(((Real)rand()/RAND_MAX) - 0.5)*pGrid->U[k][j][i].d * vflow;
        pGrid->U[k][j][i].M2 += amp*(((Real)rand()/RAND_MAX) - 0.5)*pGrid->U[k][j][i].d * vflow;
      

	
	pGrid->U[k][j][i].E = P/Gamma_1 + (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) + SQR(pGrid->U[k][j][i].M3))/(2.0*pGrid->U[k][j][i].d);
	
#ifdef MHD
	pGrid->U[k][j][i].E += (SQR(pGrid->U[k][j][i].B1c)
				+ SQR(pGrid->U[k][j][i].B2c)
				+ SQR(pGrid->U[k][j][i].B3c))*0.5;
#endif
      }
    }
  }


  
/* enroll new history variables, only once  */

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
  const Real x = time / end_time;
  const Real fact = -1.0 * floor(x * num_steps);

  return pow(4.0, fact);
}


static inline Real get_tfloor(Real cstcool)
{
  // cs*tcool = Tfloor^(2.5-alpha)/(f*Lambda_0)

  return pow(f * lambda_0 * cstcool, 2.0/(5.0-2.0*alpha));
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
  Real temp, tfloor;
  Real s = (Real) num_steps,
    x = pG->time / (n_kh*end_time);

  Real deltaE[2];
#ifdef MPI_PARALLEL
  Real deltaE_global[2];
  int ierr;
#endif  /* MPI_PARALLEL */

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  tfloor = get_tfloor(get_cstcool(pG->time));
  
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
