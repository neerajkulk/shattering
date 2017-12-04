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



/* -------------------------------------------------------------------------- */
/* global variables and prototypes to generate a random initial
   condition using FFTs
*/

/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
   For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis-nghost,gnx1))+      \
  SQR(KCOMP(j,gjs-nghost,gnx2))+      \
			    SQR(KCOMP(k,gks-nghost,gnx3))))

/* FFTW - Variables, Plan, etc. */
static struct ath_3d_fft_plan *plan;
static ath_fft_data *fd=NULL;   /* unnormalized */
static Real ***dd=NULL;

/* Cutoff wavenumbers, G&O spect peak, power law spect exponent, 2 pi/L */
static Real klow,khigh,expo,dkx;
/* Number of cells in local grid, number of cells in global grid */
static int nx1,nx2,nx3,gnx1,gnx2,gnx3;
/* Starting and ending indices for global grid */
static int gis,gie,gjs,gje,gks,gke;
/* Seed for random number generator */
long int rseed;

/* Functions appear in this file in the same order that they appear in
   the prototypes below */

/* Function prototypes for generating velocity perturbations */
static void pspect(ath_fft_data *ampl);
static inline void transform();
static inline void generate();
static void perturb(GridS *pGrid);

/* Function prototypes for initializing and interfacing with Athena */
static void initialize_fourier(GridS *pGrid, DomainS *pD);

/* Function prototypes for Numerical Recipes functions */
static double ran2(long int *idum);

/* end FFT initial condition stuff */
/* -------------------------------------------------------------------------- */



// Global variables uses

static Real mach;
static Real alpha;              // exponent for cooling function
static int num_steps;           // number of steps for cooling function
static Real lambda_0;           // normalization of the cooling curve
static Real n_kh;               // t_sim / t_kh
static Real end_time;
static Real drat;
static int cooling_flag;
static int heating_flag;
static Real f, gm;

#ifdef VISCOSITY
static Real nu_param;
#endif // VISCOSITY

#ifdef THERMAL_CONDUCTION
static Real kappa_fun(const Real d, const Real T,
                      const Real x1, const Real x2, const Real x3);

static Real f_sp;             /* normalization for the conductivity */
#endif  /* THERMAL_CONDUCTION */




// tfloor drop times. time spent in each step reduces by a factor of 2
static const Real geometric [5] = {0.51612903225806451613, 0.77419354838709677419, 0.90322580645161290323, 0.96774193548387096774, 1.0};




static int cooling_flag;
static int heating_flag;


  /*adding in mike's temperature depedant viscosity routine*/

  #ifdef VISCOSITY
/* viscosity won't work without conduction... */
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3);
#endif  /* VISCOSITY */
 


// helper functions to reduce code
//
static inline Real window(Real x, Real width, Real a);
static inline Real get_tfloor(Real cstcool);
static inline Real get_cstcool(const Real time);
static void integrate_cooling(GridS *pG);

// history outputs


static Real hst_cstcool(const GridS *pG, const int i, const int j, const int k);
static Real hst_tfloor(const GridS *pG, const int i, const int j, const int k);

static Real hst_rho_hot(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_cold(const GridS *pG, const int i, const int j, const int k);


static Real hst_rho_v_hot(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_v_cold(const GridS *pG, const int i, const int j, const int k);

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
  static int frst=1;  /* flag so new history variables enrolled only once */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  
  mach         = par_getd("problem", "mach");
  drat         = par_getd("problem", "drat");
  alpha        = par_getd("cooling", "alpha");
  end_time     = par_getd("time", "tlim");
  gm           = par_getd("problem", "gamma");
  num_steps    = par_geti("cooling", "steps");

  cooling_flag = par_geti("cooling", "cooling");
  heating_flag = par_geti("cooling", "heating");


#ifdef VISCOSITY
  NuFun_i = nu_fun;
  NuFun_a = NULL;
  nu_param = par_getd("problem", "nu_iso");
#endif  /* VISCOSITY */

#ifdef THERMAL_CONDUCTION
  KappaFun_i  = kappa_fun;
  KappaFun_a  = NULL;
  f_sp        = 0.0;
#endif  /* THERMAL_CONDUCTION */



  
  //gm and f are used in cooling routines
  f = pow(gm, 1.5)/(gm-1.0)/(2.0-alpha);
  f = 1.0/f;



  
  // set cstcool = 0.25 initially in cold gas. for a given drat, lambda_0 is determined.
  lambda_0 = 1.0/f;
  lambda_0 *= 1/pow(drat,2.5 - alpha);
  lambda_0 *= 4.0;

  
  initialize_fourier(pGrid, pDomain);
  generate();
  
  
  /* this sets the grid density field */
  perturb(pGrid);


  
  /*write initial conditions*/
  Real noise = par_getd("problem", "noise");
  Real amp   =  par_getd("problem", "amp");
  Real a     = par_getd("problem", "a");
  Real width = par_getd("problem", "width");
  Real vflow = mach * sqrt(gm);
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

	//put fourier noise in density with parameter noise.

	pGrid->U[k][j][i].d *= (noise/2.0);
	pGrid->U[k][j][i].d += 1.0 - noise;
	
	/* density perturbations are now normalized to have a mean of ~1 with drho/rho ~ noise */
	
	//set up hot cold kh beam
	pGrid->U[k][j][i].d *= 1.0 + (drat - 1.0) * window(x2, width, a);
	pGrid->U[k][j][i].M1 = (vflow * (1.0 - window(x2, width, a)))*pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

	
	/* // white noise perturbations to velocity */
	/* pGrid->U[k][j][i].M1 += amp*(((Real)rand()/RAND_MAX) - 0.5)*pGrid->U[k][j][i].d * vflow; */
        /* pGrid->U[k][j][i].M2 += amp*(((Real)rand()/RAND_MAX) - 0.5)*pGrid->U[k][j][i].d * vflow; */
      

	
	pGrid->U[k][j][i].E = P/Gamma_1 + (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) + SQR(pGrid->U[k][j][i].M3))/(2.0*pGrid->U[k][j][i].d);
	
#ifdef MHD
	pGrid->U[k][j][i].E += (SQR(pGrid->U[k][j][i].B1c)
				+ SQR(pGrid->U[k][j][i].B2c)
				+ SQR(pGrid->U[k][j][i].B3c))*0.5;
#endif
      }
    }
  }

  
  /* free arrays allocated in initialize_fourier() */
  free_3d_array(dd);
  ath_3d_fft_free(fd);
  
  

  
  /*OLD VISCOSITY ROUTINE*/
  
/* #ifdef VISCOSITY */
/*   extern Real nu_iso, nu_aniso; */
/*   nu_aniso = 0.0; */
/*   Real dx = 1.0/(par_getd("domain1","Nx1")); */
/*   nu_iso =  par_getd_def("problem","nu_iso",0.0); // in units of dx*cs */
/*   nu_iso *= dx * sqrt(gm); */
/* #endif /\* VISCOSITY *\/ */


  
/* enroll new history variables, only once  */
  dump_history_enroll(hst_cstcool, "cstcool");
  dump_history_enroll(hst_tfloor, "tfloor");
  dump_history_enroll(hst_rho_hot, "rho_hot");
  dump_history_enroll(hst_rho_v_hot, "rho_v_hot");
  dump_history_enroll(hst_rho_cold, "rho_cold");
  dump_history_enroll(hst_rho_v_cold, "rho_v_cold");
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
        //integrate_cooling(pGrid);
#endif  /* BAROTROPIC */

	/*find nans*/

  int i, j, k;
  int is, ie, js, je, ks, ke;

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	
	
	if (isnan(pGrid->U[k][j][i].E)) {
	  ath_error("bad energy of %f \n", pGrid->U[k][j][i].E);
	}
	
	if (isnan(pGrid->U[k][j][i].d)) {
	  ath_error("bad density %f \n", pGrid->U[k][j][i].d);
	}
	
	
	if (isnan(pGrid->U[k][j][i].M1)) {
	  ath_error("bad momentum %f \n", pGrid->U[k][j][i].M1);
	}
	
	
      }
    }
  }
  
  
  
  
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
  int index = 0;

  // TODO: this isn't great, but this function isn't called all-that often
  while (x > geometric[index]) {
    index += 1;
  }

  return pow(4.0, (Real)(-index-1.0));

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



/* ========================================================================== */
/* history outputs
 */

Real hst_cstcool(const GridS *pG, const int i, const int j, const int k)
{
  return get_cstcool(pG->time);
}

Real hst_tfloor(const GridS *pG, const int i, const int j, const int k)
{
  return get_tfloor(get_cstcool(pG->time));
}


Real hst_rho_hot(const GridS *pG, const int i, const int j, const int k)
{
  Real tcut = 0.5;
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




Real hst_rho_v_hot(const GridS *pG, const int i, const int j, const int k)
{
  Real tcut = 0.5;
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



Real hst_rho_cold(const GridS *pG, const int i, const int j, const int k)
{
  Real tfloor = get_tfloor(get_cstcool(pG->time));
  Real tcut = 2.0 * tfloor;
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
  Real tfloor = get_tfloor(get_cstcool(pG->time));
  Real tcut = 2.0 * tfloor;
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



#ifdef VISCOSITY
static Real nu_fun(const Real d, const Real T, const Real x1, const Real x2, const Real x3){
  return 0.0;
  // return (nu_param * pow(T,2.5));
}
#endif  /* VISCOSITY */


#ifdef THERMAL_CONDUCTION
static Real kappa_fun(const Real d, const Real T,
                      const Real x1, const Real x2, const Real x3){
  return 0.0;
}
#endif  /* THERMAL_CONDUCTION */








/* ========================================================================== */
/* functions to make a random perturbation using FFTs */

static void initialize_fourier(GridS *pGrid, DomainS *pD)
{
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  int nbuf, mpierr, nx1gh, nx2gh, nx3gh;

  /* -----------------------------------------------------------
   * Variables within this block are stored globally, and used
   * within preprocessor macros.  Don't create variables with
   * these names within your function if you are going to use
   * OFST(), KCOMP(), or KWVM() within the function! */

  /* Get local grid size */
  nx1 = (ie-is+1);
  nx2 = (je-js+1);
  nx3 = (ke-ks+1);

  /* Get global grid size */
  gnx1 = pD->Nx[0];
  gnx2 = pD->Nx[1];
  gnx3 = pD->Nx[2];

  /* Get extents of local FFT grid in global coordinates */
  gis=is+pGrid->Disp[0];  gie=ie+pGrid->Disp[0];
  gjs=js+pGrid->Disp[1];  gje=je+pGrid->Disp[1];
  gks=ks+pGrid->Disp[2];  gke=ke+pGrid->Disp[2];
  /* ----------------------------------------------------------- */

  /* Get size of arrays with ghost cells */
  nx1gh = nx1 + 2*nghost;
  nx2gh = nx2 + 2*nghost;
  nx3gh = nx3 + 2*nghost;

  expo = par_getd("problem","expo");

  /* Cutoff wavenumbers of spectrum */
  klow  = par_getd("problem","klow");  /* in integer units */
  khigh = par_getd("problem","khigh"); /* in integer units */
  dkx   = 2.0*PI/(pGrid->dx1*gnx1);    /* convert k from integer */

  /* Allocate memory for components of velocity perturbation */
  dd = (Real***) calloc_3d_array(nx3gh, nx2gh, nx1gh, sizeof(Real));

  /* Initialize the FFT plan */
  plan = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);

  fd = ath_3d_fft_malloc(plan);

  return;
}


/*  Power spectrum returned in ampl
 *  - klow   = multiple of 2 pi/L for cut-off at low  wavenumbers
 *  - khigh  = multiple of 2 pi/L for cut-off at high wavenumbers
 *  - expo   = exponent of power law
 *  - ispect = integer flag which specifies spectrum
 *
 *  Note that the fourier amplitudes are stored in an array with no
 *  ghost zones
 */
static void pspect(ath_fft_data *ampl)
{
  int i,j,k;
  double q1,q2,q3;

  /* set random amplitudes with gaussian deviation */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        q1 = ran2(&rseed);
        q2 = ran2(&rseed);
        q3 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*PI*q2);
        q1 = ran2(&rseed);
        ampl[OFST(i,j,k)][0] = q3*cos(2.0*PI*q1);
        ampl[OFST(i,j,k)][1] = q3*sin(2.0*PI*q1);
      }
    }
  }

  /* set power spectrum */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        /* compute k/dkx */
        q3 = KWVM(i,j,k);
        if ((q3 > klow) && (q3 < khigh)) { /* decreasing power law */
          q3 *= dkx; /* multiply by 2 pi/L */
	  
          ampl[OFST(i,j,k)][0] /= pow(q3,(expo+2.0));
          ampl[OFST(i,j,k)][1] /= pow(q3,(expo+2.0));
        } else {                /* introduce cut-offs at klow and khigh */
          ampl[OFST(i,j,k)][0] = 0.0;
          ampl[OFST(i,j,k)][1] = 0.0;
        }
	
      }
    }
  }
  ampl[0][0] = 0.0;
  ampl[0][1] = 0.0;
  
  return;
}


static inline void transform()
{
  /* Transform velocities from k space to physical space */
  ath_3d_fft(plan, fd);
  
  /* Transform velocities from k space to physical space */
  
  /* Should technically renormalize (divide by gnx1*gnx2*gnx3) here,
     but since we're going to renormalize to get the desired energy
     injection rate anyway, there's no point */
  
  return;
}


static inline void generate()
{
  /* Generate new perturbations following appropriate power spectrum */
  pspect(fd);
  
  khigh *= 2.0;
  expo = 0.0;

  /* Generate new perturbations following appropriate power spectrum */


  /* Transform perturbations to real space, but don't normalize until
   * just before we apply them in perturb() */
  transform();

  return;
}


static void perturb(GridS *pGrid)
{
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  int ind, mpierr;
  Real rms[2], grms[2];
  
  /* copy perturbation from fourier array to dd (= "delta density") */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ind = OFST(i-is,j-js,k-ks);
        dd[k][j][i] = fd[ind][0];
      }
    }
  }


  /* normalize to have an RMS=1 */
  rms[0] = rms[1] = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        rms[0] += SQR(dd[k][j][i]);
        rms[1] += 1;
      }
    }
  }

#ifdef MPI_PARALLEL
  mpierr = MPI_Allreduce(rms, grms, 2, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[normalize]: MPI_Allreduce error = %d\n", mpierr);
  rms[0] = grms[0];
  rms[1] = grms[1];
#endif /* MPI_PARALLEL */

  rms[0] = sqrt(rms[0]/rms[1]);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        dd[k][j][i] /= rms[0];
      }
    }
  }


  /* transform with an arbitrary nonlinear function which yields
     large-ish density fluctuations, but which never cross zero */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d = exp(1.0 + dd[k][j][i]);

        /* flatten out high density stuff */
        if (pGrid->U[k][j][i].d > 1.0) {
          pGrid->U[k][j][i].d = 1.0 + log(pGrid->U[k][j][i].d);
        }

      }
    }
  }

  return;
}
/* -------------------------------------------------------------------------- */



/* ========================================================================== */
/* random number generator used in FFT code
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
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)
#define NTAB 32

/*! \fn double ran2(long int *idum){
 *  \brief The routine ran2() is extracted from the Numerical Recipes in C
 *
 * The routine ran2() is extracted from the Numerical Recipes in C
 * (version 2) code.  I've modified it to use doubles instead of
 * floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
* thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. */

double ran2(long int *idum){
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

#undef OFST
#undef KCOMP
#undef KWVM
/* end random number generator */
/* -------------------------------------------------------------------------- */
