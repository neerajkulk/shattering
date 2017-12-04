#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"


/* ========================================================================== */
/* global variables and prototypes */

#ifdef MPI_PARALLEL
#include "mpi.h"
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else /* DOUBLE_PREC */
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */
/* -------------------------------------------------------------------------- */


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


/* -------------------------------------------------------------------------- */
/* functions for history outputs
 */
#if (NSCALARS > 1)
static Real hst_hot_mom(const GridS *pG, const int i, const int j, const int k);
#endif  /* NSCALARS */
#if (NSCALARS > 0)
static Real hst_cold_mom(const GridS *pG, const int i, const int j, const int k);
#endif  /* NSCALARS */
static Real hst_Sdye(const GridS *pG, const int i, const int j, const int k);
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/* global variables used throughout the simulation
 */


  #ifdef VISCOSITY
/* viscosity won't work without conduction... */
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3);
#endif  /* VISCOSITY */
 
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/* simple cooling integrator
 */

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
static void set_vars(DomainS *pDomain);
static Real window(Real x, Real width, Real a);

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



/* ========================================================================== */
/* main problem function -- sets initial conditions for the simulation
 */

void problem(DomainS *pDomain)
{
  GridS *pGrid = (pDomain->Grid);

  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;

  Real x1, x2, x3, r;
  Real lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Real ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  const Real P = 1.0, vy = 0.0;
  Real noise = par_getd("problem", "noise");
   
#ifdef MHD
  const Real Bx = 0.01;
#endif  /* MHD */


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
  Real  end_time     = par_getd("time", "tlim");

  Real vflow = mach * sqrt(gm);

  Real drat = 1.0 / get_tfloor(get_cstcool(0.0));

  ath_pout(0, "drat =  %f\n", drat);

  
  Prim1DS W;
  Cons1DS U1d;


  /* initialize global variables in a separate function so it can be
     done identically here and in read_restart()  */
  set_vars(pDomain);

  initialize_fourier(pGrid, pDomain);
  generate();


  /* this sets the grid density field */
  perturb(pGrid);


  /* start with a uniform grid */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	
        pGrid->U[k][j][i].d *= (noise/2.0);
        pGrid->U[k][j][i].d += 1.0 - noise;

	/* density perturbations are now normalized to have a mean of ~1 with drho/rho ~ noise */
	
	//set up hot cold kh beam
	pGrid->U[k][j][i].d *= 1.0 + (drat-1) * window(x2, width, a);
	pGrid->U[k][j][i].M1 = (vflow * (1.0 - window(x2, width, a)))*pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
	

	
	
#ifdef MHD
        pGrid->U[k][j][i].B1c = Bx;
        pGrid->U[k][j][i].B2c = 0.0;
        pGrid->U[k][j][i].B3c = 0.0;

        pGrid->B1i[k][j][i] = Bx;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = 0.0;

        if (i == ie && ie > is) pGrid->B1i[k][j][i+1] = Bx;
#endif /* MHD */
	
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


static void set_vars(DomainS *pDomain)
{
#ifdef VISCOSITY
  NuFun_i = nu_fun;
  NuFun_a = NULL;
#endif  /* VISCOSITY */


  mach         = par_getd("problem", "mach");
  alpha        = par_getd("cooling", "alpha");
  end_time     = par_getd("time", "tlim");
  n_kh         = par_getd("problem", "n_kh");
  gm           = par_getd("problem", "gamma");
  
  num_steps    = par_geti("cooling", "steps");
  
  cooling_flag = par_geti("cooling", "cooling");
  heating_flag = par_geti("cooling", "heating");
  

  
  dump_history_enroll(hst_Sdye,     "dye entropy");
#if (NSCALARS > 0)
  dump_history_enroll(hst_cold_mom, "cold momentum");
#endif  /* NSCALARS */
#if (NSCALARS > 1)
  dump_history_enroll(hst_hot_mom,  "hot momentum");
#endif  /* NSCALARS */



#ifdef REPORT_NANS
  nan_dump_count = 0;
#endif  /* REPORT_NANS */

}

/* end problem() */
/* -------------------------------------------------------------------------- */



/* ========================================================================== */
/* history outputs
 */
#if (NSCALARS > 1)
static Real hst_hot_mom(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[1] * pG->U[k][j][i].M1);
}
#endif  /* NSCALARS */


#if (NSCALARS > 0)
static Real hst_cold_mom(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M1);
}
#endif  /* NSCALARS */


static Real hst_Sdye(const GridS *pG, const int i, const int j, const int k)
{
  Real entropy;
  entropy = (-1.0)*(pG->U[k][j][i].d)*log(pG->U[k][j][i].d);
  return entropy;
}

/* end history outputs */
/* -------------------------------------------------------------------------- */



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
  int nl, nd, ntot;
  DomainS *pDomain;

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {

        pDomain = &(pM->Domain[nl][nd]);

        set_vars(pDomain);
      }
    }
  }

  return;
}


ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}


VOutFun_t get_usr_out_fun(const char *name)
{
  return NULL;
}


void Userwork_in_loop(MeshS *pM)
{
  int nl, nd, ntot;
  GridS *pGrid;

  gnx1 = pM->Nx[0];
  gnx2 = pM->Nx[1];
  gnx3 = pM->Nx[2];

  /* loop over the mesh to find the grid on this processor */
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pGrid = pM->Domain[nl][nd].Grid;

#ifndef BAROTROPIC
        //integrate_cooling(pGrid);
#endif  /* BAROTROPIC */
      }
    }
  }

  return;
}


void Userwork_after_loop(MeshS *pM)
{
  /* Don't free memory here if doing any analysis because final
   * output hasn't been written yet!! */
  return;
}


void Userwork_before_loop(MeshS *pM)
{
  /* Don't free memory here if doing any analysis because final
   * output hasn't been written yet!! */
  int nl, nd, ntot;

  /* report nans first, so we can fix them before they propagate into
     the following functions. */
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef REPORT_NANS
        ntot = report_nans(pM, &(pM->Domain[nl][nd]),1);
#endif  /* REPORT_NANS */
      }
    }
  }

  return;
}

/* end user functions */
/* -------------------------------------------------------------------------- */



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

/* end of shearflow.c */





#ifdef VISCOSITY
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3)
{				/* simple test of temperature depedant viscosity  */
  return 0.0;
  // return (nu_param * pow(T,2.5));
}
#endif  /* VISCOSITY */


#ifdef THERMAL_CONDUCTION
static Real kappa_fun(const Real d, const Real T,
                      const Real x1, const Real x2, const Real x3)
{
  return 0.0;
}
#endif  /* THERMAL_CONDUCTION */
 
