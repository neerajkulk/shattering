#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"

#define REPORT_NANS

#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix);
static Real get_velocity_shift(MeshS *pM);
static void boost_frame(DomainS *pDomain, Real dvx);
static OutputS nan_dump;
static int nan_dump_count;
static Real netboost=0.0;
static Real t_boostdump=0.0;



#endif  /* REPORT_NANS */




#ifdef MPI_PARALLEL
#include "mpi.h"
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else /* DOUBLE_PREC */
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */

/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
 * For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis-nghost,gnx1))+ \
                            SQR(KCOMP(j,gjs-nghost,gnx2))+ \
                            SQR(KCOMP(k,gks-nghost,gnx3))))

/* FFTW - Variables, Plan, etc. */
static struct ath_3d_fft_plan *plan;
static ath_fft_data *fd=NULL;   /* unnormalized */
static Real ***dd=NULL;

static ath_fft_data *fv1=NULL, *fv2=NULL, *fv3=NULL;
static Real ***dv1=NULL, ***dv2=NULL, ***dv3=NULL;

/* Cutoff wavenumbers, G&O spect peak, power law spect exponent, 2 pi/L */
static Real klow,khigh,expo,dkx;
/* Number of cells in local grid, number of cells in global grid */
static int nx1,nx2,nx3,gnx1,gnx2,gnx3;
/* Starting and ending indices for global grid */
static int gis,gie,gjs,gje,gks,gke;
/* Seed for random number generator */
long int rseed;

/* Functions appear in this file in the same order that they appear in the
 * prototypes below */

/* Function prototypes for generating velocity perturbations */
static void pspect(ath_fft_data *ampl);
static inline void transform();
static inline void generate();
static void perturb(GridS *pGrid);

/* Function prototypes for initializing and interfacing with Athena */
static void initialize(GridS *pGrid, DomainS *pD);

/* Function prototypes for Numerical Recipes functions */
static double ran2(long int *idum);


/* global definitions for the SD cooling curve using the
   Townsend (2009) exact integration scheme */
#define nfit_cool 7

static Real sdT[nfit_cool] = {
  1.0e-5,
  0.0017235,
  0.02,
  0.13,
  0.7,
  5.0,
  100.0};

static const Real sdL[nfit_cool] = {
  5.890872e-13,
  15.438249,
  66.831473,
  2.773501,
  1.195229,
  1.842056,
  6.10541
};

static const Real sdexpt[nfit_cool] = {
  6.0,
  0.6,
  -1.7,
  -0.5,
  0.22,
  0.4,
  0.4
};

static Real Yk[nfit_cool];
/* -- end piecewise power-law fit */


/* must call init_cooling() in both problem() and read_restart() */
static void init_cooling();
static void test_cooling();

static Real sdLambda(const Real T);
static Real tcool(const Real d, const Real T);

static Real Y(const Real T);
static Real Yinv(const Real Y1);

static Real newtemp_townsend(const Real d, const Real T, const Real dt_hydro);

static void integrate_cooling(GridS *pG);

/* ========================================================================== */

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

  /* Should technically renormalize (divide by gnx1*gnx2*gnx3) here, but
   * since we're going to renormalize to get the desired energy injection
   * rate anyway, there's no point */

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

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ind = OFST(i-is,j-js,k-ks);
        dd[k][j][i] = fd[ind][0];
      }
    }
  }

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
  //ath_pout(0, "rms = %f\n", rms);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        dd[k][j][i] /= rms[0];
      }
    }
  }
    
    /*play with numbers here to flatten out high density areas and change size of perturbations */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d *= exp(1.0 + dd[k][j][i]);
          /*Flatten out high density stuff*/
          if (pGrid->U[k][j][i].d > 1.0) {
               pGrid->U[k][j][i].d = 1.0 + log(pGrid->U[k][j][i].d);
          }
          
    
      }
    }
  }

  return;
}


static void initialize(GridS *pGrid, DomainS *pD)
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
  klow = par_getd("problem","klow"); /* in integer units */
  khigh = par_getd("problem","khigh"); /* in integer units */
  dkx = 2.0*PI/(pGrid->dx1*gnx1); /* convert k from integer */

  /* Allocate memory for components of velocity perturbation */
  if ((dd=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
    ath_error("[problem]: Error allocating memory for vel pert\n");
  }

  /* Initialize the FFT plan */
  plan = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);

  fd = ath_3d_fft_malloc(plan);



  return;
}

void problem(DomainS *pDomain)
{
  GridS *pGrid = (pDomain->Grid);
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  int ixs,jxs,kxs;
  Real x1,x2,x3,r,theta,rc;
  Real lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Real ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  /* Ensure a different initial random seed for each process in an MPI calc. */
  rseed = -11;
#ifdef MPI_PARALLEL
  rseed -= myID_Comm_world;
#endif
  initialize(pGrid, pDomain);

//  init_cooling();
    
    Real P=1.0,vx=0.0,vy=0.0;
    

  /* Initialize uniform density */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
          cc_pos(pGrid,i,j,0,&x1,&x2,&x3);

        pGrid->U[k][j][i].d = 1.0;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = P/Gamma_1 + (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) + SQR(pGrid->U[k][j][i].M3))/(2.0*pGrid->U[k][j][i].d);
      }
    }
  }

  /* Set the initial perturbations.  Note that we're putting in too much
   * energy this time.  This is okay since we're only interested in the
   * saturated state. */
    
  generate();
  
  perturb(pGrid);
  
 /* this results in a grid with an average density of 2 with flucions between 0~5*/
    

  ath_pout(0,"De-allocating driving memory.\n");

  /* Free Athena-style arrays */
  free_3d_array(dd);

  /* Free FFTW-style arrays */
  ath_3d_fft_free(fd);



#ifdef REPORT_NANS
  nan_dump_count = 0;
#endif


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


/* ========================================================================== */

/*
 *  Function Userwork_in_loop
 *
 *  Drive velocity field for turbulence in GMC problems
 */

void Userwork_in_loop(MeshS *pM)
{
    int nl, nd, ntot;
    GridS *pGrid;
    
    int is, ie, js, je, ks, ke;
    int i,j,k;
    Real KE,TE,Temp,dE;
    Real FloorTemp;
    Real dvx;
    for (nl=0; nl<=(pM->NLevels)-1; nl++) {
        for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                pGrid = pM->Domain[nl][nd].Grid;
                
                
                
                is = pGrid->is; ie = pGrid->ie;
                js = pGrid->js; je = pGrid->je;
                ks = pGrid->ks; ke = pGrid->ke;
                
                FloorTemp = par_getd("problem", "Floor");
                
                dvx = get_velocity_shift(pM);
                
		


		/*  recover boost routine */
                
		  netboost+=dvx;

               /*
                
                if (pGrid->time >= t_boostdump) {
                    FILE *outfile;
                    
                    outfile = fopen("boost.dat", "a");
                    fprintf(outfile, "%.20f\t%f\n",
                            pGrid->time, netboost*pow(10,16));
                    fclose(outfile);
                    
                    Real dt_boost = par_getd("problem", "boostdump");
                    t_boostdump+= dt_boost;
                }
                
                */
                
      
                
                
                
        /* maximum velocity at each timestep*/
                if (pGrid->time >= t_boostdump) {
                    FILE *outfile;
                    outfile = fopen("velocity.dat", "a");
                    
                    
                    
                    Real maxvel=0.0;
                    
                    for(k=ks; k<=ke; k++){
                        for (j=js; j<=je; j++){
                            for (i=is; i<=ie; i++){
                                
                                if (maxvel < pGrid->U[k][j][i].M1/pGrid->U[k][j][i].d) {
                                    maxvel = pGrid->U[k][j][i].M1/pGrid->U[k][j][i].d;
                                    
                                
                                }
                                
                                
                                
                            }
                        }
                    }

                    
                    
                    fprintf(outfile, "%.20f\t%f\n",
                            pGrid->time, maxvel);
                    fclose(outfile);
                    
                    Real dt_boost = par_getd("problem", "boostdump");
                    t_boostdump+= dt_boost;
                }
                
                
		
                
                
                for(k=ks; k<=ke; k++){
                    for (j=js; j<=je; j++){
                        for (i=is; i<=ie; i++){
                            
                            /*start cooling routine*/
                            
                            if (pGrid->U[k][j][i].d != pGrid->U[k][j][i].d || pGrid->U[k][j][i].d <= 0.0)
                                ath_error("bad density of %f at cell (%d,%d,%d)\n",pGrid->U[k][j][i].d, k,j,i);
                            
                            
                            KE = (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)+ SQR(pGrid->U[k][j][i].M3)) / (2.0*pGrid->U[k][j][i].d);
                            

                            
                            TE = pGrid->U[k][j][i].E - KE;
                            
                            
                            Temp = ((2.0/3.0)*TE)/(pGrid->U[k][j][i].d);
                            
                            if (Temp != Temp || Temp <= 0.0){
                                ath_error("Floor temp at (%d,%d,%d), T = %f, KE = %f, E = %f, d = %f\n",
                                          k,j,i, Temp, KE, pGrid->U[k][j][i].E, pGrid->U[k][j][i].d);
                            }
                            
                            dE = SQR(pGrid->U[k][j][i].d)*(1.0/Temp)*pGrid->dt;
                            
                            TE -= dE;
                            
                            /* temperature floor */
                            
                            
                            if (TE < 3.0*pGrid->U[k][j][i].d*FloorTemp/2.0){
                                TE = 3.0*pGrid->U[k][j][i].d*FloorTemp/2.0;
                            }
                            
                            pGrid->U[k][j][i].E = KE + TE;
                            
                            /*end cooling routine*/
                            
                            
                            
                            
                            
                            
                        }
                    }
               }

		boost_frame(&(pM->Domain[nl][nd]), dvx);
                
                

                

                
                
                
            }
        }
    }
    
   
    return;
}

/* ========================================================================== */

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
        //if(ntot > 0)
        //report_nans(pM, &(pM->Domain[nl][nd]),1);
#endif
#ifdef INSTANTCOOL
        after_cool(pM, &(pM->Domain[nl][nd]),1);
#endif
      }
    }
  }

  return;
}

void problem_write_restart(MeshS *pM, FILE *fp)
{  return;  }

void problem_read_restart(MeshS *pM, FILE *fp)
{
  /* Ensure a different initial random seed for each process in an MPI calc. */
  rseed = -11;
#ifdef MPI_PARALLEL
  rseed -= myID_Comm_world;
#endif
  int nl, nd, ntot;

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {

        initialize(pM->Domain[nl][nd].Grid, &(pM->Domain[nl][nd]));
      }
    }
  }

  init_cooling();


  /* Free Athena-style arrays */
  free_3d_array(dd);

  /* Free FFTW-style arrays */
  ath_3d_fft_free(fd);


  /* Free Athena-style arrays */
  free_3d_array(dv1);
  free_3d_array(dv2);
  free_3d_array(dv3);

  /* Free FFTW-style arrays */
  ath_3d_fft_free(fv1);
  ath_3d_fft_free(fv2);
  ath_3d_fft_free(fv3);

  return;
}

/* ========================================================================== */

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


/* ================================================================ */
/* cooling routines */

static void init_cooling()
{
  int k, n=nfit_cool-1;
  Real term;
  const Real mu = 0.62, mu_e = 1.17;

  /* convert T in the cooling function from keV to code units */
  for (k=0; k<=n; k++)
    sdT[k] /= (8.197 * mu);

  /* populate Yk following equation A6 in Townsend (2009) */
  Yk[n] = 0.0;
  for (k=n-1; k>=0; k--){
    term = (sdL[n]/sdL[k]) * (sdT[k]/sdT[n]);

    if (sdexpt[k] == 1.0)
      term *= log(sdT[k]/sdT[k+1]);
    else
      term *= ((1.0 - pow(sdT[k]/sdT[k+1], sdexpt[k]-1.0)) / (1.0-sdexpt[k]));

    Yk[k] = Yk[k+1] - term;
  }

  return;
}

/* piecewise power-law fit to the cooling curve with temperature in
   keV and L in 1e-23 erg cm^3 / s */
static Real sdLambda(const Real T)
{
  int k, n=nfit_cool-1;

  /* first find the temperature bin */
  for(k=n; k>=0; k--){
    if (T >= sdT[k])
      break;
  }

  /* piecewise power-law; see equation A4 of Townsend (2009) */
  /* (factor of 1.311e-5 takes lambda from units of 1e-23 erg cm^3 /s
     to code units.) */
  return (1.311e-5 * sdL[k] * pow(T/sdT[k], sdexpt[k]));
}

static Real tcool(const Real d, const Real T)
{
  const Real mu = 0.62, mu_e = 1.17;

  /* equation 13 of Townsend (2009) */
  return (SQR(mu_e) * T) / (Gamma_1 * d * sdLambda(T));
}

/* see sdLambda() or equation A1 of Townsend (2009) for the
   definition */
static Real Y(const Real T)
{
  int k, n=nfit_cool-1;
  Real term;

  /* first find the temperature bin */
  for(k=n; k>=0; k--){
    if (T >= sdT[k])
      break;
  }

  /* calculate Y using equation A5 in Townsend (2009) */
  term = (sdL[n]/sdL[k]) * (sdT[k]/sdT[n]);

  if (sdexpt[k] == 1.0)
    term *= log(sdT[k]/T);
  else
    term *= ((1.0 - pow(sdT[k]/T, sdexpt[k]-1.0)) / (1.0-sdexpt[k]));

  return (Yk[k] + term);
}

static Real Yinv(const Real Y1)
{
  int k, n=nfit_cool-1;
  Real term;

  /* find the bin i in which the final temperature will be */
  for(k=n; k>=0; k--){
    if (Y(sdT[k]) >= Y1)
      break;
  }


  /* calculate Yinv using equation A7 in Townsend (2009) */
  term = (sdL[k]/sdL[n]) * (sdT[n]/sdT[k]);
  term *= (Y1 - Yk[k]);

  if (sdexpt[k] == 1.0)
    term = exp(-1.0*term);
  else{
    term = pow(1.0 - (1.0-sdexpt[k])*term,
               1.0/(1.0-sdexpt[k]));
  }

  return (sdT[k] * term);
}

static Real newtemp_townsend(const Real d, const Real T, const Real dt_hydro)
{
  Real term1, Tref;
  int n=nfit_cool-1;

  Tref = sdT[n];

  term1 = (T/Tref) * (sdLambda(Tref)/sdLambda(T)) * (dt_hydro/tcool(d, T));

  return Yinv(Y(T) + term1);
}

static void integrate_cooling(GridS *pG)
{
  int i, j, k;
  int is, ie, js, je, ks, ke;

  PrimS W;
  ConsS U;
  Real temp, tcloud = Gamma_1 / 1000;

  Real deltaE;
#ifdef MPI_PARALLEL
  Real deltaE_global;
  int ierr;
#endif  /* MPI_PARALLEL */

  /* ath_pout(0, "integrating cooling using Townsend (2009) algorithm.\n"); */

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  deltaE = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

        W = Cons_to_Prim(&(pG->U[k][j][i]));

        /* find temp in keV */
        temp = W.P/W.d;
        temp = newtemp_townsend(W.d, temp, pG->dt);

        /* apply a temperature floor (nans tolerated) */
        if (isnan(temp) || temp < tcloud)
          temp = tcloud;

        W.P = W.d * temp;
        U = Prim_to_Cons(&W);

        deltaE += pG->U[k][j][i].E - U.E;
        pG->U[k][j][i].E = U.E;
      }
    }
  }


#ifdef MPI_PARALLEL
  ierr = MPI_Allreduce(&deltaE, &deltaE_global, 1, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[report_nans]: MPI_Allreduce returned error %d\n", ierr);

  deltaE = deltaE_global;
#endif  /* MPI_PARALLEL */

  deltaE = deltaE / (gnx1*gnx2*gnx3);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += deltaE;
      }
    }
  }

  return;

}

static void test_cooling()
{
  int i, npts=100;
  Real logt, temp, tc, logdt, dt;
  Real err;

  FILE *outfile;

  outfile = fopen("lambda.dat", "w");
  for(i=0; i<npts; i++){
    logt = log(1.0e-3) + (log(5.0)-log(1.0e-3))*((double) i/(npts-1));
    temp = exp(logt);

    fprintf(outfile, "%e\t%e\n", temp, sdLambda(temp));
  }
  fclose(outfile);


  temp = 10.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-10kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }

  temp = 3.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-3kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }

  temp = 1.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-1kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }

  temp = 0.3;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-0.3kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }

  temp = 0.1;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-0.1kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }


  ath_error("check cooling stuff.\n");

  return;
}
/* end cooling routines */
/* ================================================================ */


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
  Real tfloor = 1.0e-6, tceil = 100.0;
  Real rhofloor = 1.0e-3;
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
        beta = fabs(200.0);
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
        } else if (temp > tceil) {
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



/* adding bits to energy and momentum when boosting frame */
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
                pGrid->U[k][j][i].E += 0.5 * d * SQR(dvx); /*second order term*/
                pGrid->U[k][j][i].E -= dvx * pGrid->U[k][j][i].M1; /* first order term */
#endif  /* ISOTHERMAL */
                pGrid->U[k][j][i].M1 -= dvx * d;
                
            }
        }
    }
    
    return;
}



static Real get_velocity_shift(MeshS *pM)
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
                        
                        /* s is some weighting factor... maybe dx1*d(ln rho)/dx */
                        s = 1.0;
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


