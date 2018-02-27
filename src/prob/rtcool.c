#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"



/* Functions and variables for boosting frame  */
static Real netboost=0.0;
static Real t_boostdump=0.0;
static Real get_velocity_shift(MeshS *pM);
static void boost_frame(DomainS *pDomain, Real dvx);




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
  Real x1,x2,x3;
  Real lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Real ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  
  /* Ensure a different initial random seed for each process in an MPI calc. */
  rseed = -11;
#ifdef MPI_PARALLEL
  rseed -= myID_Comm_world;
#endif
  initialize(pGrid, pDomain);
    
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

  /* Set the initial perturbations.  Now apply fourier noise*/
    
  generate();
  
  perturb(pGrid);
  
  /* Now, we have fourier noise in density with an average density of ~2.0286. The range of density values values go from 0.128823<rho<4.7754. We need to tweak these density values for our problem*/



  /* Define additional parameters to write initial condition  */
 
  Real rhoh = 1.0;  //density of hot gas
  Real alpha = par_getd("problem","alpha"); // lambda = lambda_0 * T^alpha
  Real drat = par_getd("problem","drat"); // density ratio between hot and cold.
  Reah rhoc = drat; // density of cold gas

  /* Set cstcool in the cold gas to be 1 for all simulations:  cs*tcool = T^(2.5-alpha)/(f*lambda_0) 
     To change L/cstcool, we change the domain size of our simulation in the athinput file. 
  */

 
  Real f = pow(gm, 1.5)/(gm-1.0)/(2.0-alpha);
  f = 1.0/f;
  
  Real lambda_0 = 1.0/f;
  lambda_0 *= 1/pow(rhoc,2.5 - alpha);

  
  Real n = 2.0, amp=lx/50.0, thickness=lx/50.0;
  Real a = 2.0*atanh(0.9)/thickness;
  Real vb = par_getd("problem","boost");
  Real noise = par_getd("problem","noise");
 
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,0,&x1,&x2,&x3);
        
        
	/* Tweak density values to have a mean of 1 and fluctions given by the variable 'noise' (read from athinput file) */
	
	/* Step 1: Center density values around 0 (density values can temporarily be negative) */
	
	pGrid->U[k][j][i].d -= 2.0286;

	/* Step 2: Scale range of densities according to the variable 'noise' . (dmax-dmin)/2 = 2.32329*/
	
	pGrid->U[k][j][i].d = pGrid->U[k][j][i].d * (noise/2.32329);
	
	/* Step 3: Shift density profile so that mean is centered around 1.  */
	pGrid->U[k][j][i].d += 1.0;
	
	/* Now add interface between hot and cold gas */
        
	pGrid->U[k][j][i].d *= tanh(a*(x1-amp*sin(n*PI*x2/ly)))*(dbig-dsmal)/2.0 + (dbig+dsmal)/2.0;
	
        
	pGrid->U[k][j][i].E = P/Gamma_1 + (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2))/(2.0*pGrid->U[k][j][i].d);
        
        
        
	
      }
    }
  }
  
  
  
  
  
  ath_pout(0,"De-allocating driving memory.\n");
  
  /* Free Athena-style arrays */
  free_3d_array(dd);
  
  /* Free FFTW-style arrays */
  ath_3d_fft_free(fd);
  
  



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


