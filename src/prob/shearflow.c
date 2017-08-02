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
static OutputS nan_dump;
static int nan_dump_count;
#endif  /* REPORT_NANS */


static Real hotm1(const GridS *pG, const int i, const int j, const int k);
static Real coldm1(const GridS *pG, const int i, const int j, const int k);
static Real hst_Sdye(const GridS *pG, const int i, const int j, const int k);
static Real color(const GridS *pG, const int i, const int j, const int k);

static Real FloorTemp;
static Real CeilingTemp;
static Real coolon;
static Real heaton;
static Real tsim;
static int stepcooling;
static Real steps;
static Real coolinglaw;

static Real netboost=0.0;
static Real t_boostdump=0.0;

extern Real nu_iso, nu_aniso;


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
  Real x1,x2,x3;
  Real lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Real ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  Real vx;
  Real n=2.0, amp=lx/50.0, thickness=lx/50.0;
  Real a = 2.0*atanh(0.9)/thickness;
  Real vb = par_getd("problem","boost");
  Real noise = par_getd("problem","noise");
  Real width = par_getd("problem","frac_cold")*(ly/2.0);
  Real intrfc;
  Real v0 = par_getd("problem","v0");
  Real dbig = par_getd("problem","dbig");
  Real dsmal = par_getd("problem","dsmal");
  Real r; /*r in cylindrical coordinates*/

  FloorTemp = par_getd("problem", "Floor");
  CeilingTemp = par_getd("problem", "Ceiling");
  coolon = par_getd("problem", "coolon");
  heaton = par_getd("problem", "heaton");
  tsim = par_getd("time", "tlim");
  stepcooling = par_geti("problem", "stepcooling");
  steps = par_getd("problem", "steps");
  coolinglaw = par_getd("problem", "coolinglaw");
  // ath_error("temp floor is %f\n", FloorTemp);


  Prim1DS W;
  Cons1DS U1d;

  /* Ensure a different initial random seed for each process in an MPI calc. */
  rseed = -11;
#ifdef MPI_PARALLEL
  rseed -= myID_Comm_world;
#endif
  initialize(pGrid, pDomain);

  //  init_cooling();

  Real P=1.0,vy=0.0;

#ifdef MHD
  Real Bx = 0.01;
#endif


  /* Initialize uniform density */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,0,&x1,&x2,&x3);

        pGrid->U[k][j][i].d = 1.0;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

#ifdef MHD
        pGrid->U[k][j][i].B1c = 0.01;
        pGrid->U[k][j][i].B2c = 0.0;
        pGrid->U[k][j][i].B3c = 0.0;

        pGrid->B1i[k][j][i] = Bx;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = 0.0;

        if (i == ie && ie > is) pGrid->B1i[k][j][i+1] = Bx;
#endif /*MHD*/

        pGrid->U[k][j][i].E = P/Gamma_1 + (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) + SQR(pGrid->U[k][j][i].M3))/(2.0*pGrid->U[k][j][i].d);

#ifdef MHD
        pGrid->U[k][j][i].E += (SQR(pGrid->U[k][j][i].B1c)
                                + SQR(pGrid->U[k][j][i].B2c)
                                + SQR(pGrid->U[k][j][i].B3c))*0.5;
#endif
      }
    }
  }

  /* Set the initial perturbations.  Note that we're putting in too much
   * energy this time.  This is okay since we're only interested in the
   * saturated state. */

  generate();

  perturb(pGrid);

  /* this results in a grid with an average density of 2 with flucions between 0~5*/

  /* divide by 2 so avg density of cool side is 1 w/ density perturbations by a factor of 2.5 ... then add a smooth interface.*/

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

        r = sqrt(x3*x3 + x2*x2);

        intrfc = 0.5*(dbig-dsmal)*(tanh((-r+width)*a)+tanh((r+width)*a))+dsmal;
        vx = v0 - 0.5*v0*(tanh((-r+width)*a)+tanh((r+width)*a));


        pGrid->U[k][j][i].d /= 2.0;

        /* At this stage we have a density distribuiton with a mean value of 1 with a range of 0-2.5. rescale to given noise.*/

        pGrid->U[k][j][i].d /= 1.0/noise;

        pGrid->U[k][j][i].d += 1.0 - noise;

        /* add hot-cold interface */

        pGrid->U[k][j][i].d *= intrfc;

        /*shift stuff above the floor temp*/

        pGrid->U[k][j][i].d /= (1.0 + 1.8*noise) ;


        pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*vx;

        pGrid->U[k][j][i].E = P/Gamma_1 + (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2))/(2.0*pGrid->U[k][j][i].d);


#ifdef MHD
        pGrid->U[k][j][i].E += (SQR(pGrid->U[k][j][i].B1c)
                                + SQR(pGrid->U[k][j][i].B2c)
                                + SQR(pGrid->U[k][j][i].B3c))*0.5;
#endif

#if (NSCALARS > 0)
        /*cold gas dye*/
        pGrid->U[k][j][i].s[0] =  0.5*(tanh((-r+width)*a)+tanh((r+width)*a))*pGrid->U[k][j][i].d;

        /*hot gas dye*/
        pGrid->U[k][j][i].s[1] = pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0] ;

        /*gas is initialized by 2 dyes to be either hot or cold*/

#endif

        /* s = ρc. c is the specific dye (roughly 0-1 specifying what fraction of the gas is cold) the volume integral of s is conserved. */

        /*passively advected scalar with 1 in cold gas and 0 in hot gas*/

      }
    }
  }


  dump_history_enroll(hst_Sdye, "dye entropy");
  dump_history_enroll(hotm1, "hot momentum");
  dump_history_enroll(coldm1, "cold momentum");

#ifdef VISCOSITY

  Real reynolds;

  reynolds  = par_getd_def("problem","reynolds",0.0);


  nu_iso = lx * v0 / reynolds;
  nu_aniso = 0.0;

#endif // VISCOSITY





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
  DomainS *pD;
  int is, ie, js, je, ks, ke;
  int i,j,k;
  Real KE,TE,ME,temp,dE;
  Real dvx;
  Real deltaE;
  int gnx1, gnx2, gnx3;
  Real counter, tshift;
  Real tmpfloor;
#ifdef MPI_PARALLEL
  Real deltaE_global;
  int ierr;
#endif /*MPI_PARALLEL*/

  gnx1 = pM->Nx[0];
  gnx2 = pM->Nx[1];
  gnx3 = pM->Nx[2];

  /* loop over the mesh to find the grid on this processor */
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pGrid = pM->Domain[nl][nd].Grid;

        deltaE = 0;
        is = pGrid->is; ie = pGrid->ie;
        js = pGrid->js; je = pGrid->je;
        ks = pGrid->ks; ke = pGrid->ke;

        tmpfloor = FloorTemp;
        counter = 1.0;
        tshift = tsim/steps;

        for(k=ks; k<=ke; k++){
          for (j=js; j<=je; j++){
            for (i=is; i<=ie; i++){

              /* start cooling routine */
              if (pGrid->U[k][j][i].d != pGrid->U[k][j][i].d
                  || pGrid->U[k][j][i].d <= 0.0){
                ath_error("bad density of %f at cell (%d,%d,%d)\n",
                          pGrid->U[k][j][i].d, k,j,i);
              }

              KE = (SQR(pGrid->U[k][j][i].M1)
                    + SQR(pGrid->U[k][j][i].M2)
                    + SQR(pGrid->U[k][j][i].M3)) / (2.0*pGrid->U[k][j][i].d);


#ifdef MHD

              ME =  (SQR(pGrid->U[k][j][i].B1c)
                     + SQR(pGrid->U[k][j][i].B2c)
                     + SQR(pGrid->U[k][j][i].B3c))*0.5;


#endif // MHD

              TE = pGrid->U[k][j][i].E - KE;

#ifdef MHD
              TE -= ME;
#endif // MHD

              temp = ((2.0/3.0)*TE)/(pGrid->U[k][j][i].d);

              if (temp != temp || temp <= 0.0){
                ath_error("Floor temp at (%d,%d,%d), T = %f, KE = %f, E = %f, d = %f\n",k,j,i, temp, KE, pGrid->U[k][j][i].E,pGrid->U[k][j][i].d);
              }

              dE = SQR(pGrid->U[k][j][i].d)*(pow(temp,coolinglaw))*pGrid->dt;

              if (coolon == 1.0) {
                TE -= dE;
              }


              if (stepcooling ==1) {

                if (pGrid->time > tshift) {
                  tmpfloor = tmpfloor*(steps - counter)/steps;
                  tshift += tsim/steps;
                  counter += 1.0;
                }
              }




              if (stepcooling == 2) {

                if (pGrid->time > tshift) {
                  tmpfloor /= pow(4.0,(1.0/(2.5-coolinglaw)));  /*cstcool drops by a factor of 4*/
                  tshift += tsim/steps;
                  counter += 1.0;
                }

              }


              /* apply temperature floor and write to the grid */
              if (TE < 3.0*pGrid->U[k][j][i].d*tmpfloor/2.0){
                TE = 3.0*pGrid->U[k][j][i].d*tmpfloor/2.0;
              }


              /* apply temperature ceiling and write to the grid */
              if (TE > 3.0*pGrid->U[k][j][i].d*CeilingTemp/2.0){
                TE = 3.0*pGrid->U[k][j][i].d*CeilingTemp/2.0;
              }
              deltaE += pGrid->U[k][j][i].E - (KE + TE);
              // ath_pout(0, "temp floor is %f\n", FloorTemp);
#ifdef MHD

              deltaE -= ME;

#endif // MHD

              pGrid->U[k][j][i].E = KE + TE;



#ifdef MHD

              pGrid->U[k][j][i].E += ME;

#endif // MHD

            }
          }
        }

        /* heating to match cooling */

#ifdef MPI_PARALLEL
        ierr = MPI_Allreduce(&deltaE, &deltaE_global, 1, MPI_RL, MPI_SUM,
                             MPI_COMM_WORLD);
        if (ierr){
          ath_error("[Userwork_in_loop]: MPI_Allreduce returned error %d\n",
                    ierr);
        }

        deltaE = deltaE_global;
#endif

        /* MPI_PARALLEL */

        deltaE = deltaE/(gnx1*gnx2*gnx3);

        for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {

              if (heaton == 1.0) {
                pGrid->U[k][j][i].E += deltaE;

              }

            }
          }
        }


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
  GridS *pGrid;
  DomainS *pDomain;
  int is, ie, js, je, ks, ke;
  int i,j,k;
  int ixs,jxs,kxs;
  Real x1,x2,x3;
  Real lx,ly;
  Real v0;
  Real r; /*r in cylindrical coordinates*/
  Prim1DS W;
  Cons1DS U1d;
  Real reynolds;


  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {

        pDomain = &(pM->Domain[nl][nd]);

        pGrid = pM->Domain[nl][nd].Grid;
        is = pGrid->is, ie = pGrid->ie;
        js = pGrid->js, je = pGrid->je;
        ks = pGrid->ks, ke = pGrid->ke;
        lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
        ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];


        FloorTemp = par_getd("problem", "Floor");
        CeilingTemp = par_getd("problem", "Ceiling");
        coolon = par_getd("problem", "coolon");
        heaton = par_getd("problem", "heaton");
        tsim = par_getd("time", "tlim");
        stepcooling = par_geti("problem", "stepcooling");
        steps = par_getd("problem", "steps");
        coolinglaw = par_getd("problem", "coolinglaw");

        dump_history_enroll(hst_Sdye, "dye entropy");
        dump_history_enroll(hotm1, "hot momentum");
        dump_history_enroll(coldm1, "cold momentum");

#ifdef VISCOSITY
        reynolds  = par_getd_def("problem","reynolds",0.0);
        v0  = par_getd_def("problem","v0",0.0);
        nu_iso = lx * v0 / reynolds;
        nu_aniso = 0.0;
#endif // VISCOSITY

      }
    }
  }


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

  /* Real tfloor    = 1.0e-2 / drat; */
  /* Real tceil     = 100.0; */
  /* Real rhofloor  = 1.0e-2; */
#ifdef MHD
  Real betafloor = 3.0e-3;
#endif
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





static Real hotm1(const GridS *pG, const int i, const int j, const int k)
{
  return  (pG->U[k][j][i].s[1] * pG->U[k][j][i].M1)/( pG->U[k][j][i].d);
}

static Real coldm1(const GridS *pG, const int i, const int j, const int k)
{
  return  (pG->U[k][j][i].s[0] * pG->U[k][j][i].M1)/( pG->U[k][j][i].d);
}

static Real hst_Sdye(const GridS *pG, const int i, const int j, const int k)
{
  Real entropy;
  entropy = (-1.0)*(pG->U[k][j][i].d)*log(pG->U[k][j][i].d);
  return entropy;
}
