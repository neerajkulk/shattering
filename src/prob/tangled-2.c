#include "copyright.h"
/*============================================================================*/
/*! \file field_loop.c
 *  \brief Problem generator for advection of a field loop test.
 *
 * PURPOSE: Problem generator for advection of a field loop test.  Can only
 *   be run in 2D or 3D.  Input parameters are:
 *   -  problem/rad   = radius of field loop
 *   -  problem/amp   = amplitude of vector potential (and therefore B)
 *   -  problem/vflow = flow velocity
 *   -  problem/drat  = density ratio in loop.  Enables density advection and
 *                      thermal conduction tests.
 *   The flow is automatically set to run along the diagonal.
 *
 *   Various test cases are possible:
 *   - (iprob=1): field loop in x1-x2 plane (cylinder in 3D)
 *   - (iprob=2): field loop in x2-x3 plane (cylinder in 3D)
 *   - (iprob=3): field loop in x3-x1 plane (cylinder in 3D)
 *   - (iprob=4): rotated cylindrical field loop in 3D.
 *   - (iprob=5): spherical field loop in rotated plane
 *
 *   A sphere of passive scalar can be added to test advection of scalars.
 *
 *   The temperature in the loop can be changed using drat to test conduction.
 *
 * REFERENCE: T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD
 *   via constrined transport", JCP, 205, 509 (2005)			      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MPI_PARALLEL
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */

#define REPORT_NANS

typedef struct SeedParticle_s{
  int label;
  Real x;
  Real y;
  Real z;
} SeedParticle;

static int nSeedParticles = 1000;
static SeedParticle *SeedParticles;

static OutputS SeedParticle_dump;

static void dump_seed_particles(MeshS *pM, OutputS *pOut);


#ifdef REPORT_NANS
static void report_nans(MeshS *pM, DomainS *pDomain);
static OutputS nan_dump;
static int nan_dump_count;
#endif  /* REPORT_NANS */

static void check_div_b(GridS *pGrid);

/* add a force-free term to B */
void add_term(Real3Vect ***A, GridS *pG,
              Real theta, Real phi,
              Real alpha, Real beta, Real amp);

Real randomreal(Real min, Real max);
static Real RandomNormal(Real mu, Real sigma);

/* unit vectors */
Real3Vect get_e1(Real theta, Real phi);
Real3Vect get_e2(Real theta, Real phi);
Real3Vect get_e3(Real theta, Real phi);



/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  Real x1, x2, x3, r;

  Real3Vect ***A;               /* vector potential */
  Real3Vect ***J;               /* current */

  int nterms;
  Real theta, phi, alpha, beta, amp;

  Real x1min, x1max, x2min, x2max, x3min, x3max, tmp;


  Real jmag, bmag, JdB, JcBforce, norm, Brms, Bs, Bmax, ncells;
  Real3Vect JcB;

#ifdef MPI_PARALLEL
  Real scal[5], my_scal[5];
  int ierr;
#endif  /* MPI_PARALLEL */

  srand(-4);                    /* SAME seed on all procs */

#ifdef REPORT_NANS
  nan_dump_count = 0;
#endif

  x1min = par_getd("domain1", "x1min"); x1max = par_getd("domain1", "x1max");
  x2min = par_getd("domain1", "x2min"); x2max = par_getd("domain1", "x2max");
  x3min = par_getd("domain1", "x3min"); x3max = par_getd("domain1", "x3max");

  SeedParticles = (SeedParticle*) calloc_1d_array(nSeedParticles, sizeof(SeedParticle));
  for (i=0; i<nSeedParticles; i++) {
    SeedParticles[i].label = i;

    if (myID_Comm_world == 0) {

      tmp = ((Real) rand())/RAND_MAX;
      SeedParticles[i].x = x1min + (x1max - x1min) * tmp;
      tmp = ((Real) rand())/RAND_MAX;
      SeedParticles[i].y = x2min + (x2max - x2min) * tmp;
      tmp = ((Real) rand())/RAND_MAX;
      SeedParticles[i].z = x3min + (x3max - x3min) * tmp;

      if (pGrid->Nx[2]==1)
        SeedParticles[i].z = 0.0;

      r = sqrt(SQR(SeedParticles[i].x)
               + SQR(SeedParticles[i].y)
               + SQR(SeedParticles[i].z));

      while (r > 0.25) {
        tmp = ((double) rand())/RAND_MAX;
        SeedParticles[i].x = x1min + (x1max - x1min) * tmp;
        tmp = ((double) rand())/RAND_MAX;
        SeedParticles[i].y = x2min + (x2max - x2min) * tmp;
        tmp = ((double) rand())/RAND_MAX;
        SeedParticles[i].z = x3min + (x3max - x3min) * tmp;

        if (pGrid->Nx[2]==1)
          SeedParticles[i].z = 0.0;

        r = sqrt(SQR(SeedParticles[i].x)
                 + SQR(SeedParticles[i].y)
                 + SQR(SeedParticles[i].z));
      }
    } else {
      /* particle is not on the grid.  set to zero and sum
         later */
      SeedParticles[i].x = 0.0;
      SeedParticles[i].y = 0.0;
      SeedParticles[i].z = 0.0;
    }
  }

#ifdef MPI_PARALLEL
  Real *particle_pos, *my_particle_pos;
  int ierr;

  particle_pos    = (Real*) calloc_1d_array(3*nSeedParticles, sizeof(Real));
  my_particle_pos = (Real*) calloc_1d_array(3*nSeedParticles, sizeof(Real));

  for (i=0; i<nSeedParticles; i++) {
    my_particle_pos[3*i+0] = SeedParticles[i].x;
    my_particle_pos[3*i+1] = SeedParticles[i].y;
    my_particle_pos[3*i+2] = SeedParticles[i].z;
  }

  ierr = MPI_Allreduce(my_particle_pos, particle_pos,
                       3*nSeedParticles, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[report_nans]: MPI_Allreduce returned error %d\n", ierr);

  for (i=0; i<nSeedParticles; i++) {
    SeedParticles[i].x = particle_pos[3*i+0];
    SeedParticles[i].y = particle_pos[3*i+1];
    SeedParticles[i].z = particle_pos[3*i+2];
  }

  free_1d_array((void*) particle_pos);
  free_1d_array((void*) my_particle_pos);
#endif  /* MPI_PARALLEL */


  SeedParticle_dump.n = 101;
  SeedParticle_dump.t = 0.0;
  SeedParticle_dump.dt = par_getd("output2", "dt"); /* careful!
                                                       assume this is
                                                       the vtk outptu */
  SeedParticle_dump.num = 0;
  SeedParticle_dump.out = "seed";


  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d = 1.0;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E = 1.0/Gamma_1;
#endif
      }
    }
  }

  A = (Real3Vect***) calloc_3d_array(nx3, nx2, nx1, sizeof(Real3Vect));

  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        A[k][j][i].x3 = A[k][j][i].x2 = A[k][j][i].x1 = 0.0;
      }
    }
  }

  nterms = par_getd("problem", "nterms");
  for (i=0; i<nterms; i++) {
    if (myID_Comm_world == 0) {
      theta = randomreal(0.0, 2.0*PI);
      phi   = randomreal(0.0, 2.0*PI);
      beta  = randomreal(0.0, 2.0*PI);

      amp   = RandomNormal(1.0, 0.25);
      alpha = 50.0;
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

  for (k=ks; k<=ke+2; k++) {
    for (j=js; j<je+2; j++) {
      for (i=is; i<ie+2; i++) {
        A[k][j][i].x1 /= sqrt(nterms);
        A[k][j][i].x2 /= sqrt(nterms);
        A[k][j][i].x3 /= sqrt(nterms);

        cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
        r = sqrt(x1*x1 + x2*x2 + x3*x3);
        if (r > 0.25)
          A[k][j][i].x1 = A[k][j][i].x2 = A[k][j][i].x3 = 0.0;
      }
    }
  }


  Brms = Bmax = 0.0;
  ncells = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ncells += 1.0;

        Bs = (SQR(A[k][j][i].x1) +
              SQR(A[k][j][i].x2) +
              SQR(A[k][j][i].x3));

        Brms += Bs;
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

  ath_pout(0, "Arms = %f\tAmax = %f\n", Brms, sqrt(Bmax));



  /* take curl here */
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        pGrid->B1i[k][j][i] = (A[k][j+1][i].x3 - A[k][j][i].x3)/pGrid->dx2 -
          (A[k+1][j][i].x2 - A[k][j][i].x2)/pGrid->dx3;
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][j][i] = (A[k+1][j][i].x1 - A[k][j][i].x1)/pGrid->dx3 -
          (A[k][j][i+1].x3 - A[k][j][i].x3)/pGrid->dx1;
      }
    }
  }

   for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[k][j][i] = (A[k][j][i+1].x2 - A[k][j][i].x2)/pGrid->dx1 -
          (A[k][j+1][i].x1 - A[k][j][i].x1)/pGrid->dx2;
      }
    }
  }

  free_3d_array((void***) A);
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

#ifdef MHD
  check_div_b(pGrid);
#endif MHD

  J = (Real3Vect***) calloc_3d_array(nx3, nx2, nx1, sizeof(Real3Vect));

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        J[k][j][i].x1 = (pGrid->B3i[k][j+1][i] - pGrid->B3i[k][j][i])/pGrid->dx2 -
          (pGrid->B2i[k+1][j][i] - pGrid->B2i[k][j][i])/pGrid->dx3;


        J[k][j][i].x2 = (pGrid->B1i[k+1][j][i] - pGrid->B1i[k][j][i])/pGrid->dx3 -
          (pGrid->B3i[k][j][i+1] - pGrid->B3i[k][j][i])/pGrid->dx1;


        J[k][j][i].x3 = (pGrid->B2i[k][j][i+1]-pGrid->B2i[k][j][i])/pGrid->dx1 -
          (pGrid->B1i[k][j+1][i]-pGrid->B1i[k][j][i])/pGrid->dx2;
      }
    }
  }

  JdB = JcBforce = 0.0;
  ncells = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ncells += 1.0;

        bmag = sqrt(SQR(pGrid->U[k][j][i].B1c) +
                    SQR(pGrid->U[k][j][i].B2c) +
                    SQR(pGrid->U[k][j][i].B3c));

        jmag = sqrt(SQR(J[k][j][i].x1) +
                    SQR(J[k][j][i].x2) +
                    SQR(J[k][j][i].x3));

        norm = bmag * jmag + TINY_NUMBER;

        if (bmag > 1.0e-6) {
          /* cross product */
          JcB.x1 =   J[k][j][i].x2 * pGrid->U[k][j][i].B3c
                   - J[k][j][i].x3 * pGrid->U[k][j][i].B2c;

          JcB.x2 =   J[k][j][i].x3 * pGrid->U[k][j][i].B1c
                   - J[k][j][i].x1 * pGrid->U[k][j][i].B3c;

          JcB.x3 =   J[k][j][i].x1 * pGrid->U[k][j][i].B2c
                   - J[k][j][i].x2 * pGrid->U[k][j][i].B1c;

          JcBforce += (SQR(JcB.x1) + SQR(JcB.x2) + SQR(JcB.x3)) / SQR(norm);


          /* dot product */
          JdB += SQR(J[k][j][i].x1 * pGrid->U[k][j][i].B1c +
                     J[k][j][i].x2 * pGrid->U[k][j][i].B2c +
                     J[k][j][i].x3 * pGrid->U[k][j][i].B3c) / SQR(norm);
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  my_scal[0] = JcBforce;
  my_scal[1] = ncells;
  my_scal[2] = JdB;

  ierr = MPI_Allreduce(&my_scal, &scal, 3, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[problem]: MPI_Allreduce returned error %d\n", ierr);

  JcBforce = scal[0];
  ncells   = scal[1];
  JdB      = scal[2];
#endif

  JcBforce = sqrt(JcBforce/ncells);
  JdB      = sqrt(JdB/ncells);

  ath_pout(0, "JxB = %f\tJ.D = %f\n", JcBforce, JdB);

  free_3d_array((void***) J);

  Brms = Bmax = 0.0;
  ncells = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ncells += 1.0;

        Bs = (SQR(pGrid->U[k][j][i].B1c) +
              SQR(pGrid->U[k][j][i].B2c) +
              SQR(pGrid->U[k][j][i].B3c));

        Brms += Bs;
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

void problem_write_restart(__attribute__((unused))MeshS *pM, __attribute__((unused))FILE *fp)
{
  return;
}

void problem_read_restart(__attribute__((unused))MeshS *pM, __attribute__((unused))FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(__attribute__((unused))const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(__attribute__((unused))const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  if (myID_Comm_world == 0 && pM->time >= SeedParticle_dump.t) {
    dump_seed_particles(pM, &SeedParticle_dump);
    SeedParticle_dump.t += SeedParticle_dump.dt;
    SeedParticle_dump.num += 1;
  }

  int nl, nd;
  GridS *pGrid;
  int i, j, k, p;
  Real d, vx, vy, vz, ignore;

  Real Brms=0.0, ncells=0.0;
#ifdef MPI_PARALLEL
  Real scal[2], my_scal[2];
  Real *particle_pos, *my_particle_pos;
  int ierr;
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


  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      pGrid = pM->Domain[nl][nd].Grid;
      if (pGrid != NULL) {
#ifdef MHD
        check_div_b(pGrid);
#endif

        for (k=pGrid->ks; k<=pGrid->ke; k++) {
          for (j=pGrid->js; j<=pGrid->je; j++) {
            for (i=pGrid->is; i<=pGrid->ie; i++) {
              ncells += 1.0;
              Brms += SQR(pGrid->U[k][j][i].B1c)
                + SQR(pGrid->U[k][j][i].B2c)
                + SQR(pGrid->U[k][j][i].B3c);
            }
          }
        }

      } /* if grid */
    } /* loop over domains */
  } /* loop over levels */

#ifdef MPI_PARALLEL
  my_scal[0] = Brms;
  my_scal[1] = ncells;

  ierr = MPI_Allreduce(&my_scal, &scal, 2, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[problem]: MPI_Allreduce returned error %d\n", ierr);

  Brms   = scal[0];
  ncells = scal[1];
#endif

  Brms = sqrt(Brms/ncells);
  ath_pout(0, "Brms = %f\n", Brms);


  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      pGrid = pM->Domain[nl][nd].Grid;
      if (pGrid != NULL) {
        for (p=0; p<nSeedParticles; p++) {
          /* is the particle on this grid? */
          if (SeedParticles[p].x >= pGrid->MinX[0] &&
              SeedParticles[p].x <= pGrid->MaxX[0] &&
              SeedParticles[p].y >= pGrid->MinX[1] &&
              SeedParticles[p].y <= pGrid->MaxX[1] &&
              SeedParticles[p].z >= pGrid->MinX[2] &&
              SeedParticles[p].z <= pGrid->MaxX[2]) {
            /* get the cell where the particle lives */
            celli(pGrid, SeedParticles[p].x, 1.0/(pGrid->dx1), &i, &ignore);
            cellj(pGrid, SeedParticles[p].y, 1.0/(pGrid->dx2), &j, &ignore);
            cellk(pGrid, SeedParticles[p].z, 1.0/(pGrid->dx3), &k, &ignore);

            /* get the velocity */
            d = pGrid->U[k][j][i].d;
            vx = pGrid->U[k][j][i].M1 / d;
            vy = pGrid->U[k][j][i].M2 / d;
            vz = pGrid->U[k][j][i].M3 / d;

            /* push the particle */
            SeedParticles[p].x += vx * pGrid->dt;
            SeedParticles[p].y += vy * pGrid->dt;
            SeedParticles[p].z += vz * pGrid->dt;
          } else {
            /* particle is not on the grid.  set to zero and sum
               later */
            SeedParticles[p].x = 0.0;
            SeedParticles[p].y = 0.0;
            SeedParticles[p].z = 0.0;
          }
        } /* loop over particle */
      }   /* grid != NULL */
    }     /* domains */
  } /* levels */

#ifdef MPI_PARALLEL
  particle_pos    = (Real*) calloc_1d_array(3*nSeedParticles, sizeof(Real));
  my_particle_pos = (Real*) calloc_1d_array(3*nSeedParticles, sizeof(Real));

  for (p=0; p<nSeedParticles; p++) {
    my_particle_pos[3*p+0] = SeedParticles[p].x;
    my_particle_pos[3*p+1] = SeedParticles[p].y;
    my_particle_pos[3*p+2] = SeedParticles[p].z;
  }

  ierr = MPI_Allreduce(my_particle_pos, particle_pos,
                       3*nSeedParticles, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[report_nans]: MPI_Allreduce returned error %d\n", ierr);

  for (p=0; p<nSeedParticles; p++) {
    SeedParticles[p].x = particle_pos[3*p+0];
    SeedParticles[p].y = particle_pos[3*p+1];
    SeedParticles[p].z = particle_pos[3*p+2];
  }

  free_1d_array((void*) particle_pos);
  free_1d_array((void*) my_particle_pos);
#endif  /* MPI_PARALLEL */


  return;
}

void Userwork_after_loop(__attribute__((unused))MeshS *pM)
{
  free_1d_array((void*) SeedParticles);

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

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
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

Real randomreal(Real min, Real max)
{
  Real eta = ((Real)rand()/(Real)RAND_MAX);
  return min + eta * (max-min);
}

static Real RandomNormal(Real mu, Real sigma)
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
      x1 = randomreal(-1.0, 1.0);
      x2 = randomreal(-1.0, 1.0);
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }

  return (mu + y1 * sigma);
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
  divb[0] *= (pGrid->dx1 * pGrid->dx2 * pGrid->dx3);

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

  Real tfloor    = 1.0e-2;
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

static void dump_seed_particles(MeshS *pM, OutputS *pOut)
{
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]);
  GridS   *pG = pD->Grid;

  FILE *pfile;
  char *fname;

  int p;

  long nout, my_id;
  int i,is,js,ks,h,init_id = 0;
  short pos;
  Real3Vect cell1;
  Real weight[3][3][3];         /* weight function */
  Real dpar,u1,u2,u3,cs;

  if((fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,
                        pOut->num,"seed","lis")) == NULL){
    ath_error("[dump_seed_particles]: Error constructing filename\n");
  }

  /* open output file */
  if((pfile = fopen(fname,"wb")) == NULL){
    ath_error("[dump_seed_particles]: Unable to open lis file %s\n",fname);
  }

  fprintf(pfile, "# time = %f\n", pG->time);

  for (p=0; p<nSeedParticles; p++) {
    fprintf(pfile, "%d\t%f\t%f\t%f\n",
            SeedParticles[p].label,
            SeedParticles[p].x,
            SeedParticles[p].y,
            SeedParticles[p].z);
  }


  fclose(pfile);

  return;
}
