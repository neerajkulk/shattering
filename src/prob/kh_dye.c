#include "copyright.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"


static Real dye_mass(const GridS *pG, const int i, const int j, const int k);

static Real vflow;
static Real drat;


/* ========================================================================== */
/* problem(): sets up the initial conditions */
void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k;
  int is,ie,js,je,ks,ke;
  Real amp,a,sigma,x1,x2,x3,z1,z2;

  Real rho, vx, vy, P, dye;

#ifdef VISCOSITY
  Real reynolds_no;
#endif

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  Real lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Real ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Real pi = 3.14159265358979323846;
  amp = lx/50;


  /* Read problem parameters */
  vflow = par_getd("problem","vflow");
  drat  = par_getd("problem","drat");

  P=1.0;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);


        if (x2<(-ly/4 + amp*sin(2.0*pi*x1/lx) ) || x2>(ly/4+amp*sin(2.0*pi*x1/lx))) {
          rho = 1.0;
          vx = vflow;
        }

	
	else {
	  rho = drat;
	  vx = 0.0;
       
	}

        pGrid->U[k][j][i].d  = rho;
        pGrid->U[k][j][i].M1 = rho*vx;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

#ifndef BAROTROPIC
        pGrid->U[k][j][i].E = P/Gamma_1
          + 0.5*(SQR(pGrid->U[k][j][i].M1)
                 + SQR(pGrid->U[k][j][i].M2)
                 + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* BAROTROPIC */

#if (NSCALARS > 0)
	if ( x2<-ly/4 || x2>ly/4 ) {
	  dye = 0.0;
	}
	else {
	  dye = 1.0;
	}

        pGrid->U[k][j][i].s[0] = dye * rho;
#endif
      }
    }
  }

  dump_history_enroll(dye_mass, "dye entropy");

  return;
}
/* ========================================================================== */


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
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

static Real dye_mass(const GridS *pG, const int i, const int j, const int k)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;

  int ii;
  int ji;
  int ki;

  Real msum = 0.0;
  
  for (ki=ks; ki<=ke; ki++) {
    for (ji=js; ji<=je; ji++) {
      for (ii=is; ii<=ie; ii++) {

	msum += pG->U[ki][ji][ii].s[0];
	      
      }
    }
  }



  return msum;
}

