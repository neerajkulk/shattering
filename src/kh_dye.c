#include "copyright.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"


static Real hst_Sdye(const GridS *pG, const int i, const int j, const int k);

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

  Real rho, vx, vy, P;

#ifdef VISCOSITY
  Real reynolds_no;
#endif

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  /* Read problem parameters */
  vflow = par_getd("problem","vflow");
  drat  = par_getd("problem","drat");
  amp   = par_getd("problem","amp");

  a = 0.05;
  sigma = 0.2;
  z1 = 0.5;
  z2 = 1.5;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	if ( j<je/4 || j> 3*je/4 ) {
	  rho = 1.0
	}
	else {
	  rho = drat;
	}

        pGrid->U[k][j][i].d  = rho;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

#ifndef BAROTROPIC
        pGrid->U[k][j][i].E = P/Gamma_1
          + 0.5*(SQR(pGrid->U[k][j][i].M1)
                 + SQR(pGrid->U[k][j][i].M2)
                 + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* BAROTROPIC */

#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = 0.5*(tanh((x2-z2)/a) - tanh((x2-z1)/a) + 2.0) * rho;
#endif
      }
    }
  }

  dump_history_enroll(hst_Sdye, "dye entropy");

#ifdef VISCOSITY
  reynolds_no = par_getd("problem","Re");
  nu_iso   = 2.0 * vflow / reynolds_no; /* TODO: use L here */
  nu_aniso = 0.0;
#endif
#ifdef THERMAL_CONDUCTION
  prandtl_no = par_getd("problem", "Pr");
  kappa_iso = prandtl_no * nu_iso;
  kappa_aniso = 0.0;
#endif

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
  int nl, nd;
#ifdef VISCOSITY
  Real reynolds_no;
#endif
#ifdef THERMAL_CONDUCTION
  Real prandtl_no;
#endif
#if defined(THERMAL_CONDUCTION) && !defined(VISCOSITY)
  ath_error("[problem]: thermal conduction requires viscosity.\n");
#endif

  DomainS *pDomain;

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pDomain = &(pM->Domain[nl][nd]);
      }
    }
  }

  /* Read problem parameters */
  vflow = par_getd("problem","vflow");
  drat  = par_getd("problem","drat");

  dump_history_enroll(hst_Sdye, "dye entropy");

#ifdef VISCOSITY
  reynolds_no = par_getd("problem","Re");
  nu_iso   = 2.0 * vflow / reynolds_no; /* TODO: use L here */
  nu_aniso = 0.0;
#endif
#ifdef THERMAL_CONDUCTION
  prandtl_no = par_getd("problem", "Pr");
  kappa_iso = prandtl_no * nu_iso;
  kappa_aniso = 0.0;
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
