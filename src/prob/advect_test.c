#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

static Real gm;

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  Prim1DS W;
  Cons1DS U1d;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real pressure,drat,prat,rad,pa,da,x1,x2,x3;
  Real b0=0.0,Bx=0.0,rin;
  double theta;



  gm = par_getd_def("problem","gamma",1.0);

  drat = par_getd_def("problem","drat",1.0);

/* setup uniform ambient medium with spherical over-pressured region */

/*set up a blob moving at mach 0.5 moving diagonal to grid*/

  
  
  
  W.d = 1.0;
  W.Vx = 0.5*sqrt(gm/2.0);
  W.Vy = 0.5*sqrt(gm/2.0);
  W.Vz = 0.0;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	rad = sqrt(x1*x1 + x2*x2 + x3*x3);
		

	pGrid->U[k][j][i].d = 1.0;
	if (rad < 0.125) {pGrid->U[k][j][i].d = drat;
	  
	}

	pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*0.5*sqrt(gm/2.0);

	pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*0.5*sqrt(gm/2.0);

	
	pGrid->U[k][j][i].E = 1.0/(gm-1.0) + (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) + SQR(pGrid->U[k][j][i].M3))/(2.0*pGrid->U[k][j][i].d);

	
      }
    }
  }


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

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}
