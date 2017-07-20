#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

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

void problem(DomainS *pDomain)
{
    GridS *pGrid=(pDomain->Grid);
    int i, is = pGrid->is, ie = pGrid->ie;
    int j, js = pGrid->js, je = pGrid->je;
    int k, ks = pGrid->ks, ke = pGrid->ke;
    Real x1,x2,x3;
    Real lx,ly;
    
    
    /* size of the domain (in physical coordinates) */
    lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
    ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
    
    Real p;
    Real dbig = 10.0;
    Real dsmal = 0.1;
    Real n = 2.0, amp=lx/50.0, thickness=lx/50.0;
    Real a = 2.0*atanh(0.9)/thickness;
    


    
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
                

                
                Real density = tanh(a*(x1-amp*sin(n*PI*x2/ly)))*(dbig-dsmal)/2.0 + (dbig+dsmal)/2.0;
                
                p = 1.0 - 0.3*(2.0*x1/lx);
                
                
                
                pGrid->U[0][j][i].d = density;
                pGrid->U[0][j][i].M1 = 0.0;
                pGrid->U[0][j][i].M2 = 0.0;
#ifndef ISOTHERMAL
                pGrid->U[0][j][i].E = p/Gamma_1;
                
#endif
            }
        }
    }
    
    
    return;
    
}
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
    int nl, nd, ntot;
    Real dvx;
    
    dvx = get_velocity_shift(pM);
    
    for (nl=0; nl<=(pM->NLevels)-1; nl++) {
        for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                boost_frame(&(pM->Domain[nl][nd]), dvx);
            }
        }
    }
    
    return;
}

void Userwork_after_loop(MeshS *pM)
{
}

/*---------------------------------------------------------------------------------- special functions*/


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
                pGrid->U[k][j][i].E += 0.5 * d * SQR(dvx);
                pGrid->U[k][j][i].E -= dvx * pGrid->U[k][j][i].M1;
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






