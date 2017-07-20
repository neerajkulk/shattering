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

void problem(DomainS *pDomain)
{
    GridS *pGrid=(pDomain->Grid);
    int i, is = pGrid->is, ie = pGrid->ie;
    int j, js = pGrid->js, je = pGrid->je;
    int k, ks = pGrid->ks, ke = pGrid->ke;
    Real x1,x2,x3;
    Real density, pressure, vel, lx, ly;
    Real boost = par_getd("problem", "boost");
    
    /* size of the domain (in physical coordinates) */
    lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
    ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
    
    
    
    
    /* setting up background with initial conditions of 10^6 K gas and 10^4 K gas */
    
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
                
                if (x1 > 0.0) {
                    pressure = 0.01;
                } else {
                    pressure = 0.5;
                    
                }
                
                vel = boost;
                
                if (x1 > 0.0) {
                    density= 1.0;
                } else {
                    density= 1.0;
                    
                }
                
                
                
                pGrid->U[0][j][i].d = density;
                pGrid->U[0][j][i].M1 = density*vel;
                pGrid->U[0][j][i].M2 = 0.0;
#ifndef ISOTHERMAL
                pGrid->U[0][j][i].E = pressure/Gamma_1 + (SQR(pGrid->U[0][j][i].M1) + SQR(pGrid->U[0][j][i].M2))/(2.0*(pGrid->U[0][j][i].d));
                
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
    GridS *pGrid;
    int is, ie, js, je, ks, ke;
    int i, j, k;
    Real x1,x2,x3;
    Real t_next_velcocity_boost = 0.0;
    Real acceleration = par_getd("problem", "acceleration");
    Real dv;




    
    for (nl=0; nl<=(pM->NLevels)-1; nl++) {
        for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                pGrid = pM->Domain[nl][nd].Grid;
                
                dv = acceleration*(pGrid->dt);
                
                is = pGrid->is; ie = pGrid->ie;
                js = pGrid->js; je = pGrid->je;
                ks = pGrid->ks; ke = pGrid->ke;
                
                FILE *outfile;
                Real time_step_output = 0.001;
                Real jump = 0.01;
                
            
                
                
                
                if (pGrid->time >= t_next_velcocity_boost){
                    
                    /* add small amount to each velocity */
                    
                    
                    for (k=ks; k<=ke; k++) {
                        for (j=js; j<=je; j++) {
                            for (i=is; i<=ie; i++) {
                                cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
                                
#ifndef ISOTHERMAL
                                pGrid->U[k][j][i].E += (pGrid->U[k][j][i].M1)*dv + 0.5*(pGrid->U[k][j][i].d)*SQR(dv);
                                
#endif

                                
                                pGrid->U[k][j][i].M1 += (pGrid->U[k][j][i].d)*dv;
                                
                            }
                        }
                    }
                    
                    t_next_velcocity_boost += pGrid->dt;
                    
                }
                
                /* Timestep taken as a function of time */

            
                

        
            }
        }
    }
    
    return;
}


void Userwork_after_loop(MeshS *pM)
{
}




