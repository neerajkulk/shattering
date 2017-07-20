
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"


/* problem() function defines the initial condition */
void problem(DomainS *pDomain)
{
    /* define variables */
    GridS *pGrid = pDomain->Grid;
    int i,j,k,is,ie,js,je,ks,ke;
    Real x1,x2,x3,lx,ly;
    Real rho,P;
    
    /* limits for loops (integer coordinates) */
    is = pGrid->is;  ie = pGrid->ie;
    js = pGrid->js;  je = pGrid->je;
    ks = pGrid->ks;  ke = pGrid->ke;
    
    /* size of the domain (in physical coordinates) */
    lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
    ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  
       for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
            for (k=ks; k<=ke; k++) {
                
            cc_pos(pGrid,i,j,0,&x1,&x2,&x3);
            
                if (x1<0) {
                    P = 1.0;
                } else {
                    P = 0.995;
                }
                
                if (x1<0) {
                    rho = 0.1;
                } else {
                    rho = 1.0;
                }
                
                
                /*Write values into grid */
            
            pGrid->U[k][j][i].d = rho;
            pGrid->U[k][j][i].E = P/Gamma_1;
            pGrid->U[k][j][i].M1 = 0.0;
            pGrid->U[k][j][i].M2 = 0.0;
            pGrid->U[k][j][i].M3 = 0.0;
         
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
    int i,j,k;
    Real KE,TE,Temp,dE;
    Real FloorTemp;
    for (nl=0; nl<=(pM->NLevels)-1; nl++) {
        for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                pGrid = pM->Domain[nl][nd].Grid;
                
                is = pGrid->is; ie = pGrid->ie;
                js = pGrid->js; je = pGrid->je;
                ks = pGrid->ks; ke = pGrid->ke;
                
                FloorTemp = par_getd("problem", "Floor");

                
                for(k=ks; k<=ke; k++){
                    for (j=js; j<=je; j++){
                        for (i=is; i<=ie; i++){
                            
                            
                            if (pGrid->U[k][j][i].d != pGrid->U[k][j][i].d || pGrid->U[k][j][i].d <= 0.0)
                                ath_error("bad density of %f at cell (%d,%d,%d)\n",pGrid->U[k][j][i].d, k,j,i);
                            
                            
                            KE = (SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)+ SQR(pGrid->U[k][j][i].M3)) / (2.0*pGrid->U[k][j][i].d);
                            
                            if (KE != KE || KE >= 1000.0){
                                ath_error("bad KE of %f at (%d,%d,%d)\n", KE,
                                          k,j,i);
                            }
                            
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
                        }
                    }
                }                
            }
        }
    }
   
    
    return;
}

void Userwork_after_loop(MeshS *pM)
{
}


