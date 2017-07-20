
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* Additional output file */
static Real temperature(const GridS *pG, const int i, const int j, const int k);


/* problem() function defines the initial condition */
void problem(DomainS *pDomain)
{
    /* define variables */
    GridS *pGrid = pDomain->Grid;
    int i,j,k,is,ie,js,je,ks,ke;
    Real x1,x2,x3,lx,ly;
    Real rho,P,vx,vy,vb;
    Real amp,n,r,rhomax,rhomin,avg,diff,thickness,a,iseed;
    
    vb = par_getd("problem", "boost"); /* Boost velocity from input file*/
    
    /* limits for loops (integer coordinates) */
    is = pGrid->is;  ie = pGrid->ie;
    js = pGrid->js;  je = pGrid->je;
    ks = pGrid->ks;  ke = pGrid->ke;
    
    /* size of the domain (in physical coordinates) */
    lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
    ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  
    /*seeding random numbers across different proccesors*/
    iseed = -10;
    #ifdef MPI_PARALLEL
    iseed -= myID_Comm_world;
    #endif
    srand(iseed);
    
    for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
            for (k=ks; k<=ke; k++) {
                
            cc_pos(pGrid,i,j,0,&x1,&x2,&x3);
            
            vy = 0.0;
            vx = vb;
            rhomax = 1.0;
            rhomin = 0.1;
            amp = 0.1;
            n = 2.0;
            r = 0.05;
            thickness = 0.05;
            
            
            /* defines random to be a random number between -r/2 and r/2 */
            Real random = r/2.0 - rand()/(RAND_MAX/r);


          /*Setting up initial background density imbalance with a perturbed interface*/
            
            P = 1.0;
            
            avg = (rhomax + rhomin)/2.0;
            diff = rhomax - rhomin;
            
            /* analytic function to represnted a smoothed out interfce + a random perturbation. define thickness of interface when 0.9 of max value is reached*/
            a = 2.0*atanh(0.9)/thickness;
            
            rho = avg - diff*tanh(a*(-x1-amp*sin(n*PI*x2/ly)))/2.0 + random;
            
            /*Write values into grid */
            
            pGrid->U[k][j][i].d = rho;
            pGrid->U[k][j][i].E = P/Gamma_1 + 0.5*rho*(SQR(vx) + SQR(vy));
            
            pGrid->U[k][j][i].M1 = vx * rho;
            pGrid->U[k][j][i].M2 = vy * rho;
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
    int nl,nd;
    Real vb;

    
    /* set boundary conditions.  since this function takes a Mesh, you
     need to "find" the Domain */
    for (nl=0; nl<(pM->NLevels); nl++) {
        for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++) {
            if (pM->Domain[nl][nd].Disp[0] == 0)
                
                vb = par_getd("problem", "boost");
            


                
        }
    }
    
    
    /* par_getd statements copied from problem() */
    
    /* dump_history_enroll statements copied from problem() */
    
    /* gravity, etc. statements copied from problem() */
    
    /* DANGER: make sure the order here matches the order in write_restart() */

    
    return;
}

ConsFun_t get_usr_expr(const char *expr)
{
    if(strcmp(expr,"temperature")==0) return temperature;
    
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
   
    //printf("Floor temp %lf \n",random);
    
    return;
}

void Userwork_after_loop(MeshS *pM)
{
}
/*------------------------------------------------------------------------------
 Outputting temperature vtk*/

static Real temperature(const GridS *pG, const int i, const int j, const int k)
{
    Real TE,KE,Temp;
    
    KE = (SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)+ SQR(pG->U[k][j][i].M3)) / (2.0*pG->U[k][j][i].d);
    
    TE = pG->U[k][j][i].E - KE;
    
    Temp = ((2.0/3.0)*TE)/(pG->U[k][j][i].d);

    return Temp;
}

