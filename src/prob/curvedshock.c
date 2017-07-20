
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
    Real rho,P,vx,vy,vb;
    Real amp,n,r,rhomax,rhomin,rhoavg,rhodiff,thickness,a,iseed;
    Real Pavg,Pdiff,Pmax,Pmin,vxavg,vxdiff,vxmax,vxmin;
    
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
            

            amp = 0.2;
            n = 2.0;
            r = 0.05;
            thickness = 0.000005;
                
            rhomax = 1.0;
            rhomin = 0.357013;
            
            rhoavg = (rhomax + rhomin)/2.0;
            rhodiff = rhomax - rhomin;
                
                
            Pmax = 1.0;
            Pmin = 0.1175;
                
            Pavg = (Pmax + Pmin)/2.0;
            Pdiff = Pmax - Pmin;
                
            vxmax = -1.96071;
            vxmin = -0.7;
                
            vxavg = (vxmax + vxmin)/2.0;
            vxdiff = vxmax - vxmin;
                
                
            vy = 0.0;
            
            /* analytic function to represnted a smoothed out interfce + a random perturbation. define thickness of interface when 0.9 of max value is reached*/
            a = 2.0*atanh(0.9)/thickness;
            
            /*rho = rhoavg - rhodiff*tanh(a*(x1-amp*sin(n*PI*x2/ly)))/2.0;
            P = Pavg - Pdiff*tanh(a*(x1-amp*sin(n*PI*x2/ly)))/2.0;
            vx = vxavg - vxdiff*tanh(a*(-x1+amp*sin(n*PI*x2/ly)))/2.0;*/
                
                
                if (x1<0.1*sin(2*PI*x2/2.0) ) {
                    rho=1.0;
                } else {
                    rho=0.574822;
                }
               
                if (x1<0.1*sin(2*PI*x2/2.0) ) {
                    P=1.0;
                } else {
                    P=0.105;
                }
                
                if (x1<0.1*sin(2*PI*x2/2.0) ) {
                    vx=-1.1;
                } else {
                    vx=-1.91364;
                }
                
                
                
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
    return;
}

void Userwork_after_loop(MeshS *pM)
{
}

