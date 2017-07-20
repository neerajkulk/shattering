
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <fftw3.h>
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
    Real kvec,a,iseed,r;
    r = 1.0;
    

    
    /* limits for loops (integer coordinates) */
    is = pGrid->is;  ie = pGrid->ie;
    js = pGrid->js;  je = pGrid->je;
    ks = pGrid->ks;  ke = pGrid->ke;
    
    int n = je-js;
    Real fdata[n][n];
    Real data[n][n];
    
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
                cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
                Real kvec = sqrt(SQR(i-(is+ie)/2.0) + SQR(j-(is+ie)/2.0));
                Real random = rand()/(RAND_MAX/r);
                fdata[j][i] = pow(kvec,-1.0)*random; //need to add gaussian random noise. ;
                int n1=4,n2=4;
                int flag=1;

                
                fftw_complex in[4], out[4];
                fftw_plan p;
                
                
                p = fftw_create_plan(4, FFTW_FORWARD, FFTW_ESTIMATE);
                
                fftw_one(p, in, out);
                
                fftw_destroy_plan(p);

        

                
                
                
                
                
                
                P = 1.0;

                
                
                
            pGrid->U[0][j][i].d = 0.0;
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
}

void Userwork_after_loop(MeshS *pM)
{
}





