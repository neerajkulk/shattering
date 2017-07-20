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

/* Custom boundary conditions for the x1 direction*/

static void bc_ix1(GridS *pGrid);


/* additional history dump */

static Real vorticity(const GridS *pG, const int i, const int j, const int k);
static Real pleft(const GridS *pG, const int i, const int j, const int k);
static Real vyintegral(const GridS *pG, const int i, const int j, const int k);
static Real t_next_line_integral_output=0.0, dt_line_integral_output;


//may need to screw around with the static stuff...


void problem(DomainS *pDomain)
{
    GridS *pGrid=(pDomain->Grid);
    int i, is = pGrid->is, ie = pGrid->ie;
    int j, js = pGrid->js, je = pGrid->je;
    int k, ks = pGrid->ks, ke = pGrid->ke;
    Real x1,x2,x3;
    Real rho, pmax,pmin, avg, diff, prat, density, pressure, pi, vel, n, amp, lx, ly;
    
    dt_line_integral_output = par_getd("problem", "dt_line_integral");

    
    
    /* size of the domain (in physical coordinates) */
    lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
    ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
    
    
    
    pmax  = 1.0;
    pmin = 0.01;
    /* if prat=0.8, vx = -1.8965 (t&e find other vals)*/
    pi=3.14159;
    n = 2.0;                  /*Oscillations of perturbation*/
    amp = 0.05 ;/* Size of perturbation ~ 0.05 */
    avg = (pmin+pmax)/2.0;
    diff = pmax - pmin;
    
    
    /* setup uniform ambient medium with spherical over-pressured region */
    
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
                
                

                pressure = avg - diff*tanh(10000000.0*(x1-amp*sin(n*pi*x2/ly)))/2.0;
                
                if (x1 > amp*sin(n*pi*x2/ly)) {
                    vel= -1.96071;
                } else {
                    vel= -0.7;
                    
                }
                
                if (x1 > amp*sin(n*pi*x2/ly)) {
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
    
    
    
    
    /* Adding history dumps*/
    
    void dump_history_enroll(const ConsFun_t pfun, const char *label);
    
    dump_history_enroll(pleft, "<pbvals>");
    dump_history_enroll(vyintegral, "<vyintegral>");
    
    
    
    
    
    /* enroll special functions */
    bvals_mhd_fun(pDomain,left_x1, bc_ix1);
    
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
    if(strcmp(expr,"vorticity")==0) return vorticity;
    
    return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
    return NULL;
}



void Userwork_in_loop(MeshS *pM)
{
    int nl, nd, ntot;
    GridS *pGrid;
    int is, ie, js, je, ks, ke,jindex;
    
    Real line_integral = 0.0;
    
    FILE *outfile;
    
    
    for (nl=0; nl<=(pM->NLevels)-1; nl++) {
        for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                pGrid = pM->Domain[nl][nd].Grid;
                
                is = pGrid->is; ie = pGrid->ie;
                js = pGrid->js; je = pGrid->je;
                ks = pGrid->ks; ke = pGrid->ke;
                
                Real lx = pGrid->MaxX[0] - pGrid->MinX[0];
                Real ly = pGrid->MaxX[1] - pGrid->MinX[1];
                Real xres = lx/(pGrid->dx1);
                Real yres = ly/(pGrid->dx2);
                
                int iloc = (int)(xres/4);

                
                if (pGrid->time >= t_next_line_integral_output){
                    
                    /* calculate line integral */
                
                    for (jindex=js; jindex<=je; jindex++) {
                        line_integral += (pGrid->U[0][jindex][iloc].M2/pGrid->U[0][jindex][iloc].d)*(pGrid->dx2);
                        
                    }
                    
                    
                    t_next_line_integral_output += dt_line_integral_output;
                    
                    /* write data to disk */
                    outfile = fopen("integral.dat", "a");
                    fprintf(outfile, "%.20f\t%f\n",
                            pGrid->time, line_integral);
                    fclose(outfile);
                }
                
            }
        }
    }
    
    return;
}


void Userwork_after_loop(MeshS *pM)
{
}

/*---------------------------------------------------------------------------------- Constant pressure boundary condition for the inner x1 boundary*/



static void bc_ix1(GridS *pGrid)
{
    int is = pGrid->is;
    int js = pGrid->js, je = pGrid->je;
    int ks = pGrid->ks, ke = pGrid->ke;
    int i,j,k;
#ifdef MHD
    int ju, ku; /* j-upper, k-upper */
#endif
    Real P=1.0;                /* set a constant pressure */
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=1; i<=nghost; i++) {
                pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
                
                pGrid->U[k][j][is-i].E = P/Gamma_1 + (SQR(pGrid->U[0][j][i].M1) + SQR(pGrid->U[0][j][i].M2))/(2*(pGrid->U[0][j][i].d));
            }
        }
    }
    
    
    return;
}



/*----------------------------------------------------------------------------------*/
/* how well is pressure boundary condition satisfied */

static Real pleft(const GridS *pG, const int i, const int j, const int k) {
    
    return ((pG->U[0][j][nghost].E) - (SQR(pG->U[0][j][nghost].M1) + SQR(pG->U[0][j][nghost].M2))/(2*(pG->U[0][j][nghost].d))) ; /* start at 4?*/
}

/*------------------------------------------------------------------------------
 Outputting vorticity vtk*/

static Real vorticity(const GridS *pG, const int i, const int j, const int k)
{
    return (((pG->U[0][j][i+1].M2/pG->U[0][j][i+1].d) - (pG->U[0][j][i-1].M2/pG->U[0][j][i-1].d))/(2.0*pG->dx1) -
            
            ((pG->U[0][j+1][i].M1/pG->U[0][j+1][i].d) - (pG->U[0][j-1][i].M1/pG->U[0][j-1][i].d))/(2.0*pG->dx2));
}

/*------------------------------------------------------------------------------
 Line integral of vy along x1 at quarter of x domain*/



static Real vyintegral(const GridS *pG, const int i, const int j, const int k) {
    
    int jindex;
    int is = pG->is, ie = pG->ie;
    int js = pG->js, je = pG->je;
    Real lx = pG->MaxX[0] - pG->MinX[0];
    Real ly = pG->MaxX[1] - pG->MinX[1];
    Real xres = lx/(pG->dx1);
    Real yres = ly/(pG->dx2);
    
    Real sum = 0;
    
    /* x1 coordinate of line integral at res/4*/
    
    int iloc = (int)(xres/4);
    
    
    for (jindex=js; jindex<=je; jindex++) {
        sum = sum + (pG->U[0][jindex][iloc].M2/pG->U[0][jindex][iloc].d)*(pG->dx2);
        
    }
    
    
    return sum;
}




