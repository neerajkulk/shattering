#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "ath_array.h"
#include "ath_error.h"
#include "ath_vtk.h"
#include "par.h"

/* variables defined in vtk.c */
extern int Nx, Ny, Nz;          /* size of the grid in cell coordinates */
extern double dx, dy, dz;       /* size of each cell in physical coordinates */

extern Real3Vect ***vel;
extern float     ***rho;
extern float     ***press;

static float *rvals, *tvals, *kvals;
static int   *num;

void calc_profiles();
void write_data(char *outfname);


int main (int argc, char *argv[])
{
  FILE *fp;

  int i, c;

  char *vtkfile, *outfname, buf[512];

  char *definput = "input.pro";
  char *athinput = definput;


  /* parse command line options */
  for (i=1; i<argc; i++) {
    if (*(argv[i]) == '-') {
      switch(*(argv[i]+1)) {
      case 'i':                      /* -i <file>   */
        athinput = argv[++i];
        break;
      default:
        break;
      }
    }
  }

  par_open(athinput);
  par_cmdline(argc, argv);


  /* read the input file */
  vtkfile  = par_gets("files", "vtk_file");
  sprintf(buf, "%s.flines", vtkfile);
  outfname = par_gets_def("files", "out_file", buf);

  par_dump(2, stdout);
  par_close();


  /* read the VTK file */
  fp = fopen(vtkfile, "r");
  vtkread(fp);
  fclose(fp);


  c = MAX(Ny, Nz);
  c = MAX(Nx, c);

  rvals = calloc_1d_array(c, sizeof(float));
  tvals = calloc_1d_array(c, sizeof(float));
  kvals = calloc_1d_array(c, sizeof(float));

  num = calloc_1d_array(c, sizeof(int));

  calc_profiles();


  /* save the data to disk */
  write_data(outfname);

  /* Free the arrays used by read_vtk */
  cleanup_vtk();
  free_1d_array((void*)  rvals);
  free_1d_array((void*)  tvals);
  free_1d_array((void*)  kvals);
  free_1d_array((void*)  num);

  return 0;
}

void calc_profiles()
{
  int i, j, k, c, max;
  double r;

  max = MAX(Ny, Nz);
  max = MAX(Nx, max);

  for (i=0; i<max; i++) {
    num[i] = 0;
    tvals[i] = rvals[i] = kvals[i] = 0.0;
  }

  for (k=0; k<Nz; k++){
    for (j=0; j<Ny; j++){
      for (i=0; i<Nx; i++){
        r = sqrt(SQR(i - 0.5*Nx)
                 + SQR(j - 0.5*Ny)
                 + SQR(k - 0.5*Nz));
        c = (int) floor(r);

        num[c] += 1;
        tvals[c] += press[k][j][i]/rho[k][j][i];
        rvals[c] += rho[k][j][i];
        kvals[c] += (SQR(vel[k][j][i].x1)
                     + SQR(vel[k][j][i].x2)
                     + SQR(vel[k][j][i].x3));
      }
    }
  }

  for (i=0; i<max; i++) {
    if (num[i] > 0){
      tvals[i] /= num[i];
      rvals[i] /= num[i];
      kvals[i] /= num[i];
    }
  }

  return;
}


void write_data(char *outfname)
{
  int i, c;

  FILE *outfile;
  outfile = fopen(outfname, "w");

  c = MAX(Ny, Nz);
  c = MAX(Nx, c);

  fprintf(outfile, "[1] = r, [2] = T, [3] = v^2");

  for (i=0; i<c; i++) {
    fprintf(outfile, "%d\t%e\t%e\n", i, tvals[i], kvals[i]);
  }

  fclose(outfile);

  return;
}
