/* declarations to go at the top of the file */

/* global declarations for radial profile outputs */
/*  */
static OutputS profile_dump;
void dump_profile(MeshS *pM, OutputS *pOut);

static int nmax_profile;
static Real *num = NULL,
  *radius        = NULL,
  *d_profile     = NULL,
  *press_profile = NULL,
  *temp_profile  = NULL,
  *mach_profile  = NULL; /* etc. see TODO below */

#ifdef MPI_PARALLEL
static Real *num_global = NULL,
  *radius_global        = NULL,
  *d_profile_global     = NULL,
  *press_profile_global = NULL,
  *temp_profile_global  = NULL,
  *mach_profile_global  = NULL; /* etc. see TODO below */
#endif  /* MPI_PARALLEL */

static int calc_profiles(MeshS *pM,
                         int *nmax,
                         Real *radius,
                         Real *d_profile,
                         Real *press_profile,
                         Real *temp_profile,
                         Real *mach_profile) /* TODO: make sure this
                                                matches the arrays
                                                declared above! */
/*  */
/* end  global declarations for radial profile outputs */



/* in problem() */
profile_dump.n      = 100;
profile_dump.dt     = par_getd("profile", "dt");
profile_dump.t      = pM->time;
profile_dump.num    = 0;
profile_dump.out    = "prim";
profile_dump.nlevel = -1;       /* dump all levels */




/* in write_restart() -- CAREFUL WITH ORDER!! */
fwrite(&profile_dump.num, sizeof(int),   1, fp);
fwrite(&profile_dump.t,   sizeof(Real),  1, fp);
fwrite(&profile_dump.dt,  sizeof(Real),  1, fp);

/* in read_restart()  -- CAREFUL WITH ORDER!! */
profile_dump.n      = 100;
profile_dump.out    = "prim";
profile_dump.nlevel = -1;       /* dump all levels */

fread(&profile_dump.num, sizeof(int),   1, fp);
fread(&profile_dump.t,   sizeof(Real),  1, fp);
fread(&profile_dump.dt,  sizeof(Real),  1, fp);




/* in Userwork_in_loop() */
/*  */

/* check whether we need to do an output */
if(pM->time >= profile_dump.t){

  /* first, update output time */
  profile_dump.t += profile_dump.dt;

  /* next, calculate the radial profiles and store them in the global array */
  ierr = calc_profiles(pM,
                       radius,
                       d_profile,
                       press_profile,
                       temp_profile,
                       mach_profile); /* TODO: match list of arrays */

  if (ierr)
    ath_error("[Userwork_in_loop]: calc_profiles() returned error code %d\n", ierr);


  /* finally, write the data to disk, but only on the root process */
#ifdef MPI_PARALLEL
  if (myID_Comm_world == 0){
#endif /* MPI_PARALLEL */
    dump_profile(pM, &profile_dump);
#ifdef MPI_PARALLEL
  }
#endif /* MPI_PARALLEL */

  profile_dump.num += 1;
}



/* in Userwork_after_loop() */
/*  */
/* free memory if necessary */
if (num           != NULL) free_3d_array((void*) num);
if (radius        != NULL) free_3d_array((void*) radius);
if (d_profile     != NULL) free_3d_array((void*) d_profile);
if (press_profile != NULL) free_3d_array((void*) press_profile);
if (temp_profile  != NULL) free_3d_array((void*) temp_profile);
if (mach_profile  != NULL) free_3d_array((void*) mach_profile);

#ifdef MPI_PARALLEL
if (num_global           != NULL) free_3d_array((void*) num_global);
if (radius_global        != NULL) free_3d_array((void*) radius_global);
if (d_profile_global     != NULL) free_3d_array((void*) d_profile_global);
if (press_profile_global != NULL) free_3d_array((void*) press_profile_global);
if (temp_profile_global  != NULL) free_3d_array((void*) temp_profile_global);
if (mach_profile_global  != NULL) free_3d_array((void*) mach_profile_global);
#endif /* MPI_PARALLEL */
 }




void calc_profiles(MeshS *pM,
                   int   *nmax,
                   Real  *radius,
                   Real  *d_profile,
                   Real  *press_profile,
                   Real  *temp_profile,
                   Real  *mach_profile) /* TODO: make sure this
                                           matches the arrays declared
                                           above! */
{
   int nl, nd, ntot;
   GridS *pGrid;
   DomainS *pD;
   int is, ie, js, je, ks, kj;
   int i, j, k, s;
   double r;

   Real Nx = pM->Nx[0];
   Real Ny = pM->Nx[1];
   Real Nz = pM->Nx[2];

   /*Size of arrays*/
   Real c = MAX(Nx, Ny);
   c = MAX(c, Nz);

#ifdef MPI_PARALLEL
   int ierr;
#endif /*MPI_PARALLEL*/

   PrimS W;
   ConsS U;


   /* only call calloc() once */
   if (num == NULL){
     num           = calloc_1d_array(c, sizeof(Real));
     radius        = calloc_1d_array(c, sizeof(Real));
     d_profile     = calloc_1d_array(c, sizeof(Real));
     press_profile = calloc_1d_array(c, sizeof(Real));
     temp_profile  = calloc_1d_array(c, sizeof(Real));
     mach_profile  = calloc_1d_array(c, sizeof(Real));

#ifdef MPI_PARALLEL
     num_global           = calloc_1d_array(c, sizeof(Real));
     radius_global        = calloc_1d_array(c, sizeof(Real));
     d_profile_global     = calloc_1d_array(c, sizeof(Real));
     press_profile_global = calloc_1d_array(c, sizeof(Real));
     temp_profile_global  = calloc_1d_array(c, sizeof(Real));
     mach_profile_global  = calloc_1d_array(c, sizeof(Real));
#endif /* MPI_PARALLEL */
   }


   /* zero out arrays */
   for(i=0; i<c; i++){
     num[i] = d_profile[i] = press_profile[i] = temp_profile[i] \
       = mach_profile[i] = 0.0;
   }


   /* loop over mesh to find the grid */
   for(nl=0; nl<=(pM->NLevels)-1; nl++){
     for(nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++){
       if(pM->Domain[nl][nd].Grid != NULL){
         pGrid = pM->Domain[nl][nd].Grid;

         is = pGrid->is; ie=pGrid->ie;
         js = pGrid->js; je=pGrid->je;
         ks = pGrid->ks; ke=pGrid->ks;

         for(k=ks; k<=ke; k++){
           for(j=js; j<=je; j++){
             for(i=is; i<=ie; i++){
               /* calculate radius in cell coordinates and
                  corresponding index */
               r = sqrt(SQR(i-0.5*Nx) + SQR(j-0.5*Ny) + SQR(k-0.5*Nz));
               s = (int) floor(r);

               W = Cons_to_Prim(&(pG->U[k][j][i]));

               /* update arrays */
                num[s]           += 1.0;
                radius[s]        += r;
                d_profile[s]     += W.d;
                press_profile[s] += W.P;
                temp_profile[s]  += W.P/W.d;
                mach_profile[s]  += W.d * (SQR(W.V1)+SQR(W.v2)+SQR(W.V3)) / W.P;
             }
           }
         }
#ifdef MPI_PARALLEL
         /* TODO: consider packing into one 2D array and using a
            single MPI call */
         ierr = MPI_ALLreduce(&num, &num_global, c,
                              MPI_RL, MPI_SUM, MPI_COMM_WORLD);
         if(ierr)
           ath_error("[calc_profiles]: MPI_Allreduce returned error %d\n", ierr);

         ierr = MPI_ALLreduce(&radius, &radius_global, c,
                              MPI_RL, MPI_SUM, MPI_COMM_WORLD);
         if(ierr)
           ath_error("[calc_profiles]: MPI_Allreduce returned error %d\n", ierr);

         ierr = MPI_ALLreduce(&d_profile, &d_profile_global, c,
                              MPI_RL, MPI_SUM, MPI_COMM_WORLD);
         if(ierr)
           ath_error("[calc_profiles]: MPI_Allreduce returned error %d\n", ierr);

         ierr = MPI_ALLreduce(&press_profile, &press_profile_global, c,
                              MPI_RL, MPI_SUM, MPI_COMM_WORLD);
         if(ierr)
           ath_error("[calc_profiles]: MPI_Allreduce returned error %d\n", ierr);

         ierr = MPI_ALLreduce(&temp_profile, &temp_profile_global, c,
                              MPI_RL, MPI_SUM, MPI_COMM_WORLD);
         if(ierr)
           ath_error("[calc_profiles]: MPI_Allreduce returned error %d\n", ierr);

         ierr = MPI_ALLreduce(&mach_profile, &mach_profile_global, c,
                              MPI_RL, MPI_SUM, MPI_COMM_WORLD);
         if(ierr)
           ath_error("[calc_profiles]: MPI_Allreduce returned error %d\n", ierr);

         for(i=0; i<c; i++){
           num[i]           = num_global[i];
           radius[i]        = radius_global[i];
           d_profile[i]     = d_profile_global[i];
           press_profile[i] = press_profile_global[i];
           temp_profile[i]  = temp_profile_global[i];
           mach_profile[i]  = mach_profile_global[i];
         }
#endif /*MPI_PARALLEL*/

         for(i=0; i<c; i++){
           /* Divide by number of cells to get proper average */
           tvals[i] /= num[i];
           rvals[i] /= num[i];
           kvals[i] /= num[i];
         }


       } /* if(pM->Domain[nl][nd].Grid != NULL) */
     }   /* loop over domains */
   }     /* loop over SMR levels */

   return;
}





void dump_profile(MeshS *pM, OutputS *pOut)
{
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
  Real x1,x2,x3;
  char zone_fmt[20], fmt[80];
  int col_cnt;

  /* Anoop, add whatever variables you need here */
  int c;

  /* Anoop, mostly ignore this stuff until you see things flagged TODO */

  sprintf(fmt," %%12.8e"); /* Use a default format */

  col_cnt = 1;

  /* construct output filename. */
  if((fname = ath_fname(NULL,
                        pM->outfilename,
                        NULL,
                        NULL,
                        num_digit,
                        pOut->num,
                        NULL,
                        "pro")) == NULL){
    ath_error("[dump_profile]: Error constructing filename\n");
  }

  /* open output file */
  if((pfile = fopen(fname,"w")) == NULL){
    ath_error("[dump_profile]: Unable to open text file %s\n",fname);
  }
  free(fname);

  /* Upper and Lower bounds on i,j,k for data dump */
  sprintf(zone_fmt,"%%%dd", (int)(2+log10((double)(nmax_profile))));

  /* Write out some header information */
  if (pG->Nx[0] > 1) {
    fprintf(pfile,"# N = %d\n", nmax_profile);
  }
  fprintf(pfile,"# RADIAL PROFILE at Time= %g\n", pM->time);

  /* write out column headers.  Note column number is embedded in header */
  fprintf(pfile,"# [%d]=i",col_cnt);
  col_cnt++;

  if (pG->Nx[0] > 1) {
    fprintf(pfile," [%d]=r",col_cnt);
    col_cnt++;
  }

  /* TODO: change to the vars you want */
  /* MM recommends: density, pressure, temperature, entropy, mach
     number, maybe radial velocity.  try to keep this code flexible so
     you can change it */
  /* later, let's add magnetic energy and plasma beta */

  fprintf(pfile," [%d]=rho",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=P",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=T",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=M",col_cnt);
  col_cnt++;

  fprintf(pfile,"\n");

  /* Write out data */

  for(c=0; c<nmax_profile; c++){
    fprintf(pfile, zone_fmt, c);
    fprintf(pfile, fmt, radius[c]);

    /* Dump all variables */
    fprintf(pfile, fmt, d_profile[c]);
    fprintf(pfile, fmt, press_profile[c]);
    fprintf(pfile, fmt, temp_profile[c]);
    fprintf(pfile, fmt, mach_profile[c]); /* TODO: match definitions */


    fprintf(pfile,"\n");
  }

  fclose(pfile);

  return;
}
