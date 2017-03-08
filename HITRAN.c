#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define prs 950.0
#define dp 100.0
#define gra 9.81
#define h2o_mol_wgt 1.8e-2
#define h2o_mx_rat 0.01
#define avo 6.022140e23
#define pi 4*atan(1.0)

long double *arange(long double start, long double end, long double step)
{
   // 'arange' routine.
   int i; 
   int arr_size = ((end - start) / step) + 1;
   long double *output = malloc(arr_size * sizeof(long double));

   for(i=0;i<arr_size;i++)
   {
      output[i] = start + (step * i);
   }

   return output;
}

int main()
{
   long double nu_ini = 100.0, nu_end = 2000.0, nu_step = 1.0e-3;
   long double *delnu = arange(nu_ini, nu_end, nu_step);
   long double *nu, *inten, A, *gam_air, gam_self, E_pprime, *n_air, *del_air;
   long double *gamma, *f, *f_sum;
   
   int i, j, dum, lines=0, ID, delnu_size = (((nu_end - nu_ini)/nu_step) + 1);
   FILE *fp = fopen("h2o_HITRAN.par","r");
   char string[320];

   while(!feof(fp))
   {
     dum = fgetc(fp);
     if(dum == '\n')
     {
       lines++;
     }
   }
  
   rewind(fp);

   nu       = malloc(lines * sizeof(long double));
   inten    = malloc(lines * sizeof(long double));
   gam_air  = malloc(lines * sizeof(long double));
   n_air    = malloc(lines * sizeof(long double));
   del_air  = malloc(lines * sizeof(long double));
   gamma    = malloc(lines * sizeof(long double));
   f        = malloc(delnu_size * sizeof(long double));
//   f_sum    = malloc(delnu_size * sizeof(long double));
  

   i=0;
   while(fgets(string, 320, fp))
   {
      sscanf(string, "%2d %12Lf %10Le %10Le %5Lf %5Lf %10Lf %4Lf %8Lf", &ID, &nu[i], &inten[i], &A, &gam_air[i], &gam_self, &E_pprime, &n_air[i], &del_air[i]);

      i++;
   }

   
   // gamma calculation
   for(i=0;i<lines;i++)
   {
      gamma[i] = pow((296.0/300.0),n_air[i]) * (gam_air[i]*(prs/1013.0));
   }

   // Initialize data container
   for(i=0;i<delnu_size;i++)
   {
//      f_sum[i] = 0.0;
      f[i] = 0.0;
   }
  
   // 'f': absorption coefficient calculation
   for(i=0;i<lines;i++)
      for(j=0;j<delnu_size;j++)
      {
         f[j] += inten[i] * ((1.0/pi) * (gamma[i] / (pow(gamma[i],2.0) + pow(delnu[j] - (nu[i] + del_air[i]*prs/1013.0),2.0))));
//         f_sum[j] = f_sum[j] + f[j];
      }

   FILE *file = fopen("plottest","w");
   
   for(i=0;i<delnu_size;i++)
   {
      fprintf(file, "%Le\n", f[i]);
   }

//   printf("\n%Le, %Le\n", f[0], f[1899999]);
   free(nu);
   free(inten);
   free(gam_air);
   free(n_air);
   free(del_air);
   free(delnu);
   free(gamma);
   free(f);
   fclose(fp);
   fclose(file);
   return 0;
}
