#include <stdio.h>
#include <math.h>

#define prs = 950
#define dp = 100
#define gra = 9.81
#define h2o_mol_wgt = 1.8e-2
#define h2o_mx_rat = 0.01
#define avo = 6.022140e23

int main()
{
   int i, dum, lines, ID;
   FILE *fp = fopen("text.txt","r");

   long double nu, inten, A, gam_air, gam_self, E_pprime, n_air, del_air;
//   long double nu_ini = 100.0, nu_end = 2000.0, nu_step = 1.0e-3, *delnu;

   char string[161];


   while(fgets(string, 161, fp))
   {
      sscanf(string, "%2d %12Lf %10Le %10Le %5Lf %5Lf %10Lf %4Lf %8Lf", &ID, &nu, &inten, &A, &gam_air, &gam_self, &E_pprime, &n_air, &del_air);

      printf("%12.6Lf, %10.3Le, %5.4Lf, %4.2Lf, %8.6Lf\n", nu, inten, gam_air, n_air, del_air);
   }

   fclose(fp);
   return 0;
}
