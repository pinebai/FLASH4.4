/* sedov.c                                                                 */
/*                                                                         */
/* Author: Aamer Haque                                                     */
/* Email: ahaque@flash.uchicago.edu                                        */
/*                                                                         */
/* Program that computes the analytical solution to the Sedov problem.     */
/* Reference: L.I. Sedov, Similarity and Dimensional Methods in Mechanics, */
/*            Academic Press, New York, 1959.                              */
/* License: Non-commerical use ok. Premission to redistribute provided     */
/*          that the author is credited.                                   */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
	
#define INT int	
#define REAL double	

#define ABS(a)      fabs(a)
#define SQR(a)     ((a)*(a))
#define CUBE(a)    ((a)*(a)*(a))
#define SQRT(a)    (sqrt(a))
#define POW(a,b)   (pow(a,b))
#define LOG(a)     (log10(a))
#define LN(a)      (log(a))
#define EXP(a)     (exp(a))

#define PI   3.14159265358979

#define MAX_POINTS      10000
#define MAX_VV_POINTS    1000
#define MAX_V_ITS         200
#define MAX_ALPHA_ITS      25

#define E_TOL       1.0e-7
#define V_MIN_TOL   1.0e-3

#define TRUE    (1==1)
#define FALSE   (!TRUE)

typedef int BOOLEAN;


/* Variables read from parameter file */
INT geo;                     /* Integer-valued geometry     */
REAL E0;                     /* Energy source               */
REAL g;                      /* Gamma                       */
REAL rho0;                   /* Initial density             */
REAL t;                      /* Time                        */
REAL r_max;                  /* Maximum size of domain      */
	
REAL nu;                     /* Real-valued geometry        */
INT n;	                     /* Number of spatial points    */
BOOLEAN rmode;               /* Use dimensional coordinates */
	
/* Look-up table to determine V and r */
REAL VV[MAX_VV_POINTS];      /* Similarity variable V       */
REAL rr[MAX_VV_POINTS];      /* Distance from origin        */

/* Solution for each variable */
REAL r[MAX_POINTS];          /* Distance from origin        */
REAL Vs[MAX_POINTS];         /* Similarity variable V       */
REAL rho[MAX_POINTS];        /* Density                     */
REAL p[MAX_POINTS];          /* Pressure                    */
REAL u[MAX_POINTS];          /* Velocity                    */

/* Energy variables */
REAL KE[MAX_POINTS];         /* Kinetic energy              */
REAL IE[MAX_POINTS];         /* Internal energy             */
REAL Etot[MAX_POINTS];       /* Total energy                */
REAL KE_cum[MAX_POINTS];     /* Cumulative kinetic energy   */
REAL IE_cum[MAX_POINTS];     /* Cumulative internal energy  */
REAL Etot_cum[MAX_POINTS];   /* Cumulative total energy     */

/* Additional parameters needed for the solution */
REAL a[8];                   /* Alpha array                 */
REAL V_min;                  /* Minimum value of V          */
REAL V_max;                  /* Maximum value of V          */
REAL k1;                     /* k1 parameter for solution   */
REAL k2;                     /* k2 parameter for solution   */

/* Values for the post-shock state */
REAL shock_speed;            /* Shock speed                 */
REAL r_shock;                /* Location of shock           */
REAL rho_shock;              /* Post-shock density          */
REAL p_shock;                /* Post-shock pressure         */
REAL u_shock;                /* Post-shock velocity         */



void read_params(char filename[])
{
   FILE *infile;
   char junk[80];

   infile = fopen(filename,"r");
   if (infile==NULL) {
      printf("\nCannot open file: %s\n",filename);
      exit(1);
   }
   
   printf("\nReading parameter file: %s\n",filename);

   fscanf(infile,"%d %s",&geo,junk);
   printf("d=%d\n",geo);

   fscanf(infile,"%lg %s",&E0,junk);
   printf("E0=%g\n",E0);
   if (E0<=0.0) {
      printf("Error: rho0 must be greater than 0\n");
      exit(1);
   }

   fscanf(infile,"%lg %s",&g,junk);
   printf("g=%g\n",g);
   if ((g<1.1)||(g>1.9)) {
      printf("Error: g must be in [1.1,1.9]\n");
      exit(1);
   }

   fscanf(infile,"%lg %s",&rho0,junk);
   printf("rho0=%g\n",rho0);
   if (rho0<=0.0) {
      printf("Error: rho0 must be greater than 0\n");
      exit(1);
   }

   fscanf(infile,"%lg %s",&t,junk);
   printf("t=%g\n",t);
   if (t<=0.0) {
      printf("Error: t must be greater than 0\n");
      exit(1);
   }

   fscanf(infile,"%lg %s",&r_max,junk);
   printf("r_max=%g\n\n",r_max);
   if (r_max<=0.0) {
      printf("Error: r_max must be greater than 0\n");
      exit(1);
   }

   fclose(infile);
}



void gen_r_uniform() 
{
   INT i;
   REAL dx;
   
   n = MAX_POINTS;
   dx = r_max / ((REAL) n);
   
   for (i=0; i<n; i++) 
      r[i] = (i+0.5)*dx;
}	



void read_r(char filename[])
{
   FILE *infile;
   char junk[80];
   INT i;

   n = 0;

   infile = fopen(filename,"r");
   if (infile==NULL) {
      printf("Error: Cannot read file %s\n",filename);
      exit(1);
   }

   fscanf(infile,"%s %d",junk,&n);

   for (i=0; i<n; i++)
      fscanf(infile,"%lg",&r[i]);

   fclose(infile);
}



void calc_a()
{
   switch (geo) {
   case 1:
      nu = 1.0;
      break;
   case 2:
      nu = 2.0;
      break;
   case 3:
      nu = 3.0;
      break;
   default:
      printf("Invalid geometry\n");
      exit(1);
      break;
   }

   a[2] = (1.0-g) / (2.0*(g-1) + nu);
   a[3] = nu / (2.0*(g-1.0) + nu);
   a[5] = 2.0 / (g-2.0);
   a[6] = g / (2.0*(g-1.0) + nu);
   
   a[1] = ( ((nu+2.0)*g)/(2.0+nu*(g-1.0)) ) * ( (2.0*nu*(2.0-g))/(g*SQR(nu+2.0)) - a[2]);
   a[4] = (a[1]*(nu+2.0)) / (2.0-g);
   a[7] = (2.0 + nu*(g-1))*a[1]/(nu*(2.0-g));
}



void calc_k()
{
   switch (geo) {
   case 1:
      k1 =  ( (5.0*g-4.0) / (3.0*(g-1.0)*(2.0-g)) ) * LN(2.0*g-1.0);
      k1 += ( (4.0+g-3.0*g*g) / (3.0*(g-1.0)*(2.0-g)) ) * LN(g+1.0);
      k1 -= ( 2.0 / (3.0*(g-1.0)) ) * LN(2.0);
      k1 -= ( 1.0 / (g-1.0) ) * LN(g);
      k1 -= LN(g-1.0);
      k1 =  EXP(k1);
      
      k2 =  (7.0/3.0) * LN(2.0);
      k2 += ( (5.0*g-4.0) / (3.0*(2.0-g)) ) * LN(2.0*g-1.0);
      k2 -= LN(9.0);
      k2 -= ( 2.0*(g+1.0) / (3.0*(2.0-g)) ) * LN(g+1.0);
      k2 =  EXP(k2);
      break;
   case 2:
      k1 =  ( (g+1.0)/(g-1.0) ) * LN(g+1.0);
      k1 += ( (3.0*g-4.0)/((g-1.0)*(2.0-g)) ) * LN(g);
      k1 -= LN(g-1.0);
      k1 -= ( 2.0/((g-1.0)*(2.0-g)) ) * LN(2.0);
      k1 =  EXP(k1);
      
      k2 =  ( 2.0*(g-1.0)/(2.0-g) ) * LN(g);
      k2 -= ( (4.0-g)/(2.0-g) ) * LN(2.0);
      k2 =  EXP(k2);
      break;
   case 3:	 	
      k1 =  ( (3.0*g*(g+1.0)) / ((g-1.0)*(3.0*g-1.0)) ) * LN(g+1.0);
      k1 -= LN(g-1.0);
      k1 -= ( 6.0/(5*(g-1.0)) ) * LN(2);
      k1 -= ( (7.0*g-1.0) / ((g-1.0)*(3.0*g-1.0)) ) * LN(g);
      k1 += ( (13*g*g - 7.0*g + 12.0) / (5.0*(g-1.0)*(2.0-g)*(3*g-1.0)) ) 
	 * LN((2.0*g+1.0)/(7.0-g));
      k1 =  EXP(k1);
	 
      k2 =  LN(0.32);
      k2 += ( (g+1.0)/(3.0*g-1.0) ) * LN(g+1.0);
      k2 -= (6.0/5.0) * LN(2.0);
      k2 -= ( 4.0*g/(3.0*g-1.0) ) * LN(g);      
      k2 += ( (13*g*g - 7.0*g + 12.0) / (5.0*(2.0-g)*(3*g-1.0)) ) 
	 * LN((2.0*g+1.0)/(7.0-g));
      k2 =  EXP(k2);
      break;
   }
}



void calc_V_min_max()
{
   if ((geo>=3) && (g > 7.0)) {
      V_min = 4.0 / (5.0*(g+1.0));
      V_max = 0.4;
   }
   else {
      V_min = 2.0 / ((nu+2.0)*g);
      V_max = 4.0 / ((nu+2.0)*(g+1.0));
   }   
}



void calc_shock(REAL alpha)
{
   REAL Es;

   Es = E0 / alpha;

   switch (geo) {
   case 1:
      shock_speed = (2.0/3.0)*POW(Es/rho0,1.0/3.0)*POW(t,-1.0/3.0);
      break;
   case 2:
      shock_speed = 0.5*POW(Es/rho0,0.25)/SQRT(t);
      break;
   case 3:
      shock_speed = 0.4*POW(Es/rho0,0.2)*POW(t,-0.6);
      break;
   }

   r_shock = POW(Es/rho0,1.0/(2.0+nu)) * POW(t,2.0/(2.0+nu));
   rho_shock = ((g+1.0)/(g-1.0))*rho0;
   p_shock   = (2.0*rho0*SQR(shock_speed)) / (g+1.0);
   u_shock   = (2.0*shock_speed) / (g+1.0);
}



void calc_para_shock(REAL Es)
{
   switch (geo) {
   case 1:
      shock_speed = (2.0/3.0)*POW(Es/rho0,1.0/3.0)*POW(t,-1.0/3.0);
      break;
   case 2:
      shock_speed = 0.5*POW(Es/rho0,0.25)/SQRT(t);
      break;
   case 3:
      shock_speed = 0.4*POW(Es/rho0,0.2)*POW(t,-0.6);
      break;
   }

   r_shock = POW(Es/rho0,1.0/(2.0+nu)) * POW(t,2.0/(2.0+nu));
   rho_shock = ((g+1.0)/(g-1.0))*rho0;
   p_shock   = (2.0*rho0*SQR(shock_speed)) / (g+1.0);
   u_shock   = (2.0*shock_speed) / (g+1.0);
}




void calc_beta(REAL V, REAL beta[])
{
   beta[0] = (nu+2.0)*(g+1.0)*V/4.0;
   beta[1] = ((g+1.0)/(g-1.0)) * ( (nu+2.0)*g*V/2.0 - 1.0 );
   beta[2] = ((nu+2.0)*(g+1.0)) / ( (nu+2.0)*(g+1.0) -2.0*(2.0 + nu*(g-1.0)) );
   beta[2] *= 1.0 - (2.0 + nu*(g-1.0))*V/2.0;
   beta[3] = ((g+1.0)/(g-1.0)) * (1.0 - (nu+2.0)*V/2.0);
}



void gen_V()
{
   REAL dV;
   INT i;

   n = MAX_POINTS;

   dV = (V_max - V_min) / (MAX_POINTS-1);

   for (i=0; i<(MAX_POINTS-1); i++)
      Vs[i] = (i+0.5)*dV + V_min;
}



void gen_r_V()
{
   REAL beta[4];
   REAL fval;
   INT i;

   for (i=0; i<n; i++) {
      calc_beta(Vs[i],beta);
      fval =  (-2.0/(2.0+nu))*LN(beta[0]);
      fval += -a[2]*LN(beta[1]);
      fval += -a[1]*LN(beta[2]);
      fval =  EXP(fval);
 
      fval *= r_shock;
      r[i] = fval;
   }

   r[n-1] = r_shock+1.0e-7;
}



void gen_VV()
{
   REAL dVV;
   INT i;

   dVV = (V_max - V_min) / MAX_VV_POINTS;

   for (i=0; i<MAX_VV_POINTS; i++)
      VV[i] = (i+0.5)*dVV + V_min;
}


void gen_rr()
{
   REAL beta[4];
   REAL fval;
   INT i;

   for (i=0; i<MAX_VV_POINTS; i++) {
      calc_beta(VV[i],beta);
      fval =  (-2.0/(2.0+nu))*LN(beta[0]);
      fval += -a[2]*LN(beta[1]);
      fval += -a[1]*LN(beta[2]);
      fval =  EXP(fval);
 
      fval *= r_shock;
      rr[i] = fval;
   }
}


REAL find_VV(REAL rval, REAL *V1, REAL *V2)
{
   REAL lambda, VV_val;
   INT i;

   VV_val = V_max;

   if (rval <= rr[0]) {
      lambda = rval/rr[0];
      VV_val = lambda*VV[0] + (1.0 - lambda)*V_min;      

      *V1 = V_min;
      *V2 = VV[0];

      return(VV_val);
   }

   if (rval >= rr[MAX_VV_POINTS-1]) {
      lambda = (rval - rr[MAX_VV_POINTS-1])/(r_shock - rr[MAX_VV_POINTS-1]);
      VV_val = lambda*V_max + (1.0 - lambda)*rr[MAX_VV_POINTS-1];      

      *V1 = VV[MAX_VV_POINTS-1];
      *V2 = V_max;

      return(VV_val);
   }


   for (i=1; i<MAX_VV_POINTS; i++) 
      if (rval < rr[i]) {
	 lambda = (rval - rr[i-1])/(rr[i] - rr[i-1]);
	 VV_val = lambda*VV[i] + (1.0-lambda)*VV[i-1];

	 *V1 = VV[i-1];
	 *V2 = VV[i];
	 break;
      }

   return(VV_val);
}



REAL f(REAL V, REAL rval, REAL rs)
{
   REAL beta[4];
   REAL fval;

   calc_beta(V,beta);
   fval =  (-2.0/(2.0+nu))*LN(beta[0]);
   fval += -a[2]*LN(beta[1]);
   fval += -a[1]*LN(beta[2]);
   fval =  EXP(fval);
   fval =  (rval/rs) - fval;
   
   return(fval);
}



REAL calc_V_val(REAL rval, REAL rs)
{
   REAL V1, V2, Vm, V;
   REAL f1, f2, fm;
   INT i;
      
   if (ABS(rval)<1.0e-10)
      return(V_min);
   if (ABS(rval-rs)/rs < 1.0e-10)
      return(V_max);
  
   Vm = find_VV(rval,&V1,&V2);

   f1 = f(V1,rval,r_shock);
   f2 = f(V2,rval,r_shock);
   fm = f(Vm,rval,r_shock);

   for (i=0; i<MAX_V_ITS; i++) {
      if (f1*fm<0.0) {
	 V2 = Vm;
	 f2 = fm;
      }
      else {
	 V1 = Vm;
	 f1 = fm;
      }

      Vm = 0.5*(V1 + V2);
      fm = f(Vm,rval,r_shock);

      if ((ABS(V2-V1)/Vm) < 1.0e-9)
	 break;
   }

   V = Vm;

   return(V);
}



void calc_Vs()
{
   INT i;
   
   for (i=0; i<n; i++) {
      if (r[i] > r_shock)
	 Vs[i] = -1.0;
      else 
	 Vs[i] = calc_V_val(r[i],r_shock);
   }
}	
	


REAL calc_rho_val(REAL beta[])
{
   REAL rho_val;
 
   rho_val =  a[3]*LN(beta[1]);
   rho_val += a[5]*LN(beta[3]);
   rho_val += a[4]*LN(beta[2]);
   rho_val =  EXP(rho_val);
   rho_val *= rho_shock;

   return(rho_val);
}



REAL calc_p_val(REAL beta[])
{
   REAL p_val;

   p_val =  (2.0*nu/(2.0+nu))*LN(beta[0]);
   p_val += (a[5]+1.0)*LN(beta[3]);
   p_val += (a[4]-2.0*a[1])*LN(beta[2]);
   p_val =  EXP(p_val);
   p_val *= p_shock;

   return(p_val);
}



REAL calc_u_val(REAL beta[])
{
   REAL u_val;

   u_val = beta[0];
   u_val *= u_shock/r_shock;

   return(u_val);
}



void calc_sol_V_min(REAL rval, REAL Eval, REAL *rho_min, REAL *p_min, REAL *u_min)
{
   REAL val;

   val = POW(rval,nu/(g-1.0));
   val *= POW(t,-2.0*nu/((nu+2.0)*(g-1.0)));
   val *= POW(Eval/rho0,-nu/((nu+2.0)*(g-1.0)));
   
   val *= k1*rho0;
   *rho_min = val;

   val =  POW(t,-2.0*nu/(nu+2.0));
   val *= POW(Eval/rho0,2.0/(nu+2.0));
   val *= k2*rho0;
   *p_min = val;

   val = 2.0*rval / ((2.0+nu)*g*t);
   *u_min = val;
}



void calc_sol(REAL Es)
{
   REAL beta[4];
   INT i;
   
   for (i=0; i<n; i++) {
      if (r[i] >= r_shock) {
	 rho[i] = rho0;
	 p[i] = 0.0;
	 u[i] = 0.0;
      }
      else if (ABS((Vs[i]-V_min)/V_min) < V_MIN_TOL)
	 calc_sol_V_min(r[i],Es,&rho[i],&p[i],&u[i]);
      else {
	 calc_beta(Vs[i],beta);
	 rho[i]= calc_rho_val(beta);
	 p[i]  = calc_p_val(beta);
	 u[i]  = r[i]*calc_u_val(beta);
      }
   }
}



REAL calc_E_total()
{
   REAL r1, r2, dVol;
   REAL KE_sum, IE_sum, Etot_sum;
   INT i;

   KE_sum   = 0.0;
   IE_sum   = 0.0;
   Etot_sum = 0.0;

   r1 = 0.0;

   switch (geo) {
   case 1:
      for (i=0; i<(n-1); i++) {
	 r2 = 0.5*(r[i+1]+r[i]);
	 dVol = r2 - r1;

	 KE[i] = 2.0 * (0.5 * rho[i] * SQR(u[i])) * dVol;
	 IE[i] = 2.0 * (p[i] / (g-1.0)) * dVol;
	 Etot[i] = KE[i] + IE[i];

	 KE_sum   += KE[i];
	 IE_sum   += IE[i];
	 Etot_sum += Etot[i];

	 KE_cum[i] = KE_sum;
	 IE_cum[i] = IE_sum;
	 Etot_cum[i] = Etot_sum;

	 r1 = r2;
      }
      break;
   case 2:
      for (i=0; i<(n-1); i++) {
	 r2 = 0.5*(r[i+1]+r[i]);
	 dVol = PI*(SQR(r2)-SQR(r1));

	 KE[i] = (0.5 * rho[i] * SQR(u[i])) * dVol;
	 IE[i] = (p[i] / (g-1.0)) * dVol;
	 Etot[i] = KE[i] + IE[i];

	 KE_sum   += KE[i];
	 IE_sum   += IE[i];
	 Etot_sum += Etot[i];

	 KE_cum[i] = KE_sum;
	 IE_cum[i] = IE_sum;
	 Etot_cum[i] = Etot_sum;
	
	 r1 = r2;
      }
      break;
   case 3:
      for (i=0; i<(n-1); i++) {
	 r2 = 0.5*(r[i+1]+r[i]);
	 dVol = 4.0*PI*(CUBE(r2)-CUBE(r1))/3.0;

	 KE[i] = (0.5 * rho[i] * SQR(u[i])) * dVol;
	 IE[i] = (p[i] / (g-1.0)) * dVol;
	 Etot[i] = KE[i] + IE[i];
	 
	 KE_sum   += KE[i];
	 IE_sum   += IE[i];
	 Etot_sum += Etot[i];

	 KE_cum[i] = KE_sum;
	 IE_cum[i] = IE_sum;
	 Etot_cum[i] = Etot_sum;

	 r1 = r2;
      }
      break;
   }

   return(Etot_sum);
}



REAL calc_sedov(REAL alpha)
{
   REAL E_diff;

   calc_shock(alpha);
   if (rmode) {
      gen_rr();
      calc_Vs();
   }
   else 
      gen_r_V();
   calc_sol(E0/alpha);

   E_diff = E0 - calc_E_total();
   return(E_diff);
}



REAL calc_para_sedov(REAL Es)
{
   REAL E_total, alpha;

   calc_para_shock(Es);
   if (rmode) {
      gen_rr();
      calc_Vs();
   }
   else
      gen_r_V();
   calc_sol(Es);

   E_total = calc_E_total();
   alpha = E_total / Es;

   return(alpha);
}



REAL calc_alpha0()
{
   REAL z, alpha0;

   z = -(1.1409 + 0.11737*LOG(g-1.0));
   alpha0 = 0.31246*POW(g-1.0,z);

   if (geo<3)
      alpha0 = 1.2*alpha0;

   return(alpha0);
}



REAL calc_bi_alpha(REAL alpha_min, REAL alpha_max)
{
   REAL alpha1, alpha2, alpham;
   REAL f1, f2, fm;
   REAL alpha;
   INT i;

   alpha1 = alpha_min;
   alpha2 = alpha_max;
   alpham = 0.5*(alpha_min + alpha_max);
   alpha  = alpham;

   f1 = calc_sedov(alpha1);
   f2 = calc_sedov(alpha2);
   fm = calc_sedov(alpham);
   
   for (i=0; i<(MAX_ALPHA_ITS); i++) {
      if (f1*fm < 0.0) { 
	 alpha2 = alpham;
	 f2 = calc_sedov(alpha2);
      }
      else {
	 alpha1 = alpham;
	 f1 = calc_sedov(alpha1);
      }
	 
      alpham = 0.5*(alpha1 + alpha2);
      fm = calc_sedov(alpham);

      if ((ABS(fm)/E0)<E_TOL)
	 break;

      alpha = alpham;

      printf("bi alpha=%g, E_diff=%g\n",alpha,fm);
   }

   return(alpha);
}



REAL calc_sec_alpha(REAL alpha0, INT *its)
{
   REAL alpha1, alpha2, alpha;
   REAL f_min, alpha_min;
   REAL f1, f2;
   INT i;

   alpha1 = alpha0;
   alpha2 = 1.01*alpha1;
   f1 = calc_sedov(alpha1);
   f2 = calc_sedov(alpha2);

   f_min = ABS(f1);
   alpha_min = alpha1;  

   for (i=0; i<(MAX_ALPHA_ITS); i++) {
      alpha = alpha2 - f2*(alpha2 - alpha1)/(f2 - f1);

      if (alpha < 0.1) {
	 printf("alpha min trouble\n");
	 alpha = 0.1;
	 break;
      }
      if (alpha > 10.0) {
	 printf("alpha max trouble\n");
	 alpha = 10.0;
	 break;
      }

      alpha1 = alpha2;
      alpha2 = alpha;

      f1 = f2;
      f2 = calc_sedov(alpha2);
   
      if (ABS(f2) < f_min) {
	 f_min = ABS(f2);
	 alpha_min = alpha;
      }

      if ((ABS(f2)/E0)<E_TOL)
	 break;
        
      printf("sec alpha=%g, E_diff=%g\n",alpha2,f2);
   }
   
   *its = i;

   return(alpha_min);   
}



REAL calc_alpha(REAL alpha0)
{
   REAL alpha;
   INT its;

   printf("Iteration history:\n");
   printf("alpha0=%g\n",alpha0);

   if (alpha0<0.33)
      alpha = calc_bi_alpha(0.8*alpha0,1.2*alpha0);
   else
      alpha = calc_sec_alpha(alpha0,&its);

   return(alpha);
}



void dump_sol(REAL alpha) 
{
   FILE *outfile;
   INT i;
   
   outfile = fopen("sedov.dat","w");
   fprintf(outfile,"# t=%g rs=%g s=%g alpha=%g\n",t,r_shock,shock_speed,alpha);
   fprintf(outfile,"# r V-V_min rho p u\n");
   for (i=0; i<n; i++) {
      fprintf(outfile,"%g %g %g %g %g\n",r[i],Vs[i]-V_min,rho[i],p[i],u[i]);
   }
   fclose(outfile);

   outfile = fopen("sedov_nd.dat","w");
   fprintf(outfile,"# t=%g rs=%g s=%g alpha=%g\n",t,r_shock,shock_speed,alpha);
   fprintf(outfile,"# r V rho p u\n");
   for (i=0; i<n; i++) {
      fprintf(outfile,"%g %g %g %g %g\n",r[i]/r_shock,(Vs[i]-V_min)/(V_max-V_min),
	      rho[i]/rho_shock,p[i]/p_shock,u[i]/u_shock);
   }
   fclose(outfile);

   outfile = fopen("E.dat","w");
   fprintf(outfile,"# t=%g rs=%g s=%g alpha=%g\n",t,r_shock,shock_speed,alpha);
   fprintf(outfile,"# r V-V_min KE IE E KE_cum IE_cum E_cum\n");
   for (i=0; i<n; i++) {
      fprintf(outfile,"%g %g %g %g %g %g %g %g\n",r[i],Vs[i]-V_min,
	      KE[i],IE[i],Etot[i],
	      KE_cum[i],IE_cum[i],Etot_cum[i]);
   }
   fclose(outfile);

   outfile = fopen("E_nd.dat","w");
   fprintf(outfile,"# t=%g rs=%g s=%g alpha=%g\n",t,r_shock,shock_speed,alpha);
   fprintf(outfile,"# r V KE IE E KE_cum IE_cum E_cum\n");
   for (i=0; i<n; i++) {
      fprintf(outfile,"%g %g %g %g %g %g %g %g\n",r[i]/r_shock,(Vs[i]-V_min)/(V_max-V_min),
	      KE[i]/E0,IE[i]/E0,Etot[i]/E0,
	      KE_cum[i]/E0,IE_cum[i]/E0,Etot_cum[i]/E0);
   }
   fclose(outfile);
}


void dump_alpha(REAL alpha)
{
   FILE *outfile;

   printf("\nalpha=%g\n",alpha);

   outfile = fopen("alpha.dat","a");
   fprintf(outfile,"%g %g\n",g,alpha);
   fclose(outfile);
}


void dump_E_diff(REAL alpha, REAL E_diff)
{
   FILE *outfile;

   printf("\nE_diff=%g\n",E_diff);

   outfile = fopen("E_diff.dat","a");
   fprintf(outfile,"%g %g %g\n",g,alpha,E_diff);
   fclose(outfile);
}



void print_help()
{
   printf("\nsedov\n");
   printf("============================================================\n\n");
   printf("./sedov <param file> [mode] [alpha] {optfile}\n\n");
   printf("<param file>    : Parmater file (required)\n\n");
   printf("[mode] must be exactly one of the following:\n");
   printf("-r              : Use dimensional coordinates\n");
   printf("-v              : Use dimensionless coordinates\n\n");
   printf("[alpha] must be exactly one of the following:\n");
   printf("-para           : Parameteric mode\n");
   printf("-auto           : Choose initial alpha automatically\n");  
   printf("-alpha=<a>      : Set alpha to <a>\n");
   printf("-alpha0=<a0>    : Set initial value of alpha to <a0>\n\n");   
   printf("{optfile} is an optional data file:\n");
   printf("-rfile=<file>   : Dump solution for r's listed in <file>\n\n");
   printf("Example:\n\n");
   printf("./sedov sedov.param -alpha=0.851 -rfile=r.dat\n\n");
   printf("============================================================\n\n");
}



BOOLEAN get_argv_sparam(char ins[], const char form[], char outs[])
{
   INT status;

   status = sscanf(ins,form,outs);

   if ((status==0)||(status==EOF)) 
      return(FALSE);
   else
      return(TRUE);
}



BOOLEAN get_argv_gparam(char ins[], const char form[], REAL *out)
{
   INT status;

   status = sscanf(ins,form,out);

   if ((status==0)||(status==EOF)) 
      return(FALSE);
   else
      return(TRUE);
}



int main(int argc, char *argv[])
{
   char rfilename[160];
   REAL alpha0, alpha;
   REAL E_diff;
   INT i;
   
   printf("Running:\n");
   for (i=0; i<argc; i++)
      printf("%s ",argv[i]);
   printf("\n\n");

   if (argc < 4) {
      print_help();
      printf("Error: Not enough command line parameters\n");
      exit(1);
   }

   read_params(argv[1]);
   
   calc_a();
   calc_k();
   calc_V_min_max();
   gen_VV();

   if (strncmp(argv[2],"-r",80)==0) {
      gen_r_uniform();
      rmode = TRUE;
   }
   else if (strncmp(argv[2],"-v",80)==0) {
      gen_V();
      rmode = FALSE;
   }
   else {
      printf("Error: Unknown parameter %s\n",argv[2]);
      exit(1);
   }

   if (strncmp(argv[3],"-para",80)==0)
      alpha = calc_para_sedov(E0);
   else if (strncmp(argv[3],"-auto",80)==0) {
      alpha0 = calc_alpha0();
      alpha = calc_alpha(alpha0);
   }
   else {
      if (get_argv_gparam(argv[3],"-alpha=%lg",&alpha)) {
	 E_diff = calc_sedov(alpha);
	 dump_E_diff(alpha,E_diff);
      }
      else if (get_argv_gparam(argv[3],"-alpha0=%lg",&alpha0))
	 alpha = calc_alpha(alpha0);
      else {
	 printf("Error: Could not parse %s\n",argv[3]);
	 exit(1);
      }
   }

   if (argc>4) {
      if (!get_argv_sparam(argv[4],"-rfile=%s",rfilename)) {
	 printf("Error: Could not parse %s\n",argv[4]);
	 exit(1);
      }
      rmode = TRUE;
      read_r(rfilename);
      E_diff = calc_sedov(alpha);  
      dump_E_diff(alpha,E_diff);
   }

   dump_sol(alpha);
   dump_alpha(alpha);

   return(0);
}
