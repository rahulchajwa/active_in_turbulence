//Pseudo-spectral turbulence with particles code and laplacian in orientation

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define Nx 512
#define Ny 512
#define N 100000
#define PI 3.1415926536

#define beta 0.175 //motility strength
#define ell 1.0   //orientation laplacian coupling

#define mu 0.01
#define famp 0.1
#define k0 3


double pi;
double dx;
double dy;
double dt;
double t_max;
double Re;
double nu;
double gammad;
double nut;
double rho;

//********************* Gaussian random number generator ********************************
double AWGN_generator()
{/* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */
 
  double temp1;
  double temp2;
  double result;
  int p;

  p = 1;

  while( p > 0 )
  {
	temp2 = ( rand() / ( (double)RAND_MAX ) ); /*  rand() function generates an
                                                       integer between 0 and  RAND_MAX,
                                                       which is defined in stdlib.h.
                                                   */

    if ( temp2 == 0 )
    {// temp2 is >= (RAND_MAX / 2)
      p = 1;
    }// end if
    else
    {// temp2 is < (RAND_MAX / 2)
       p = -1;
    }// end else

  }// end while()

  temp1 = cos( ( 2.0 * (double)PI ) * rand() / ( (double)RAND_MAX ) );
  result = sqrt( -2.0 * log( temp2 ) ) * temp1;

  return result;	// return the generated random sample to the caller

}// end AWGN_generator()


//**************** Particle trajectory *************************************** 
double F1(double t1, double x1, double y1, double p1, double p2, double U1, double U2, double nU1, double nU2, double gU1[4])
{
return(beta*p1 + U1);
}

double F2(double t1, double x1, double y1, double p1, double p2, double U1, double U2, double nU1, double nU2, double gU1[4])
{
return(beta*p2 + U2);
}

double F3(double t1, double x1, double y1, double p1, double p2, double U1, double U2, double nU1, double nU2, double gU1[4])
{
double alpha;
alpha = 1;
return(-1.0*(p1*p1 + p2*p2 -1.0)*p1 + alpha*0.5*2*gU1[0]*p1 + alpha*0.5*(gU1[1] + gU1[2])*p2 + 0.5*(gU1[1] - gU1[2])*p2 + ell*ell*nU1); //note the sign difference in the laplacian term
}

double F4(double t1, double x1, double y1, double p1, double p2, double U1, double U2, double nU1, double nU2, double gU1[4])
{
double alpha;
alpha = 1;
return(-1.0*(p1*p1 + p2*p2 -1.0)*p2 + alpha*0.5*2*gU1[3]*p2 + alpha*0.5*(gU1[1] + gU1[2])*p1 + 0.5*(gU1[2] - gU1[1])*p1 + ell*ell*nU2); //note the sign difference in the laplacian term
}

//*****************  Function for calculating velocity field ****************************
//***************************************************************************************

void velocity(double ***w_ca, double ***velo, double ***gvelo, double ***lapvelo)
{
    double l;
    int i,j, kx[Nx], K;
    for(i = 0; i < Nx; i++)
    {

        if(i < Nx/2)
            kx[i] = i;
        else
            kx[i] = -Nx+i;

    }
    kx[Nx/2] = 0;


    //double psi_cap[Nx][Ny][2], psi_x[Nx][Ny][2], psi_y[Nx][Ny][2], u_x[Nx][Ny][2], u_y[Nx][Ny][2];

      double ***psi_cap, ***psi_x, ***psi_y, ***psi_xx, ***psi_yy,***psi_xy, ***psi_yx, ***u_x, ***u_y, ***u_xx, ***u_xy, ***u_yx, ***u_yy;

      psi_cap = (double ***)malloc(Nx*sizeof(double**));
      psi_x = (double ***)malloc(Nx*sizeof(double**));
      psi_y = (double ***)malloc(Nx*sizeof(double**));
      psi_xx = (double ***)malloc(Nx*sizeof(double**));
      psi_xy = (double ***)malloc(Nx*sizeof(double**));
      psi_yx = (double ***)malloc(Nx*sizeof(double**));
      psi_yy = (double ***)malloc(Nx*sizeof(double**));
      u_x = (double ***)malloc(Nx*sizeof(double**));
      u_y = (double ***)malloc(Nx*sizeof(double**));
      u_xx = (double ***)malloc(Nx*sizeof(double**));
      u_xy = (double ***)malloc(Nx*sizeof(double**));
      u_yx = (double ***)malloc(Nx*sizeof(double**));
      u_yy = (double ***)malloc(Nx*sizeof(double**));
      

 /*     if (psi_cap == NULL || psi_x == NULL || psi_y == NULL || u_x == NULL || u_y == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	} */


for (i =0 ; i < Nx ; i++)
{
     psi_cap[i] = (double **)malloc(Ny * sizeof(double*));
     psi_x[i] = (double **)malloc(Ny * sizeof(double*));
     psi_y[i] = (double **)malloc(Ny * sizeof(double*));
     psi_xx[i] = (double **)malloc(Ny * sizeof(double*));
     psi_xy[i] = (double **)malloc(Ny * sizeof(double*));
     psi_yx[i] = (double **)malloc(Ny * sizeof(double*));
     psi_yy[i] = (double **)malloc(Ny * sizeof(double*));
     u_x[i] = (double **)malloc(Ny * sizeof(double*));
     u_y[i] = (double **)malloc(Ny * sizeof(double*));
     u_xx[i] = (double **)malloc(Ny * sizeof(double*));
     u_xy[i] = (double **)malloc(Ny * sizeof(double*));
     u_yx[i] = (double **)malloc(Ny * sizeof(double*));
     u_yy[i] = (double **)malloc(Ny * sizeof(double*));
     
 /*    if (psi_cap[i] == NULL || psi_x[i] == NULL || psi_y[i] == NULL || u_x[i] == NULL || u_y[i] == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	} */
     

      for (j =0 ; j < Ny ; j++)
     {
      psi_cap[i][j] = (double *)malloc(2 * sizeof(double));
      psi_x[i][j] = (double *)malloc(2 * sizeof(double));
      psi_y[i][j] = (double *)malloc(2 * sizeof(double));
      psi_xx[i][j] = (double *)malloc(2 * sizeof(double));
      psi_xy[i][j] = (double *)malloc(2 * sizeof(double));
      psi_yx[i][j] = (double *)malloc(2 * sizeof(double));
      psi_yy[i][j] = (double *)malloc(2 * sizeof(double));
      u_x[i][j] = (double *)malloc(2 * sizeof(double));
      u_y[i][j]= (double *)malloc(2 * sizeof(double));
      u_xx[i][j] = (double *)malloc(2 * sizeof(double));
      u_xy[i][j] = (double *)malloc(2 * sizeof(double));
      u_yx[i][j] = (double *)malloc(2 * sizeof(double));
      u_yy[i][j] = (double *)malloc(2 * sizeof(double));
      


/*	if (psi_cap[i][j] == NULL || psi_x[i][j] == NULL || psi_y[i][j] == NULL || u_x[i][j] == NULL || u_y[i][j] == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	} */

     }
}


      

	for( i = 0; i < Nx; i++)
   	 {
        	for( j = 0; j < Ny; j++)
        	{
            	        l = kx[i]*kx[i] + kx[j]*kx[j];
            		if(l != 0)
            		{
                		psi_cap[i][j][0] = -w_ca[i][j][0]/l;
                		psi_cap[i][j][1] = -w_ca[i][j][1]/l;
            		}
           		 else
            		{
                		psi_cap[i][j][0] = 0.0;
                		psi_cap[i][j][1] = 0.0;
            		}
        	}
    	}

///********** de - aliasing step **********************************************************

    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
	    K = kx[i]*kx[i] + kx[j]*kx[j]; 

            if(K>= 2*Nx*Ny/9)
            {
                
		psi_x[i][j][0] = 0.0;
                psi_x[i][j][1] = 0.0;


                psi_y[i][j][0] = 0.0;
                psi_y[i][j][1] = 0.0;

		psi_xx[i][j][0] = 0.0;
                psi_xx[i][j][1] = 0.0;


                psi_xy[i][j][0] = 0.0;
                psi_xy[i][j][1] = 0.0;

		psi_yx[i][j][0] = 0.0;
                psi_yx[i][j][1] = 0.0;


                psi_yy[i][j][0] = 0.0;
                psi_yy[i][j][1] = 0.0;


            }
            else
            {
                

		psi_x[i][j][0] = -kx[i]*psi_cap[i][j][1];
                psi_x[i][j][1] = kx[i]*psi_cap[i][j][0];


                psi_y[i][j][0] = -kx[j]*psi_cap[i][j][1];
                psi_y[i][j][1] = kx[j]*psi_cap[i][j][0];

		psi_xx[i][j][0] = -kx[i]*psi_x[i][j][1];		//for velocity strain tensor
		psi_xx[i][j][1] = kx[i]*psi_x[i][j][0];

		psi_xy[i][j][0] = -kx[j]*psi_x[i][j][1];
		psi_xy[i][j][1] = kx[j]*psi_x[i][j][0];

		psi_yx[i][j][0] = -kx[i]*psi_y[i][j][1];
		psi_yx[i][j][1] = kx[i]*psi_y[i][j][0];

		psi_yy[i][j][0] = -kx[j]*psi_y[i][j][1];
		psi_yy[i][j][1] = kx[j]*psi_y[i][j][0];


            }
        }
    }

//******************************************************************************************

	fftw_plan ifftpsix, ifftpsiy, ifftpsixx, ifftpsixy, ifftpsiyx, ifftpsiyy;

    fftw_complex *psix_r, *psix_c, *psiy_r, *psiy_c, *psixx_r, *psixx_c, *psixy_r, *psixy_c, *psiyx_r, *psiyx_c, *psiyy_r, *psiyy_c ;

    psix_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    psix_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    psiy_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    psiy_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    psixx_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    psixx_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  
    psixy_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    psixy_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    psiyx_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    psiyx_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    psiyy_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    psiyy_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    ifftpsix = fftw_plan_dft_2d(Nx, Ny, psix_c, psix_r, FFTW_BACKWARD, FFTW_MEASURE);
    ifftpsiy = fftw_plan_dft_2d(Nx, Ny, psiy_c, psiy_r, FFTW_BACKWARD, FFTW_MEASURE);

    ifftpsixx = fftw_plan_dft_2d(Nx, Ny, psixx_c, psixx_r, FFTW_BACKWARD, FFTW_MEASURE);
    ifftpsixy = fftw_plan_dft_2d(Nx, Ny, psixy_c, psixy_r, FFTW_BACKWARD, FFTW_MEASURE);

    ifftpsiyx = fftw_plan_dft_2d(Nx, Ny, psiyx_c, psiyx_r, FFTW_BACKWARD, FFTW_MEASURE);
    ifftpsiyy = fftw_plan_dft_2d(Nx, Ny, psiyy_c, psiyy_r, FFTW_BACKWARD, FFTW_MEASURE);

for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            psix_c[j + Ny*i][0] = psi_x[i][j][0];
            psix_c[j + Ny*i][1] = psi_x[i][j][1];
            
            psiy_c[j + Ny*i][0] = psi_y[i][j][0];
            psiy_c[j + Ny*i][1] = psi_y[i][j][1];
            
            psixx_c[j + Ny*i][0] = psi_xx[i][j][0];
            psixx_c[j + Ny*i][1] = psi_xx[i][j][1];
            
	    psixy_c[j + Ny*i][0] = psi_xy[i][j][0];
            psixy_c[j + Ny*i][1] = psi_xy[i][j][1];
	    
            psiyx_c[j + Ny*i][0] = psi_yx[i][j][0];
            psiyx_c[j + Ny*i][1] = psi_yx[i][j][1];
            
            psiyy_c[j + Ny*i][0] = psi_yy[i][j][0];
            psiyy_c[j + Ny*i][1] = psi_yy[i][j][1];

        }
    }

    fftw_execute(ifftpsix);
    fftw_execute(ifftpsiy);
    fftw_execute(ifftpsixx);
    fftw_execute(ifftpsixy);
    fftw_execute(ifftpsiyx);
    fftw_execute(ifftpsiyy);

    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            u_y[i][j][0] = -psix_r[j + Ny*i][0];
            u_y[i][j][1] = -psix_r[j + Ny*i][1];

            u_x[i][j][0] = psiy_r[j + Ny*i][0];
            u_x[i][j][1] = psiy_r[j + Ny*i][1];

	    u_xx[i][j][0] = psiyx_r[j + Ny*i][0];
            u_xx[i][j][1] = psiyx_r[j + Ny*i][1];

            u_xy[i][j][0] = psiyy_r[j + Ny*i][0];
            u_xy[i][j][1] = psiyy_r[j + Ny*i][1];

	    u_yx[i][j][0] = -psixx_r[j + Ny*i][0];
            u_yx[i][j][1] = -psixx_r[j + Ny*i][1];

            u_yy[i][j][0] = -psixy_r[j + Ny*i][0];
            u_yy[i][j][1] = -psixy_r[j + Ny*i][1];        

        }
    }

    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {

	    velo[i][j][0] = u_x[i][j][0];
	
            velo[i][j][1] = u_y[i][j][0];

	    gvelo[i][j][0] = u_xx[i][j][0];
	    gvelo[i][j][1] = u_xy[i][j][0];
	    gvelo[i][j][2] = u_yx[i][j][0];
            gvelo[i][j][3] = u_yy[i][j][0];

        }
    }
    
    
   //using same psi pointers to get the laplacian of the velocity
   
   

///********** de - aliasing step **********************************************************

    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
	    K = kx[i]*kx[i] + kx[j]*kx[j]; 

            if(K>= 2*Nx*Ny/9)
            {
                
		psi_x[i][j][0] = 0.0;
                psi_x[i][j][1] = 0.0;


                psi_y[i][j][0] = 0.0;
                psi_y[i][j][1] = 0.0;


            }
            else
            {
                

		psi_x[i][j][0] = K*kx[i]*psi_cap[i][j][1];    //note the reversed sign since laplacian acquires an overall negative sign in fourier space
                psi_x[i][j][1] = -K*kx[i]*psi_cap[i][j][0];


                psi_y[i][j][0] = K*kx[j]*psi_cap[i][j][1];
                psi_y[i][j][1] = -K*kx[j]*psi_cap[i][j][0];


            }
        }
    }
    
    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            psix_c[j + Ny*i][0] = psi_x[i][j][0];
            psix_c[j + Ny*i][1] = psi_x[i][j][1];
            
            psiy_c[j + Ny*i][0] = psi_y[i][j][0];
            psiy_c[j + Ny*i][1] = psi_y[i][j][1];

        }
    }

    fftw_execute(ifftpsix);
    fftw_execute(ifftpsiy);
    
    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            u_y[i][j][0] = -psix_r[j + Ny*i][0];
            u_y[i][j][1] = -psix_r[j + Ny*i][1];

            u_x[i][j][0] = psiy_r[j + Ny*i][0];
            u_x[i][j][1] = psiy_r[j + Ny*i][1];


        }
    }

    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {

	    lapvelo[i][j][0] = u_x[i][j][0];
	
            lapvelo[i][j][1] = u_y[i][j][0];

        }
    }
    
   
   
   


    fftw_destroy_plan(ifftpsix);
    fftw_destroy_plan(ifftpsiy);
    fftw_destroy_plan(ifftpsixx);
    fftw_destroy_plan(ifftpsixy);
    fftw_destroy_plan(ifftpsiyx);
    fftw_destroy_plan(ifftpsiyy);
    fftw_free(psix_c);
    fftw_free(psix_r);
    fftw_free(psiy_c);
    fftw_free(psiy_r);
    fftw_free(psixx_c);
    fftw_free(psixx_r);
    fftw_free(psixy_c);
    fftw_free(psixy_r);
    fftw_free(psiyx_c);
    fftw_free(psiyx_r);
    fftw_free(psiyy_c);
    fftw_free(psiyy_r);

// deallocate triple pointer memory

	for (i = 0; i < Nx; i++) 
	{
		for (j = 0; j < Ny; j++)
		{
			free(psi_cap[i][j]);
			free(psi_x[i][j]);
    			free(psi_y[i][j]);
    			free(u_x[i][j]);
    			free(u_y[i][j]);
			free(u_xx[i][j]);
    			free(u_xy[i][j]);
			free(u_yx[i][j]);
    			free(u_yy[i][j]);
			free(psi_xx[i][j]);
    			free(psi_xy[i][j]);
			free(psi_yx[i][j]);
    			free(psi_yy[i][j]);
		}

			free(psi_cap[i]);
			free(psi_x[i]);
    			free(psi_y[i]);
    			free(u_x[i]);
    			free(u_y[i]);
			free(u_xx[i]);
    			free(u_xy[i]);
			free(u_yx[i]);
    			free(u_yy[i]);
			free(psi_xx[i]);
    			free(psi_xy[i]);
			free(psi_yx[i]);
    			free(psi_yy[i]);
	}

    free(psi_cap);
    free(psi_x);
    free(psi_y);
    free(u_x);
    free(u_y);
    free(u_xx);
    free(u_xy);
    free(u_yx);
    free(u_yy);
    free(psi_xx);
    free(psi_xy);
    free(psi_yx);
    free(psi_yy);

}

//**********************************************************************************
//*************** Function for convection step *************************************

double convection(double ***psi_x,double ***psi_y,double ***w_x,double ***w_y,double ***c1,double ***c2)
{

    int i, j;
	
    //double psixx[Nx][Ny][2], psiyy[Nx][Ny][2], w_xx[Nx][Ny][2], w_yy[Nx][Ny][2];

    double ***psixx, ***psiyy, ***w_xx, ***w_yy;

      psixx = (double ***)malloc(Nx*sizeof(double**));
      psiyy = (double ***)malloc(Nx*sizeof(double**));
      w_xx = (double ***)malloc(Nx*sizeof(double**));
      w_yy = (double ***)malloc(Nx*sizeof(double**));

/*     if (psixx == NULL || psiyy == NULL || w_xx == NULL || w_yy == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	} */


  for (i =0 ; i < Nx ; i++)
  {
     psixx[i] = (double **)malloc(Ny * sizeof(double*));
     psiyy[i] = (double **)malloc(Ny * sizeof(double*));
     w_xx[i] = (double **)malloc(Ny * sizeof(double*));
     w_yy[i] = (double **)malloc(Ny * sizeof(double*));

/*     if (psixx[i] == NULL || psiyy[i] == NULL || w_xx[i] == NULL || w_yy[i] == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	} */
     

      for (j =0 ; j < Ny ; j++)
     {
      psixx[i][j] = (double *)malloc(2 * sizeof(double));
      psiyy[i][j] = (double *)malloc(2 * sizeof(double));
      w_xx[i][j] = (double *)malloc(2 * sizeof(double));
      w_yy[i][j]= (double *)malloc(2 * sizeof(double));

/*        if (psixx[i][j] == NULL || psiyy[i][j] == NULL || w_xx[i][j] == NULL || w_yy[i][j] == NULL) 
	  {
		fprintf(stderr, "Out of memory");
		exit(0);
	  } */
     }
  }
   

    fftw_plan ifftpsix, ifftpsiy, ifftwx, ifftwy, fftc1, fftc2;

    fftw_complex *psix_c, *psix_r, *psiy_c, *psiy_r, *wx_c, *wx_r, *wy_c, *wy_r, *c1f, *c2f, *c1r, *c2r;

    psix_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    psiy_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    psix_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    psiy_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    wx_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    wy_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    wx_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    wy_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    c1r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    c2r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    c1f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    c2f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    ifftpsix = fftw_plan_dft_2d(Nx, Ny, psix_c, psix_r, FFTW_BACKWARD, FFTW_MEASURE);
    ifftpsiy = fftw_plan_dft_2d(Nx, Ny, psiy_c, psiy_r, FFTW_BACKWARD, FFTW_MEASURE);
    ifftwx = fftw_plan_dft_2d(Nx, Ny, wx_c, wx_r, FFTW_BACKWARD, FFTW_MEASURE);
    ifftwy = fftw_plan_dft_2d(Nx, Ny, wy_c, wy_r, FFTW_BACKWARD, FFTW_MEASURE);

    fftc1 = fftw_plan_dft_2d(Nx, Ny, c1r, c1f, FFTW_FORWARD, FFTW_MEASURE);
    fftc2 = fftw_plan_dft_2d(Nx, Ny, c2r, c2f, FFTW_FORWARD, FFTW_MEASURE);

    for( i = 0; i < Nx; i++)
    {
        for( j = 0; j < Ny; j++)
        {
            psix_c[j + Ny*i][0] = psi_x[i][j][0];
            psix_c[j + Ny*i][1] = psi_x[i][j][1];
            psiy_c[j + Ny*i][0] = psi_y[i][j][0];
            psiy_c[j + Ny*i][1] = psi_y[i][j][1];

            wx_c[j + Ny*i][0] = w_x[i][j][0];
            wx_c[j + Ny*i][1] = w_x[i][j][1];
            wy_c[j + Ny*i][0] = w_y[i][j][0];
            wy_c[j + Ny*i][1] = w_y[i][j][1];
        }
    }


    fftw_execute(ifftpsix);
    fftw_execute(ifftpsiy);
    fftw_execute(ifftwx);
    fftw_execute(ifftwy);

    for( i = 0; i < Nx; i++)
    {
        for( j = 0; j < Ny; j++)
        {
            psixx[i][j][0] = psix_r[j + Ny*i][0];
            psixx[i][j][1] = psix_r[j + Ny*i][1];

            psiyy[i][j][0] = psiy_r[j + Ny*i][0];
            psiyy[i][j][1] = psiy_r[j + Ny*i][1];

            w_xx[i][j][0] = wx_r[j + Ny*i][0];
            w_xx[i][j][1] = wx_r[j + Ny*i][1];

            w_yy[i][j][0] = wy_r[j + Ny*i][0];
            w_yy[i][j][1] = wy_r[j + Ny*i][1];
        }
    }


   for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            c1r[j + Ny*i][0] = psiyy[i][j][0]*w_xx[i][j][0]/(Nx*Ny);
            c2r[j + Ny*i][0] = psixx[i][j][0]*w_yy[i][j][0]/(Nx*Ny);
            c1r[j + Ny*i][1] = 0.0;
            c2r[j + Ny*i][1] = 0.0;
        }
    }


    fftw_execute(fftc1);
    fftw_execute(fftc2);

    for( i = 0; i < Nx; i++)
    {
        for( j = 0; j < Ny; j++)
        {
            c1[i][j][0] = c1f[j + Ny*i][0];
            c1[i][j][1] = c1f[j + Ny*i][1];

            c2[i][j][0] = c2f[j + Ny*i][0];
            c2[i][j][1] = c2f[j + Ny*i][1];
        }
    }






    fftw_destroy_plan(ifftpsix);
    fftw_destroy_plan(ifftpsiy);
    fftw_destroy_plan(ifftwx);
    fftw_destroy_plan(ifftwy);
    fftw_destroy_plan(fftc1);
    fftw_destroy_plan(fftc2);

    fftw_free(psix_c);
    fftw_free(psix_r);
    fftw_free(psiy_c);
    fftw_free(psiy_r);
    fftw_free(wx_c);
    fftw_free(wx_r);
    fftw_free(wy_c);
    fftw_free(wy_r);
    fftw_free(c1f);
    fftw_free(c2f);
    fftw_free(c1r);
    fftw_free(c2r);

// deallocate triple pointer memory

	for (i = 0; i < Nx; i++) 
	{
		for (j = 0; j < Ny; j++)
		{
			free(psixx[i][j]);
    			free(psiyy[i][j]);
    			free(w_xx[i][j]);
    			free(w_yy[i][j]);
		}

			free(psixx[i]);
    			free(psiyy[i]);
    			free(w_xx[i]);
    			free(w_yy[i]);
	}


    free(psixx);
    free(psiyy);
    free(w_xx);
    free(w_yy);


}






//**************************************************************************************
//*************** Function for calculating Runge kutta step ****************************

void runge(double ***w_ca, double ***rhs)
{

    double l;
    int i, j, kx[Nx], K;
    
     for(i = 0; i < Nx; i++)
    {

        if(i < Nx/2)
            kx[i] = i;
        else
            kx[i] = -Nx+i;

    }
    kx[Nx/2] = 0;

    //double psi_cap[Nx][Ny][2], psi_x[Nx][Ny][2], psi_y[Nx][Ny][2], w_x[Nx][Ny][2], w_y[Nx][Ny][2], c1[Nx][Ny][2], c2[Nx][Ny][2];

    double ***psi_cap, ***psi_x, ***psi_y, ***w_x, ***w_y, ***c1, ***c2;
  
      psi_cap = (double ***)malloc(Nx*sizeof(double**));
      psi_x = (double ***)malloc(Nx*sizeof(double**));
      psi_y = (double ***)malloc(Nx*sizeof(double**));
      w_x = (double ***)malloc(Nx*sizeof(double**));
      w_y = (double ***)malloc(Nx*sizeof(double**));
      c1 = (double ***)malloc(Nx*sizeof(double**));
      c2 = (double ***)malloc(Nx*sizeof(double**));

/*      if (psi_cap == NULL || psi_x == NULL || psi_y == NULL || w_x == NULL || w_y == NULL || c1 == NULL || c2 == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	} */


  for (i =0 ; i < Nx ; i++)
  {

     psi_cap[i] = (double **)malloc(Ny * sizeof(double*));
     psi_x[i] = (double **)malloc(Ny * sizeof(double*));
     psi_y[i] = (double **)malloc(Ny * sizeof(double*));
     w_x[i] = (double **)malloc(Ny * sizeof(double*));
     w_y[i] = (double **)malloc(Ny * sizeof(double*));
     c1[i] = (double **)malloc(Ny * sizeof(double*));
     c2[i] = (double **)malloc(Ny * sizeof(double*));

/*     if (psi_cap[i] == NULL || psi_x[i] == NULL || psi_y[i] == NULL || w_x[i] == NULL || w_y[i] == NULL || c1[i] == NULL || c2[i] == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	} */
     

      for (j =0 ; j < Ny ; j++)
     {

      psi_cap[i][j] = (double *)malloc(2 * sizeof(double));
      psi_x[i][j] = (double *)malloc(2 * sizeof(double));
      psi_y[i][j] = (double *)malloc(2 * sizeof(double));
      w_x[i][j] = (double *)malloc(2 * sizeof(double));
      w_y[i][j]= (double *)malloc(2 * sizeof(double));
      c1[i][j]= (double *)malloc(2 * sizeof(double));
      c2[i][j]= (double *)malloc(2 * sizeof(double));

/*        if (psi_cap[i][j] == NULL || psi_x[i][j] == NULL || psi_y[i][j] == NULL || w_x[i][j] == NULL || w_y[i][j] == NULL || c1[i][j] == NULL || c2[i][j] == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	} */

     }
  }
   


    //double w_n[Nx][Ny][2];

    //double c1[Nx][Ny][2], c2[Nx][Ny][2];

    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            l = kx[i]*kx[i] + kx[j]*kx[j];
            if(l != 0)
            {
                psi_cap[i][j][0] = -w_ca[i][j][0]/l;
                psi_cap[i][j][1] = -w_ca[i][j][1]/l;
            }
            else
            {
                psi_cap[i][j][0] = 0.0;
                psi_cap[i][j][1] = 0.0;
            }
        }
    }

//*************** de-aliasing step **************************************************


    for( i = 0; i < Nx; i++)
    {
        for( j = 0; j < Ny; j++)
        {
		psi_x[i][j][0] = -kx[i]*psi_cap[i][j][1];
                psi_x[i][j][1] = kx[i]*psi_cap[i][j][0];

                w_x[i][j][0] = -kx[i]*w_ca[i][j][1];
                w_x[i][j][1] = kx[i]*w_ca[i][j][0];

                psi_y[i][j][0] = -kx[j]*psi_cap[i][j][1];
                psi_y[i][j][1] = kx[j]*psi_cap[i][j][0];

                w_y[i][j][0] = -kx[j]*w_ca[i][j][1];
                w_y[i][j][1] = kx[j]*w_ca[i][j][0];
	}
    }
	
    

    for( i = 0; i < Nx; i++)
    {
        for( j = 0; j < Ny; j++)
        {
	    K = kx[i]*kx[i] + kx[j]*kx[j];
            if(K >= 2*Nx*Ny/9)
            {

		psi_x[i][j][0] = 0.0;
                psi_x[i][j][1] = 0.0;

                w_x[i][j][0] = 0.0;
                w_x[i][j][1] = 0.0;

                psi_y[i][j][0] = 0.0;
                psi_y[i][j][1] = 0.0;

                w_y[i][j][0] = 0.0;
                w_y[i][j][1] = 0.0;
                
            }
        }
    }

//***********************************************************************************

    for( i = 0; i < Nx; i++)
    {
        for( j = 0; j < Ny; j++)
        {
            c1[i][j][0] = 0;
            c1[i][j][1] = 0;

            c2[i][j][0] = 0;
            c2[i][j][1] = 0;
        }
    }


   convection(psi_x,psi_y,w_x,w_y,c1,c2);                   //calling convection function to calculate RHS of the time stepping

  for( i = 0; i < Nx; i++)
    {
        for( j = 0; j < Ny; j++)
        {
            l = kx[i]*kx[i] + kx[j]*kx[j];
            rhs[i][j][0] = (-(nu*l + mu)*w_ca[i][j][0] - c1[i][j][0] + c2[i][j][0]);
            rhs[i][j][1] = (-(nu*l + mu)*w_ca[i][j][1] - c1[i][j][1] + c2[i][j][1]);
        }
    }


/********* de-aliasing ******************************/

for( i = 0; i < Nx; i++)
    {
        for( j = 0; j < Ny; j++)
        {
	    K = kx[i]*kx[i] + kx[j]*kx[j];
            if(K >= 2*Nx*Ny/9)
            {

		rhs[i][j][0] = 0.0;
                rhs[i][j][1] = 0.0;

                
            }
        }
    }

 

// deallocate triple pointer memory

	for (i = 0; i < Nx; i++) 
	{
		for (j = 0; j < Ny; j++)
		{
			free(psi_cap[i][j]);
			free(psi_x[i][j]);
    			free(psi_y[i][j]);
    			free(w_x[i][j]);
    			free(w_y[i][j]);
			free(c1[i][j]);
			free(c2[i][j]);
		}

			free(psi_cap[i]);
			free(psi_x[i]);
    			free(psi_y[i]);
    			free(w_x[i]);
    			free(w_y[i]);
			free(c1[i]);
			free(c2[i]);
	}


free(psi_cap);
free(psi_x);
free(psi_y);
free(w_x);
free(w_y);
free(c1);
free(c2);


}




//**********************************************************************************************
//********************** Main code *************************************************************


int main()
{

clock_t tic = clock();

pi = 4*atan(1);
dx = 2*pi/Nx;  // box size 2*pi
dy = 2*pi/Ny;  // box size 2*pi 
dt = 0.001;
t_max = 100.0;
//t_max = 1500.0;
Re = 10000;
nu = 0.00005;
gammad = 1;
nut = 0.1;
rho = 1;

int count,i,j,k, m1, n1,m2,n2, N_tot;
char file[64];
FILE *fp; // for fluid flow

N_tot = Nx*Ny;

char filename[64], enstropyfile[64], spectrumfile[64];
FILE *gp; // for particles
FILE *hp; // for enstropy
FILE *np; // for energy and enstropy spectrum

double x, y, l, t, fun, phi;

double buf[N_tot]; //for reading initial condition

double a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,d1,d2,d3,d4;

double X[N][2], U[2], lapU[2], gradU[4], P[N][2],r,theta;

//initializing particle positions and orientations
for(i =0; i<N;i++)
{
//r = 2*pi*pi*((double)rand()/((double)RAND_MAX));
//theta = 2*pi*((double)rand()/((double)RAND_MAX));
//X[i][0] = sqrt(r)*cos(theta);
//X[i][1] = sqrt(r)*sin(theta);
X[i][0] = 2*pi*((double)rand()/((double)RAND_MAX));
X[i][1] = 2*pi*((double)rand()/((double)RAND_MAX));
theta = 2*pi*((double)rand()/((double)RAND_MAX));
P[i][0] = cos(theta);
P[i][1] = sin(theta);
// X[i][0] = pi + 0.5; X[i][1] = pi + 0.5;

// P[i][0] = 0.0; P[i][1] = -0.05;

}

//double V[Nx][Ny][2]; 

//double w[Nx][Ny], psi[Nx][Ny];

//double psi[Nx][Ny];

//double w_o[Nx][Ny]; //initial vorticity

double ener_spec[Nx], ens_spec[Nx], enstropy; 		// energy and enstropy spectrum 

int kx[Nx], kr;
    for( i = 0; i < Nx; i++)
    {

        if(i < Nx/2)
            kx[i] = i;
        else
            kx[i] = -Nx+i;

    }
    kx[Nx/2] = 0;

double **w, **w_o, **Force_r;

w = (double **)malloc(Nx * sizeof(double*));
w_o = (double **)malloc(Nx * sizeof(double*));
Force_r = (double **)malloc(Nx * sizeof(double*));

if (w == NULL || w_o == NULL || Force_r == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	}

for (i =0 ; i < Nx ; i++)
 {
 w[i] = (double *)malloc(Ny * sizeof(double));
 w_o[i] = (double *)malloc(Ny * sizeof(double));
Force_r[i] = (double *)malloc(Ny * sizeof(double));

        if (w[i] == NULL || w_o[i] == NULL ) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	}

 }




//double wx[Nx][Ny], wy[Nx][Ny];

//double w_2[Nx][Ny][2],w_3[Nx][Ny][2],w_4[Nx][Ny][2], w_n[Nx][Ny][2];

//double k1[Nx][Ny][2],k2[Nx][Ny][2],k3[Nx][Ny][2],k4[Nx][Ny][2];

//double w_cap[Nx][Ny][2];

//double psi_cap[Nx][Ny][2];

double ***w_2 , ***w_3, ***w_4 , ***w_n, ***w_cap, ***Force, ***V, ***lapV, ***k1, ***k2, ***k3, ***k4, ***gradV;

w_2 = (double ***)malloc(Nx*sizeof(double**));
w_3 = (double ***)malloc(Nx*sizeof(double**));
w_4 = (double ***)malloc(Nx*sizeof(double**));
w_n = (double ***)malloc(Nx*sizeof(double**));
w_cap = (double ***)malloc(Nx*sizeof(double**));
Force = (double ***)malloc(Nx*sizeof(double**));
V = (double ***)malloc(Nx*sizeof(double**));
lapV = (double ***)malloc(Nx*sizeof(double**));
k1 = (double ***)malloc(Nx*sizeof(double**));
k2 = (double ***)malloc(Nx*sizeof(double**));
k3 = (double ***)malloc(Nx*sizeof(double**));
k4 = (double ***)malloc(Nx*sizeof(double**));
gradV = (double ***)malloc(Nx*sizeof(double**));

if (w_2 == NULL || w_3 == NULL || w_4 == NULL || w_n == NULL || w_cap == NULL || V == NULL || lapV == NULL || k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	}

for (i =0 ; i < Nx ; i++)
{
     w_2[i] = (double **)malloc(Ny * sizeof(double*));
     w_3[i] = (double **)malloc(Ny * sizeof(double*));
     w_4[i] = (double **)malloc(Ny * sizeof(double*));
     w_n[i] = (double **)malloc(Ny * sizeof(double*));
     w_cap[i] = (double **)malloc(Ny * sizeof(double*));
     Force[i] = (double **)malloc(Ny * sizeof(double*));
     V[i] = (double **)malloc(Ny * sizeof(double*));
     lapV[i] = (double **)malloc(Ny * sizeof(double*));
     k1[i] = (double **)malloc(Ny * sizeof(double*));
     k2[i] = (double **)malloc(Ny * sizeof(double*));
     k3[i] = (double **)malloc(Ny * sizeof(double*));
     k4[i] = (double **)malloc(Ny * sizeof(double*));
     gradV[i] = (double **)malloc(Ny * sizeof(double*));

        if (w_2[i] == NULL || w_3[i] == NULL || w_4[i] == NULL || w_n[i] == NULL || w_cap[i] == NULL || V[i] == NULL || lapV[i] == NULL || k1[i] == NULL || k2[i] == NULL || k3[i] == NULL || k4[i] == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	}
     

      for (j =0 ; j < Ny ; j++)
     {
      w_2[i][j] = (double *)malloc(2 * sizeof(double));
      w_3[i][j] = (double *)malloc(2 * sizeof(double));
      w_4[i][j] = (double *)malloc(2 * sizeof(double));
      w_n[i][j] = (double *)malloc(2 * sizeof(double));
      w_cap[i][j]= (double *)malloc(2 * sizeof(double));
      Force[i][j]= (double *)malloc(2 * sizeof(double));
      V[i][j] = (double *)malloc(2 * sizeof(double));
      lapV[i][j] = (double *)malloc(2 * sizeof(double));
      k1[i][j] = (double *)malloc(2 * sizeof(double));
      k2[i][j] = (double *)malloc(2 * sizeof(double));
      k3[i][j] = (double *)malloc(2 * sizeof(double));
      k4[i][j] = (double *)malloc(2 * sizeof(double));
      gradV[i][j] = (double *)malloc(4 * sizeof(double));       // velocity gradient pointer ( Vx_x, Vx_y, Vy_x, Vy_y) by formal definition of strain tensor

if (w_2[i][j] == NULL || w_3[i][j] == NULL || w_4[i][j] == NULL || w_n[i][j] == NULL || w_cap[i][j] == NULL || V[i][j] == NULL || lapV[i][j] == NULL || k1[i][j] == NULL || k2[i][j] == NULL || k3[i][j] == NULL || k4[i][j] == NULL) 
	{
		fprintf(stderr, "Out of memory");
		exit(0);
	}

     }
}


//creating fourier transform plan for vorticity

fftw_plan fftw, ifftw;

    fftw_complex *w_r, *w_c;

    w_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
    w_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

    fftw = fftw_plan_dft_2d(Nx, Ny, w_r, w_c, FFTW_FORWARD, FFTW_MEASURE);
    ifftw = fftw_plan_dft_2d(Nx, Ny, w_c, w_r, FFTW_BACKWARD, FFTW_MEASURE);

//reading initial condition from initial.dat

FILE* ptr = fopen("initial.dat","r"); 
    if (ptr==NULL) 
    { 
        printf("no such file."); 
        return 0; 
    } 
    
i=0;
    while (fscanf(ptr,"%*f %*f %*f %*f %lf %*f %*f %*f %*f",&buf[i])==1)
{ 
i=i+1;
}  


//initialising vorticity
count=0;
for(i=0; i<Nx; i++)
    {
        for(j=0; j<Ny; j++)
        {
	    w_o[i][j] = buf[count];
	    count ++;
        }
    }
//creating column for fft
    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            w_r[j + Ny*i][0] = w_o[i][j];
        }
    }


    fftw_execute(fftw);

   for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            w_cap[i][j][0] = w_c[j + Ny*i][0]/(Nx*Ny);
            w_cap[i][j][1] = w_c[j + Ny*i][1]/(Nx*Ny);
        }
    }

    //  ABC forcing in real-space
    for( i=0; i<Nx; i++)
    {
        for( j=0; j<Ny; j++)
        {
            x = dx*i;
            y = dy*j;
            Force_r[i][j] = k0*cos(k0*x);
        }
    }
    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            w_r[j + Ny*i][0] = Force_r[i][j];
        }
    }

    fftw_execute(fftw);

    //  ABC forcing in fourier-space
    for(i = 0; i < Nx; i++)
    {
        for(j = 0; j < Ny; j++)
        {
            Force[i][j][0] = w_c[j + Ny*i][0]/(Nx*Ny);
            Force[i][j][1] = w_c[j + Ny*i][1]/(Nx*Ny);
        }
    }
	

    t = 0;
    count = 0;

           sprintf(filename, "trajectory%d.dat", count);
	   gp = fopen( filename, "w+");

	for(i=0;i<N;i++)
        {
        fprintf(gp, "%lf\t%lf\t%lf\t%lf\t%lf\n", t,X[i][0],X[i][1],P[i][0],P[i][1]);
        }
        fclose(gp);

	   sprintf(enstropyfile, "Enstropy.dat");
	   hp = fopen( enstropyfile, "w+");



    do													//do-while loop begins here
    {


//******************* Runge kutta march for vorticity field  **********************************
//*********************************************************************************************
	runge(w_cap,k1);
        for(i = 0; i < Nx; i++)
        {
            for(j = 0; j < Ny; j++)
            {
                w_2[i][j][0] = w_cap[i][j][0]+ (dt/2)*(k1[i][j][0] -famp*Force[i][j][0]);
                w_2[i][j][1] = w_cap[i][j][1]+ (dt/2)*(k1[i][j][1] -famp*Force[i][j][1]);
            }
        }

        runge(w_2,k2);
        for(i = 0; i < Nx; i++)
        {
            for(j = 0; j < Ny; j++)
            {
                w_3[i][j][0] = w_2[i][j][0]+ (dt/2)*(k2[i][j][0] -famp*Force[i][j][0]);
                w_3[i][j][1] = w_2[i][j][1]+ (dt/2)*(k2[i][j][1] -famp*Force[i][j][1]);
            }
        }
        runge(w_3,k3);
        for(i = 0; i < Nx; i++)
        {
            for(j = 0; j < Ny; j++)
            {
                w_4[i][j][0] = w_3[i][j][0]+ (dt/2)*(k3[i][j][0] -famp*Force[i][j][0]);
                w_4[i][j][1] = w_3[i][j][1]+ (dt/2)*(k3[i][j][1] -famp*Force[i][j][1]);
            }
        }
                  
        runge(w_4,k4);
        for(i = 0; i < Nx; i++)
        {
            for(j = 0; j < Ny; j++)
            {
               w_n[i][j][0] = w_cap[i][j][0]+ (dt/6)*(k1[i][j][0] + 2*k2[i][j][0] + 2*k3[i][j][0] + k4[i][j][0] -famp*Force[i][j][0] );
               w_n[i][j][1] = w_cap[i][j][1]+ (dt/6)*(k1[i][j][1] + 2*k2[i][j][1] + 2*k3[i][j][1] + k4[i][j][1] -famp*Force[i][j][1] );
            }
        }


        for(i = 0; i < Nx; i++)
        {
            for(j = 0; j < Ny; j++)
            {
                w_cap[i][j][0] = w_n[i][j][0];		//feeding back for next step
                w_cap[i][j][1] = w_n[i][j][1];
            }
        }

	//if(count%100 == 0)
        //{
            for(i = 0; i < Nx; i++)
            {
                for(j = 0; j < Ny; j++)
                {
                    w_c[j + Ny*i][0] = w_cap[i][j][0];
                    w_c[j + Ny*i][1] = w_cap[i][j][1];
                }
            }

            fftw_execute(ifftw);

            for( i = 0; i < Nx; i++)
            {
                for(j = 0; j < Ny; j++)
                {
                    w[i][j] = w_r[j + Ny*i][0];
                }
            }



//*************************************************************************************************
//*********************	Calculate Velocity feild **************************************************

	velocity(w_cap,V, gradV, lapV);


//***********************************************************************************************

//************************* Particle evolution **************************************************

for(k = 0;k<N;k++)							//particle loop begins here
{

	//************	Periodic ********************
	if(X[k][0]>2*pi)
	{
	X[k][0] = X[k][0] - 2*pi;
	}
	if(X[k][0]<0)
	{
	X[k][0] = X[k][0] + 2*pi;
	}
	if(X[k][1]>2*pi)
	{
	X[k][1] = X[k][1] - 2*pi;
	}
	if(X[k][1]<0)
	{
	X[k][1] = X[k][1] + 2*pi;
	}

       //*******************************************

	
	for(i=1;i<=Nx;i++)
	{
		
 		   if(i*2*pi/Nx > X[k][0])
		   {
		   break;					// finding x grid point 
		   }
	}

	for(j = 1;j<=Ny;j++)
	{
		   if(j*2*pi/Ny > X[k][1] )
		   {
		   break;					// finding y grid point 
		   }
	}
	
	
        m1 = i-1;
	n1 = j-1;
	
	if(m1<Nx -1)
	{
	m2 = m1+1;
	}
	if(m1 == Nx -1)
	{
	m2 = 0;
	}
	if(n1<Nx -1)
	{
	n2 = n1+1;
	}
	if(n1 == Nx-1)
	{
	n2 = 0;
	}



 //     printf("%d\t%d\n", m,n);

      /* X[0] = X[0] + dt*V[m][n][0] + 0.5*dt*P[0];
       X[1] = X[1] + dt*V[m][n][1] + 0.5*dt*P[1];
       P[0] = P[0] - dt*0.5*(P[0]*P[0] + P[1]*P[1] -0.01)*P[0] + 2*dt* gradV[m][n][0]*P[0] + 2*dt*0.5* (gradV[m][n][1] + gradV[m][n][2])*P[1] + dt*0.5* (gradV[m][n][1] - gradV[m][n][2])*P[1];
       P[1] = P[1] - dt*0.5*(P[0]*P[0] + P[1]*P[1] -0.01)*P[1] + 2*dt* gradV[m][n][3]*P[1] + 2*dt*0.5* (gradV[m][n][1] + gradV[m][n][2])*P[0] + dt*0.5* (gradV[m][n][2] - gradV[m][n][1])*P[0]; */

//************* bilinear interpolation ********************************


	U[0] = ( V[m1][n1][0]*((m1+1)*dx - X[k][0])*((n1+1)*dy - X[k][1]) + V[m2][n1][0]*(X[k][0] - m1*dx)*((n1+1)*dy - X[k][1]) + V[m1][n2][0]*((m1+1)*dx - X[k][0])*(X[k][1] - n1*dy) + V[m2][n2][0]*(X[k][0] - m1*dx)*(X[k][1] - n1*dy) )/(dx*dy);

	U[1] = ( V[m1][n1][1]*((m1+1)*dx - X[k][0])*((n1+1)*dy - X[k][1]) + V[m2][n1][1]*(X[k][0] - m1*dx)*((n1+1)*dy - X[k][1]) + V[m1][n2][1]*((m1+1)*dx - X[k][0])*(X[k][1] - n1*dy) + V[m2][n2][1]*(X[k][0] - m1*dx)*(X[k][1] - n1*dy) )/(dx*dy);

	gradU[0] = ( gradV[m1][n1][0]*((m1+1)*dx - X[k][0])*((n1+1)*dy - X[k][1]) + gradV[m2][n1][0]*(X[k][0] - m1*dx)*((n1+1)*dy - X[k][1]) + gradV[m1][n2][0]*((m1+1)*dx - X[k][0])*(X[k][1] - n1*dy) + gradV[m2][n2][0]*(X[k][0] - m1*dx)*(X[k][1] - n1*dy) )/(dx*dy);

	gradU[1] = ( gradV[m1][n1][1]*((m1+1)*dx - X[k][0])*((n1+1)*dy - X[k][1]) + gradV[m2][n1][1]*(X[k][0] - m1*dx)*((n1+1)*dy - X[k][1]) + gradV[m1][n2][1]*((m1+1)*dx - X[k][0])*(X[k][1] - n1*dy) + gradV[m2][n2][1]*(X[k][0] - m1*dx)*(X[k][1] - n1*dy) )/(dx*dy);

	gradU[2] = ( gradV[m1][n1][2]*((m1+1)*dx - X[k][0])*((n1+1)*dy - X[k][1]) + gradV[m2][n1][2]*(X[k][0] - m1*dx)*((n1+1)*dy - X[k][1]) + gradV[m1][n2][2]*((m1+1)*dx - X[k][0])*(X[k][1] - n1*dy) + gradV[m2][n2][2]*(X[k][0] - m1*dx)*(X[k][1] - n1*dy) )/(dx*dy);

	gradU[3] = ( gradV[m1][n1][3]*((m1+1)*dx - X[k][0])*((n1+1)*dy - X[k][1]) + gradV[m2][n1][3]*(X[k][0] - m1*dx)*((n1+1)*dy - X[k][1]) + gradV[m1][n2][3]*((m1+1)*dx - X[k][0])*(X[k][1] - n1*dy) + gradV[m2][n2][3]*(X[k][0] - m1*dx)*(X[k][1] - n1*dy) )/(dx*dy);
	
	lapU[0] = ( lapV[m1][n1][0]*((m1+1)*dx - X[k][0])*((n1+1)*dy - X[k][1]) + lapV[m2][n1][0]*(X[k][0] - m1*dx)*((n1+1)*dy - X[k][1]) + lapV[m1][n2][0]*((m1+1)*dx - X[k][0])*(X[k][1] - n1*dy) + lapV[m2][n2][0]*(X[k][0] - m1*dx)*(X[k][1] - n1*dy) )/(dx*dy);

	lapU[1] = ( lapV[m1][n1][1]*((m1+1)*dx - X[k][0])*((n1+1)*dy - X[k][1]) + lapV[m2][n1][1]*(X[k][0] - m1*dx)*((n1+1)*dy - X[k][1]) + lapV[m1][n2][1]*((m1+1)*dx - X[k][0])*(X[k][1] - n1*dy) + lapV[m2][n2][1]*(X[k][0] - m1*dx)*(X[k][1] - n1*dy) )/(dx*dy);





//********** Runge Kutta march for particles

	a1=F1(t,X[k][0],X[k][1],P[k][0],P[k][1],U[0],U[1], lapU[0], lapU[1], gradU);
        b1=F2(t,X[k][0],X[k][1],P[k][0],P[k][1],U[0], U[1], lapU[0], lapU[1], gradU);
	c1=F3(t,X[k][0],X[k][1],P[k][0],P[k][1],U[0], U[1], lapU[0], lapU[1], gradU);
	d1=F4(t,X[k][0],X[k][1],P[k][0],P[k][1],U[0], U[1], lapU[0], lapU[1], gradU);

	a2=F1(t+dt/2, X[k][0] + dt*a1/2, X[k][1] + dt*b1/2, P[k][0] + dt*c1/2, P[k][1] + dt*d1/2, U[0], U[1], lapU[0], lapU[1], gradU);
	b2=F2(t+dt/2, X[k][0] + dt*a1/2, X[k][1] + dt*b1/2, P[k][0] + dt*c1/2, P[k][1] + dt*d1/2, U[0], U[1], lapU[0], lapU[1], gradU);
	c2=F3(t+dt/2, X[k][0] + dt*a1/2, X[k][1] + dt*b1/2, P[k][0] + dt*c1/2, P[k][1] + dt*d1/2, U[0], U[1], lapU[0], lapU[1], gradU);
	d2=F4(t+dt/2, X[k][0] + dt*a1/2, X[k][1] + dt*b1/2, P[k][0] + dt*c1/2, P[k][1] + dt*d1/2, U[0], U[1], lapU[0], lapU[1], gradU);

	a3=F1(t+dt/2, X[k][0] + dt*a2/2, X[k][1] + dt*b2/2, P[k][0] + dt*c2/2, P[k][1] + dt*d2/2, U[0], U[1], lapU[0], lapU[1], gradU );
	b3=F2(t+dt/2, X[k][0] + dt*a2/2, X[k][1] + dt*b2/2, P[k][0] + dt*c2/2, P[k][1] + dt*d2/2, U[0], U[1], lapU[0], lapU[1], gradU );
	c3=F3(t+dt/2, X[k][0] + dt*a2/2, X[k][1] + dt*b2/2, P[k][0] + dt*c2/2, P[k][1] + dt*d2/2, U[0], U[1], lapU[0], lapU[1], gradU );
	d3=F4(t+dt/2, X[k][0] + dt*a2/2, X[k][1] + dt*b2/2, P[k][0] + dt*c2/2, P[k][1] + dt*d2/2, U[0], U[1], lapU[0], lapU[1], gradU );

	a4=F1(t+dt, X[k][0] + dt*a3 , X[k][1] + dt*b3, P[k][0] + dt*c3, P[k][1] + dt*d3, U[0], U[1], lapU[0], lapU[1], gradU );
	b4=F2(t+dt, X[k][0] + dt*a3 , X[k][1] + dt*b3, P[k][0] + dt*c3, P[k][1] + dt*d3, U[0], U[1], lapU[0], lapU[1], gradU );
	c4=F3(t+dt, X[k][0] + dt*a3 , X[k][1] + dt*b3, P[k][0] + dt*c3, P[k][1] + dt*d3, U[0], U[1], lapU[0], lapU[1], gradU );
	d4=F4(t+dt, X[k][0] + dt*a3 , X[k][1] + dt*b3, P[k][0] + dt*c3, P[k][1] + dt*d3, U[0], U[1], lapU[0], lapU[1], gradU );

	X[k][0] = X[k][0] + dt*( a1 + 2*a2 + 2*a3 + a4)/6 ;
	X[k][1] = X[k][1] + dt*( b1 + 2*b2 + 2*b3 + b4)/6 ;
	P[k][0] = P[k][0] + dt*( c1 + 2*c2 + 2*c3 + c4)/6 ;
	P[k][1] = P[k][1] + dt*( d1 + 2*d2 + 2*d3 + d4)/6 ;

	//************	Periodic ********************
	if(X[k][0]>2*pi)
	{
	X[k][0] = X[k][0] - 2*pi;
	}
	if(X[k][0]<0)
	{
	X[k][0] = X[k][0] + 2*pi;
	}
	if(X[k][1]>2*pi)
	{
	X[k][1] = X[k][1] - 2*pi;
	}
	if(X[k][1]<0)
	{
	X[k][1] = X[k][1] + 2*pi;
	}


}												//particle loop ends here

	

	//if(count%2000 == 0)
	if(count%200 == 0)
	{
	
	    clock_t toc = clock();

            printf("Elapsed: %f minutes\n", (double)(toc - tic) / (60* CLOCKS_PER_SEC));

	    //sprintf(file, "%d.dat", count/2000);
	    sprintf(file, "%d.dat", count/200);
	    fp = fopen( file, "w");

	    //sprintf(filename,"trajectory%d.dat",count/2000);
	    sprintf(filename,"trajectory%d.dat",count/200);
	    gp = fopen(filename,"w");
	
	
	    for(i = 0; i < Nx; i++)
            {
                for(j = 0; j < Ny; j++)
		{
		fprintf(fp, "%lf\t%lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n",2*pi*i/Nx, 2*pi*j/Ny, V[i][j][0],V[i][j][1], w[i][j], gradV[i][j][0], gradV[i][j][1], gradV[i][j][2], gradV[i][j][3] );
		//fprintf(fp, "%lf\t%lf\t%lf\n",2*pi*i/Nx, 2*pi*j/Ny, w[i][j]);
		}
		fprintf(fp,"\n");
                    
            }

            /****** calculating and printing enstropy **********************/

	    enstropy = 0.0;

            for(i = 0; i < Nx; i++)
            {
                for(j = 0; j < Ny; j++)
                {
                    enstropy += 0.5*(w[i][j]*w[i][j])/(4*pi*pi);
                }
            }
	   fprintf(hp,"%lf\t%lf\n",t,enstropy);	

           /****************************************************************/	




		for(k=0;k<N;k++)
		{
		fprintf(gp,"%lf\t%lf\t%lf\t%lf\n",X[k][0],X[k][1],P[k][0],P[k][1]);
		}

		
	  fclose(fp);
	  fclose(gp);

	}

        //if(count%20000 == 0)
        if(count%2000 == 0)
	{

	    //sprintf(spectrumfile, "spectrum%d.dat", count/20000);
	    sprintf(spectrumfile, "spectrum%d.dat", count/2000);
	    np = fopen( spectrumfile, "w");
	    
	    for( i = 0; i < Nx; i++)
            {
		ener_spec[i] =0.0;
		ens_spec[i]=0.0;
	    }

            for( i = 0; i < Nx; i++)
            {
                for( j = 0; j < Ny; j++)
                {
                    kr = floor(sqrt(kx[i]*kx[i]+kx[j]*kx[j]));
                    ener_spec[kr] += 0.5*(w_cap[i][j][0]*w_cap[i][j][0] + w_cap[i][j][1]*w_cap[i][j][1])/(kx[i]*kx[i]+kx[j]*kx[j]);
                    ens_spec[kr] += 0.5*(w_cap[i][j][0]*w_cap[i][j][0] + w_cap[i][j][1]*w_cap[i][j][1]);
                }
            }	
	    
	    for( i = 0; i < Nx; i++)
            {
		fprintf(np,"%d\t%.15lf\t%.15lf\n",i, ener_spec[i],ens_spec[i]);	
	    }

	    fclose(np);     

	}



        t = t + dt;

        count++;

    }while(t <= t_max);	

fclose(hp);											//do-while loop ends here


//****************** deallocate fftw plan memory ****************************************

fftw_destroy_plan(fftw);
fftw_destroy_plan(ifftw);
fftw_free(w_r);
fftw_free(w_c);


//******************* deallocate triple pointer memory ************************************

for (i = 0; i < Nx; i++) 
	{
		free(w[i]);
		free(w_o[i]);
		free(Force_r[i]);
	}

free(w);
free(w_o);
free(Force_r);




for (i = 0; i < Nx; i++) 
	{
		for (j = 0; j < Ny; j++)
		{
			free(V[i][j]);
			free(lapV[i][j]);
			free(gradV[i][j]);
			free(w_2[i][j]);
			free(w_3[i][j]);
			free(w_4[i][j]);
			free(w_n[i][j]);
			free(w_cap[i][j]);
			free(Force[i][j]);
			free(k1[i][j]);
			free(k2[i][j]);
			free(k3[i][j]);
			free(k4[i][j]);
			
		}

			free(V[i]);
			free(lapV[i]);
			free(gradV[i]);
			free(w_2[i]);
			free(w_3[i]);
			free(w_4[i]);
			free(w_n[i]);
			free(w_cap[i]);
			free(Force[i]);
			free(k1[i]);
			free(k2[i]);
			free(k3[i]);
			free(k4[i]);
	}

free(V);
free(lapV);
free(gradV);
free(w_2);
free(w_3);
free(w_4);
free(w_n);
free(w_cap);
free(Force);
free(k1);
free(k2);
free(k3);
free(k4);

clock_t toc = clock();

printf("Total Elapsed: %f hours\n", (double)(toc - tic) / (3600* CLOCKS_PER_SEC));

return(0);
}

