/*	CCP5 Mesoscopic simulation option July 2011;
    Square lid-driven cavity simulations

*/
#include                <stdio.h>
#include                <math.h>
#include                <stdlib.h>

#define 		XLENGTH		101       /* System dimensions   */
#define 		YLENGTH 	111
#define			YLENGTH_SYS 101
#define         Q           9
#define 		RHO0		2.0
#define			T		    4000
#define			UX    		0.02     /*	Lid x velocity	      */

/*	Function declarations */
void	tables                   ( void );
void 	initial_values 		     ( void );
void 	propagate 		         ( void );
void 	Calc_obs 		         ( void );
void 	collide 		         ( void );
void 	DirichletBoundaryApply1  ( void );
int     Report			         ( void );
double 	f_equ_EILBGK             ( int i, double ux, double uy, double density);
void 	stream_function		     ( void );

/*	External variables : universal scope */
int		next_x	[XLENGTH] 	            [Q];						            /* Propagation look-up tables                  */
int		next_y	           	[YLENGTH] 	[Q];
int    	Cix 				            [Q] 	= { 0,-1, 0, 1, 1, 1, 0,-1, -1};
int    	Ciy 				            [Q] 	= { 0, 1, 1, 1, 0,-1,-1,-1, 0 };
double 	t  		 		                [Q] 	= { 4.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0,1.0/9.0, 1.0/36.0, 1.0/9.0};
double 	f        [XLENGTH]	[YLENGTH]	[Q];                                   /* primary quantity; particle dist. fn          */
double 	b        [XLENGTH]	[YLENGTH];	    
double 	ux       [XLENGTH]	[YLENGTH];                                         /* velocity                                     */
double 	uy       [XLENGTH]	[YLENGTH];
double 	rho      [XLENGTH]	[YLENGTH];                                         /* density                                      */
double	psi      [XLENGTH]	[YLENGTH];                                         /* stream function                              */
double 	OMEGA = 1.00;

int main ( void )
{	int	count;
	
	initial_values 	();
	tables		    ();
	for (count = 0;count < T; count++)
	{	printf("%d\n",count);
		DirichletBoundaryApply1();
		propagate		();
		Calc_obs		();
        collide			();
	}
    propagate	         ();
	Calc_obs		     ();
	stream_function	     ();
	Report		         ();
	return	0;
} /* end of main */

/* definition of functions */

void initial_values (void)
{	int	x, y, i;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
			for (i = 0; i < Q; i++)
			f[x][y][i] = b[x][y] = t[i] * RHO0;
	return;
}

void	tables	(void)
/*  Sets propagation look-up tables. 
    Implicit periodic boundary conditions "unscrambled" by "DirichletBoundaryApply" 
*/
{   int     x, y, i;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
			for (i = 0; i < Q; i++)
			{   next_x [x][i] = x + Cix[i];
		 	    if (next_x[x][i] > XLENGTH - 1)
			    	next_x [x][i] = 0;
	 		    if (next_x[x][i] < 0)
				    next_x [x][i] = XLENGTH - 1;
 			    next_y [y][i] = y + Ciy[i];
			    if (next_y[y][i] > YLENGTH - 1)
				    next_y [y][i] = 0;
	 		    if (next_y[y][i] < 0)
				    next_y [y][i] = YLENGTH - 1;
			}
    return;
}

void 	propagate ( void )
/*	Invokes propagation look-up tables. 
	Propagates each link in turn to a background lattice, then overwrites */
{	int     x, y, i;
	for (i = 0; i < Q; i++)
	{	for (x = 0; x < XLENGTH; x++)
			for (y = 0; y < YLENGTH; y++)
				b[ next_x[x][i] ] [ next_y[y][i] ] = f[x][y][i];
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
				f[x][y][i] = b [x][y];
	}
	return;
}

void 	Calc_obs ( void )
/* Calculates macroscopic observables */
{	int     x, y, i;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
		{	ux[x][y] = uy[x][y] = rho[x][y] = 0.0;
			for (i = 0; i < Q; i++)
			{	rho[x][y] += f[x][y][i];
				ux [x][y] += f[x][y][i] * Cix[i];
				uy [x][y] += f[x][y][i] * Ciy[i];
			}
			ux[x][y] /= RHO0;
            uy[x][y] /= RHO0;
		}
	/* Impose lid motion */
	for (x = 1; x < XLENGTH-2; x++)
    {   ux[x][YLENGTH_SYS] = UX;
        uy[x][YLENGTH_SYS] = 0.0;
    }
    
	return;
}

void 	collide ( void )
/* LBGK collision or over-relaxation step for all bulk fluid sites */
{	int     x, y, i;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
			for (i = 0; i < Q; i++)	
			    f[x][y][i] += OMEGA * ( f_equ_EILBGK (i, ux[x][y], uy[x][y], rho[x][y]) - f[x][y][i] );			           
   return;
}

int	Report ( void )
{	int	x, y;
	FILE	*f1;
	/* Data handling for WINDOWS / Excel */
    f1 = fopen ("data_stfn.csv","w");
	for (y = 0; y < YLENGTH_SYS; y++)
	{	for (x = 0; x < XLENGTH; x++)
       		fprintf(f1,"%f,",psi[x][y]);
      	fprintf(f1, "\n");
    }
	fflush (f1); fclose (f1);
	f1 = fopen ("data_pressure.csv","w");
	for (y = 0; y < YLENGTH_SYS; y++)
	{	for (x = 0; x < XLENGTH; x++)
       		fprintf(f1,"%f,",rho[x][y]);
       	fprintf(f1, "\n");
    }
	fflush (f1); fclose (f1); 
    /* Data handling for GNUplot */
    f1 = fopen ("data_stfn.dat","w");
	for (y = 0; y < YLENGTH_SYS; y++)
	{	for (x = 0; x < XLENGTH; x++)
           	fprintf(f1,"%d\t%d\t%f\n",x,y,psi[x][y]);
        fprintf(f1, "\n");
    } 
	fflush (f1); fclose (f1);
	f1 = fopen ("data_pressure.dat","w");
	for (y = 0; y < YLENGTH_SYS; y++)
	{	for (x = 0; x < XLENGTH; x++)
       		fprintf(f1,"%d\t%d\t%f\n",x,y,rho[x][y]);
       	fprintf(f1, "\n");
    }
	fflush (f1); fclose (f1);
	return	0;	
}

double 	f_equ_EILBGK	(int i, double ux, double uy, double density)
/* Equilibrium distribution function for D2Q9 LBGK model */
{	double uDOTC, uDOTu;
	uDOTC = (ux * Cix[i]+uy * Ciy[i]);
	uDOTu = (ux*ux+uy*uy);
	return (t[i]*density + t[i]*RHO0*(3.0*uDOTC-1.5*uDOTu+4.5*uDOTC*uDOTC));
}
	
void DirichletBoundaryApply1 ( void )
/*    Designed to pre-condition lattice, so underlying periodic bondary conditions in subsequent call to "propagate" will redistribute the fs    
*/
{	int     x, y, xp, yp, i;
	double  temp;
    
	/*  o(1.5) accurate "mid-link" distribution function bounce-back in e.g. line x = -0.5
    /* lhs and rhs */
	xp = XLENGTH-1;
    for (y=0; y<YLENGTH;y++)
    {   yp		     = next_y[y][1];
		temp         = f[0 ][y ][1];
		f[0 ][y ][1] = f[xp][yp][5];
		f[xp][yp][5] = temp; 
		
		yp		     = next_y[y][8];
		temp         = f[0 ][y ][8];
		f[0 ][y ][8] = f[xp][yp][4];
		f[xp][yp][4] = temp; 
		      
		yp		     = next_y[y][7];
		temp         = f[0 ][ y][7];
		f[0 ][y ][7] = f[xp][yp][3];
		f[xp][yp][3] = temp;
    }
	/* top and bottom */
	yp = YLENGTH-1;
    for (x=0; x<XLENGTH;x++)
    {   xp           = next_x[x][5];
		temp         = f[x ][0 ][5];
		f[x ][0 ][5] = f[xp][yp][1];
		f[xp][yp][1] = temp;
		
		xp           = next_x[x][6];
		temp         = f[x ][0 ][6];
		f[x ][0 ][6] = f[xp][yp][2];
		f[xp][yp][2] = temp;
		 
		xp           = next_x[x][7];
		temp         = f[x ][0 ][7];
		f[x ][0 ][7] = f[xp][yp][3];
        f[xp][yp][3] = temp;   
    }
    
     /* o(1) accurate "on-node" distribution function bounce-back in e.g. lines x = 0. */
	 /* lhs and rhs 
	for (y=0; y<YLENGTH;y++)
    {   temp       = f[0][y][1];
        f[0][y][1] = f[0][y][5];
        f[0][y][5] = temp;
        temp       = f[0][y][4];
        f[0][y][4] = f[0][y][8];
        f[0][y][8] = temp;
        temp       = f[0][y][7];
        f[0][y][7] = f[0][y][3];
        f[0][y][3] = temp;   
    }
    */
    /* top and bottom 
    for (x=0; x<XLENGTH;x++)
    {   temp       = f[x][0][1];
        f[x][0][1] = f[x][0][5];
        f[x][0][5] = temp;
        temp       = f[x][0][2];
        f[x][0][2] = f[x][0][6];
        f[x][0][6] = temp;
        temp       = f[x][0][3];
        f[x][0][3] = f[x][0][7];
        f[x][0][7] = temp;   
    }
    */
    return;
}

void	stream_function	( void )
/*   Rectangular stream function for sub-domain y<YLENGTH_SYS. 
	 Psi is defined 0.0 on bottom boundary.
*/
{   int     x, y;
    for ( x=0; x < XLENGTH-1; x++)
    {   psi[x][0] = 0.0;
        for ( y = 1; y < YLENGTH_SYS; y++)
            psi[x][y] = psi[x][y-1] + 0.5 * ( ux[x][y] + ux[x][y-1]);
    }
    return;
}

