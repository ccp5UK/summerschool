#include                <stdio.h>
#include                <math.h>
#include                <stdlib.h>

#define		DELTA2      	0.25							/* 0.01 Lishchuk interface 0.0001 Gunstensen interface */
#define		XLENGTH     	131
#define		YLENGTH     	131
#define		Q           	9
#define		w0          	4.0 / 9.0
#define		w1          	1.0 / 9.0
#define		w2          	1.0 / 36.0
#define		RHO_0		    2.0
#define		EVAP_LIM	    1.0e-8
#define		DROP_RAD    	35
#define		T           	10000
#define     BETA_LKR    	0.69							/* Recolour parameter   */

void		tables				( void );
void		propagate			( void );
void		initialise			( void );
double		calc_obs			( void );				    /* Returns red mass as a check            */
void		calc_obs2			( void );        			/* Generalised collision note "trigger" for algorithm extensions.   */
void	    collide             ( void );
void 		ReColour            ( void );
void		density_data        ( void );

double		f		[XLENGTH]	[YLENGTH]  [Q];
double		r		[XLENGTH]	[YLENGTH]  [Q];
double		b		[XLENGTH]	[YLENGTH]  [Q];

int         next_x  [XLENGTH]              [Q];                                    /* Propagation look up table; periodic bc   */
int         next_y              [YLENGTH]  [Q];

int		    cx 					            [Q] = { 0,-1, 0, 1, 1, 1, 0,-1,-1};     /* Lattice basis                            */
int		    cy					            [Q] = { 0, 1, 1, 1, 0,-1,-1,-1, 0};
double		t 			              		[Q] = {w0,w2,w1,w2,w1,w2,w1,w2,w1};     /* Link weights (D2Q9 specific)             */

double		rhoN	[XLENGTH]	[YLENGTH];                                           	/* Phase field                             */
double		ux		[XLENGTH]	[YLENGTH];                                          	/* Velocity                                */
double		uy		[XLENGTH]	[YLENGTH];
double		rho		[XLENGTH]	[YLENGTH];                                            	/* Density                                 */
double		psi		[XLENGTH]	[YLENGTH];                                             	/* Stream function                         */
double		rho_r	[XLENGTH]	[YLENGTH];                                           	/* Coloured density                        */
double		rho_b	[XLENGTH]	[YLENGTH];
double		f_x		[XLENGTH]	[YLENGTH];                                          	/* Colour field / interface normal         */
double		f_y		[XLENGTH]	[YLENGTH];
double		gradx	[XLENGTH]	[YLENGTH];                                           	/* Phase field gradient                    */
double		grady	[XLENGTH]	[YLENGTH];  
double		FF_x	[XLENGTH]	[YLENGTH];                                          	/* NS level interface body force           */
double		FF_y	[XLENGTH]	[YLENGTH];   
double		K    	[XLENGTH]	[YLENGTH];   
double		back	[XLENGTH]	[YLENGTH];
int		    flag	[XLENGTH]	[YLENGTH];                                           	/* Interfacial node marker                 */

double		OMEGA_B = 1.0, OMEGA_R = 1.0;

int	main	(void)
{ 	int    t, count;
	double Red_mass;
	
	tables    ( );
    initialise( );
    for (t = 0; t < T; t++)
	{	calc_obs  (        );
		calc_obs2 (        );
        collide   (        ); 
		ReColour  (        );
		propagate (        );
		printf("time = %d\n",t);
	}
	density_data ( );
	return  0;
}

/*  Function definitions
*/

void    density_data (  )
{	int     x, y;
	FILE    *f1;
	/* Data handling for WINDOWS */
    f1 = fopen("density.csv", "w");  
    fprintf(f1,"\n\ndensity\n");
	for (y = YLENGTH - 1; y > -1; y--)
	{	for (x = 0; x < XLENGTH; x++)
            	fprintf(f1,"%lf,",rho[x][y]);
        	fprintf(f1, "\n");
    	}
	fprintf(f1,"\n\nrhoN\n");
	for (y = YLENGTH - 1; y > -1; y--)
	{	for (x = 0; x < XLENGTH; x++)
			fprintf(f1,"%lf,",rhoN[x][y]);
		fprintf(f1, "\n");
	}
	fprintf(f1,"\n\nK\n");
	for (y = YLENGTH - 1; y > -1; y--)
	{	for (x = 0; x < XLENGTH; x++)
			fprintf(f1,"%lf,", K[x][y] );
		fprintf(f1, "\n");
	}
    fprintf(f1,"\n\npsi\n");
	for (y = YLENGTH - 1; y > -1; y--)
	{	for (x = 0; x < XLENGTH; x++)
			fprintf(f1,"%lf,", /* psi[x][y] */ pow( ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y], 0.5) * 1000.0 );
		fprintf(f1, "\n");
	}    
	fflush(f1); fclose(f1);
	
    /* Data handling for GNUplot */
	f1 = fopen ("data_density.dat","w");
	for (y = 0; y < YLENGTH; y++)
	{	for (x = 0; x < XLENGTH; x++)
       		fprintf(f1,"%d\t%d\t%f\n",x,y,rho[x][y]);
       	fprintf(f1, "\n");
    }
	fflush (f1); fclose (f1);
	f1 = fopen ("data_rhoN.dat","w");
	for (y = 0; y < YLENGTH; y++)
	{	for (x = 0; x < XLENGTH; x++)
       		fprintf(f1,"%d\t%d\t%f\n",x,y,rhoN[x][y]);
       	fprintf(f1, "\n");
    }
	fflush (f1); fclose (f1);
	f1 = fopen ("data_nK.dat","w");
	for (y = 0; y < YLENGTH; y++)
	{	for (x = 0; x < XLENGTH; x++)
       		fprintf(f1,"%d\t%d\t%f\n",x,y,K[x][y]);
       	fprintf(f1, "\n");
    }
	fflush (f1); fclose (f1);
	f1 = fopen ("data_psi.dat","w");
	for (y = 0; y < YLENGTH; y++)
	{	for (x = 0; x < XLENGTH; x++)
       		fprintf(f1,"%d\t%d\t%f\n",x,y,/* psi[x][y] */ pow( ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y], 0.5) * 1000.0 );
       	fprintf(f1, "\n");
    }
	fflush (f1); fclose (f1);
    	
	return;
}

void	tables	(void)
{	int     x, y, i;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
			for (i = 0; i < Q; i++)
			{	next_x [x][i] = x + cx[i];
		 		if (next_x[x][i] > XLENGTH - 1)
					next_x [x][i] = 0;
	 			if (next_x[x][i] < 0)
					next_x [x][i] = XLENGTH - 1;
 				next_y [y][i] = y + cy[i];
			 	if (next_y[y][i] > YLENGTH - 1)
					next_y [y][i] = 0;
	 			if (next_y[y][i] < 0)
					next_y [y][i] = YLENGTH - 1;
			}
    return;
}

void	initialise (void)
{	int        i, x, y, x0 = XLENGTH / 2, y0 = YLENGTH / 2;
	double     dist_centre;
    for (x = 0; x < XLENGTH; x++)
       for (y = 0; y < YLENGTH; y++)
           for (i = 0; i < Q; i++)
               f[x][y][i] = r[x][y][i] = b[x][y][i] = back[x][y] = 0.0;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
		{	dist_centre = sqrt ( (x -x0) * (x - x0) + (y - y0) * (y - y0) );
			if (dist_centre < DROP_RAD)
				for (i = 0; i < Q; i++)
					r [x][y][i] = RHO_0 * t[i];
		}
	for (x = 0; x < XLENGTH; x++)
        for (y = 0; y < YLENGTH; y++)
            for (i = 0; i < Q; i++)
                b [x][y][i] = RHO_0 * t[i] -  r [x][y][i];
    return;
}

double	calc_obs ( void )
{	int         x, X, y, i;
	double      redmass = 0.0, totmass = 0.0;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH;y++)
		{	rho[x][y] = rho_r[x][y] = rho_b[x][y] = ux[x][y] = uy[x][y] = 0.0;
 			for (i = 0; i < Q; i++)
			{	f    [x][y][i]	 = r[x][y][i] + b[x][y][i];
				rho  [x][y]	+= f[x][y][i];
				rho_r[x][y]	+= r[x][y][i];
				rho_b[x][y]	+= b[x][y][i];
				ux   [x][y]	+= f[x][y][i] * cx[i];
				uy   [x][y]	+= f[x][y][i] * cy[i];
			}
			rhoN [x][y]  = (rho_r[x][y] - rho_b[x][y])/rho[x][y];
			ux   [x][y] /= rho  [x][y];
			uy   [x][y] /= rho  [x][y];
			redmass     += rho_r[x][y];
			totmass     += rho  [x][y];
		}		
    return redmass;
}

void	calc_obs2 ( void )
/*  Calculates observables associated with external, immersed boundary force, F */
{	int		x, y, i, n_x, n_y;
	double	rho_N_nn, dxnx, dxny, dynx, dyny, Kurv, fnorm;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
		{	gradx[x][y] = grady[x][y]  = f_x[x][y] = f_y[x][y] = FF_x [x][y] = FF_y [x][y] = K[x][y] = 0.0;
			flag[x][y] = -1;
		}
	/* \rho^N gradients and interface unit normal using "classic" stencil, which is o(4) accurate */ 
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
		{	for (i = 1; i < Q; i++)
			{	n_x = next_x[x][i]; n_y   = next_y[y][i];
				rho_N_nn = rhoN[ n_x  ][n_y  ];
				gradx[x][y] += t[i] * rho_N_nn * cx[i];
				grady[x][y] += t[i] * rho_N_nn  * cy[i];
 			}
			gradx[x][y] *= 3.0; grady[x][y] *= 3.0;
 			if ( ( fnorm = sqrt ( gradx[x][y] * gradx[x][y] + grady[x][y] * grady[x][y] ) ) > 2.5e-8 )
			{	f_x [x][y] = gradx[x][y] / fnorm;
				f_y [x][y] = grady[x][y] / fnorm;
				flag[x][y] = 1;
			}
			else
           		f_x[x][y] = f_y[x][y] = 0.0;
		}
	/* Calculates surface curvature to O(4). Note use of second numnerical derivatives  */
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
			if ( flag[x][y] == 1 )
			{	dxnx = dxny = dynx = dyny = 0.0;
				for (i = 1; i < Q; i++)
				{	n_x = next_x[x][i];  n_y = next_y[y][i];
					dxnx  += t[i] * f_x[n_x][n_y] * cx[i];
					dxny  += t[i] * f_y[n_x][n_y] * cx[i];
					dynx  += t[i] * f_x[n_x][n_y] * cy[i];
					dyny  += t[i] * f_y[n_x][n_y] * cy[i];
				}
				K   [x][y] = 3.0 * f_x[x][y]*f_y[x][y]*(dynx + dxny) - 3.0 * f_x[x][y]*f_x[x][y]*dyny - 3.0 * f_y[x][y]*f_y[x][y]*dxnx;
				FF_x[x][y] = 0.5 * DELTA2 * K[x][y] * gradx[x][y];  /* NS level interface force components */
				FF_y[x][y] = 0.5 * DELTA2 * K[x][y] * grady[x][y];                  
				ux  [x][y] +=0.5 * FF_x[x][y] / rho[x][y];         /* Guo velocity correction */
				uy  [x][y] +=0.5 * FF_y[x][y] / rho[x][y];
			}	
	return;
}

void	collide ( )
/*  Collide the single fluid, ignoring colour, using a graded or interpolated lattice viscosity.
    Uses Guo's source term ie. applies interface forces / stresses on "flagged" nodes.
*/	
{   int		x, y, i;
	double	omega, udotc, u2, f0, uu, vv;
	for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
		{	/* Interpolate separated fluids' viscosities. More important than you might think! */
 			if ( rho_r[x][y] <= EVAP_LIM )
    			omega = OMEGA_B;
			else
				if ( rho_b[x][y] > EVAP_LIM)
					omega = 2.0 * OMEGA_R * OMEGA_B / ( OMEGA_B * (1.0 + tanh (2.5 * rhoN[x][y])) + OMEGA_R * (1.0 - tanh(2.5 * rhoN[x][y])));
 				else
					omega = OMEGA_R;
			uu = ux[x][y]; 
			vv = uy[x][y];
			u2 = ux[x][y] * uu + uy[x][y] * vv;
	 		for (i = 0; i < Q; i++)
			{	udotc = uu * cx[i] + vv * cy[i];
				f0    = t[i] * rho[x][y] * (1.0 + 3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2);
				/* LBGK collision */
				f[x][y][i] += omega * ( f0-f[x][y][i] );
				/* Add source term */
				f[x][y][i] += ( 1.0 - omega / 2.0 ) * ( 3.0 * t[i] * (( cx[i] - uu) * FF_x[x][y] + ( cy[i] - vv) * FF_y[x][y])
					                                  + 9.0 * t[i] * (  cx[i] * uu               +   cy[i] * vv) * ( cx[i] * FF_x[x][y] + cy[i] * FF_y[x][y]) );
			}
		}
	return;
}

void    ReColour (void)
/* d'Ortona's segregation method */
{	int x, y, i;
    for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
        if ( flag[x][y] == 1 )
        {  r [x][y][0] = rho_r[x][y] / rho[x][y] * f[x][y][0];
           for (i = 1; i < Q; i++)
	           r [x][y][i] = rho_r[x][y] / rho[x][y] * ( f[x][y][i] + BETA_LKR * t[i] * ( rho[x][y] - rho_r[x][y]) * ( f_x[x][y] * cx[i] + f_y[x][y] * cy[i] ));
           for (i = 0; i < Q; i++)
               b[x][y][i] = f[x][y][i] - r[x][y][i];
        } 	
    /* Clean-up */
    for (x = 0; x < XLENGTH; x++)
		for (y = 0; y < YLENGTH; y++)
		{	if (rho_r[x][y] <= EVAP_LIM)
			for (i = 0; i < Q; i++)
			{	b[x][y][i] = f[x][y][i]; r[x][y][i] = 0.0;
			}
			if (rho_b[x][y] <= EVAP_LIM)
				for (i = 0; i < Q; i++)
				{	r[x][y][i] = f[x][y][i]; b[x][y][i] = 0.0;
				}
		}
    return;
}

void	propagate (void)
{	int	x, y, i, new_x, new_y;
	/* propagate red distribution function */
	for (i = 1; i < Q; i++)
 	{	for (x = 0; x < XLENGTH; x++)
			for (y = 0; y < YLENGTH; y++)
			{	new_x			  = next_x[x][i];
				new_y			  = next_y[y][i];
				back[new_x][new_y]= r     [x][y][i];
			}
 		for (x = 0; x < XLENGTH; x++)
			for (y = 0; y < YLENGTH; y++)
				r[x][y][i] = back[x][y];
	}
	/* propagate blue distribution function */
	for (i = 1; i < Q; i++)
 	{	for (x = 0; x < XLENGTH; x++)
			for (y = 0; y < YLENGTH; y++)
			{	new_x			   = next_x[x][i];
				new_y			   = next_y[y][i];
				back[new_x][new_y] = b     [x][y][i];
			}
 		for (x = 0; x < XLENGTH; x++)
			for (y = 0; y < YLENGTH; y++)
				b[x][y][i] = back[x][y];
	}
	return;
}
