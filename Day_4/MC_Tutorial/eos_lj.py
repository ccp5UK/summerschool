#!/usr/bin/env python3
# eos_lj.py

# DISCLAIMER
# (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
# The author makes no warranties about, and disclaims liability for all uses of, this software.
# The author does not recommend use of this software for any purpose.

# This Python code uses the fitting function described and parametrized in
# M Thol, G Rutkai, R Span, J Vrabec, R Lustig, Int J Thermophys 36, 25 (2015)
# Those authors also supply C++ codes (in the supplementary information of those papers)
# They are NOT responsible for this Python code, which was written independently by Michael P Allen
# A similar notation, consistent with the papers, is retained for clarity.
# Formulae for P, E/N etc in terms of the scaled free energy derivatives a_res(0,1) etc
# may be found in the above paper and in the paper dealing with the full LJ potential:
# M Thol, G Rutkai, A Koester, R Lustig, R Span, J Vrabec, J Phys Chem Ref Data 45, 023101 (2016)

"""Equation of State for Lennard-Jones pair potential with rcut=2.5."""

def power ( tau, delta, c ):
    """The power basis function."""

    import numpy as np
    
    # f[0,0] is n*(tau**t)*(delta**d)
    # f[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # f[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j

    f = np.full ( (3,3), c['n'] * (tau**c['t']) * (delta**c['d']), dtype=np.float64 )
    f[1,:] = f[1,:] * c['t']
    f[2,:] = f[2,:] * c['t'] * ( c['t'] - 1.0 )
    f[:,1] = f[:,1] * c['d']
    f[:,2] = f[:,2] * c['d'] * ( c['d'] - 1.0 )

    return f

def expon ( tau, delta, c ):
    """The exponential basis function."""

    import numpy as np
    
    # f[0,0] is n*(tau**t)*(delta**d)*exp(-delta**l)
    # f[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # f[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j

    f = np.full ( (3,3), c['n'] * (tau**c['t']) * (delta**c['d']) * np.exp(-delta**c['l']), dtype=np.float64 )
    f[1,:] = f[1,:] * c['t']
    f[2,:] = f[2,:] * c['t'] * ( c['t'] - 1.0 )
    f[:,1] = f[:,1] * (c['d']-c['l']*delta**c['l'])
    f[:,2] = f[:,2] * ( (c['d']-c['l']*delta**c['l']) * (c['d']-1.0-c['l']*delta**c['l']) - (c['l']**2)*delta**c['l'] )

    return f

def gauss ( tau, delta, c ):
    """The gaussian basis function."""

    import numpy as np

    # f[0,0] is n*(tau**t)*exp(-beta*(tau-gamma)**2)*(delta**d)*exp(-eta*(delta-epsilon)**2)
    # f[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # f[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j

    f = np.full ( (3,3), c['n']*(tau**c['t'])*np.exp(-c['beta']*(tau-c['gamma'])**2)
                      *(delta**c['d'])*np.exp(-c['eta']*(delta-c['epsilon'])**2), dtype=np.float64)
    f[1,:] = f[1,:] * ( c['t'] - 2.0*c['beta']*tau*(tau-c['gamma']) )
    f[2,:] = f[2,:] * ( ( c['t'] - 2.0*c['beta']*tau*(tau-c['gamma']) )**2 - c['t'] - 2*c['beta']*tau**2 )
    f[:,1] = f[:,1] * ( c['d'] - 2.0*c['eta']*delta*(delta-c['epsilon']) )
    f[:,2] = f[:,2] * ( ( c['d'] - 2.0*c['eta']*delta*(delta-c['epsilon']) )**2 - c['d'] - 2*c['eta']*delta**2 )

    return f

def a_res_cutshift ( temp, rho ):
    """Reduced residual free energy and scaled derivatives for Lennard-Jones potential cut-and-shifted at 2.5 sigma."""

    import numpy as np
    
    # Given temperature and density in LJ units.
    # In a[i,j], index i refers to the tau-derivative and index j to the delta-derivative
    # The derivatives are multiplied by the corresponding powers of tau and delta
    # a[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # a[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j

    temp_crit = 1.086 # Critical temperature
    rho_crit  = 0.319 # Critical density

    tau   = temp_crit / temp # Reduced inverse temperature
    delta = rho / rho_crit   # Reduced density

    a = np.zeros ( (3,3), dtype=np.float64 )

    # Coefficients taken from Table 1 of 
    # M Thol, G Rutkai, R Span, J Vrabec, R Lustig, Int J Thermophys 36, 25 (2015) 

    cp = np.empty ( 6, dtype = [ ('n',np.float64), ('t',np.float64), ('d',np.float64) ] )
    cp['n'] = [ 0.015606084, 1.7917527, -1.9613228, 1.3045604, -1.8117673, 0.15483997 ]
    cp['t'] = [ 1.000, 0.304, 0.583, 0.662, 0.870, 0.870 ]
    cp['d'] = [ 4.0,   1.0,   1.0,   2.0,   2.0,   3.0   ]
    for c in cp:
        a = a + power ( tau, delta, c )

    ce = np.empty ( 6, dtype = [ ('n',np.float64), ('t',np.float64), ('d',np.float64), ('l',np.float64) ] )
    ce['n'] = [ -0.094885204, -0.20092412,  0.11639644, -0.50607364, -0.58422807, -0.47510982 ]
    ce['t'] = [  1.250, 3.000, 1.700, 2.400, 1.960, 1.286 ]
    ce['d'] = [  5.0,   2.0,   2.0,   3.0,   1.0,   1.0   ]
    ce['l'] = [  1.0,   2.0,   1.0,   2.0,   2.0,   1.0   ]
    for c in ce:
        a = a + expon ( tau, delta, c )

    cg = np.empty ( 9, dtype = [ ('n',np.float64), ('t',np.float64), ('d',np.float64), ('eta',np.float64),
                                 ('beta',np.float64), ('gamma',np.float64), ('epsilon',np.float64)] )
    cg['n']       = [  0.0094333106, 0.30444628,  -0.0010820946, -0.099693391, 0.0091193522,
                       0.12970543,   0.023036030, -0.082671073,  -2.2497821 ]
    cg['t']       = [  3.600, 2.080, 5.240, 0.960, 1.360, 1.655, 0.900, 0.860,  3.950 ]
    cg['d']       = [  1.0,   1.0,   2.0,   3.0,   3.0,   2.0,   1.0,   2.0,    3.0   ]
    cg['eta']     = [  4.70,  1.92,  2.70,  1.49,  0.65,  1.73,  3.70,  1.90,  13.2   ]
    cg['beta']    = [ 20.0,   0.77,  0.5,   0.8,   0.4,   0.43,  8.0,   3.3,  114.0   ]
    cg['gamma']   = [  1.0,   0.5,   0.8,   1.5,   0.7,   1.6,   1.3,   0.6,    1.3   ]
    cg['epsilon'] = [  0.55,  0.7,   2.0,   1.14,  1.2,   1.31,  1.14,  0.53,   0.96  ]
    for c in cg:
        a = a + gauss ( tau, delta, c )

    return a

def eos ( temperature, density ):
    """Function to be imported, returns results as a dictionary"""
    
    import numpy as np
    
    a_res = a_res_cutshift ( temperature, density )

    # Results to be returned as a dictionary
    results={}
    # Temperature
    results['temperature'] = temperature
    # Density
    results['density'] = density
    # Pressure
    results['P'] = density * temperature * ( 1.0 + a_res[0,1] )
    # Potential energy U/N
    results['u'] = temperature * a_res[1,0]
    # Energy E/N
    results['e'] = temperature * ( 1.5 + a_res[1,0] )
    # Constant-volume heat capacity Cv/NkB
    results['c_V'] = 1.5 - a_res[2,0]
    # Constant-pressure heat capacity
    results['c_P'] = 2.5 - a_res[2,0]+(1.0+a_res[0,1]-a_res[1,1])*(1.0+a_res[0,1]-a_res[1,1])/(1.0+2.0*a_res[0,1]+a_res[0,2]) - 1.0
    # Isothermal compressibility
    results['kappa_T'] = 1.0 / ( density * temperature * ( 1.0 + 2.0 * a_res[0,1] +  a_res[0,2] ) )
    # Reduced chemical potential beta*mu
    results['betamu'] = np.log(density) + a_res[0,0] + a_res[0,1]
    # Chemical potential
    results['mu'] = temperature * ( np.log(density) + a_res[0,0] + a_res[0,1] )
    # Activity exp(betamu)
    results['z'] = density * np.exp ( a_res[0,0] + a_res[0,1] )

    return results

# The following to be executed if this file is invoked as a script

if __name__ == "__main__":

    import sys
    import numpy as np
  
    temperature = float(input('Enter temperature: '))
    density     = float(input('Enter density: '))

    # Write out parameters
    print ( "{:40}{:15.6f}".format('Temperature T', temperature    ) )
    print ( "{:40}{:15.6f}".format('Density rho',   density) )

    # Results for cut-and-shifted potential from Thol et al (2015) fitting formula
    print('')
    print('Lennard-Jones potential cut-and-shifted at 2.5 sigma')
    print('')

    a_res = a_res_cutshift ( temperature, density )

    for (i,j), aij in np.ndenumerate ( a_res ):
        if i+j > 2: # Only interested in some of the results
            continue
        print ( "{:4}{:1d}{:<35d}{:15.6f}".format('Ares', i, j, aij ) )

    P       = density * temperature * ( 1.0 + a_res[0,1] )
    kappa_T = 1.0 / ( density * temperature * ( 1.0 + 2.0 * a_res[0,1] +  a_res[0,2] ) )
    u       = temperature * a_res[1,0]
    e       = temperature * ( 1.5 + a_res[1,0] )
    c_V     = 1.5 - a_res[2,0]
    c_P     = 2.5 - a_res[2,0]+(1.0+a_res[0,1]-a_res[1,1])*(1.0+a_res[0,1]-a_res[1,1])/(1.0+2.0*a_res[0,1]+a_res[0,2]) - 1.0
    betamu  = np.log(density) + a_res[0,0] + a_res[0,1]
    mu      = temperature * ( np.log(density) + a_res[0,0] + a_res[0,1] )
    z       = density * np.exp ( a_res[0,0] + a_res[0,1] )
    print('')
    print ( "{:40}{:15.6f}".format('Pressure P',                         P       ) )
    print ( "{:40}{:15.6f}".format('Energy per atom E/N',                e       ) )
    print ( "{:40}{:15.6f}".format('Potential energy per atom PE/N',     u       ) )
    print ( "{:40}{:15.6f}".format('Heat capacity C_V/Nk_B',             c_V     ) )
    print ( "{:40}{:15.6f}".format('Heat capacity C_P/Nk_B',             c_P     ) )
    print ( "{:40}{:15.6f}".format('Isothermal compressibility kappa_T', kappa_T ) )
    print ( "{:40}{:15.6f}".format('Chemical potential mu',              mu      ) )
    print ( "{:40}{:15.6f}".format('Reduced chemical potential beta*mu', betamu  ) )
    print ( "{:40}{:15.6f}".format('Activity z = exp(beta*mu)',          z       ) )
