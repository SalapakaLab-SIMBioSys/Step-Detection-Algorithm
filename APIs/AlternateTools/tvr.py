# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:01:15 2019

@author: siva_
"""

# =============================================================================
# Description
# =============================================================================

# Total variation reconstruction 
# Code source: https://cvxopt.org/examples/book/tv.html
# from Convex Optimization by Boyd & Vandenberghe

# =============================================================================


# =============================================================================
# Imports
# ============================================================================= 
from math import pi
from cvxopt import blas, lapack, solvers
from cvxopt import matrix, spmatrix, sin, mul, div, normal
solvers.options['show_progress'] = 0

import sys
sys.path.append('../../')
# ============================================================================= 



#--------------------------------------------------------------------------
# Total variation smoothing.
#
# minimize (1/2) * ||x-corr||_2^2 + delta * || D*x ||_1
#
# minimize    (1/2) * ||x-corr||_2^2 + delta * 1'*y
# subject to  -y <= D*x <= y
#
# Variables x (n), y (n-1).
#--------------------------------------------------------------------------
def tv(delta, corr):
    """
        minimize    (1/2) * ||x-corr||_2^2 + delta * sum(y)
        subject to  -y <= D*x <= y

    Variables x (n), y (n-1).
    """
    #-- For compatibility with SDA code and numpy --#
    corr = matrix(corr) 
    n = corr.size[0]
    #-----------------------------------------------#
    
    q = matrix(0.0, (2*n-1,1))
    q[:n] = -corr
    q[n:] = delta

    def P(u, v, alpha = 1.0, beta = 0.0):
        """
            v := alpha*u + beta*v
        """

        v *= beta
        v[:n] += alpha*u[:n]


    def G(u, v, alpha = 1.0, beta = 0.0, trans = 'N'):
        """
           v := alpha*[D, -I;  -D, -I] * u + beta * v  (trans = 'N')
           v := alpha*[D, -I;  -D, -I]' * u + beta * v  (trans = 'T')

        For an n-vector z, D*z = z[1:] - z[:-1].
        For an (n-1)-vector z, D'*z = [-z;0] + [0; z].
        """

        v *= beta
        if trans == 'N':
            y = u[1:n] - u[:n-1]
            v[:n-1] += alpha*(y - u[n:])
            v[n-1:] += alpha*(-y - u[n:])
        else:
            y = u[:n-1] - u[n-1:]
            v[:n-1] -= alpha * y
            v[1:n] += alpha * y
            v[n:] -= alpha * (u[:n-1] + u[n-1:])

    h = matrix(0.0, (2*(n-1),1))


    # Customized solver for KKT system with coefficient
    #
    #     [  I    0    D'   -D' ]
    #     [  0    0   -I    -I  ]
    #     [  D   -I   -D1    0  ]
    #     [ -D   -I    0    -D2 ].

    # Diagonal and subdiagonal.
    Sd = matrix(0.0, (n,1))
    Se = matrix(0.0, (n-1,1))

    def Fkkt(W):
        """
        Factor the tridiagonal matrix

             S = I + 4.0 * D' * diag( d1.*d2./(d1+d2) ) * D

        with d1 = W['di'][:n-1]**2 = diag(D1^-1)
        d2 = W['di'][n-1:]**2 = diag(D2^-1).
        """

        d1 = W['di'][:n-1]**2
        d2 = W['di'][n-1:]**2
        d = 4.0*div( mul(d1,d2), d1+d2)
        Sd[:] = 1.0
        Sd[:n-1] += d
        Sd[1:] += d
        Se[:] = -d
        lapack.pttrf(Sd, Se)

        def g(x, y, z):

            """
            Solve

                [  I   0   D'  -D' ] [x[:n]   ]    [bx[:n]   ]
                [  0   0  -I   -I  ] [x[n:]   ] =  [bx[n:]   ]
                [  D  -I  -D1   0  ] [z[:n-1] ]    [bz[:n-1] ]
                [ -D  -I   0   -D2 ] [z[n-1:] ]    [bz[n-1:] ].

            First solve

                S*x[:n] = bx[:n] + D' * ( (d1-d2) ./ (d1+d2) .* bx[n:]
                    + 2*d1.*d2./(d1+d2) .* (bz[:n-1] - bz[n-1:]) ).

            Then take

                x[n:] = (d1+d2)^-1 .* ( bx[n:] - d1.*bz[:n-1]
                         - d2.*bz[n-1:]  + (d1-d2) .* D*x[:n] )
                z[:n-1] = d1 .* (D*x[:n] - x[n:] - bz[:n-1])
                z[n-1:] = d2 .* (-D*x[:n] - x[n:] - bz[n-1:]).
            """

            # y = (d1-d2) ./ (d1+d2) .* bx[n:] +
            #     2*d1.*d2./(d1+d2) .* (bz[:n-1] - bz[n-1:])
            y = mul( div(d1-d2, d1+d2), x[n:]) + \
                mul( 0.5*d, z[:n-1]-z[n-1:] )

            # x[:n] += D*y
            x[:n-1] -= y
            x[1:n] += y

            # x[:n] := S^-1 * x[:n]
            lapack.pttrs(Sd, Se, x)

            # u = D*x[:n]
            u = x[1:n] - x[0:n-1]

            # x[n:] = (d1+d2)^-1 .* ( bx[n:] - d1.*bz[:n-1]
            #     - d2.*bz[n-1:]  + (d1-d2) .* u)
            x[n:] = div( x[n:] - mul(d1, z[:n-1]) -
                mul(d2, z[n-1:]) + mul(d1-d2, u), d1+d2 )

            # z[:n-1] = d1 .* (D*x[:n] - x[n:] - bz[:n-1])
            # z[n-1:] = d2 .* (-D*x[:n] - x[n:] - bz[n-1:])
            z[:n-1] = mul(W['di'][:n-1], u - x[n:] - z[:n-1])
            z[n-1:] = mul(W['di'][n-1:], -u - x[n:] - z[n-1:])

        return g

    return solvers.coneqp(P, q, G, h, kktsolver = Fkkt)['x'][:n]
#--------------------------------------------------------------------------