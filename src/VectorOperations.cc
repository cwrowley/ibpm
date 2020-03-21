// VectorOperations.cc
//
// Description:
// Definition of VectorOperations functions.
//
// Author(s):
// Clancy Rowley
// Zhanhua Ma
//
// Date: 15 Jul 2008
//
// $Revision$
// $LastChangedDate$
// $LastChangedBy$
// $HeadURL$

#include "Array.h"
#include "BC.h"
#include "Grid.h"
#include "Scalar.h"
#include "Flux.h"
#include "BoundaryVector.h"
#include "VectorOperations.h"
#include "NavierStokesModel.h"
#include <fftw3.h>
#include <iostream>

using Array::Array2;

namespace ibpm {

// Compute the curl of Flux q, as a Scalar object f
void Curl(const Flux& q, Scalar& f ) {
    assert( q.Nx() == f.Nx() );
    assert( q.Ny() == f.Ny() );
    assert( q.Ngrid() == f.Ngrid() );
    int nx = q.Nx();
    int ny = q.Ny();
    
    // Curl (u,v) = v_x - u_y

    // Start with finest grid, to coarsest grid
    for (int lev=0; lev<q.Ngrid(); ++lev ) {
        // compute curl at all nodes
        double dx = q.Dx(lev);
        double bydx2 = 1. / (dx * dx);
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny; ++j) {
                f(lev,i,j) = ( q(lev,Y,i,j) - q(lev,Y,i-1,j) 
                    + q(lev,X,i,j-1) - q(lev,X,i,j) ) * bydx2;
            }
        }
    }
}

Scalar Curl(const Flux& q) {
    Scalar omega( q.getGrid() );
    Curl( q, omega );
    return omega;
}

// Return the curl of Scalar f, as a Flux object q.
void Curl(const Scalar& f, Flux& q) {
    assert( f.Nx() == q.Nx() );
    assert( f.Ny() == q.Ny() );
    assert( f.Ngrid() == q.Ngrid() );
    int nx = f.Nx();
    int ny = f.Ny();

    // boundary condition object for use in computing curl on finer grids
    BC bc(nx,ny);
    
    // From coarsest grid to finest
    for (int lev=f.Ngrid()-1; lev >= 0; --lev) {
        // For outermost grid, all boundaries are zero
        if (lev == f.Ngrid()-1) {
            bc = 0.;
        }
        // Otherwise, get bc from next coarser grid
        else {
            f.getBC( lev, bc );
        }

        // X direction: u = df/dy
        
        // Compute all points except boundaries
        for (int j=1; j<ny-1; ++j) {
            for (int i=1; i<nx; ++i) {
                q(lev,X,i,j) = f(lev,i,j+1) - f(lev,i,j);
            }
        }
        // top and bottom boundaries
        for (int i=1; i<nx; ++i) {
            q(lev,X,i,0) = f(lev,i,1) - bc.bottom(i);
            q(lev,X,i,ny-1) = bc.top(i) - f(lev,i,ny-1);
        }
        // left and right boundaries
        for (int j=0; j<ny; ++j) {
            q(lev,X,0,j) = bc.left(j+1) - bc.left(j);
            q(lev,X,nx,j) = bc.right(j+1) - bc.right(j);
        }

        // Y direction: v = -df/dx
        
        // Compute all points except boundaries
        for (int j=1; j<ny; ++j) {
            for (int i=1; i<nx-1; ++i) {
                q(lev,Y,i,j) = f(lev,i,j) - f(lev,i+1,j);
            }
        }
        
        // left and right boundaries
        for (int j=1; j<ny; ++j) {
            q(lev,Y,0,j) = bc.left(j) - f(lev,1,j);
            q(lev,Y,nx-1,j) = f(lev,nx-1,j) - bc.right(j);
        }
        // top and bottom boundaries
        for (int i=0; i<nx; ++i) {
            q(lev,Y,i,0) = bc.bottom(i) - bc.bottom(i+1);
            q(lev,Y,i,ny) = bc.top(i) - bc.top(i+1);
        }
    }
}

Flux Curl(const Scalar& f) {
    Flux q( f.getGrid() );
    Curl( f, q );
    return q;
}

// Return g = L f where L is the discrete Laplacian.
// Assumes boundary values of f are zero
void Laplacian(const Scalar& f, Scalar& g) {
    assert( f.Nx() == g.Nx() );
    assert( f.Ny() == g.Ny() );
    assert( f.Ngrid() == g.Ngrid() );

    // Laplacian = - Curl( Curl( ) )
    Flux q = Curl( f );
    Curl( q, g );
    g *= -1;
    
   // // Single-grid version
//    assert( f.Ngrid() == 1 );
//    double bydx2 = 1. / (f.Dx() * f.Dx());
//    int nx = f.Nx();
//    int ny = f.Ny();
//    
//    for (int i=2; i<nx-1; ++i) {
//        for (int j=2; j<ny-1; ++j) {
//            g(0,i,j) = f(0,i-1,j) + f(0,i+1,j) + f(0,i,j-1) + f(0,i,j+1) -4*f(0,i,j);
//            g(0,i,j) *= bydx2;
//        }
//    }
//    // left and right edges
//    for (int j=2; j<ny-1; ++j) {
//        g(0,1,j) = f(0,2,j) + f(0,1,j+1) + f(0,1,j-1) - 4*f(0,1,j);
//        g(0,1,j) *= bydx2;
//        g(0,nx-1,j) = f(0,nx-2,j) + f(0,nx-1,j+1) + f(0,nx-1,j-1) - 4*f(0,nx-1,j);
//        g(0,nx-1,j) *= bydx2;
//    }
//    // top and bottom edges
//    for (int i=2; i<nx-1; ++i) {
//        g(0,i,1) = f(0,i,2) + f(0,i+1,1) + f(0,i-1,1) - 4*f(0,i,1);
//        g(0,i,1) *= bydx2;
//        g(0,i,ny-1) = f(0,i,ny-2) + f(0,i+1,ny-1) + f(0,i-1,ny-1) - 4*f(0,i,ny-1);
//        g(0,i,ny-1) *= bydx2;
//    }
//    // corners
//    g(0,1,1) = ( f(0,1,2) + f(0,2,1) - 4 * f(0,1,1) ) * bydx2;
//    g(0,nx-1,1) = ( f(0,nx-1,2) + f(0,nx-2,1) - 4 * f(0,nx-1,1) ) * bydx2;
//    g(0,1,ny-1) = ( f(0,1,ny-2) + f(0,2,ny-1) - 4 * f(0,1,ny-1) ) * bydx2;
//    g(0,nx-1,ny-1) = ( f(0,nx-2,ny-1) + f(0,nx-1,ny-2) - 4 * f(0,nx-1,ny-1) ) * bydx2;
}

void Laplacian( const Array2<double>& f,
                double dx,
                const BC& bc,
                Array2<double>& g ) {
    assert( f.Nx() == g.Nx() );
    assert( f.Ny() == g.Ny() );
    assert( dx > 0. );
    assert( bc.Nx() == f.Nx() + 1 );
    assert( bc.Ny() == f.Ny() + 1 );
    double bydx2 = 1. / (dx * dx);
    int nx = bc.Nx();
    int ny = bc.Ny();
    
    for (int i=2; i<nx-1; ++i) {
        for (int j=2; j<ny-1; ++j) {
            g(i,j) = f(i-1,j) + f(i+1,j) + f(i,j-1) + f(i,j+1) -4*f(i,j);
            g(i,j) *= bydx2;
        }
    }
    // left and right edges
    for (int j=2; j<ny-1; ++j) {
        g(1,j) = bc.left(j) + f(2,j) + f(1,j+1) + f(1,j-1) - 4*f(1,j);
        g(1,j) *= bydx2;
        g(nx-1,j) = bc.right(j) + f(nx-2,j) + f(nx-1,j+1) + f(nx-1,j-1) - 4*f(nx-1,j);
        g(nx-1,j) *= bydx2;
    }
    // top and bottom edges
    for (int i=2; i<nx-1; ++i) {
        g(i,1) = bc.bottom(i) + f(i,2) + f(i+1,1) + f(i-1,1) - 4*f(i,1);
        g(i,1) *= bydx2;
        g(i,ny-1) = bc.top(i) + f(i,ny-2) + f(i+1,ny-1) + f(i-1,ny-1) - 4*f(i,ny-1);
        g(i,ny-1) *= bydx2;
    }
    // corners
    g(1,1) = ( bc.left(1) + bc.bottom(1) + f(1,2) + f(2,1) - 4 * f(1,1) ) * bydx2;
    g(nx-1,1) = ( bc.right(1) + bc.bottom(nx-1) + f(nx-1,2) + f(nx-2,1) - 4 * f(nx-1,1) ) * bydx2;
    g(1,ny-1) = ( bc.left(ny-1) + bc.top(1) + f(1,ny-2) + f(2,ny-1) - 4 * f(1,ny-1) ) * bydx2;
    g(nx-1,ny-1) = ( bc.right(ny-1) + bc.top(nx-1) + f(nx-2,ny-1) + f(nx-1,ny-2) - 4 * f(nx-1,ny-1) ) * bydx2;        
}
    
Scalar Laplacian( const Scalar& f ) {
    Scalar g( f.getGrid() );
    Laplacian( f, g );
    return g;
}
    
    
// ~~~~~~~~~~~~~~~~~~~~~~
// Inner product of two Scalars, taken over the finest domain only
double FineGridInnerProduct( const Scalar& f, const Scalar& g ) {
    assert( f.Ngrid() == g.Ngrid() );
    assert( f.Nx() == g.Nx() );
    assert( f.Ny() == g.Ny() );
    
    double ip = 0.;
    int nx = f.Nx();
    int ny = f.Ny();
    
    // Finest grid interior points
    double dx2 = f.Dx() * f.Dx();
    for (int i = 1; i < nx; ++i) {
        for ( int j = 1; j < ny; ++j) {
            ip += f( 0 ,i, j ) * g( 0 ,i, j ) * dx2;
        }
    }
    
    return ip;
}


// ~~~~~~~~~~~~~~~~~~~~~~

// Inner product of two Scalars. 
double InnerProduct (const Scalar& f, const Scalar& g){
    assert( f.Ngrid() == g.Ngrid() );
    assert( f.Nx() == g.Nx() );
    assert( f.Ny() == g.Ny() );
    
    int nx = f.Nx();
    int ny = f.Ny();
    int nx2 = f.NxExt();  // # coarse cells outside each fine domain
    int ny2 = f.NyExt();
    double dx2 = f.Dx() * f.Dx();
    
    double ip = FineGridInnerProduct( f, g );
   
    // Coarser grids
    for (int lev=1; lev < f.Ngrid(); ++lev) {
        dx2 = f.Dx(lev) * f.Dx(lev);        
        // Interface points
        // corners
        ip += f(lev,nx2,ny2) * g(lev,nx2,ny2) * dx2 * 15./16;
        ip += f(lev,nx/2+nx2,ny2) * g(lev,nx/2+nx2,ny2) * dx2 * 15./16;
        ip += f(lev,nx2,ny/2+ny2) * g(lev,nx2,ny/2+ny2) * dx2 * 15./16;
        ip += f(lev,nx/2+nx2,ny/2+ny2) * g(lev,nx/2+nx2,ny/2+ny2) * dx2 * 15./16;
        // edges
        for (int j=ny2+1; j < ny/2 + ny2; ++j) {
            // left & right
            int i = nx2;
            ip += f(lev,i,j) * g(lev, i,j) * dx2 * 0.75;
            i = nx/2 + nx2;
            ip += f(lev,i,j) * g(lev,i,j) * dx2 * 0.75;
        }
        for (int i=nx2+1; i< nx/2 + nx2; ++i) {
            // top & bottom
            int j=ny2;
            ip += f(lev,i,j) * g(lev, i,j) * dx2 * 0.75;
            j = ny/2+ny2;
            ip += f(lev,i,j) * g(lev,i,j) * dx2 * 0.75;
        }
        // Left border
        for (int i = 1; i < nx2; ++i) {
            for ( int j = 1; j < ny; ++j) {
                ip += f(lev, i, j) * g(lev, i, j) * dx2;
            }
        }
        // Right border
        for (int i = nx/2 + nx2 + 1; i < nx; ++i ) {
            for (int j = 1; j < ny; ++j) {
                ip += f(lev, i, j) * g(lev, i, j) * dx2;
            }
        }
        for (int i = nx2; i < nx/2 + nx2 + 1; ++i ) {
            // Bottom border
            for (int j=1; j < ny2; ++ j ) {
                ip += f(lev, i, j) * g(lev, i, j) * dx2;
            }
            // Top border
            for (int j = ny/2 + ny2 + 1; j < ny; ++j) {
                ip += f(lev, i, j) * g(lev, i, j) * dx2;
            }
        }
    }
    return ip;
}

// Inner product of Flux p and Flux q.
// Note that this is not multiplied by dx * dx, since the Fluxes are already
// multiplied by these (i.e., inner product is really over *velocities*).

// Finest grid, all interior points
double FineGridInnerProduct( const Flux& p, const Flux& q ) { 
    assert( p.Ngrid() == q.Ngrid() );
    assert( p.Nx() == q.Nx() );
    assert( p.Ny() == q.Ny() );
    
    int nx = p.Nx();
    int ny = p.Ny();

    double ip = 0.;
    
    // X-fluxes
    for (int j=0; j<ny; ++j) {
        for (int i=1; i<nx; ++i){
            ip += p(0,X,i,j) * q(0,X,i,j);
        }
    }

    // Y-fluxes    
    for (int i=0; i<nx; ++i) {
        for (int j=1; j<ny; ++j){
            ip += p(0,Y,i,j) * q(0,Y,i,j);
        }
    }

    return ip;
}
    
// All grids
double InnerProduct (const Flux& p, const Flux& q){
    assert( p.Ngrid() == q.Ngrid() );
    assert( p.Nx() == q.Nx() );
    assert( p.Ny() == q.Ny() );
    int nx = p.Nx();
    int ny = p.Ny();
    int nx2 = p.NxExt();
    int ny2 = p.NyExt();
    
    double ip = FineGridInnerProduct( p, q );
    
    // X-fluxes, coarser grids
    for (int lev=1; lev < p.Ngrid(); ++lev) {
        // left and right interfaces (edges)
        for (int j=ny2; j<ny/2+ny2; ++j) {
            ip += p(lev,X,nx2,j) * q(lev,X,nx2,j) * 0.75;
            ip += p(lev,X,nx/2+nx2,j) * q(lev,X,nx/2+nx2,j) * 0.75;
        }
        // left and right coarse points
        for (int j=0; j<ny; ++j) {
            for (int i=1; i<nx2; ++i) {
                ip += p(lev,X,i,j) * q(lev,X,i,j);
            }
            for (int i=nx/2+nx2+1; i<nx; ++i) {
                ip += p(lev,X,i,j) * q(lev,X,i,j);
            }
        }
        // top and bottom coarse points
        for (int i=nx2; i<nx/2+nx2+1; ++i) {
            for (int j=0; j<ny2; ++j) {
                ip += p(lev,X,i,j) * q(lev,X,i,j);
            }
            for (int j=ny/2+ny2; j<ny; ++j) {
                ip += p(lev,X,i,j) * q(lev,X,i,j);
            }
        }
    }
    
    
    // Y-fluxes, coarser grids
    for (int lev=1; lev < p.Ngrid(); ++lev) {
        // left and right interfaces (edges)
        for (int i=nx2; i<nx/2+nx2; ++i) {
            ip += p(lev,Y,i,ny2) * q(lev,Y,i,ny2) * 0.75;
            ip += p(lev,Y,i,ny/2+ny2) * q(lev,Y,i,ny/2+ny2) * 0.75;
        }
        // left and right coarse points
        for (int i=0; i<nx; ++i) {
            for (int j=1; j<ny2; ++j) {
                ip += p(lev,Y,i,j) * q(lev,Y,i,j);
            }
            for (int j=ny/2+ny2+1; j<ny; ++j) {
                ip += p(lev,Y,i,j) * q(lev,Y,i,j);
            }
        }
        // top and bottom coarse points
        for (int j=ny2; j<ny/2+ny2+1; ++j) {
            for (int i=0; i<nx2; ++i) {
                ip += p(lev,Y,i,j) * q(lev,Y,i,j);
            }
            for (int i=nx/2+nx2; i<nx; ++i) {
                ip += p(lev,Y,i,j) * q(lev,Y,i,j);
            }
        }
    }
    
    return ip;
}

/*  Take the inner product of two vorticity fields by
    taking the L2 inner product of the first vorticity field
    with the streamfunction of the second.  This is equiv-
    alent to taking the L2 inner product of the two fluxes.
 */
// Fine grid
double FineGridVorticityInnerProduct( const Scalar& omega1, const Scalar& omega2, const NavierStokesModel& model ) { 
    Scalar psi2 = model.vorticityToStreamfunction( omega2 );
    return FineGridInnerProduct( omega1, psi2 );
}
    
// All grids
double VorticityInnerProduct( const Scalar& omega1, const Scalar& omega2, const NavierStokesModel& model ) { 
    Scalar psi2 = model.vorticityToStreamfunction( omega2 );
    return InnerProduct( omega1, psi2 );
}
    
// Return cross product of a Flux q and a Scalar f, as a Flux object.
//   q x f = ( f v, -f u )
//
// Step 1: Convert flux to velocities (u,v) at vertices
// Step 2: Compute f * v and -f * u
// Step 3: Convert the velocities to fluxes through edges

Flux CrossProduct(const Flux& q, const Scalar& f){
    assert( q.Nx() == f.Nx());
    assert( q.Ny() == f.Ny());
    assert( q.Ngrid() == f.Ngrid() );

    Scalar u( f.getGrid() );
    Scalar v( f.getGrid() );
    	
    FluxToXVelocity( q, u );
    u *= f;
    u *= -1;
    
    FluxToYVelocity( q, v );
    v *= f;
    
    Flux cross( q.getGrid() );
    VelocityToFlux( v, u, cross );      // cross = ( f v, -f u )

    return cross;
}

// Return cross product of two Flux objects q1, q2, as a Scalar object.
//   q1 x q2 = u1 v2 - u2 v1
//
// Step 1: Convert the fluxes to velocities at nodes
// Step 2: Compute u1 * v2 - u2 * v1

Scalar CrossProduct(const Flux& q1, const Flux& q2){
    assert( q1.Nx() == q2.Nx() );
    assert( q1.Ny() == q2.Ny() );
    assert( q1.Ngrid() == q2.Ngrid() );
    
    const Grid& grid = q1.getGrid();
    Scalar u( grid );
    Scalar v( grid );
    
    FluxToXVelocity( q1, u );
    FluxToYVelocity( q2, v );

    Scalar f = u * v;  // f is now u1 * v2
    
    FluxToXVelocity( q2, u );
    FluxToYVelocity( q1, v );
    
    f -= u * v;  // f is now u1 * v2 - u2 * v1
    
    f.coarsify();   // fill in overlapping grid regions
    return f;
};

void FluxToXVelocity(const Flux& q, Scalar& u) {
    assert( q.Nx() == u.Nx() );
    assert( q.Ny() == u.Ny() );
    assert( q.Ngrid() == u.Ngrid() );
    int nx = q.Nx();
    int ny = q.Ny();
    int nx2 = q.NxExt();
    int ny2 = q.NyExt();
    double oneOver2Delta = 1./ (2 * q.Dx());
    
    // Compute interior points (A)
    for (int i=1; i < nx; ++i ){
        for (int j=1; j < ny; ++j ) {
            u(0,i,j) = ( q(0,X,i,j) + q(0,X,i,j-1) ) * oneOver2Delta;
        }
    }
    
    // Compute border points for each coarse grid
    //
    // These points are accessed by loops below labeled B, C, D, E, F
    // which, for a grid with nx = ny = 8, correspond to the following:
    //
    //     1 2 3 4 5 6 7
    //  7  B C C C C C B
    //  6  B F E E E F B
    //  5  B D 0 0 0 D B
    //  4  B D 0 0 0 D B
    //  3  B D 0 0 0 D B
    //  2  B F E E E F B
    //  1  B C C C C C B
	
    for (int lev=1; lev < q.Ngrid(); ++lev) {
        double bydx = 1. / q.Dx(lev);
        // left and right borders (excluding interface) (B)
        for (int j=1; j<ny; ++j) {
            for (int i=1; i<nx2; ++i) {
                u(lev,i,j) = ( q(lev,X,i,j) + q(lev,X,i,j-1) ) * 0.5 * bydx;
            }
            for (int i=nx/2+nx2+1; i<nx; ++i) {
                u(lev,i,j) = ( q(lev,X,i,j) + q(lev,X,i,j-1) ) * 0.5 * bydx;
            }
        }
        // top and bottom borders (excluding interfaces) (C)
        for (int i=nx2; i<=nx/2+nx2; ++i) {
            for (int j = 1; j<ny2; ++j ) {
                u(lev,i,j) = ( q(lev,X,i,j) + q(lev,X,i,j-1) ) * 0.5 * bydx;
            }
            for (int j = ny/2+ny2+1; j<ny; ++j) {
                u(lev,i,j) = ( q(lev,X,i,j) + q(lev,X,i,j-1) ) * 0.5 * bydx;                
            }
        }
        
        // left and right interfaces, excluding corners (D)
        for ( int j=ny2+1; j<ny/2+ny2; ++j ) {
            u(lev,nx2,j) = ( q(lev,X,nx2,j) + q(lev,X,nx2,j-1) ) * 0.5 * bydx;
            u(lev,nx/2+nx2,j) = ( q(lev,X,nx/2+nx2,j) + q(lev,X,nx/2+nx2,j-1) ) * 0.5 * bydx;
        }
        // top and bottom interfaces, excluding corners (E)
        for ( int i=nx2+1; i<nx/2+nx2; ++i ) {
            int ii, jj;  // fine coords
            u.getGrid().c2f(i,ny2,ii,jj);
            u(lev,i,ny2) = ( q(lev,X,i,ny2-1) * 2./3 + q(lev-1,X,ii,jj) * 1./3 +
                    ( q(lev-1,X,ii-1,jj) + q(lev-1,X,ii+1,jj) ) * 1./6 ) * bydx;
            u.getGrid().c2f(i,ny/2+ny2,ii,jj);
            u(lev,i,ny/2+ny2) = ( q(lev,X,i,ny/2+ny2) * 2./3 + q(lev-1,X,ii,jj-1) * 1./3 +
                    ( q(lev-1,X,ii-1,jj-1) + q(lev-1,X,ii+1,jj-1) ) * 1./6 ) * bydx;
        }
        // corners (F)
        // lower left
        int i=nx2;
        int j=ny2;
        u(lev,i,j) = ( q(lev,X,i,j-1) * 8./15 + q(lev,X,i,j) * 6./15 +
                      q(lev-1,X,1,0) * 2./15 ) * bydx;

        // lower right
        i = nx/2 + nx2;
        u(lev,i,j) = ( q(lev,X,i,j-1) * 8./15 + q(lev,X,i,j) * 6./15 +
                      q(lev-1,X,nx-1,0) * 2./15 ) * bydx;

        // upper left
        i = nx2; j = ny/2 + ny2;
        u(lev,i,j) = ( q(lev,X,i,j) * 8./15 + q(lev,X,i,j-1) * 6./15 +
                      q(lev-1,X,1,ny-1) * 2./15 ) * bydx;

        // upper right
        i = nx/2 + nx2;
        u(lev,i,j) = ( q(lev,X,i,j) * 8./15 + q(lev,X,i,j-1) * 6./15 +
                      q(lev-1,X,nx-1,ny-1) * 2./15 ) * bydx;
    }
}

void FluxToYVelocity(const Flux& q, Scalar& v) {
    assert( q.Nx() == v.Nx() );
    assert( q.Ny() == v.Ny() );
    assert( q.Ngrid() == v.Ngrid() );
    int nx = q.Nx();
    int ny = q.Ny();
    int nx2 = q.NxExt();
    int ny2 = q.NyExt();
    double oneOver2Delta = 1./ (2 * q.Dx());
    
    // Compute interior points (A)
    for (int j=1; j < ny; ++j ){
        for (int i=1; i < nx; ++i ) {
            v(0,i,j) = ( q(0,Y,i-1,j) + q(0,Y,i,j) ) * oneOver2Delta;
        }
    }
    
    // Compute border points for each coarse grid
    //
    // These points are accessed by loops below labeled B, C, D, E, F
    // which, for a grid with nx = ny = 8, correspond to the following:
    //
    // FIXME: diagram below should be transposed
    //     1 2 3 4 5 6 7
    //  7  B C C C C C B
    //  6  B F E E E F B
    //  5  B D 0 0 0 D B
    //  4  B D 0 0 0 D B
    //  3  B D 0 0 0 D B
    //  2  B F E E E F B
    //  1  B C C C C C B
    
    for (int lev=1; lev < q.Ngrid(); ++lev) {
        double bydx = 1. / q.Dx(lev);
        // top and bottom borders (excluding interface) (B)
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny2; ++j) {
                v(lev,i,j) = ( q(lev,Y,i,j) + q(lev,Y,i-1,j) ) * 0.5 * bydx;
            }
            for (int j=ny/2+ny2+1; j<ny; ++j) {
                v(lev,i,j) = ( q(lev,Y,i,j) + q(lev,Y,i-1,j) ) * 0.5 * bydx;
            }
        }
        // left and right borders (excluding interfaces) (C)
        for (int j=ny2; j<=ny/2+ny2; ++j) {
            for (int i = 1; i<nx2; ++i ) {
                v(lev,i,j) = ( q(lev,Y,i,j) + q(lev,Y,i-1,j) ) * 0.5 * bydx;
            }
            for (int i = nx/2+nx2+1; i<nx; ++i) {
                v(lev,i,j) = ( q(lev,Y,i,j) + q(lev,Y,i-1,j) ) * 0.5 * bydx;                
            }
        }
        
        // top and bottom interfaces, excluding corners (D)
        for ( int i=nx2+1; i<nx/2+nx2; ++i ) {
            v(lev,i,ny2) = ( q(lev,Y,i,ny2) + q(lev,Y,i-1,ny2) ) * 0.5 * bydx;
            v(lev,i,ny/2+ny2) = ( q(lev,Y,i,ny/2+ny2) + q(lev,Y,i-1,ny/2+ny2) ) * 0.5 * bydx;
        }
        // left and right interfaces, excluding corners (E)
        for ( int j=ny2+1; j<ny/2+ny2; ++j ) {
            int ii, jj;  // fine coords
            v.getGrid().c2f(nx2,j,ii,jj);
            v(lev,nx2,j) = ( q(lev,Y,nx2-1,j) * 2./3 + q(lev-1,Y,ii,jj) * 1./3 +
                            ( q(lev-1,Y,ii,jj-1) + q(lev-1,Y,ii,jj+1) ) * 1./6 ) * bydx;
            v.getGrid().c2f(nx/2+nx2,j,ii,jj);
            v(lev,nx/2+nx2,j) = ( q(lev,Y,nx/2+nx2,j) * 2./3 + q(lev-1,Y,ii-1,jj) * 1./3 +
                               ( q(lev-1,Y,ii-1,jj-1) + q(lev-1,Y,ii-1,jj+1) ) * 1./6 ) * bydx;
        }
        // corners (F)
        int j=ny2;
        int i=nx2;
        v(lev,i,j) = ( q(lev,Y,i-1,j) * 8./15 + q(lev,Y,i,j) * 6./15 +
                      q(lev-1,Y,0,1) * 2./15 ) * bydx;
        
        j = ny/2 + ny2;
        v(lev,i,j) = ( q(lev,Y,i-1,j) * 8./15 + q(lev,Y,i,j) * 6./15 +
                      q(lev-1,Y,0,ny-1) * 2./15 ) * bydx;
        
        j = ny2; i = nx/2 + nx2;
        v(lev,i,j) = ( q(lev,Y,i,j) * 8./15 + q(lev,Y,i-1,j) * 6./15 +
                      q(lev-1,Y,nx-1,1) * 2./15 ) * bydx;
        
        j = ny/2 + ny2;
        v(lev,i,j) = ( q(lev,Y,i,j) * 8./15 + q(lev,Y,i-1,j) * 6./15 +
                      q(lev-1,Y,nx-1,ny-1) * 2./15 ) * bydx;
    }
}

// Convert u-velocities at vertices to x-fluxes through edges.
// Does not touch the y-component of the Flux q passed in.
void XVelocityToFlux(const Scalar& u, Flux& q) {
    assert( u.Nx() == q.Nx() );
    assert( u.Ny() == q.Ny() );
    assert( u.Ngrid() == q.Ngrid() );
    int nx = u.Nx();
    int ny = u.Ny();
    int nx2 = u.NxExt();
    int ny2 = u.NyExt();
    const Grid& g = q.getGrid();

    // Finest grid, fluxes: (nx = ny = 8)
    //     0 1 2 3 4 5 6 7 8
    //  7  A C C C C C C C A
    //  6  A D D D D D D D A
    //  5  A D D D D D D D A
    //  4  A D D D D D D D A
    //  3  A D D D D D D D A
    //  2  A D D D D D D D A
    //  1  A D D D D D D D A
    //  0  A C C C C C C C A
    //
    // Coarser grids, fluxes: (nx = ny = 8)
    //     0 1 2 3 4 5 6 7 8
    //  7  A C C C C C C C A
    //  6  A B B B B B B B A
    //  5  A D D G G G D D A
    //  4  A D D G G G D D A
    //  3  A D D G G G D D A
    //  2  A D D G G G D D A
    //  1  A B B B B B B B A
    //  0  A C C C C C C C A
        
    for (int lev=0; lev < u.Ngrid(); ++lev) {
        double dx = g.Dx(lev);
        // Interior points
        // If on fine grid, compute all interior points
        if (lev == 0) {
            // compute interior points on finest grid, minus top and bottom rows (D)
            for (int i=1; i<nx; ++i) {
                for (int j=1; j<ny-1; ++j) {
                    q(0,X,i,j) = ( u(0,i,j) + u(0,i,j+1) ) * 0.5 * dx;
                }
            }            
        }
        else {  // not the finest grid
            for (int i=1; i<nx; ++i) {
                // top and bottom portions of coarse grid, excluding outer interface (B)
                for (int j=1; j<ny2; ++j) {
                    q(lev,X,i,j) = ( u(lev,i,j) + u(lev,i,j+1) ) * 0.5 * dx;
                }
                for (int j=ny/2+ny2; j<ny-1; ++j) {
                    q(lev,X,i,j) = ( u(lev,i,j) + u(lev,i,j+1) ) * 0.5 * dx;
                }
            }
            // left and right portions of coarse grid (D)
            for (int j=ny2; j<ny/2+ny2; ++j) {
                for (int i=1; i<=nx2; ++i) {
                    q(lev,X,i,j) = ( u(lev,i,j) + u(lev,i,j+1) ) * 0.5 * dx;
                }
                for (int i=nx/2+nx2; i<nx; ++i) {
                    q(lev,X,i,j) = ( u(lev,i,j) + u(lev,i,j+1) ) * 0.5 * dx;
                }
            }
            // get interior portion of coarse grid from fine grid (G)
            for (int i=nx2+1; i<nx/2+nx2; ++i) {
                for (int j=ny2; j<ny/2+ny2; ++j) {
                    int ii,jj; // fine gridpoints
                    g.c2f(i,j,ii,jj);
                    q(lev,X,i,j) = q(lev-1,X,ii,jj) + q(lev-1,X,ii,jj+1);
                }
            }            
        }
        // Boundary points
        // left and right boundaries of coarsest grid are zero (A)
        if (lev == g.Ngrid()-1) {
            for (int j=0; j<ny; ++j) {
                q(lev,X,0,j) = 0;
                q(lev,X,nx,j) = 0;
            }
        }
        // left and right boundaries of finer grids take values from coarser grid (A)
        else {
            for (int j=0; j<ny-1; j+=2) {
                int ii,jj;  // coarse indices
                g.f2c(0,j,ii,jj);
                q(lev,X,0,j) = ( 0.75 * u(lev+1,nx2,jj) + 0.25 * u(lev+1,nx2,jj+1) ) * dx;
                q(lev,X,nx,j) = ( 0.75 * u(lev+1,nx/2+nx2,jj) + 0.25 * u(lev+1,nx/2+nx2,jj+1) ) * dx;
                q(lev,X,0,j+1) = ( 0.25 * u(lev+1,nx2,jj) + 0.75 * u(lev+1,nx2,jj+1) ) * dx;
                q(lev,X,nx,j+1) = ( 0.25 * u(lev+1,nx/2+nx2,jj) + 0.75 * u(lev+1,nx/2+nx2,jj+1) ) * dx;
            }
        }
        for (int i=1; i<nx; ++i) {
            // outer interface, get values from coarser grid
            // (or zero, for coarsest) (C)
            if (lev == u.Ngrid()-1) {
                // on coarsest grid: zero bcs
                q(lev,X,i,0) =  u(lev,i,1) * 0.5 * dx;
                q(lev,X,i,ny-1) =  u(lev,i,ny-1) * 0.5 * dx;
            }
            else {
                // on intermediate grid: get bcs from coarser grid
                for (int i=2; i<nx; i += 2) {
                    // points that correspond to coarse points
                    int ii,jj; // coarse points
                    g.f2c(i,0,ii,jj);
                    q(lev,X,i,0) = ( u(lev,i,1) + u(lev+1,ii,ny2) ) * 0.5 * dx;
                    q(lev,X,i,ny-1) = ( u(lev,i,ny-1) + u(lev+1,ii,ny/2+ny2) ) * 0.5 * dx;
                }
                for (int i=1; i<nx; i += 2) {
                    // points that do not correspond to coarse points
                    int ii,jj; // coarse points
                    g.f2c(i,0,ii,jj);
                    q(lev,X,i,0) = ( 0.5 * u(lev,i,1) +
                                  0.25 * u(lev+1,ii,ny2) + 0.25 * u(lev+1,ii+1,ny2) ) * dx;
                    q(lev,X,i,ny-1) = ( 0.5 * u(lev,i,ny-1) +
                                     0.25* u(lev+1,ii,ny/2+ny2) + 0.25 * u(lev+1,ii+1,ny/2+ny2) ) * dx;
                }
            }
        }
    }
}

// Convert v-velocities at vertices to y-fluxes through edges.
// Does not touch the x-component of the Flux q passed in.
void YVelocityToFlux(const Scalar& v, Flux& q) {
    assert( v.Nx() == q.Nx() );
    assert( v.Ny() == q.Ny() );
    assert( v.Ngrid() == q.Ngrid() );
    int nx = v.Nx();
    int ny = v.Ny();
    int nx2 = v.NxExt();
    int ny2 = v.NyExt();
    const Grid& g = v.getGrid();
    
    // Finest grid, fluxes: (nx = ny = 8)
    //     0 1 2 3 4 5 6 7
    //  8  A A A A A A A A
    //  7  C D D D D D D C
    //  6  C D D D D D D C
    //  5  C D D D D D D C
    //  4  C D D D D D D C
    //  3  C D D D D D D C
    //  2  C D D D D D D C
    //  1  C D D D D D D C
    //  0  A A A A A A A A
    //
    // Coarser grids, fluxes: (nx = ny = 8)
    //     0 1 2 3 4 5 6 7
    //  8  A A A A A A A A
    //  7  C B D D D D B C
    //  6  C B D D D D B C
    //  5  C B G G G G B C
    //  4  C B G G G G B C
    //  3  C B G G G G B C
    //  2  C B D D D D B C
    //  1  C B D D D D B C
    //  0  A A A A A A A A

    for (int lev=0; lev < g.Ngrid(); ++lev) {
        double dx = g.Dx(lev);
        // Interior points
        // If on fine grid, compute all interior points
        if (lev == 0) {
            // compute interior points on finest grid, minus top and bottom rows (D)
            for (int j=1; j<ny; ++j) {
                for (int i=1; i<nx-1; ++i) {
                    q(0,Y,i,j) = ( v(0,i,j) + v(0,i+1,j) ) * 0.5 * dx;
                }
            }            
        }
        else {  // not the finest grid
            for (int j=1; j<ny; ++j) {
                // left and right portions of coarse grid, excluding outer interface (B)
                for (int i=1; i<nx2; ++i) {
                    q(lev,Y,i,j) = ( v(lev,i,j) + v(lev,i+1,j) ) * 0.5 * dx;
                }
                for (int i=nx/2+nx2; i<nx-1; ++i) {
                    q(lev,Y,i,j) = ( v(lev,i,j) + v(lev,i+1,j) ) * 0.5 * dx;
                }
            }
            // top and bottom portions of coarse grid (D)
            for (int i=nx2; i<nx/2+nx2; ++i) {
                for (int j=1; j<=ny2; ++j) {
                    q(lev,Y,i,j) = ( v(lev,i,j) + v(lev,i+1,j) ) * 0.5 * dx;
                }
                for (int j=ny/2+ny2; j<ny; ++j) {
                    q(lev,Y,i,j) = ( v(lev,i,j) + v(lev,i+1,j) ) * 0.5 * dx;
                }
            }
            // get interior portion of coarse grid from fine grid (G)
            for (int j=ny2+1; j<ny/2+ny2; ++j) {
                for (int i=nx2; i<nx/2+nx2; ++i) {
                    int ii,jj; // fine gridpoints
                    g.c2f(i,j,ii,jj);
                    q(lev,Y,i,j) = q(lev-1,Y,ii,jj) + q(lev-1,Y,ii+1,jj);
                }
            }            
        }
        // Boundary points
        // top and bottom boundaries of coarsest grid are zero (A)
        if (lev == g.Ngrid()-1) {
            for (int i=0; i<nx; ++i) {
                q(lev,Y,i,0) = 0;
                q(lev,Y,i,ny) = 0;
            }
        }
        // top and bottom boundaries of finer grids take values from coarser grid (A)
        else {
            for (int i=0; i<nx-1; i+=2) {
                int ii,jj;  // coarse indices
                g.f2c(i,0,ii,jj);
                q(lev,Y,i,0) = ( 0.75 * v(lev+1,ii,ny2) + 0.25 * v(lev+1,ii+1,ny2) ) * dx;
                q(lev,Y,i,ny) = ( 0.75 * v(lev+1,ii,ny/2+ny2) + 0.25 * v(lev+1,ii+1,ny/2+ny2) ) * dx;
                q(lev,Y,i+1,0) = ( 0.25 * v(lev+1,ii,ny2) + 0.75 * v(lev+1,ii+1,ny2) ) * dx;
                q(lev,Y,i+1,ny) = ( 0.25 * v(lev+1,ii,ny/2+ny2) + 0.75 * v(lev+1,ii+1,ny/2+ny2) ) * dx;
            }
        }
        for (int j=1; j<ny; ++j) {
            // outer interface, get values from coarser grid
            // (or zero, for coarsest) (C)
            if (lev == g.Ngrid()-1) {
                // on coarsest grid: zero bcs
                q(lev,Y,0,j) =  v(lev,1,j) * 0.5 * dx;
                q(lev,Y,nx-1,j) =  v(lev,nx-1,j) * 0.5 * dx;
            }
            else {
                // on intermediate grid: get bcs from coarser grid
                for (int j=2; j<ny; j += 2) {
                    // points that correspond to coarse points
                    int ii,jj; // coarse points
                    g.f2c(0,j,ii,jj);
                    q(lev,Y,0,j) = ( v(lev,1,j) + v(lev+1,nx2,jj) ) * 0.5 * dx;
                    q(lev,Y,nx-1,j) = ( v(lev,nx-1,j) + v(lev+1,nx/2+nx2,jj) ) * 0.5 * dx;
                }
                for (int j=1; j<ny; j += 2) {
                    // points that do not correspond to coarse points
                    int ii,jj; // coarse points
                    g.f2c(0,j,ii,jj);
                    q(lev,Y,0,j) = ( 0.5 * v(lev,1,j) +
                                    0.25 * v(lev+1,nx2,jj) + 0.25 * v(lev+1,nx2,jj+1) ) * dx;
                    q(lev,Y,nx-1,j) = ( 0.5 * v(lev,nx-1,j) +
                                       0.25* v(lev+1,nx/2+nx2,jj) + 0.25 * v(lev+1,nx/2+nx2,jj+1) ) * dx;
                }
            }
        }
    }
}

// Convert u- and v-velocities at vertices to fluxes through edges
void VelocityToFlux(const Scalar& u, const Scalar& v, Flux& q) {
    XVelocityToFlux( u, q );
    YVelocityToFlux( v, q );
}

// Convert fluxes through edges to u- and v-velocities at vertices
void FluxToVelocity(const Flux& q, Scalar& u, Scalar& v) {
    FluxToXVelocity( q, u );
    FluxToYVelocity( q, v );
}

} // namespace ibpm
