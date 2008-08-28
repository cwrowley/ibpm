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
// $Revision: $
// $LastChangedDate: $
// $LastChangedBy: zma $
// $HeadURL: $

#include "Grid.h"
#include "Scalar.h"
#include "Flux.h"
#include "BoundaryVector.h"
#include "VectorOperations.h"
#include <fftw3.h>

// Return the curl of Flux q, as a Scalar object
Scalar Curl(const Flux& q) { 
    const Grid& grid = q.getGrid();
    Scalar f(grid);
    int nx = grid.getNx();
    int ny = grid.getNy();
    double byDeltaSquared = 1. / ( grid.getDx() * grid.getDx() );

	// Compute curl at interior nodes
	for (int i = 1; i < nx; ++i) {
		for (int j=1; j<ny; ++j) {
			f(i,j) = ( q(Y,i,j) - q(Y,i-1,j) - q(X,i,j) + q(X,i,j-1) ) *
			    byDeltaSquared;
		}
	}
	
	// Return zero values at the boundary nodes.
	for (int j = 0; j < ny+1; ++j) {
		f(0,j) = 0;
		f(nx,j) = 0;
	}
	for (int i = 0; i < nx+1; ++i) {
		f(i,0) = 0;
		f(i,ny) = 0;
	}		
	return f;
}

// Return the curl of Scalar f, as a Flux object.
Flux Curl(const Scalar& f) {
    const Grid& grid = f.getGrid();
    Flux q(grid);
    
    int nx = grid.getNx();
    int ny = grid.getNy();
    
    // X direction
    for (int i=0; i<=nx; ++i) {
        for (int j=0; j<ny; ++j) {
            q(X,i,j) = f(i,j+1) - f(i,j);
        }
    }
    
    // Y direction
    for (int i=0; i<nx; ++i) {
        for (int j=0; j<=ny; ++j) {
            q(Y,i,j) = f(i,j) - f(i+1,j);
        }
    }
    return q;
}

// Inner product of two Scalars. 
double InnerProduct (const Scalar& f, const Scalar& g){
    const Grid& grid = f.getGrid();
	int nx = grid.getNx();
	int ny = grid.getNy();
	assert(nx == g.getGrid().getNx());
	assert(ny == g.getGrid().getNy());
	
	double ip = 0;
	for (int i = 1; i < nx; ++i) {
		for ( int j = 1; j < ny; ++j) {
			ip += f(i, j) * g(i, j);
		}
		ip += f(i, 0) * g(i, 0) / 2;
		ip += f(i, ny) * g(i, ny) / 2;			
	}
	for ( int j = 1; j < ny; ++j) {
		ip += f(0, j) * g(0, j)/2;
		ip += f(nx, j) * g(nx, j)/2;
	}
	ip += ( f(0,0) * g(0,0) + f(nx,0) * g(nx,0) 
		  + f(0, ny) * g(0, ny) + f(nx, ny) * g(nx,ny) )/4;

    double dx = grid.getDx();
	return ip * dx * dx;
}

// Inner product of Flux p and Flux q.
// Note that this is not multiplied by dx * dx, since the Fluxes are already
// multiplied by these (i.e., inner product is really over *velocities*).
double InnerProduct (const Flux& p, const Flux& q){
    const Grid& grid = p.getGrid();
	int nx = grid.getNx();
	int ny = grid.getNy();
	assert(nx == q.getGrid().getNx());
	assert(ny == q.getGrid().getNy());

    double ip = 0;
    // Add x-fluxes
    for (int j=0; j<ny; ++j) {
        for (int i=1; i<nx; ++i){
            ip += p(X,i,j) * q(X,i,j);
        }
		ip += p(X,0,j) * q(X,0,j) / 2;
		ip += p(X,nx,j) * q(X,nx,j) / 2;		
    }
    
    // Add y-fluxes
    for (int i=0; i<nx; ++i) {
        for (int j=1; j<ny; ++j) {
            ip += p(Y,i,j) * q(Y,i,j);
        }
		ip += p(Y,i,0) * q(Y,i,0) / 2;
		ip += p(Y,i,ny) * q(Y,i,ny) / 2;	
    }
    return ip;
}

// Sum of X-component of Flux q
double XSum( const Flux& q ) {
   const Grid& grid = q.getGrid();
   int nx = grid.getNx();
   int ny = grid.getNy();
   double qsumx = 0;
   // Sum x-fluxes 
   for (int j=0; j<ny; ++j) {
       for (int i=0; i<=nx; ++i){  		
	   qsumx += q(X,i,j);
       }
   }
   return qsumx;
}

// Sum of Y-component of Flux q
double YSum( const Flux& q) {
   const Grid& grid = q.getGrid();
   int nx = grid.getNx();
   int ny = grid.getNy();
   double qsumy = 0;
   // Sum y-fluxes
   for (int i=0; i<nx; ++i) {
       for (int j=0; j<=ny; ++j) {
           qsumy += q(Y,i,j);
       } 
   } 
   return qsumy;
}	

// Routine for computing X & Y forces
void computeNetForce( BoundaryVector& f, double& xforce, double& yforce) {
   xforce = 0;
   yforce = 0;
   for( int i=0; i<f.getNumPoints(); i++ ) {
      xforce += f(X,i);
      yforce += f(Y,i);
   }  
}

// sine transform of a Scalar object using fft (type DST-I).
// (fftw library is used(kind: RODFT00); Only interior nodes are considered.)
Scalar SinTransform(const Scalar& f) {
	Scalar f_fft(f.getGrid());
	int nx = f.getGrid().getNx();
	int ny = f.getGrid().getNy();
	
	double *in;
	in = (double*) fftw_malloc(sizeof(double) * (nx-1) * (ny-1));
	
	fftw_plan p;
	p = fftw_plan_r2r_2d(nx-1, ny-1, in, in,
						 FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE); 	
	
	for (int i = 0; i < nx-1; ++i){
		for (int j = 0; j<ny-1; ++j){
			*(in+i*(ny-1)+j) = f(i+1,j+1);
		}
	}
	
	fftw_execute(p); 

	for (int i = 0; i < nx-1; ++i){
		for (int j = 0; j<ny-1; ++j){
			f_fft(i+1,j+1) = *(in+i*(ny-1)+j);	     
		}	  
	}
	
	fftw_destroy_plan(p);
	fftw_free(in);
	
	//boundary elements are set zero:
	for (int i=0; i<nx+1; ++i) {
		f_fft(i, 0) = 0;
		f_fft(i,ny) = 0;
	}
	for (int j=0; j<ny+1; ++j) {       
		f_fft(0,j)=0;
		f_fft(nx,j)=0 ;  
	}           
		
	return f_fft;
 }

// Return cross product of a Flux q and a Scalar f, as a Flux object.
//   q x f = ( f v, -f u )
//
// Step 1: Convert flux to velocities (u,v) at vertices
// Step 2: Compute f * v and -f * u
// Step 3: Convert the velocities to fluxes through edges

Flux CrossProduct(const Flux& q, const Scalar& f){
    const Grid& grid = f.getGrid();
	int nx = grid.getNx();
	int ny = grid.getNy();
	assert(nx == q.getGrid().getNx());
	assert(ny == q.getGrid().getNy());                                                                                                                                                                 

    Scalar u(grid);
    Scalar v(grid);
    
    FluxToXVelocity( q, u );
    u *= f;
    u *= -1;
    
    FluxToYVelocity( q, v );
    v *= f;
    
    Flux cross(grid);
    VelocityToFlux( v, u, cross );      // cross = ( f v, -f u )

    return cross;
}

// Return cross product of two Flux objects q1, q2, as a Scalar object.
//   q1 x q2 = u1 v2 - u2 v1
//
// Step 1: Convert the fluxes to velocities at nodes
// Step 2: Compute u1 * v2 - u2 * v1

Scalar CrossProduct(const Flux& q1, const Flux& q2){
	assert(q1.getGrid().getNx() == q2.getGrid().getNx());
	assert(q1.getGrid().getNy() == q2.getGrid().getNy());
	
    const Grid& grid = q1.getGrid();
    Scalar u(grid);
    Scalar v(grid);
    
    FluxToXVelocity( q1, u );
    FluxToYVelocity( q2, v );

    Scalar f = u * v;  // f is now u1 * v2
    
    FluxToXVelocity( q2, u );
    FluxToYVelocity( q1, v );
    
    f -= u * v;  // f is now u1 * v2 - u2 * v1
	
	return f;
};

void FluxToXVelocity(const Flux& q, Scalar& u) {
    const Grid& grid = u.getGrid();
	int nx = grid.getNx();
	int ny = grid.getNy();
    double oneOver2Delta = 1./ (2 * grid.getDx());

    for (int i=0; i <= nx; ++i ){
        for (int j=1; j < ny; ++j ) {
            u(i,j) = ( q(X,i,j) + q(X,i,j-1) ) * oneOver2Delta;
        }
        u(i,0) = q(X,i,0) * 2 * oneOver2Delta;
        u(i,ny) = q(X,i,ny-1) * 2 * oneOver2Delta;
	}
}

void FluxToYVelocity(const Flux& q, Scalar& v) {
    const Grid& grid = v.getGrid();
	int nx = grid.getNx();
	int ny = grid.getNy();
    double oneOver2Delta = 1./ (2 * grid.getDx());

    for (int j=0; j <= ny; ++j ){
        for (int i=1; i < nx; ++i ) {
            v(i,j) = ( q(Y,i,j) + q(Y,i-1,j) ) * oneOver2Delta;
        }
        v(0,j) = q(Y,0,j) * 2 * oneOver2Delta;
        v(nx,j) = q(Y,nx-1,j) * 2 * oneOver2Delta;
	}
}

// Convert u-velocities at vertices to x-fluxes through edges.
// Does not touch the y-component of the Flux q passed in.
void XVelocityToFlux(const Scalar& u, Flux& q) {
    const Grid& grid = u.getGrid();
    int nx = grid.getNx();
    int ny = grid.getNy();
    double deltaOver2 = grid.getDx() / 2.;
    
    for (int i=0; i <= nx; ++i ) {
        for (int j=0; j < ny; ++j ){
            q(X,i,j) = ( u(i,j) + u(i,j+1) ) * deltaOver2;
        }
    }
}

// Convert v-velocities at vertices to y-fluxes through edges.
// Does not touch the x-component of the Flux q passed in.
void YVelocityToFlux(const Scalar& v, Flux& q) {
    const Grid& grid = v.getGrid();
    int nx = grid.getNx();
    int ny = grid.getNy();
    double deltaOver2 = grid.getDx() / 2.;
    
    for (int i=0; i < nx; ++i ) {
        for (int j=0; j <= ny; ++j ){
            q(Y,i,j) = ( v(i,j) + v(i+1,j) ) * deltaOver2;
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
