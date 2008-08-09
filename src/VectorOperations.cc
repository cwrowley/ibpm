// VectorOperations.cc
//
// Description:
// Definition of VectorOperations functions.
//
// Author(s):
// Clancy Rowley, Zhanhua Ma
//
// Date: 15 Jul 2008
//
// $Revision: $
// $LastChangedDate: $
// $LastChangedBy: zma $
// $HeadURL:// $Header$

#include "Grid.h"
#include "Scalar.h"
#include "Flux.h"
#include "BoundaryVector.h"

// Return the curl of Flux q, as a Scalar object
Scalar curl(const Flux& q) { 
	Scalar f(q.getGrid());
	int nx = f.getGrid().getNx();
	int ny = f.getGrid().getNy();
	
	// Compute curl at interior nodes
	for (int i = 1; i < nx; ++i) {
		for (int j=1; j<ny; ++j) {
			f(i,j) = q(Y,i,j) - q(Y,i-1,j) - q(X,i,j) + q(X,i,j-1);
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
	double k = f(1,1);		 
	return f;
}

// Return the curl of Scalar f, as a Flux object.
Flux curl(const Scalar& f) {
	Flux q(f.getGrid());
	int nx = q.getGrid().getNx();
	int ny = q.getGrid().getNy();
	
	// X direction.
	for (int i = 0; i < nx+1; ++i) {
		for ( int j = 0; j < ny; ++j) {
			q(X,i,j) = f(i,j+1) - f(i,j);
		}					
	}		
		
	// Y direction.
	for (int i=0; i<nx; ++i) {
		for (int j=0; j<ny+1; ++j) {
			q(Y,i,j) = f(i,j) - f(i+1,j);
		}									
	}			
	return q;
}

//dot product of two Scalar products. 
double InnerProduct (const Scalar& f, const Scalar& g){
	int nx = f.getGrid().getNx();
	int ny = f.getGrid().getNy();
	assert(nx == g.getGrid().getNx());
	assert(ny == g.getGrid().getNy());
	
	double dp = 0;
	for (int i = 1; i < nx; ++i) {
		for ( int j = 1; j < ny; ++j) {
			dp += f(i, j) * g(i, j);
		}
		dp += f(i, 0) * g(i, 0) / 2;
		dp += f(i, ny) * g(i, ny) / 2;			
	}
	for ( int j = 1; j < ny; ++j) {
		dp += f(0, j) * g(0, j)/2;
		dp += f(nx, j) * g(nx, j)/2;
	}
	dp += (f(0,0) * g(0,0) + f(nx,0) * g(nx,0) 
		+ f(0, ny) * g(0, ny) + f(nx, ny) * g(nx,ny))/4;
	dp *= pow(f.getGrid().getDx(), 2);					
	return dp;
}

// dot product of Flux p and Flux q. 
double InnerProduct (const Flux& p, const Flux& q){
	int nx = p.getGrid().getNx();
	int ny = p.getGrid().getNy();
	assert(nx == q.getGrid().getNx());
	assert(ny == q.getGrid().getNy());

    double dp = 0;
    // Add x-fluxes
    for (int j=0; j<ny; ++j) {
        for (int i=1; i<nx; ++i){
            dp += p(X,i,j) * q(X,i,j);
        }
		dp += p(X,0,j) * q(X,0,j) / 2;
		dp += p(X,nx,j) * q(X,nx,j) / 2;		
    }
    
    // Add y-fluxes
    for (int i=0; i<nx; ++i) {
        for (int j=1; j<ny; ++j) {
            dp += p(Y,i,j) * q(Y,i,j);
        }
		dp += p(Y,i,0) * q(Y,i,0) / 2;
		dp += p(Y,i,ny) * q(Y,i,ny) / 2;	
    }
    double dx = p.getGrid().getDx();
    return dp * dx * dx;
}

// Return the inner product of BoundaryVectors x and y.
// NOTE: The implementation below uses only the public interface.
//       A version using private data is contained in BoundaryVector.h
// double InnerProduct(
//     const BoundaryVector& x,
//     const BoundaryVector& y
//     ) {
//     BoundaryVector::index ind;
//     double ip = 0;
// 
//     for (ind = x.begin(); ind != x.end(); ++ind) {
//         ip += x(ind) * y(ind);
//     }
//     return ip;
// }