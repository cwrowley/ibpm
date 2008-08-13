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
#include <fftw3.h>
//#include <blitz/array.h>
//
//BZ_USING_NAMESPACE(blitz)

// Declarations of functions used only in this file
Scalar fluxXAverage(const Flux& q);
Scalar fluxYAverage(const Flux& q);

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

// sine transform of a Scalar object using fft (type DST-I).
// (fftw library is used(kind: RODFT00); Only interior nodes are considered.)
Scalar sinTransform(const Scalar& f) {
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

// inverse sine transform of a Scalar object using fft.(fftw library is used.)
//Scalar sinTransformInv(const Scalar& f) {
//	// WORK HERE
//	//return *this;
//}

// Return 'cross product' of a Flux and a Scalar object, as a Flux object.
Flux crossproduct(const Flux& q, const Scalar& f){
	Flux p(q.getGrid());
	int nx = q.getGrid().getNx();
	int ny = q.getGrid().getNy();
	assert(nx == f.getGrid().getNx());
	assert(ny == f.getGrid().getNy());                                                                                                                                                                 
	
	Scalar tempY = fluxYAverage(q);//Y direction average of flux.

	// X direction flux (the p.X at i = 0 and nx are set zero)
	tempY *= f;
	for (int j= 0; j<ny; ++j)  {
		for (int i = 1; i < nx; ++i){
			p(X,i,j) = 0.5 * (tempY(i,j) + tempY(i,j+1));
		}
		p(X,0,j) = 0;
		p(X,nx,j) = 0; 
	}	
	
	Scalar tempX = fluxXAverage(q);// X direction average of flux.
	
	// Y direction flux (the p.Y at j = 0 and ny are set zero)
	tempX *= f;
	for  (int i = 0; i < nx; ++i) {
		for (int j= 1; j<ny; ++j){
			p(Y,i,j) = -0.5 * (tempX(i,j) + tempX(i+1,j));
		}
		p(Y,i,0) = 0;
		p(Y,i,ny) = 0;
	}
	
	return p;  
}

// Return q cross p, 'cross product'  of two Flux objects q, p, as a Scalar object.
Scalar crossproduct(const Flux& q, const Flux& p){
	assert(q.getGrid().getNx() == p.getGrid().getNx());
	assert(q.getGrid().getNy() == p.getGrid().getNy());
	
	Scalar f = fluxXAverage(q);
	f *= fluxYAverage(p);
	Scalar tmp = fluxYAverage(q);
	tmp *= fluxXAverage(p);
	f -= tmp;
	
	return f;
};

/// Return the average of a flux in x direction, taken at nodes, as a Scalar object. 
/// BCs: By supposing zero 'ghost flux.x' outside of the domain.  
Scalar fluxXAverage(const Flux& q){
	Scalar f(q.getGrid());
	int nx = q.getGrid().getNx();
	int ny = q.getGrid().getNy();
	
	for (int i = 0; i < nx+1; ++i) {
		for (int j= 1; j<ny; ++j) {
			f(i, j) = 0.5 * (q(X,i,j) + q(X,i,j-1)); // at 'inner nodes'
		}
		// at 'boundary nodes' j = 0 and ny:
		f(i, 0) = 0.5 * q(X, i, 0); 
		f(i, ny) = 0.5 * q(X, i, ny-1);
	}
	return f;		
}

/// Return the average of a flux in y direction, taken at nodes, as a Scalar object. 
/// BCs: By supposing zero 'ghost flux.y' outside of the domain.  
Scalar fluxYAverage(const Flux& q){
	Scalar f(q.getGrid());
	int nx = q.getGrid().getNx();
	int ny = q.getGrid().getNy();
	
	for (int j = 0; j < ny+1; ++j) {
		for (int i= 1; i<nx; ++i) {
			f(i, j) = 0.5 * (q(Y,i,j) + q(Y,i-1,j)); // at 'inner nodes'
		}
		// at 'boundary nodes' i = 0 and nx:
		f(0, j) = 0.5 * q(Y, 0, j); 
		f(nx, j) = 0.5 * q(Y, nx-1, j);
	}
	return f;		
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
