#include <string>
#include <fstream>
#include <fftw3.h>
#include "Field.h"
#include "Utils.h"
#include "QGEqn.h"

using namespace std;
using std::string;

QGEqn::QGEqn( Field* _phi, double _beta, double _f, double _nu, double _dt ) {
	phi	= _phi;
	beta	= _beta;
	f	= _f;
	nu	= _nu;
	dt	= _dt;

	F	= NULL;
}

QGEqn::~QGEqn() {
}

void QGEqn::Laplacian( Field* psi, Field* omega ) {
	int	k[2];

	psi->Forward();
	for( int i = 0; i < psi->nf; i++ ) {
		psi->IndexToMode( i, k );
		omega->kVals[i][0] = -(k[0]*k[0] + k[1]*k[1])*psi->kVals[i][0];
		omega->kVals[i][1] = -(k[0]*k[0] + k[1]*k[1])*psi->kVals[i][1];
	}
	omega->Backward();
}

void QGEqn::Derivative( Field* psi, Field* dPsi, int dim ) {
	int	k[2];

	psi->Forward();
	for( int i = 0; i < psi->nf; i++ ) {
		psi->IndexToMode( i, k );
		dPsi->kVals[i][0] = -k[dim]*psi->kVals[i][1];
		dPsi->kVals[i][1] = +k[dim]*psi->kVals[i][0];
	}
	dPsi->Backward();
}

void QGEqn::MapToBuffer( Field* psi, Field* psi_b ) {
	int	start1, start2;
	int	diff = psi_b->nx - psi->nx;

	psi->Forward();
	for( int i = 0; i < psi->nx/2+1; i++ ) {
		start1 = i*(psi->ny/2+1);
		start2 = i*(psi_b->ny/2+1);
		for( int j = 0; j < psi->ny/2+1; j++ ) {
			psi_b->kVals[start2+j][0] = psi->kVals[start1+j][0];
			psi_b->kVals[start2+j][1] = psi->kVals[start1+j][1];
		}
	}
	for( int i = psi_b->nx-1; i > psi_b->nx - psi->nx/2; i-- ) {
		start1 = (i-diff)*(psi->ny/2+1);
		start2 = i*(psi_b->ny/2+1);
		for( int j = 0; j < psi->ny/2+1; j++ ) {
			psi_b->kVals[start2+j][0] = psi->kVals[start1+j][0];
			psi_b->kVals[start2+j][1] = psi->kVals[start1+j][1];
		}
	}

	psi_b->Backward();
}

void QGEqn::MapFromBuffer( Field* psi_b, Field* psi ) {
	int 	start1, start2;
	int	diff = psi_b->nx - psi->nx;

	psi_b->Forward();
	for( int i = 0; i < psi->nx/2+1; i++ ) {
		start1 = i*(psi->ny/2+1);
		start2 = i*(psi_b->ny/2+1);
		for( int j = 0; j < psi->ny/2+1; j++ ) {
			psi->kVals[start1+j][0] = psi_b->kVals[start2+j][0];
			psi->kVals[start1+j][1] = psi_b->kVals[start2+j][1];
		}
	}
	for( int i = psi_b->nx-1; i > psi_b->nx - psi->nx/2; i-- ) {
		start1 = (i-diff)*(psi->ny/2+1);
		start2 = i*(psi_b->ny/2+1);
		for( int j = 0; j < psi->ny/2+1; j++ ) {
			psi->kVals[start1+j][0] = psi_b->kVals[start2+j][0];
			psi->kVals[start1+j][1] = psi_b->kVals[start2+j][1];
		}
	}
	psi->Backward();
}

void QGEqn::Convolve( Field* psi, Field* jacobian ) {
	int	nx2		= psi->nx + psi->nx/2;
	int	ny2		= psi->ny + psi->ny/2;
	Field*	psiBuf		= new Field( "psiBuf", nx2, ny2 );
	Field*	omega		= new Field( "omega", nx2, ny2 );
	Field*	dPsiDx		= new Field( "dPsiDx", nx2, ny2 );
	Field*	dPsiDy		= new Field( "dPsiDy", nx2, ny2 );
	Field*	dOmegaDx	= new Field( "dOmegaDx", nx2, ny2 );
	Field*	dOmegaDy	= new Field( "dOmegaDy", nx2, ny2 );
	Field*	jacBuf		= new Field( "jacBuf", nx2, ny2 );

	MapToBuffer( psi, psiBuf );

	Laplacian( psiBuf, omega );
	Derivative( psiBuf, dPsiDx, 0 );
	Derivative( psiBuf, dPsiDy, 1 );
	Derivative( omega, dOmegaDx, 0 );
	Derivative( omega, dOmegaDy, 1 );

	for( int i = 0; i < jacBuf->nr; i++ ) {
		jacBuf->xVals[i] = dPsiDx->xVals[i]*dOmegaDy->xVals[i] - dPsiDy->xVals[i]*dOmegaDx->xVals[i];
	}
	MapFromBuffer( jacBuf, jacobian );

	delete psiBuf;
	delete omega;
	delete dPsiDx;
	delete dPsiDy;
	delete dOmegaDx;
	delete dOmegaDy;
	delete jacBuf;
}

void QGEqn::Solve( Field* phiPrev ) {
	Field*	phiTemp	= NULL;
	if( phiPrev ) {
		phiTemp = new Field( "phi-temp", phi->nx, phi->ny );
		phiTemp->Copy( phi );
		SecondOrder( phiPrev );
		phiPrev->Copy( phiTemp );
		delete phiTemp;
	}
	else {
		FirstOrder();
	}
}

void QGEqn::FirstOrder() {
	int	k[2], k2;
	double	rhs[2], a11, a12, a21, a22, detInv;
	Field*	jacobian	= new Field( "jacobian", phi->nx, phi->ny );

	Convolve( phi, jacobian );

	/* update phi */
	for( int i = 1; i < phi->nf; i++ ) {
		phi->IndexToMode( i, k );
		k2 = -(k[0]*k[0] + k[1]*k[1]);
		rhs[0] = (k2 - f)*phi->kVals[i][0] - dt*jacobian->kVals[i][0];
		rhs[1] = (k2 - f)*phi->kVals[i][1] - dt*jacobian->kVals[i][1];

		a11 = k2 - f + dt*nu*k2*k2*k2; //biharmonic viscosity of opposite sign to regular viscosity
		a12 = -dt*beta*k[0];
		a21 = +dt*beta*k[0];
		a22 = a11;
		detInv = 1.0/(a11*a22 - a12*a21);
		
		phi->kVals[i][0] = detInv*(a22*rhs[0] - a12*rhs[1]);
		phi->kVals[i][1] = detInv*(a11*rhs[1] - a21*rhs[0]);
	}
	phi->Backward();

	delete jacobian;
}

void QGEqn::SecondOrder( Field* phiPrev ) {
	int	k[2], k2;
	double	rhs[2], a11, a12, a21, a22, detInv;
	Field*	jac1	= new Field( "jac-1", phi->nx, phi->ny );
	Field*	jac2	= new Field( "jac-2", phi->nx, phi->ny );

	Convolve( phi, jac1 );
	Convolve( phiPrev, jac2 );

	/* update phi */
	for( int i = 1; i < phi->nf; i++ ) {
		phi->IndexToMode( i, k );
		k2 = -(k[0]*k[0] + k[1]*k[1]);
		rhs[0] = 2.0*(k2 - f)*phi->kVals[i][0] - 0.5*(k2 - f)*phiPrev->kVals[i][0] - 2.0*dt*jac1->kVals[i][0] + dt*jac2->kVals[i][0];
		rhs[1] = 2.0*(k2 - f)*phi->kVals[i][1] - 0.5*(k2 - f)*phiPrev->kVals[i][1] - 2.0*dt*jac1->kVals[i][1] + dt*jac2->kVals[i][1];
		if( F ) {
			rhs[0] += dt*F->kVals[i][0];
			rhs[1] += dt*F->kVals[i][1];
		}

		a11 = 1.5*(k2 - f) + dt*nu*k2*k2*k2; //biharmonic viscosity of opposite sign to regular viscosity
		a12 = -dt*beta*k[0];
		a21 = +dt*beta*k[0];
		a22 = a11;
		detInv = 1.0/(a11*a22 - a12*a21);
		
		phi->kVals[i][0] = detInv*(a22*rhs[0] - a12*rhs[1]);
		phi->kVals[i][1] = detInv*(a11*rhs[1] - a21*rhs[0]);
	}
	phi->Backward();

	delete jac1;
	delete jac2;
}

void QGEqn::RKF( Field* psi, Field* rhs ) {
	Field*	jac	= new Field( "jac", psi->nx, psi->ny );
	int	k[2], k2;

	Convolve( psi, jac );

	for( int i = 1; i < psi->nf; i++ ) {
		psi->IndexToMode( i, k );
		k2 = -(k[0]*k[0] + k[1]*k[1]);
		rhs->kVals[i][0] = (-nu*k2*k2*k2 + beta*k[0] - jac->kVals[i][0])/(k2 - f);
		rhs->kVals[i][1] = (-nu*k2*k2*k2 - beta*k[0] - jac->kVals[i][1])/(k2 - f);
		if( F ) {
			rhs->kVals[i][0] += dt*F->kVals[i][0];
			rhs->kVals[i][1] += dt*F->kVals[i][1];
		}
	}
	rhs->Backward();

	delete jac;
}

void QGEqn::SolveRK4() {
	Field*  psi	= new Field( "psi", phi->nx, phi->ny );
	Field* 	k1	= new Field( "k1", phi->nx, phi->ny );
	Field* 	k2	= new Field( "k2", phi->nx, phi->ny );
	Field* 	k3	= new Field( "k3", phi->nx, phi->ny );
	Field* 	k4	= new Field( "k4", phi->nx, phi->ny );
	Field*	in	= new Field( "in", phi->nx, phi->ny );

	psi->Copy( phi );

	for( int i = 0; i < phi->nf; i++ ) {
		in->kVals[i][0] = psi->kVals[i][0];
		in->kVals[i][1] = psi->kVals[i][1];
	}
	in->Backward();
	RKF( in, k1 );

	for( int i = 0; i < phi->nf; i++ ) {
		in->kVals[i][0] = psi->kVals[i][0] + 0.5*dt*k1->kVals[i][0];
		in->kVals[i][1] = psi->kVals[i][1] + 0.5*dt*k1->kVals[i][1];
	}
	in->Backward();
	RKF( in, k2 );

	for( int i = 0; i < phi->nf; i++ ) {
		in->kVals[i][0] = psi->kVals[i][0] + 0.5*dt*k2->kVals[i][0];
		in->kVals[i][1] = psi->kVals[i][1] + 0.5*dt*k2->kVals[i][1];
	}
	in->Backward();
	RKF( in, k3 );

	for( int i = 0; i < phi->nf; i++ ) {
		in->kVals[i][0] = psi->kVals[i][0] + dt*k3->kVals[i][0];
		in->kVals[i][1] = psi->kVals[i][1] + dt*k3->kVals[i][1];
	}
	in->Backward();
	RKF( in, k4 );

	for( int i = 0; i < phi->nf; i++ ) {
		phi->kVals[i][0] = psi->kVals[i][0] + dt*(k1->kVals[i][0] + 2.0*k2->kVals[i][0] + 2.0*k3->kVals[i][0] + k4->kVals[i][0])/6.0;
		phi->kVals[i][1] = psi->kVals[i][1] + dt*(k1->kVals[i][1] + 2.0*k2->kVals[i][1] + 2.0*k3->kVals[i][1] + k4->kVals[i][1])/6.0;
	}
	phi->Backward();

	delete psi;
	delete k1;
	delete k2;
	delete k3;
	delete k4;
	delete in;
}

void QGEqn::SolveRK3() {
	Field*	psi	= new Field( "psi", phi->nx, phi->ny );
	Field*	k1	= new Field( "k1", phi->nx, phi->ny );
	Field*	k2	= new Field( "k2", phi->nx, phi->ny );
	Field*	k3	= new Field( "k3", phi->nx, phi->ny );
	Field*	in	= new Field( "in", phi->nx, phi->ny );

	psi->Copy( phi );

	for( int i = 0; i < phi->nf; i++ ) {
		in->kVals[i][0] = psi->kVals[i][0];
		in->kVals[i][1] = psi->kVals[i][1];
	}
	RKF( in, k1 );

	for( int i = 0; i < phi->nf; i++ ) {
		in->kVals[i][0] = psi->kVals[i][0] + 0.5*dt*k1->kVals[i][0];
		in->kVals[i][0] = psi->kVals[i][0] + 0.5*dt*k1->kVals[i][0];
	}
	RKF( in, k2 );

	for( int i = 0; i < phi->nf; i++ ) {
		in->kVals[i][0] = psi->kVals[i][0] - dt*k1->kVals[i][0] + 2.0*dt*k2->kVals[i][0];
		in->kVals[i][1] = psi->kVals[i][1] - dt*k1->kVals[i][1] + 2.0*dt*k2->kVals[i][1];
	}
	RKF( in, k3 );

	for( int i = 0; i < phi->nf; i++ ) {
		phi->kVals[i][0] = psi->kVals[i][0] + dt*(k1->kVals[i][0] + 4.0*k2->kVals[i][0] + k3->kVals[i][0])/6.0;
		phi->kVals[i][1] = psi->kVals[i][1] + dt*(k1->kVals[i][1] + 4.0*k2->kVals[i][1] + k3->kVals[i][1])/6.0;
	}
	phi->Backward();

	delete psi;
	delete k1;
	delete k2;
	delete k3;
	delete in;
}
