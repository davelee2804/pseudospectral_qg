#include <string>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <hdf5.h>
#include <fftw3.h>

#include "Field.h"
#include "Utils.h"
#include "Force.h"
#include "QGEqn.h"

#define NX 256
#define NY 256

using namespace std;
using std::string;

int main( int argc, char** argv ) {
	Field*		phi		= new Field( "phi", NX, NY );
	Field*		phiPrev		= new Field( "phi-prev", NX, NY );
	Field*		omega		= new Field( "omega", NX, NY );
	Field*		velx		= new Field( "vel-x", NX, NY );
	Field*		vely		= new Field( "vel-y", NX, NY );
	Field*		fields[5];
	double		f		= 0.0;
	double		beta		= 3.0;
	double		dt		= 0.125*M_PI/NX;
	int		dumpEvery	= 40;
	double		nu		= 20.0/pow(NX,4);
	double		time		= 0.0;
	QGEqn*		qg		= new QGEqn( phi, 0.0, 0.0, 0.0, dt );
	Field*		fField		= new Field( "force", NX, NY );
	Force*		force;
	double*		A;
	int		na		= 0;
	int		kmin		= 10;
	int		kmax		= 14;
	int		k[2];

	cout << "Time step (dt):\t" << dt << endl;
	cout << "Viscosity (nu):\t" << nu << endl;

	fields[0] = phi;
	fields[1] = omega;
	fields[2] = velx;
	fields[3] = vely;
	fields[4] = fField;

	MeshSave( NX, NY );

	/* setup the forcing */
	qg->F = fField;
	for( int i = 0; i < fField->nf; i++ ) {
		fField->IndexToMode( i, k );
		if( abs(k[0])*abs(k[1]) >= kmin && abs(k[0])*abs(k[1]) <= kmax ) {
			na++;
		}
	}
	A = new double[na];
	for( int i = 0; i < na; i++ ) {
		A[i] = 1.0;
	}
	force = new Force( fField, 0.5, A, kmin, kmax );

	time += dt;
	force->Eval();
	phiPrev->Copy( phi );
	qg->SolveRK4();

	qg->Laplacian( phi, omega );
	qg->Derivative( phi, velx, 1 );
	for( int i = 0; i < phi->nr; i++ ) { velx->xVals[i] *= -1.0; }
	qg->Derivative( phi, vely, 0 );
	WriteXDMFHeader( 1 );
	WriteXDMF( fields, 5, 1, time, dt );
	WriteXDMFFooter( 1 );

	qg->nu   = nu;
	qg->beta = beta;
	qg->f    = f;

	for( int s = 2; s <= 50000; s++ ) {
		cout << "step: " << s << endl;
		time += dt;
		force->Eval();
		qg->Solve( phiPrev );
		if( s%dumpEvery == 0 ) {
			qg->Laplacian( phi, omega );
			qg->Derivative( phi, velx, 1 );
			for( int i = 0; i < phi->nr; i++ ) { velx->xVals[i] *= -1.0; }
			qg->Derivative( phi, vely, 0 );
			WriteXDMFHeader( s );
			WriteXDMF( fields, 5, s, time, dt );
			WriteXDMFFooter( s );
			WriteXDMFTemporal( s, dumpEvery );
		}
	}

	delete phi;
	delete phiPrev;
	delete omega;
	delete velx;
	delete vely;
	delete qg;
	delete fField;
	delete force;

	return EXIT_SUCCESS;
}
