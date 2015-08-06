#include <string>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <hdf5.h>
#include <fftw3.h>

#include "Field.h"
#include "Utils.h"
#include "QGEqn.h"

#define NX 128
#define NY 128

using namespace std;
using std::string;

void InitPhi( Field* phi ) {
	int kx, ky;

	srand(101);
	for( int l = 0; l < 2; l++ ) {
		do {
			kx = rand()%8;
			//ky = rand()%8 - 4;
			//if( ky < 0 ) {
			//	ky = NY/2-1 - ky;
			//}
			ky = rand()%8;
		} while( kx == 0 || ky == 0 );
		phi->kVals[kx*(NY/2+1)+ky][0] = 0.1;
		cout << "kx[" << l << "]: " << kx << "\tky[" << l << "]: " << ky << endl;
	}
	phi->Backward();
}

int main( int argc, char** argv ) {
	Field*		phi		= new Field( "phi", NX, NY );
	Field*		phiPrev		= new Field( "phi-prev", NX, NY );
	Field*		omega		= new Field( "omega", NX, NY );
	Field*		velx		= new Field( "vel-x", NX, NY );
	Field*		vely		= new Field( "vel-y", NX, NY );
	Field*		fields[4];
	double		f		= 0.0;
	double		beta		= 0.5;
	double		dt		= 0.125*M_PI/NX;
	int		dumpEvery	= 10;
	double		nu		= 20.0/pow(NX,4);
	double		time		= 0.0;
	QGEqn*		qg		= new QGEqn( phi, beta, f, 0.0, dt );

	cout << "Time step (dt):\t" << dt << endl;
	cout << "Viscosity (nu):\t" << nu << endl;

	fields[0] = phi;
	fields[1] = omega;
	fields[2] = velx;
	fields[3] = vely;

	MeshSave( NX, NY );

	InitPhi( phi );
	qg->Laplacian( phi, omega );
	qg->Derivative( phi, velx, 1 );
	for( int i = 0; i < phi->nr; i++ ) { velx->xVals[i] *= -1.0; }
	qg->Derivative( phi, vely, 0 );
	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 4, 0, time, dt );
	WriteXDMFFooter( 0 );

	time += dt;
	phiPrev->Copy( phi );
	qg->SolveRK4();

	qg->Laplacian( phi, omega );
	qg->Derivative( phi, velx, 1 );
	for( int i = 0; i < phi->nr; i++ ) { velx->xVals[i] *= -1.0; }
	qg->Derivative( phi, vely, 0 );
	WriteXDMFHeader( 1 );
	WriteXDMF( fields, 4, 1, time, dt );
	WriteXDMFFooter( 1 );

	qg->nu = 20.0/pow(NX,4);

	for( int s = 2; s <= 25000; s++ ) {
		cout << "step: " << s << endl;
		time += dt;
		qg->Solve( phiPrev );
		if( s%dumpEvery == 0 ) {
			qg->Laplacian( phi, omega );
			qg->Derivative( phi, velx, 1 );
			for( int i = 0; i < phi->nr; i++ ) { velx->xVals[i] *= -1.0; }
			qg->Derivative( phi, vely, 0 );
			WriteXDMFHeader( s );
			WriteXDMF( fields, 4, s, time, dt );
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

	return EXIT_SUCCESS;
}
