#include <string>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <hdf5.h>
#include <fftw3.h>

#include "Field.h"
#include "Utils.h"
#include "QGEqn.h"

#define NX 256
#define NY 256

using namespace std;
using std::string;

double Sech( double x, double a ) { return 1.0/cosh(a*(x-M_PI)); }

int main( int argc, char** argv ) {
	Field*		phi		= new Field( "phi", NX, NY );
	Field*		phiSqSpectral	= new Field( "phi-sq-spectral", NX, NY );
	Field*		phiSqAnalytic	= new Field( "phi-sq-analytic", NX, NY );
	Field*		phiBuffer	= new Field( "phi-buffer", NX+NX/2, NY+NY/2 );
	Field*		phiSqBuffer	= new Field( "phi-sq-buffer", NX+NX/2, NY+NY/2 );
	Field*		fields[3];
	QGEqn*		qg		= new QGEqn( phi, 0.0, 0.0, 0.0, 0.0 );
	double		x, y;
	double		a = 1.5, b = 3.0;

	for( int i = 0; i < NX*NY; i++ ) {
		x = 2.0*M_PI*((double)(i/NY))/NX;
		y = 2.0*M_PI*((double)(i%NY))/NY;
		phi->xVals[i] = Sech(x,a)*Sech(y,b);
		phiSqAnalytic->xVals[i] = phi->xVals[i]*phi->xVals[i];
	}

	fields[0] = phi;
	fields[1] = phiSqSpectral;
	fields[2] = phiSqAnalytic;

	MeshSave( NX, NY );

	qg->MapToBuffer( phi, phiBuffer );
	for( int i = 0; i < phiBuffer->nr; i++ ) {
		phiSqBuffer->xVals[i] = phiBuffer->xVals[i]*phiBuffer->xVals[i];
	}
	qg->MapFromBuffer( phiSqBuffer, phiSqSpectral );

	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 3, 0, 0, 0 );
	WriteXDMFFooter( 0 );

/*
	MeshSave(NX+NX/2,NY+NY/2);
	WriteXDMFHeader( 0 );
	WriteXDMF( &phiBuffer, 1, 0, 0, 0 );
	WriteXDMFFooter( 0 );
*/

	delete phi;
	delete phiSqSpectral;
	delete phiSqAnalytic;
	delete phiBuffer;
	delete phiSqBuffer;
	delete qg;

	return EXIT_SUCCESS;
}
