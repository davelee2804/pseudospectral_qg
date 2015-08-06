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
	Field*		dPhiDxSpectral	= new Field( "dPhiDX-spectral", NX, NY );
	Field*		dPhiDxAnalytic	= new Field( "dPhiDX-analytic", NX, NY );
	Field*		dPhiDySpectral	= new Field( "dPhiDY-spectral", NX, NY );
	Field*		dPhiDyAnalytic	= new Field( "dPhiDY-analytic", NX, NY );
	Field*		fields[5];
	QGEqn*		qg		= new QGEqn( phi, 0.0, 0.0, 0.0, 0.0 );
	double		x, y;
	double		a = 2.0, b = 4.0;

	for( int i = 0; i < NX*NY; i++ ) {
		x = 2.0*M_PI*((double)(i/NY))/NX;
		y = 2.0*M_PI*((double)(i%NY))/NY;
		phi->xVals[i] = Sech(x,a)*Sech(y,b);
		dPhiDxAnalytic->xVals[i] = -a*tanh(a*(x-M_PI))*Sech(x,a)*Sech(y,b);
		dPhiDyAnalytic->xVals[i] = -b*tanh(b*(y-M_PI))*Sech(x,a)*Sech(y,b);
	}

	fields[0] = phi;
	fields[1] = dPhiDxSpectral;
	fields[2] = dPhiDxAnalytic;
	fields[3] = dPhiDySpectral;
	fields[4] = dPhiDyAnalytic;

	MeshSave( NX, NY );

	qg->Derivative( phi, dPhiDxSpectral, 0 );
	qg->Derivative( phi, dPhiDySpectral, 1 );
	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 5, 0, 0, 0 );
	WriteXDMFFooter( 0 );

	delete phi;
	delete dPhiDxSpectral;
	delete dPhiDxAnalytic;
	delete dPhiDySpectral;
	delete dPhiDyAnalytic;
	delete qg;

	return EXIT_SUCCESS;
}
