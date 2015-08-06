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
	Field*		omegaSpectral	= new Field( "omega-spectral", NX, NY );
	Field*		omegaAnalytic	= new Field( "omega-analytic", NX, NY );
	Field*		phiBuffer	= new Field( "phi-buffer", NX+NX/2, NY+NY/2 );
	Field*		omegaBuffer	= new Field( "omega-buffer", NX+NX/2, NY+NY/2 );
	Field*		fields[3];
	QGEqn*		qg		= new QGEqn( phi, 0.0, 0.0, 0.0, 0.0 );
	double		x, y;
	double		a = 2.0, b = 3.0;

	for( int i = 0; i < NX*NY; i++ ) {
		x = 2.0*M_PI*((double)(i/NY))/NX;
		y = 2.0*M_PI*((double)(i%NY))/NY;
		phi->xVals[i] = Sech(x,a)*Sech(y,b);
		omegaAnalytic->xVals[i] = a*a*Sech(x,a)*(1.0 - 2.0*Sech(x,a)*Sech(x,a))*Sech(y,b) + 
					  b*b*Sech(x,a)*Sech(y,b)*(1.0 - 2.0*Sech(y,b)*Sech(y,b));
	}

	fields[0] = phi;
	fields[1] = omegaSpectral;
	fields[2] = omegaAnalytic;

	MeshSave( NX, NY );

	qg->MapToBuffer( phi, phiBuffer );
	qg->Laplacian( phiBuffer, omegaBuffer );
	qg->MapFromBuffer( omegaBuffer, omegaSpectral );

	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 3, 0, 0, 0 );
	WriteXDMFFooter( 0 );

	delete phi;
	delete omegaSpectral;
	delete omegaAnalytic;
	delete phiBuffer;
	delete omegaBuffer;
	delete qg;

	return EXIT_SUCCESS;
}
