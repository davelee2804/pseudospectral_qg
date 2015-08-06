#include <string>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <hdf5.h>
#include <fftw3.h>

#include "Field.h"
#include "Utils.h"
#include "QGEqn.h"

#define NX 512
#define NY 512

using namespace std;
using std::string;

double Sech( double x, double a ) { return 1.0/cosh(a*(x-M_PI)); }
double Tanh( double x, double a ) { return tanh(a*(x-M_PI)); }

int main( int argc, char** argv ) {
	Field*		phi		= new Field( "phi", NX, NY );
	Field*		omegaSpectral	= new Field( "omega-spectral", NX, NY );
	Field*		omegaAnalytic	= new Field( "omega-analytic", NX, NY );
	Field*		omegaAnalytic2	= new Field( "omega-analytic-2", NX, NY );
	Field*		fields[4];
	double		x, y;
	double		a = 6.0, b = 6.0;
	double		d2x, d2y, d4x, d4y;
	int		k[2], k2;

	for( int i = 0; i < NX*NY; i++ ) {
		x = 2.0*M_PI*((double)(i/NY))/NX;
		y = 2.0*M_PI*((double)(i%NY))/NY;
		phi->xVals[i] = Sech(x,a)*Sech(y,b);
		
		d2x = a*a*(Tanh(x,a)*Tanh(x,a) - Sech(x,a)*Sech(x,a))*Sech(x,a);
		d2y = b*b*(Tanh(y,b)*Tanh(y,b) - Sech(y,b)*Sech(y,b))*Sech(y,b);
		d4x = a*a*a*a*(5*Sech(x,a)*Sech(x,a)*Sech(x,a)*Sech(x,a) - 18*Tanh(x,a)*Tanh(x,a)*Sech(x,a)*Sech(x,a) + 
			       Tanh(x,a)*Tanh(x,a)*Tanh(x,a)*Tanh(x,a))*Sech(x,a);
		d4y = b*b*b*b*(5*Sech(y,b)*Sech(y,b)*Sech(y,b)*Sech(y,b) - 18*Tanh(y,b)*Tanh(y,b)*Sech(y,b)*Sech(y,b) + 
			       Tanh(y,b)*Tanh(y,b)*Tanh(y,b)*Tanh(y,b))*Sech(y,b);

		omegaAnalytic->xVals[i] = d4x*Sech(y,b) + d2x*d2y + d4y*Sech(x,a);

		d2x = 0.5*a*a*(cosh(2*a*(x-M_PI)) - 3)*pow(Sech(x,a),3);
		d2y = 0.5*b*b*(cosh(2*b*(y-M_PI)) - 3)*pow(Sech(y,b),3);
		d4x = 0.125*a*a*a*a*(cosh(4*a*(x-M_PI)) - 76*cosh(2*a*(x-M_PI)) + 115)*pow(Sech(x,a),5);
		d4y = 0.125*b*b*b*b*(cosh(4*b*(y-M_PI)) - 76*cosh(2*b*(y-M_PI)) + 115)*pow(Sech(y,b),5);

		omegaAnalytic2->xVals[i] = d4x*Sech(y,b) + d2x*d2y + d4y*Sech(x,a);
	}

	fields[0] = phi;
	fields[1] = omegaSpectral;
	fields[2] = omegaAnalytic;
	fields[3] = omegaAnalytic2;

	MeshSave( NX, NY );

	phi->Forward();
	for( int i = 1; i < phi->nf; i++ ) {
		phi->IndexToMode( i, k );
		k2 = k[0]*k[0] + k[1]*k[1];
		omegaSpectral->kVals[i][0] = k2*k2*phi->kVals[i][0];
		omegaSpectral->kVals[i][1] = k2*k2*phi->kVals[i][1];
	}
	omegaSpectral->Backward();
	
	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 4, 0, 0, 0 );
	WriteXDMFFooter( 0 );

	delete phi;
	delete omegaSpectral;
	delete omegaAnalytic;
	delete omegaAnalytic2;

	return EXIT_SUCCESS;
}
