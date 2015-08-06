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
	Field*		omega		= new Field( "omega", NX, NY );
	Field*		pv		= new Field( "pv", NX, NY );
	Field*		fields[3];
	double		beta		= 3.0;
	double		dt		= 0.125*M_PI/NX;
	int		dumpEvery	= 40;
	double		time		= 0.0;
	double		y;

	fields[0] = phi;
	fields[1] = omega;
	fields[2] = pv;

	for( int s = dumpEvery; s <= 50000; s += dumpEvery ) {
		cout << "step: " << s << endl;
		time = s*dt;

		phi->Read( s );
		omega->Read( s );
		for( int i = 0; i < NX*NY; i++ ) {
			y = 2.0*M_PI*((double)(i%NY))/NY - M_PI;
			pv->xVals[i] = omega->xVals[i] + beta*y;
		}

		WriteXDMFHeader( s );
		WriteXDMF( fields, 3, s, time, dt );
		WriteXDMFFooter( s );
		WriteXDMFTemporal( s, dumpEvery );
	}

	delete phi;
	delete omega;
	delete pv;

	return EXIT_SUCCESS;
}
