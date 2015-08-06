#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include "fftw3.h"

#include "Field.h"
#include "Utils.h"

#define NX 256
#define NY 256

using namespace std;
using std::string;

int main( int argc, char** argv ) {
	int		step		= atoi( argv[1] );
	int		maxStep		= atoi( argv[2] );
	double		pr, pi, wr, wi;
	int		k[2];
	double*		Ex		= new double[NY/2+1];
	double*		Ey		= new double[NX/2+1];
	double*		Wx		= new double[NX/2+1];
	double*		Wy		= new double[NY/2+1];
	Field*		phi		= new Field( "phi", NX, NY );
	Field*		energy		= new Field( "energy", NX, NY );
	char		filename[20];
	ofstream	file;

	for( int s = step; s <= maxStep; s += step ) {
		cout << "step: " << s << endl;
		phi->Read( s );
		phi->Forward();

		/* kinetic energy for the x-wavenumbers */
		for( int ky = 0; ky < NY/2+1; ky++ ) {
			pr = phi->kVals[ky][0];
			pi = phi->kVals[ky][1];
			Ex[ky] = ky*ky*( pr*pr + pi*pi );
			wr = ky*ky*pr;
			wi = ky*ky*pi;
			Wy[ky] = wr*wr + wi*wi;
		}
		/* kinetic energy of the y-wavenumbers */
		for( int kx = 0; kx < NX/2+1; kx ++ ) {
			pr = phi->kVals[kx*(NY/2+1)][0];
			pi = phi->kVals[kx*(NY/2+1)][1];
			Ey[kx] = kx*kx*( pr*pr + pi*pi );
			wr = kx*kx*pr;
			wi = kx*kx*pi;
			Wx[kx] = wr*wr + wi*wi;
		}

		sprintf( filename, "phi.%.5u.en", s );
		file.open( filename );
		for( int kx = 0; kx < NX/2+1; kx++ ) {
			file << kx << "\t" << Ex[kx] << "\t" << Ey[kx] << "\t" << Wx[kx] << "\t" << Wy[kx] << endl;
		}
		file.close();

		/* kinetic energy of the full field in fourier space */
		/*
		for( int i = 0; i < phi->nf; i++ ) {
			phi->IndexToMode( i, k );
			pr = phi->kVals[i][0];
			pi = phi->kVals[i][1];
			energy->xVals[(k[0]+NX/2)*NY+k[1]] = (k[0]*k[0] + k[1]*k[1])*(pr*pr + pi*pi);
		}
		WriteXDMFHeader( 99999 );
		WriteXDMF( &energy, 1, 99999, 0.0, 0.0 );
		WriteXDMFFooter( 99999 );
		*/
	}

	delete[] Ex;
	delete[] Ey;
	delete[] Wx;
	delete[] Wy;
	delete phi;

	return EXIT_SUCCESS;
}
