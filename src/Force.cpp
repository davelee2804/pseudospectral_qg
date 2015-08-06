#include <cmath>
#include <cstdlib>
#include <string>
#include <hdf5.h>
#include <fftw3.h>

#include "Field.h"
#include "Force.h"

using namespace std;
using std::string;

/* Stochastic markovian forcing of wave numbers k_min => k_max
   Reference: Maltrud & Vallis, 1991, JFM vol. 228, 321-342 */

Force::Force( Field* _field, double _R, double* _A, int _kmin, int _kmax ) {
	field	= _field;
	R	= _R;
	A	= _A;
	kmin	= _kmin;
	kmax	= _kmax;

	srand( time( NULL ) );

	fPrev 	= new Field( "force-prev", field->nx, field->ny );
}

Force::~Force() {
	delete[] A;
	delete fPrev;
}

void Force::Eval() {
	double	theta;
	int	k[2], ki = 0;

	fPrev->Copy( field );

	for( int i = 0; i < field->nf; i++ ) {
		field->IndexToMode( i, k );
		if( k[0]*k[0] + k[1]*k[1] >= kmin*kmin && k[0]*k[0] + k[1]*k[1] <= kmax*kmax ) {
			theta = 2.0*M_PI*(((double)rand())/RAND_MAX);
			field->kVals[i][0] = /*A[ki]*/sqrt(1.0 - R*R)*cos(theta) + R*fPrev->kVals[i][0];
			field->kVals[i][1] = /*A[ki]*/sqrt(1.0 - R*R)*sin(theta) + R*fPrev->kVals[i][1];
			ki++;
		}
	}
	field->Backward();
}
