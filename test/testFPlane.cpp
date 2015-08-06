#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <hdf5.h>
#include <fftw3.h>

#include <Field.h>
#include <QGEqn.h>
#include <Utils.h>

using namespace std;
using std::string;

#define NX 256
#define NY 256

void LoadFields( Field* psi, Field* omega ) {
	ifstream 	file;
	int		nx, ny;
	double		x[2];

	file.open( "psi.txt" );
	file >> nx;
	file >> ny;
	if( nx != psi->nx || ny != psi->ny ) {
		cerr << "ERROR: file size does not match field size." << endl;
		abort();
	}
	for( int i = 0; i < psi->nr; i++ ) {
		file >> x[0];
		file >> x[1];
		file >> psi->xVals[i];
	}
	file.close();

	file.open( "psi.txt" );
	file >> nx;
	file >> ny;
	if( nx != psi->nx || ny != psi->ny ) {
		cerr << "ERROR: file size does not match field size." << endl;
		abort();
	}
	for( int i = 0; i < psi->nr; i++ ) {
		file >> x[0];
		file >> x[1];
		file >> psi->xVals[i];
	}
	file.close();
}

int main( int argc, char** argv ) {
	Field*	psi		= new Field( "psi", NX, NY );
	Field*	omega		= new Field( "omega", NX, NY );
	Field*	psiAnal		= new Field( "psi-anal", NX, NY );
	Field*	omegaAnal	= new Field( "omega-anal", NX, NY );
	Field*	fields[4];
	double	dt		= 0.25*2*M_PI/NX;
	QGEqn*	qg		= new QGEqn( psi, 0.0, 0.0, 0.0, dt );

	fields[0] = psi;
	fields[1] = psiAnal;
	fields[2] = omega;
	fields[3] = omegaAnal;

	LoadFields( psi, omega );
	psiAnal->Copy( psi );
	omegaAnal->Copy( omega );

	MeshSave( NX, NY );

	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 4, 0.0, 0.0, 0.0 );
	WriteXDMFFooter( 0 );

	delete psi;
	delete omega;
	delete psiAnal;
	delete omegaAnal;
	delete qg;

	return EXIT_SUCCESS;
}
