#include <string>
#include <hdf5.h>
#include <fftw3.h>

#include "Field.h"

using namespace std;
using std::string;

Field::Field( string _name, int _nx, int _ny ) {
	name	= _name;
	nx	= _nx;
	ny	= _ny;

	nr	= nx*ny;
	nf	= nx*(ny/2+1);

	xVals 	= (double*)fftw_malloc(nr*sizeof(double));
	kVals	= (fftw_complex*)fftw_malloc(nf*sizeof(fftw_complex));
	kBuff	= (fftw_complex*)fftw_malloc(nf*sizeof(fftw_complex));

	forward  = fftw_plan_dft_r2c_2d( nx, ny, xVals, kVals, FFTW_MEASURE );
	backward = fftw_plan_dft_c2r_2d( nx, ny, kVals, xVals, FFTW_MEASURE );

	for( int i = 0; i < nr; i++ ) {
		xVals[i] = 0.0;
	}
	for( int i = 0; i < nf; i++ ) {
		kVals[i][0] = 0.0;
		kVals[i][1] = 0.0;
	}
}

Field::~Field() {
	fftw_destroy_plan(forward);
	fftw_destroy_plan(backward);
	fftw_free(xVals);
	fftw_free(kVals);
	fftw_free(kBuff);
}

void Field::IndexToNode( int i, int* n ) {
	n[0] = i/nx;
	n[1] = i%nx;
}

void Field::IndexToMode( int i, int* m ) {
	m[0] = i/(ny/2+1);
	m[1] = i%(ny/2+1);
	if( m[0] > nx/2+1 ) {
		m[0] = m[0] - nx;
	}
}

void Field::Copy( Field* field ) {
	for( int i = 0; i < nr; i++ ) {
		xVals[i] = field->xVals[i];
	}
	for( int i = 0; i < nf; i++ ) {
		kVals[i][0] = field->kVals[i][0];
		kVals[i][1] = field->kVals[i][1];
	}
}

void Field::Forward() {
	double scale = 2.0/nr;

	fftw_execute( forward );
	for( int i = 0; i < nf; i++ ) {
		kVals[i][0] *= scale;
		kVals[i][1] *= scale;
	}
}

void Field::Backward() {
	/* preserve the input by copying it to a buffer and then copying back */
	for( int i = 0; i < nf; i++ ) {
		kBuff[i][0] = kVals[i][0];
		kBuff[i][1] = kVals[i][1];
	}

	fftw_execute( backward );
	for( int i = 0; i < nr; i++ ) {
		xVals[i] *= 0.5;
	}
	for( int i = 0; i < nf; i++ ) {
		kVals[i][0] = kBuff[i][0];
		kVals[i][1] = kBuff[i][1];
	}
}

//#define HDF5_OLD

void Field::Read( int timeStep ) {
	hid_t		file, fileData, fileSpace, memSpace;
	hsize_t		count[2], start[2];
	double*		buf		= new double[1];
	char		filename[50];

	sprintf( filename, "%s.%.5u.h5", name.c_str(), timeStep );

	file = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
#ifdef HDF5_OLD
	fileData = H5Dopen( file, "/data" );
#else
	fileData = H5Dopen( file, "/data", H5P_DEFAULT );
#endif
	fileSpace = H5Dget_space( fileData );

	count[0] = 1;
	count[1] = 1; /* no. degrees of freedom */
	memSpace = H5Screate_simple( 2, count, NULL );
	start[1] = 0;

	for( int node_i = 0; node_i < nr; node_i++ ) {
		start[0] = node_i;
		H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Sselect_all( memSpace );
		H5Dread( fileData, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, buf );
		xVals[node_i] = buf[0];
	}
	delete[] buf;

	H5Dclose( fileData );
	H5Sclose( memSpace );
	H5Sclose( fileSpace );
	H5Fclose( file );
}

void Field::Save( int timeStep ) {
	char		filename[50];
	hid_t		file, attribData_id, attrib_id, group_id, fileSpace, fileData, memSpace;
	hsize_t		a_dims, start[2], count[2], size[2];
	int		attribData;
	int		nLinearEls[2];
	double		buf[1];	/* 1 degree of freedom */

	sprintf( filename, "%s.%.5u.h5", name.c_str(), timeStep );
	file = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

	/* write the dimensionality */
	a_dims = 1;
	attribData = 2; /* no. dimensions */
	attribData_id = H5Screate_simple( 1, &a_dims, NULL );
#ifdef HDF5_OLD
	group_id  = H5Gopen( file, "/" );
	attrib_id = H5Acreate( group_id, "dimensions", H5T_STD_I32BE, attribData_id, H5P_DEFAULT );
#else
	group_id  = H5Gopen( file, "/", H5P_DEFAULT );
	attrib_id = H5Acreate( group_id, "dimensions", H5T_STD_I32BE, attribData_id, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Awrite( attrib_id, H5T_NATIVE_INT, &attribData );
	H5Aclose( attrib_id );
	H5Gclose( group_id );
	H5Sclose( attribData_id );
	
	a_dims = 2; /* dimensions of the field */
	nLinearEls[0] = nx - 1;
	nLinearEls[1] = ny - 1;
	attribData_id = H5Screate_simple( 1, &a_dims, NULL );
#ifdef HDF5_OLD 
	group_id  = H5Gopen(file, "/" );
	attrib_id = H5Acreate(group_id, "mesh resolution", H5T_STD_I32BE, attribData_id, H5P_DEFAULT );
#else
	group_id  = H5Gopen(file, "/", H5P_DEFAULT );
	attrib_id = H5Acreate(group_id, "mesh resolution", H5T_STD_I32BE, attribData_id, H5P_DEFAULT, H5P_DEFAULT );
#endif
        H5Awrite( attrib_id, H5T_NATIVE_INT, nLinearEls );
	H5Aclose( attrib_id );
	H5Gclose( group_id );
	H5Sclose( attribData_id );

	size[0] = nr;
	size[1] = 1; /* no. degrees of freedom */
	fileSpace = H5Screate_simple( 2, size, NULL );
#ifdef HDF5_OLD
	fileData  = H5Dcreate( file, "/data", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
#else
	fileData  = H5Dcreate( file, "/data", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif

	count[0] = 1;
	count[1] = 1; /* no. degrees of freedom */
	memSpace = H5Screate_simple( 2, count, NULL );
	H5Sselect_all( memSpace );
	for( int node_i = 0; node_i < nr; node_i++ ) {
		buf[0] = xVals[node_i];
		start[0] = node_i;
		start[1] = 0;
		H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Dwrite( fileData, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, buf );
	}

	H5Dclose( fileData );
	H5Sclose( memSpace );
	H5Sclose( fileSpace );
	H5Fclose( file );
}
