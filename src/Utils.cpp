#include <iostream>
#include <fstream>
#include <cmath>
#include <hdf5.h>
#include <fftw3.h>

#include "Field.h"
#include "Utils.h"

using namespace std;
using std::string;

/*
Calculate the viscosity for the corresponding minimum wave length of two grid points
Referece:	
	Davidson, pp. 582
*/
double CalcViscosity( double l, double u, int n ) {
	double	eta	= 2.0*l/n;	/* 2D Kolmogorov microscale */
	double	nu	= u*eta*eta/l;	/* viscosity */

	return nu;
}

//#define HDF5_OLD

void MeshSave( int nx, int ny ) {
	hid_t 		file, attrib_id, group_id, attribData_id, fileSpace, fileData, fileSpace2, fileData2, memSpace;
	hsize_t 	a_dims, start[2], count[2], size[2];
	int		attribData;
	char 		filename[40];
	int 		nEls[2];
	int 		nDims		= 2;
	int		nodesPerEl	= 4;
	double		min[2]		= { 0.0, 0.0 };
	double		max[2]		= { 2.0*M_PI, 2.0*M_PI };

	sprintf( filename, "mesh.h5" );
	file = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

	/* mesh dimensionality */
	a_dims = 1;
	attribData = 2; /* mesh dimensionality */
	attribData_id = H5Screate_simple( 1, &a_dims, NULL );
#ifdef HDF5_OLD
	group_id = H5Gopen( file, "/" );
	attrib_id = H5Acreate( group_id, "dimensions", H5T_STD_I32BE, attribData_id, H5P_DEFAULT );
#else
	group_id = H5Gopen( file, "/", H5P_DEFAULT );
	attrib_id = H5Acreate( group_id, "dimensions", H5T_STD_I32BE, attribData_id, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Awrite( attrib_id, H5T_NATIVE_INT, &attribData );
	H5Aclose( attrib_id );
	H5Gclose( group_id );
	H5Sclose( attribData_id );
	
	/* mesh resolution */
	nEls[0] = nx - 1;
	nEls[1] = ny - 1;
	a_dims = 2;
	attribData_id = H5Screate_simple( 1, &a_dims, NULL );
#ifdef HDF5_OLD
	group_id = H5Gopen( file, "/" );
	attrib_id = H5Acreate( group_id, "mesh resolution", H5T_STD_I32BE, attribData_id, H5P_DEFAULT );
#else
	group_id = H5Gopen( file, "/", H5P_DEFAULT );
	attrib_id = H5Acreate( group_id, "mesh resolution", H5T_STD_I32BE, attribData_id, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Awrite( attrib_id, H5T_NATIVE_INT, nEls );
	H5Aclose( attrib_id );
	H5Gclose( group_id );
	H5Sclose( attribData_id );

	/* max and min coords of mesh */
	count[0] = (hsize_t)nDims;
	fileSpace = H5Screate_simple( 1, count, NULL );
#ifdef HDF5_OLD
	fileData = H5Dcreate( file, "/min", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
#else
	fileData = H5Dcreate( file, "/min", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Dwrite( fileData, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, min );
	H5Dclose( fileData );
	H5Sclose( fileSpace );
	fileSpace = H5Screate_simple( 1, count, NULL );
#ifdef HDF5_OLD
	fileData = H5Dcreate( file, "/max", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
#else
	fileData = H5Dcreate( file, "/max", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Dwrite( fileData, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, max );
	H5Dclose( fileData );
	H5Sclose( fileSpace );

	/* write the vertices */
 	size[0] = (hsize_t)(nx*ny);
	size[1] = (hsize_t)nDims;
	fileSpace = H5Screate_simple( 2, size, NULL );
#ifdef HDF5_OLD
	fileData = H5Dcreate( file, "/vertices", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
#else
	fileData = H5Dcreate( file, "/vertices", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif
	size[0] = (hsize_t)nEls[0]*nEls[1]; /* use linear quads for writing file */
	size[1] = (hsize_t)nodesPerEl;
	fileSpace2 = H5Screate_simple( 2, size, NULL );
#ifdef HDF5_OLD
	fileData2 = H5Dcreate( file, "/connectivity", H5T_NATIVE_INT, fileSpace2, H5P_DEFAULT );
#else
	fileData2 = H5Dcreate( file, "/connectivity", H5T_NATIVE_INT, fileSpace2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif

	count[0] = 1;
	count[1] = nDims;
	memSpace = H5Screate_simple( 2, count, NULL );
	H5Sselect_all( memSpace );

	/* assume a domain of {0,2*pi} in both dimensions */
	for( int node_i = 0; node_i < nx*ny; node_i++ ) {
		double vert[2];
		vert[0] = 2*M_PI*((double)(node_i/ny))/nx;
		vert[1] = 2*M_PI*((double)(node_i%ny))/ny;
		start[1] = 0;
		start[0] = node_i;
		H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Dwrite( fileData, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, vert );
	}
	H5Sclose( memSpace );
	H5Dclose( fileData );
	H5Sclose( fileSpace );

	H5Sget_simple_extent_dims( fileSpace2, size, NULL );
	count[0] = 1;
	count[1] = size[1];
	memSpace = H5Screate_simple( 2, count, NULL );
	H5Sselect_all( memSpace );
	
	for( int el_i = 0; el_i < nEls[0]*nEls[1]; el_i++ ) {
		int elNodes[4];
		int j = el_i/nEls[1];
		elNodes[0] = el_i + j;
		elNodes[1] = elNodes[0] + 1;
		elNodes[2] = elNodes[1] + ny;
		elNodes[3] = elNodes[2] - 1;	
		start[1] = 0;
		start[0] = el_i;
		H5Sselect_hyperslab( fileSpace2, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Dwrite( fileData2, H5T_NATIVE_INT, memSpace, fileSpace2, H5P_DEFAULT, elNodes );
	}

	H5Sclose( memSpace );
	H5Dclose( fileData2 );
	H5Sclose( fileSpace2 );
	H5Fclose( file );
}

void WriteXDMFHeader( int timeStep ) {
	ofstream 	file;
	char 		filename[20];
	char		ts[6] = "00000";

	cout << "writing fields to file at time step: " << timeStep << endl;

	sprintf( ts, "%.5u", timeStep );
	sprintf( filename, "XDMF.%s.xmf", ts );
	file.open( filename );

	file << "<?xml version=\"1.0\" ?>\n";
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
	file << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
	file << "\n";
	file << "<Domain>\n";
	file << "\n";

	file.close();
}

void WriteXDMFFooter( int timeStep ) {
	ofstream 	file;
	char 		filename[20];
	char		ts[6] = "00000";

	sprintf( ts, "%.5u", timeStep );
	sprintf( filename, "XDMF.%s.xmf", ts );
	file.open( filename, ios::app );

	file << "</Domain>\n";
	file << "\n";
	file << "</Xdmf>\n";
	file << "\n";

	file.close();
}

void WriteXDMF( Field** fields, int nFields, int timeStep, double time, double dt ) {
	ofstream 	file;
	char		filename[20];
	char		varType[40];
	Field*		field;
	char		ts[6] = "00000";
	int		nEls	= (fields[0]->nx-1)*(fields[0]->ny-1);

	sprintf( ts, "%.5u", timeStep );

	/* write the timestepping info to file */
	sprintf( filename, "time.%s.ts", ts );
	file.open( filename );
	file << "time step:\t" << timeStep << "\ttime:\t" << time << "\tdt:\t" << dt << "\n";
	file.close();

	for( int field_i = 0; field_i < nFields; field_i++ ) {
		fields[field_i]->Save( timeStep );
	}

	sprintf( varType, "NumberType=\"Float\" Precision=\"8\"" );

	sprintf( filename, "XDMF.%s.xmf", ts );
	file.open( filename, ios::app );

	/* mesh */
        file << "   <Grid Name=\"FEM_Grid_" << "mesh" << "\">\n\n";
        file << "      <Time Value=\"" << time << "\" />\n\n";
        file << "         <Topology Type=\"Quadrilateral\" NumberOfElements=\"" << nEls << "\"> \n";
        file << "            <DataItem Format=\"HDF\" DataType=\"Int\"  Dimensions=\"" << nEls << " 4\">" << "mesh" << ".h5:/connectivity</DataItem>\n";
	file << "         </Topology>\n\n";
	file << "         <Geometry Type=\"XYZ\">\n";
	file << "            <DataItem ItemType=\"Function\"  Dimensions=\"" << fields[0]->nr << " 3\" Function=\"JOIN($0, $1, 0*$1)\">\n";
	file << "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << fields[0]->nr << " 1\" Name=\"XCoords\">\n";
	file << "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 " << fields[0]->nr << " 1 </DataItem>\n";
	file << "                  <DataItem Format=\"HDF\" " << varType << " Dimensions=\"" << fields[0]->nr << " 2\">" << "mesh" << ".h5:/vertices</DataItem>\n";
	file << "               </DataItem>\n";
	file << "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << fields[0]->nr << " 1\" Name=\"YCoords\">\n";
	file << "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 1 1 1 " << fields[0]->nr << " 1 </DataItem>\n";
	file << "                  <DataItem Format=\"HDF\" " << varType << " Dimensions=\"" << fields[0]->nr << " 2\">" << "mesh" << ".h5:/vertices</DataItem>\n";
	file << "               </DataItem>\n";
	file << "            </DataItem>\n";
	file << "         </Geometry>\n\n";

	/* fields */
	for( int field_i = 0; field_i < nFields; field_i++ ) {
		field = fields[field_i];
		file << "         <Attribute Type=\"Scalar\" Center=\"Node\" Name=\"" << field->name << "\">\n";
		file << "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << fields[0]->nr << " 1\" >\n";
		file << "               <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 " << fields[0]->nr << " 1 </DataItem>\n";
		file << "               <DataItem Format=\"HDF\" " << varType << " Dimensions=\"" << fields[0]->nr << " 1\">" << field->name << "." << ts << ".h5:/data</DataItem>\n";
		file << "            </DataItem>\n";
		file << "         </Attribute>\n\n";
	}
	file << "   </Grid>\n\n";

	file.close();
}

void WriteXDMFTemporal( int nTimeSteps, int dumpEvery ) {
	ofstream	file;
	char		filename[40];

	sprintf( filename, "XDMF.temporalAll.xmf" );
	file.open( filename );
	/* header info */
        file << "<?xml version=\"1.0\" ?>\n";
        file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        file << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
        file << "\n";
        file << "<Domain>\n";
        file << "\n";
	/* fields info */
	file << "   <xi:include href=\"XDMF.FilesField.xdmf\" xpointer=\"xpointer(//Xdmf/Grid)\"/>\n\n";
	/* footer info */
	file << "</Domain>\n";
        file << "\n";
        file << "</Xdmf>\n";
        file << "\n";
	file.close();

	sprintf( filename, "XDMF.temporalFields.xmf" );
	file.open( filename );
	/* header info */
        file << "<?xml version=\"1.0\" ?>\n";
        file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        file << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
        file << "\n";
        file << "<Domain>\n";
        file << "\n";
	/* fields info */
	file << "   <xi:include href=\"XDMF.FilesField.xdmf\" xpointer=\"xpointer(//Xdmf/Grid)\"/>\n\n";
	/* footer info */
	file << "</Domain>\n";
        file << "\n";
        file << "</Xdmf>\n";
        file << "\n";
	file.close();

	sprintf( filename, "XDMF.FilesField.xdmf" );
	file.open( filename );
	file << "<?xml version=\"1.0\" ?>\n";
	file << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
	file << "<Grid GridType=\"Collection\" CollectionType=\"Temporal\" Name=\"FEM_Mesh_Fields\">\n";

	for( int step_i = 0; step_i <= nTimeSteps; step_i += dumpEvery ) {
		sprintf( filename, "XDMF.%.5u.xmf", step_i );
		file << "    <xi:include href=\"" << filename << "\" xpointer=\"xpointer(//Xdmf/Domain/Grid[1])\"/>\n";
	}

	file << "</Grid>\n";
        file << "</Xdmf>\n";

	file.close();
}
