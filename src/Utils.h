double CalcViscosity( double l, double u, int n );
void   MeshSave( int nx, int ny );
void   WriteXDMFHeader( int timeStep );
void   WriteXDMFFooter( int timeStep );
void   WriteXDMF( Field** fields, int nFields, int timeStep, double time, double dt );
void   WriteXDMFTemporal( int nTimeSteps, int dumpEvery );
