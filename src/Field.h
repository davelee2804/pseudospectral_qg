class Field {
	public:
		Field( std::string _name, int _nx, int _ny );
		~Field();
		std::string 	name;
		int		nx;				/* no. nodes in x */
		int		ny;				/* no. nodes in y */
		int		nr;				/* total nodes in real space */
		int		nf;				/* total modes in fourier space */
		double* 	xVals;				/* real space data */
		fftw_complex*	kVals;				/* fourier space data */
		fftw_plan	forward;			/* real to fourier transform */
		fftw_plan	backward;			/* fourier to real transform */
		void		IndexToNode( int i, int* n );	/* node index to {x_i, x_j} coordinate */
		void		IndexToMode( int i, int* m );	/* mode index to {k_i, k_j} coordinate */
		void		Forward();
		void		Backward();
		void		Copy( Field* field );
		void		Read( int timeStep );
		void		Save( int timeStep );
	private:
		fftw_complex*	kBuff;
};
