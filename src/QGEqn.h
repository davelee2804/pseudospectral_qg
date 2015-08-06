class QGEqn {
	public:
		QGEqn( Field* _phi, double _beta, double _f, double _nu, double _dt );
		~QGEqn();
		Field* 	phi;
		double	beta;
		double 	f;
		double	nu;
		double	dt;
		Field*	F;
		void	Laplacian( Field* psi, Field* omega );
		void	Derivative( Field* psi, Field* dPsi, int dim );
		void	MapToBuffer( Field* psi, Field* psi_b );
		void	MapFromBuffer( Field* psi_b, Field* psi );
		void	Convolve( Field* psi, Field* jacobian );
		void	Solve( Field* phiPrev );
		void	SolveRK3();
		void	SolveRK4();
	private:
		void	RKF( Field* psi, Field* rhs );
		void	FirstOrder();
		void	SecondOrder( Field* phiPrev );
};
