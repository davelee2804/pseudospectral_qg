class Force {
	public:
		Force( Field* _field, double _R, double* _A, int _kmin, int _kmax );
		~Force();
		double	R;
		double*	A;
		int	kmin;
		int	kmax;
		Field* 	field;
		void	Eval();
	private:
		Field*	fPrev;
};
