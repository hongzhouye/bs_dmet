#include <Eigen/Dense>
#include <cstdio>
#include "hf.h"
#include "hubbard.h"
#include "schmidt.h"

// Reduced Hamiltonian
class HRED
{
	public:
		int Ni;			// 2 * Nimp
		MatrixXd h;
		double *V;
		void _xform_ (HUBBARD&, SCHMIDT&);
};

void HRED::_xform_ (HUBBARD& hub, SCHMIDT& sm)
{
	int i, j, k, l, mu;
	MatrixXd T = sm.T;

	// Environment's contribution to himp
	MatrixXd hc = ((sm.TE * sm.TE.transpose ()).diagonal () * hub.U).asDiagonal ();

	// himp
	h = T.transpose () * (hub.h + hc) * T;

	cout << "himp :\n" << h << "\n\n";

	// Vimp
	Ni = 2 * sm.Nimp;	
	V = _darray_gen_ (Ni * Ni * Ni * Ni);
	for (i = 0; i < Ni; i++)
		for (j = i; j < Ni; j++)
			for (k = j; k < Ni; k++)
				for (l = k; l < Ni; l++)
				{
					for (mu = 0; mu < hub.K; mu++)
						V[index4(i,j,k,l,Ni)] += T (mu, i) * T (mu, j) * 
							T (mu, k) * T (mu, l);
					V[index4(i,j,k,l,Ni)] *= hub.U;
					V[index4(i,j,l,k,Ni)] = V[index4(i,k,j,l,Ni)] = V[index4(i,k,l,j,Ni)] = 
						V[index4(i,l,j,k,Ni)] = V[index4(i,l,k,j,Ni)] = V[index4(j,i,k,l,Ni)] =
						V[index4(j,i,l,k,Ni)] = V[index4(j,k,i,l,Ni)] = V[index4(j,k,l,i,Ni)] =
						V[index4(j,l,i,k,Ni)] = V[index4(j,l,k,i,Ni)] = V[index4(k,i,j,l,Ni)] =
						V[index4(k,i,l,j,Ni)] = V[index4(k,j,i,l,Ni)] = V[index4(k,j,l,i,Ni)] =
						V[index4(k,l,i,j,Ni)] = V[index4(k,l,j,i,Ni)] = V[index4(l,i,j,k,Ni)] =
						V[index4(l,i,k,j,Ni)] = V[index4(l,j,i,k,Ni)] = V[index4(l,j,k,i,Ni)] =
						V[index4(l,k,i,j,Ni)] = V[index4(l,k,j,i,Ni)] = V[index4(i,j,k,l,Ni)];
				}

	// check Vimp
	/*for (i = 0; i < Ni; i++)
		for (j = 0; j < Ni; j++)
		{
			printf ("The (%d, %d) block:\n", i, j);
			for (k = 0; k < Ni; k++)
			{
				for (l = 0; l < Ni; l++)
					printf ("%.4f  ", V[index4(i,j,k,l,Ni)]);
				printf ("\n");
			}
			printf ("\n");
		}
	*/
}
