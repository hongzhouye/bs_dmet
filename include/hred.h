#include <Eigen/Dense>
#include "hf.h"
#include "hubbard.h"
#include "schmidt.h"

// Reduced Hamiltonian
class HRED
{
	public:
		int Ni;			// 2 * Nimp
		MatrixXd hup, hdn;
		double *Vup, *Vdn;
		void _xform_ (HUBBARD&, SCHMIDT&);
};

void HRED::_xform_ (HUBBARD& hub, SCHMIDT& sm)
{
	int i, j, k, l, mu;
	MatrixXd Tup = sm.Tup, Tdn = sm.Tdn;

	// Environment's contribution to himp
	MatrixXd hcup = ((sm.TEup * sm.TEup.transpose ()).diagonal () * hub.U).asDiagonal ();
	MatrixXd hcdn = ((sm.TEdn * sm.TEdn.transpose ()).diagonal () * hub.U).asDiagonal ();

	// himp
	hup = Tup.transpose () * (hub.hup + hcup) * Tup;
	hdn = Tdn.transpose () * (hub.hdn + hcdn) * Tdn;

	cout << "himp up:\n" << hup << "\n\n";

	// Vimp
	Ni = 2 * sm.Nimp;	
	Vup = _darray_gen_ (Ni * Ni * Ni * Ni);
	Vdn = _darray_gen_ (Ni * Ni * Ni * Ni);
	for (i = 0; i < Ni; i++)
		for (j = i; j < Ni; j++)
			for (k = j; k < Ni; k++)
				for (l = k; l < Ni; l++)
				{
					for (mu = 0; mu < hub.K; mu++)
					{
						Vup[index4(i,j,k,l,Ni)] += Tup (mu, i) * Tup (mu, j) * 
							Tup (mu, k) * Tup (mu, l);
						Vdn[index4(i,j,k,l,Ni)] += Tdn (mu, i) * Tdn (mu, j) * 
							Tdn (mu, k) * Tdn (mu, l);
					}
					Vup[index4(i,j,k,l,Ni)] *= hub.U;
					Vdn[index4(i,j,k,l,Ni)] *= hub.U;
					Vup[index4(i,j,l,k,Ni)] = Vup[index4(i,k,j,l,Ni)] = Vup[index4(i,k,l,j,Ni)] = 
						Vup[index4(i,l,j,k,Ni)] = Vup[index4(i,l,k,j,Ni)] = Vup[index4(j,i,k,l,Ni)] =
						Vup[index4(j,i,l,k,Ni)] = Vup[index4(j,k,i,l,Ni)] = Vup[index4(j,k,l,i,Ni)] =
						Vup[index4(j,l,i,k,Ni)] = Vup[index4(j,l,k,i,Ni)] = Vup[index4(k,i,j,l,Ni)] =
						Vup[index4(k,i,l,j,Ni)] = Vup[index4(k,j,i,l,Ni)] = Vup[index4(k,j,l,i,Ni)] =
						Vup[index4(k,l,i,j,Ni)] = Vup[index4(k,l,j,i,Ni)] = Vup[index4(l,i,j,k,Ni)] =
						Vup[index4(l,i,k,j,Ni)] = Vup[index4(l,j,i,k,Ni)] = Vup[index4(l,j,k,i,Ni)] =
						Vup[index4(l,k,i,j,Ni)] = Vup[index4(l,k,j,i,Ni)] = Vup[index4(i,j,k,l,Ni)];
					Vdn[index4(i,j,l,k,Ni)] = Vdn[index4(i,k,j,l,Ni)] = Vdn[index4(i,k,l,j,Ni)] = 
						Vdn[index4(i,l,j,k,Ni)] = Vdn[index4(i,l,k,j,Ni)] = Vdn[index4(j,i,k,l,Ni)] =
						Vdn[index4(j,i,l,k,Ni)] = Vdn[index4(j,k,i,l,Ni)] = Vdn[index4(j,k,l,i,Ni)] =
						Vdn[index4(j,l,i,k,Ni)] = Vdn[index4(j,l,k,i,Ni)] = Vdn[index4(k,i,j,l,Ni)] =
						Vdn[index4(k,i,l,j,Ni)] = Vdn[index4(k,j,i,l,Ni)] = Vdn[index4(k,j,l,i,Ni)] =
						Vdn[index4(k,l,i,j,Ni)] = Vdn[index4(k,l,j,i,Ni)] = Vdn[index4(l,i,j,k,Ni)] =
						Vdn[index4(l,i,k,j,Ni)] = Vdn[index4(l,j,i,k,Ni)] = Vdn[index4(l,j,k,i,Ni)] =
						Vdn[index4(l,k,i,j,Ni)] = Vdn[index4(l,k,j,i,Ni)] = Vdn[index4(i,j,k,l,Ni)];
				}

}
