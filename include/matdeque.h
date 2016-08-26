#ifndef _MATDEQUE_H_INCLUDED_
#define _MATDEQUE_H_INCLUDED_

#include <iostream>
#include "hf.h"

using namespace std;
using namespace Eigen;

class Matdeq
{
  public:
	int max_size, now, full_flag;
	MatrixXd *element;
	void _init_ (int);
	void allocate (int);
	void append (MatrixXd&);
	void print (void);
};

void Matdeq::_init_ (int n)
{
	max_size = n;
	now = 0;
	full_flag = 0;
}

void Matdeq::allocate (int K)
{
	element = new MatrixXd [max_size];
    int i = 0;
    for (i = 0; i < max_size; i++)  element[i].setZero (K, K);
}

void Matdeq::append (MatrixXd& a)
{
	element[now] = a;
	if (full_flag == 0 && now == max_size - 1)
	{
		full_flag = 1;
		now = 0;
		return;
	}
	now = (now + 1) % max_size;
}

void Matdeq::print (void)
{
	int i = 0, max = (full_flag) ? (max_size) : (now);
	for (i = 0; i < max; i++)
		cout << element[i] << endl << endl;
}

#endif
