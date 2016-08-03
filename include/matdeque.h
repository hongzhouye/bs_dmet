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
	Matdeq ();
	Matdeq (int);
	MatrixKd * element;
	void allocate (void);
	void append (MatrixKd &);
	void print (void);
};

Matdeq::Matdeq (void)
{
	max_size = 5;
	now = 0;
	full_flag = 0;
}

Matdeq::Matdeq (int n)
{
	max_size = n;
	now = 0;
	full_flag = 0;
}

void Matdeq::allocate (void)
{
	element = new MatrixKd [max_size];
}

void Matdeq::append (MatrixKd& a)
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
