#pragma once
#include <sstream>
#include "Complex.h"

class Triangle
{
public:
	complex edge[3][3][4];

	Triangle()
	{
		for (unsigned int k1 = 0; k1 < 3; k1++)
		{
			for (unsigned int k2 = 0; k2 < 3; k2++)
			{
				for (unsigned int cp = 0; cp < 4; cp++)
				{
					edge[k1][k2][cp] = 0.0;
				}
			}
		}
	}

	~Triangle(){
	}

	void setNull(){
		for (unsigned int k1 = 0; k1 < 3; k1++)
		{
			for (unsigned int k2 = 0; k2 < 3; k2++)
			{
				for (unsigned int cp = 0; cp < 4; cp++)
				{
					edge[k1][k2][cp] = complex(0.0,0.0);
				}
			}
		}
	}
};