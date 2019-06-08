#pragma once
#include <sstream>
#include "Complex.h"

class Triangle
{
public:
	complex edge[3][2];

	Triangle()
	{
		for (unsigned int k = 0; k < 3; k++)
		{
				for (unsigned int cp = 0; cp < 2; cp++)
				{
					edge[k][cp] = 0.0;
				}
		}
	}

	~Triangle() {
	}

	void setNull() {

		for (unsigned int k = 0; k < 3; k++)
		{
			for (unsigned int cp = 0; cp < 2; cp++)
			{
				edge[k][cp] = 0.0;
			}
		}

	}
};