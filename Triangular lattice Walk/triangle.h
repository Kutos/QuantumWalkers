#pragma once
#include <sstream>
#include "Complex.h"

typedef struct { complex e0[2]; complex e1[2]; complex e2[2];} triangle;

void setNulltriangle(triangle *T) {
	(*T).e0[0] = 0.0;
	(*T).e0[1] = 0.0;
	(*T).e1[0] = 0.0;
	(*T).e1[1] = 0.0;
	(*T).e2[0] = 0.0;
	(*T).e2[1] = 0.0;
}