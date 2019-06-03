#pragma once
#include<math.h>
#include<string>

class complex {

private:
	/* Coefficients */
	double r;
	double i;

	/* Operations */
	void add(complex c) {
		r += c.real();
		i += c.imaginary();
	}
	void add(double d) {
		r += d;
	}

	void subs(complex c) {
		r -= c.real();
		i -= c.imaginary();
	}
	void subs(double d) {
		r -= d;
	}

	void multiply(complex c) {
		double r1 = r, i1 = i;
		double r2 = c.real(), i2 = c.imaginary();
		r = r1 * r2 - i1 * i2;
		i = r1 * i2 + r2 * i1;
	}
	void multiply(double d) {
		r *= d;
		i *= d;
	}

	void divide(complex c) {
		double r1 = r, i1 = i;
		double r2 = c.real(), i2 = c.imaginary();
		double d = r2 * r2 + i2 * i2;

		r = (r1 * r2 + i1 * i2) / d;
		i = (r2 * i1 - r1 * i2) / d;
	}
	void divide(double d) {
		r *= 1 / d;
		i *= 1 / d;
	}

public:
	/* Construteurs */
	complex(double real, double imag) {
		r = real;
		i = imag;
	}
	complex(double real) {
		r = real;
		i = 0.0;
	}
	complex(int real, int imag) {
		r = double(real);
		i = double(imag);
	}
	complex(int real) {
		r = double(real);
		i = 0.0;
	}
	complex(const char* cs) {

		std::string str(cs);
		int numr = str.find('r');
		int numi = str.find('i');

		char* real = new char[numr];
		for (int i = 0; i < numr; i++) {
			real[i] = cs[i];
		}

		char* imag = new char[numi - numr - 1];
		for (int i = numr + 1; i < numi; i++) {
			imag[i - (numr + 1)] = cs[i];
		}

		r = atof(real);
		i = atof(imag);
		delete[] real;
		delete[] imag;
	}
	complex() {
		r = 0.0;
		i = 0.0;
	}
	/* Valeurs */
	double real() {
		return r;
	}

	double imaginary() {
		return i;
	}

	double abs2() {
		return r * r + i * i;
	}

	complex conjugate() {
		return complex(r, -i);
	}

	/* Setters */
	void setReal(double real) {
		r = real;
	}
	void setImaginary(double imaginary) {
		i = imaginary;
	}

	/* Operators */
	complex operator +=(complex const& c) {
		add(c);
		return *this;
	}
	complex operator +=(double const& d) {
		add(d);
		return *this;
	}
	complex operator +(complex const& c) {
		complex result(r, i);
		result += c;
		return result;
	}
	complex operator +(double const& d) {
		complex result(r, i);
		result += d;
		return result;
	}

	complex operator -=(complex const& c) {
		subs(c);
		return *this;
	}
	complex operator -=(double const& d) {
		subs(d);
		return *this;
	}
	complex operator -(complex const& c) {
		complex result(r, i);
		result -= c;
		return result;
	}
	complex operator -(double const& d) {
		complex result(r, i);
		result -= d;
		return result;
	}

	complex operator *=(complex const& c) {
		multiply(c);
		return *this;
	}
	complex operator *=(double const& d) {
		multiply(d);
		return *this;
	}
	complex operator *(complex const& c) {
		complex result(r, i);
		result *= c;
		return result;
	}
	complex operator *(double const& d) {
		complex result(r, i);
		result *= d;
		return result;
	}

	complex operator /=(complex const& c) {
		divide(c);
		return *this;
	}
	complex operator /=(double const& d) {
		divide(d);
		return *this;
	}
	complex operator /(complex const& c) {
		complex result(r, i);
		result /= c;
		return result;
	}
	complex operator /(double const& d) {
		complex result(r, i);
		result /= d;
		return result;
	}

	friend complex operator +(double const&, complex const&);
	friend complex operator -(double const&, complex const&);
	friend complex operator *(double const&, complex const&);
	friend complex operator /(double const&, complex const&);
};
complex operator +(double const& d, complex const& c) {
	complex result(c.r, c.i);
	result += d;
	return result;
}
complex operator -(double const& d, complex const& c) {
	complex result(c.r, c.i);
	result -= d;
	return result;
}
complex operator *(double const& d, complex const& c) {
	complex result(c.r, c.i);
	result *= d;
	return result;
}
complex operator /(double const& d, complex const& c) {
	complex result(d);
	result /= c;
	return result;
}
complex exp(complex c) {
	return complex(exp(c.real()))* complex(cos(c.imaginary()), sin(c.imaginary()));
}
const char* c_complex(complex c) {
	std::string s = "";
	if (c.real() > 0) s += "+";
	s += std::to_string(c.real());
	s += "r";
	if (c.imaginary() > 0) s += "+";
	s += std::to_string(c.imaginary());
	s += "i";
	return (s.c_str());
}
std::string str_complex(complex c) {
	std::string s = "";
	if (c.real() >= 0) s += "+";
	s += std::to_string(c.real());
	s += "r";
	if (c.imaginary() >= 0) s += "+";
	s += std::to_string(c.imaginary());
	s += "i";
	return s;
}