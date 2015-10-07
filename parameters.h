/*
 *  parameters.h
 *  
 *
 *  Created by Alec Rasmussen on 8/24/15.
 *  
 *
 */
#include <complex>
// using std::complex;

#ifndef PHDATA
#define PHDATA
typedef std::complex<double> phdata;
#endif

#ifndef DCOMPLEX_H_
#define DCOMPLEX_H_
const phdata I = (0.0,1.0); // sqrt(-1)
#endif

#ifndef TWOPI
#define TWOPI
const double twoPI = 2*M_PI;
#endif

#ifndef SOL
#define SOL
const double c = 299792458; // Speed of light in vacuum
#endif

#ifndef RXDATA
#define RXDATA
#include <iostream>
#include <cmath>
#include <vector>
// using std::vector;
//#include <Array>

// Define receive and transmit data structures
class rxdata{
public:
	double minaz;
	double maxaz;
	double azstep;
	double freqmin;
	double freqmax;
	double freqbins;
	double grazing; //Degrees
	double range;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	std::vector<double> az;
	std::vector<double> freq;
	std::vector<std::vector<phdata> > rxphase;
	
	rxdata();
	
	void printAz()
	{
		for ( unsigned int i = 0; i < az.size(); i++)
		{
			std::cout << az[i] << ", ";
		}
		std::cout << std::endl;
	}
};

rxdata::rxdata(){
	minaz = -5*M_PI/180;
	maxaz = 5*M_PI/180;
	azstep = 1*M_PI/180;
	freqmin = 9 * pow(10,9);
	freqmax = 11 * pow(10,9);
	freqbins = 64;
	grazing = 30;
	range = 100;
	rxphase.assign(freqbins, std::vector<phdata> (3001, {1,0}) );
}

#endif

#ifndef IMGDATA
#define IMGDATA
// #include <complex>
// Define the image grid
class imgdata {
public:
	double x0;
	double y0;
	double z0;
	double x_extent;
	double y_extent;
	double z_extent;
	int Nfft;
	int Nx;
	int Ny;
	int Nz;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	std::vector<double> r_vec;
	std::vector<phdata> rangeprofile;
	//double x_mat[128][128];
	//double y_mat[128][128];
	//double z_mat[128][128];
	phdata img_final[128][128][1];
	
	imgdata();
};

imgdata::imgdata(){
	x0 = 0;
	y0 = 0;
	z0 = 0;
	x_extent = 10;
	y_extent = 10;
	z_extent = 1;
	Nfft = pow(2,8);
	Nx = 128;
	Ny = 128;
	Nz = 1;
	for (int i = 0; i<Nx; i++)
	{
		x.push_back(-x_extent/2+i*((x_extent/Nx)));
	}
	for (int i = 0; i<Ny; i++)
	{
		y.push_back(-y_extent/2+i*((y_extent/Ny)));
	}
	if (z_extent)
	{
		for (int i = 0; i<Nz; i++)
		{
			z.push_back(-z_extent/2+i*((z_extent/Nz)));
		}
	}
}
#endif

#ifndef BACKPROJECT
#define BACKPROJECT

#include <cmath>
#include <complex>
#include <fftw3.h>

void backproject(rxdata& rx, rxdata& tx, imgdata& image);//{/*...*/};

#endif

#ifndef SWAP
#define SWAP

void swap(std::complex<double> *var1, std::complex<double> *var2);//{/*...*/};

#endif

#ifndef FFTSHIFT
#define FFTSHIFT

#include <cmath>

void fftshift(std::vector<phdata> *argin, std::size_t count);//{/*...*/};

#endif

#ifndef IFFTSHIFT
#define IFFTSHIFT

#include <cmath>

void ifftshift(std::vector<phdata> *argin, std::size_t count);//{/*...*/};

#endif

#ifndef DIFF_RANGE
#define DIFF_RANGE

double dR(double txazim, double txgraz, double rxazim, double rxgraz,
	  double x, double y, double z); //{/*...*/};

#endif
