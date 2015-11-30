/*
 *  parameters.h
 *  
 *
 *  Created by Alec Rasmussen on 8/24/15.
 *  
 *
 */
#include <complex>
#include <cmath>

#ifndef PHDATA
#define PHDATA
typedef std::complex<double> phdata;
#endif

#ifndef DCOMPLEX_H_
#define DCOMPLEX_H_
//#define J phdata(0.0,1.0)
const phdata J(0.0,1.0); // sqrt(-1)
#endif

#ifndef TWOPI
#define TWOPI
const double twoPI = 2*M_PI;
#endif

#ifndef SOL
#define SOL
const double C = 299792458.0; // Speed of light in vacuum
#endif

#ifndef RXDATA
#define RXDATA
#include <iostream>
#include <vector>

// Define receive and transmit data structures
class rxdata{
public:
	double minaz;
	double maxaz;
	int azpos;
	double cFreq;
	double bWidth;
	int freqbins;
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
		
		std::cout << az.size() << std::endl;
	}

};

rxdata::rxdata(){
	//minaz = -15*M_PI/180;
	//maxaz = 15*M_PI/180;
	//azstep = 0.01*M_PI/180;
	//freqmin = 9 * pow(10,9);
	//freqmax = 11 * pow(10,9);
	//freqbins = 1024;
	//grazing = 30*M_PI/180;
	range = 100;
	//rxphase.assign(freqbins, std::vector<phdata> (3001, {1,0}) );
	az.clear();
	freq.clear();
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
	std::vector<std::vector<phdata> > img_final;
	
	imgdata();
	
	void makeX()
	{
		for (int i = 0; i<Nx; i++)
		{
			x.push_back(x0+(-x_extent/2+i*(x_extent/(Nx-1))));
		}
	}
	
	void printX()
	{
		std::vector<double>::iterator xit;
		for ( xit = x.begin(); xit != x.end(); xit++)
		{
			std::cout << *xit << ", ";
		}
		std::cout << std::endl;
	}
	
	void printY()
	{
		std::vector<double>::iterator yit;
		for ( yit = y.begin(); yit != y.end(); yit++)
		{
			std::cout << *yit << ", ";
		}
		std::cout << std::endl;
	}
	
	void makeY()
	{
		for (int i = 0; i<Ny; i++)
		{
			y.push_back(y0+(-y_extent/2+i*(y_extent/(Ny-1))));
		}
	}
	
	void makeZ()
	{
		if (z_extent)
		{
			for (int i = 0; i<Nz; i++)
			{
				z.push_back(z0+(-z_extent/2+i*(z_extent/(Nz-1))));
			}
		}
		else 
		{
			z.push_back(0);
		} 
	}
};

imgdata::imgdata(){
	x0 = 0;
	y0 = 0;
	z0 = 0;
	x_extent = 20;
	y_extent = 20;
	z_extent = 0;
	Nfft = pow(2,14); //16,384
	Nx = 256;
	Ny = 256;
	Nz = 1;


}
#endif

#ifndef BACKPROJECT
#define BACKPROJECT

#include <fftw3.h>

void backproject(rxdata &rx, rxdata& tx, imgdata &image);//{/*...*/};

#endif

#ifndef FFTSHIFT
#define FFTSHIFT

void fftshift(std::vector<phdata> *argin, std::size_t count);//{/*...*/};

#endif

#ifndef IFFTSHIFT
#define IFFTSHIFT

void ifftshift(std::vector<phdata> *argin, std::size_t count);//{/*...*/};

#endif

#ifndef DIFF_RANGE
#define DIFF_RANGE

double dR(double txazim, double txgraz, double rxazim, double rxgraz,
	  double x, double y, double z); //{/*...*/};

#endif

#ifndef RANGE_VEC
#define RANGE_VEC

void rvec_Create(rxdata& rx, imgdata& image); //{/*...*/};

#endif

#ifndef FACTORIZE
#define FACTORIZE

void factorize(rxdata& rx, rxdata& tx, imgdata& image, int numSteps, int currStep); //{/*...*/};

#endif
