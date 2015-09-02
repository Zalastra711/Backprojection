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
	minaz = -15;
	maxaz = 15;
	azstep = 0.01;
	freqmin = 9 * pow(10,9);
	freqmax = 11 * pow(10,9);
	freqbins = 512;
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
	double x_mat[128][128];
	double y_mat[128][128];
	double z_mat[128][128];
	std::vector<std::vector<phdata> > img_final;
	
	imgdata();
};

imgdata::imgdata(){
	x0 = 0;
	y0 = 0;
	z0 = 0;
	x_extent = 10;
	y_extent = 10;
	z_extent = 0;
	Nfft = pow(2,14);
	Nx = 128;
	Ny = 128;
	Nz = 128;
}
#endif

#ifndef BACKPROJECT
#define BACKPROJECT

#include <cmath>
#include <complex>
#include <fftw3.h>

void backproject(rxdata& rx, rxdata& tx, imgdata& image);//{/*...*/};

#endif
