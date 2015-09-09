/*
 *  BPmain.cpp
 *  
 *
 *  Created by Alec Rasmussen on 8/24/15.
 *  
 *
 */
#include <iostream>
// #include <stdio.h>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <numeric>
#include "parameters.h"

// Compile with:
// g++ -I .. -fopenmp <filename.cc> <path to library> -lfftw2 -lfftw3_omp

//using namespace std;
//using namespace fftwpp;

// This program is designed to take input phase history data from a file 
// and use "fast backprojection" to create an image of a scene.


//using namespace std;
//using namespace fftw;

void swap(complex<double> *var1, complex<double> *var2)
{
	complex<double> tmp = *var1;
	*var1 = *var2;
	*var2 = tmp;
}

void fftshift(complex<double> *argin, int count)
{
	int k = 0;
	int c = (int) floor((float)count/2);
	
	if (count % 2 == 0) //define case for even length vector
	{
		for (k=0; k < c; k++)
		{	
			swap(&argin[k], &argin[k+c]);
		}
	}
	else
	{
		complex<double> tmp = argin[0];
		for (k=0; k<c;k++)
		{
			argin[k] = argin[c+k+1];
			argin[c+k+1] = argin[k+1];
		}
		argin[c]=tmp;
	}
}

void ifftshift(complex<double> *argin, int count)
{
	int k = 0;
	int c = (int)floor((float)count/2);
	if (count % 2 == 0)
	{
		for (k=0; k<c;k++)
		{
			swap(&argin[k], &argin[k+c]);
		}
	}
	else
	{
		complex<double> tmp = argin[count-1];
		for (k=c-1; k>=0; k--)
		{
			argin[c+k+1] = argin[k];
			argin[k] = argin[c+k];
		}
		argin[c] = tmp;
	}
}

void backproject(rxdata& rx, rxdata& tx, imgdata& image){
	
	std::vector<phdata> filteredData;
	
	fftw_complex *in, *out;
	fftw_plan plan;
	double deltaFreq, diffFreq, sampleRange;
	
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image.Nfft );
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image.Nfft );
	plan = fftw_plan_dft_1d(image.Nfft, in, out, FFTW_BACKWARD, FFTW_MEASURE);
	
	double x,y,z;
	std::vector<double> fx,fy,fz,radonFilter;
	
	int pulses = rx.az.size();
	double fcenter = std::accumulate(rx.freq.begin(),rx.freq.end(),0.0)/rx.freq.size();
	
	for ( unsigned int i=0; i < tx.freq.size(); i++ )
	{
		diffFreq = tx.freq[i+1] - tx.freq[i];
		deltaFreq += diffFreq;
	};
	
	// Assumes uniformly sampled frequency vector
	deltaFreq /= tx.freq.size(); 
	// Range sample spacing corresponding to collected frequency samples
	sampleRange = c/(image.Nfft*deltaFreq); 
	
	// Define the range at each sample
	for ( int i = 0; i < image.Nfft; i++ )
	{
		image.r_vec.push_back((((-image.Nfft+(image.Nfft%2))/2)+i)/sampleRange);
	};
	
	// Define the inverse Radon Transform frequency ramp filter
	for ( unsigned int i=0; i< tx.freq.size(); i++ )
	{
		radonFilter.push_back(fabs(tx.freq[i])); 
		// Assumes the frequencies are the same for every pulse
	};
	
	// Determine size of desired zero-padding
	int L = image.Nfft - rx.freqbins;
	
	int zeros = (L+L%2)/2;
	std::vector<phdata> zeroPad(zeros);
	//vector<phdata>::iterator it;
	
	// Determine the Amplitude scaling factor based on zero padding
	double ampScaleFactor = image.Nfft/rx.freqbins/fcenter;
	
	// iterate through pulses to backproject across the aperture
	for ( int i=0; i<pulses; i++ )
	{
		x = (cos(tx.az[i]*M_PI/180)*cos(tx.grazing*M_PI/180) 
		     + cos(rx.az[i]*M_PI/180)*cos(rx.grazing*M_PI/180) )/2;
		y = (sin(tx.az[i]*M_PI/180)*cos(tx.grazing*M_PI/180) 
		     + sin(rx.az[i]*M_PI/180)*cos(rx.grazing*M_PI/180) )/2;
		z = (sin(tx.grazing*M_PI/180) + sin(rx.grazing*M_PI/180) )/2;
		// Need to figure out how to handle the Bistatic look angle
		for ( unsigned int n=0; i< tx.freq.size(); i++ )
		{
			fx.push_back(tx.freq[n] * x);
			fy.push_back(tx.freq[n] * y);
			fz.push_back(tx.freq[n] * z);
			//elBiLook.pushback(atan(fz(j),sqrt(fx(j)^2,fy(j)^2)));
		}
		
		for ( int m=0; m<radonFilter.size(); m++ )
		{
			filteredData.push_back(radonFilter[m]*rx.rxphase[m][i]);
		}
		
		//zero pad the data for the FFT
		//it=filteredData.begin();
		filteredData.insert( filteredData.begin(), zeroPad.begin(), zeroPad.end()); 
		filteredData.insert( filteredData.end(), zeroPad.begin(), zeroPad.end());
		
		// Recast the vector to the type needed for FFTW and compute FFT
		in = reinterpret_cast<fftw_complex*>(&filteredData);
		void fftw_execute(const fftw_plan plan);
		//std::cout << out << std::endl;
		
		// Interpolation step next...
	}
	
	fftw_destroy_plan(plan);
	//fftw_free(in);
	//fftw_free(out);
}

int main()
{
	int imgdim = 2;
	
	rxdata rx,tx;
	
	// Define the azimuth vectors
	for (int i = 0; i<((rx.maxaz-rx.minaz)/rx.azstep); i++)
	{
		rx.az.push_back(rx.minaz + (i*rx.azstep));
		tx.az.push_back(tx.minaz + (i*tx.azstep));
		
	}
	
	rx.printAz();
	tx.printAz();
	
	
	// Define the frequency vector
	for  (int i = 0; i<tx.freqbins; i++)
	{
		tx.freq.push_back(tx.freqmin + (((tx.freqmax-tx.freqmin)/tx.freqbins)*i));
	}
	
	imgdata image;
	
	// Define phase history data for single point scatterer at the origin
	// rx.rxphase(size_type 512, vector<complex<double> > vec(3001,1);
	
	backproject(rx, tx, image);
		
	// Define the file to write the image data to
	// Note that the file will be overwritten through this function
	//ofstream imgdata ("stuff1.csv", ios::trunc | ios::out);
	//if (imgdata.isopen())
	//{
	// Write a header to the columns of the csv
	//imgdata << "'Index', 'Time', 'Frequency', 'Pulse'" << endl;
	//}
	//else cout << "Unable to open file";

	// Define the phase history data file to be read
	// This is a matlab file in this case, and so requires the inclusion
	// of the matlab mat.h header
	//ifstream phdata ("Sentra_el30.0000.mat",ios::ate | ios::in);
	//if (phdata.isopen())
	//{
		
	

	return 0;
}

