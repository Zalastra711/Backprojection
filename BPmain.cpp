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
#include <cstring>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <numeric>
#include <algorithm>
#include "parameters.h"

// Compile with:
// g++ -I .. -fopenmp <filename.cc> <path to library> -lfftw2 -lfftw3_omp

//using namespace std;
//using namespace fftwpp;

// This program is designed to take input phase history data from a file 
// and use "fast backprojection" to create an image of a scene.


//using namespace std;
//using namespace fftw;

/*
void swap(std::vector<phdata> *var1, std::vector<phdata> *var2){
	std::vector<phdata> tmp = *var1;
	*var1 = *var2;
	std::cout << &var1 << ", " << &var2 << ", " << &tmp << std::endl;
	*var2 = tmp;
};
*/

// Create FFTSHIFT and IFFTSHIFT functions
void fftshift(std::vector<phdata> &argin, std::size_t count){
	int k = 0;
	int c = (int) floor((float)count/2);
	
	if (count % 2 == 0) //define case for even length vector
	{
		for (k=0; k < c; k++)
		{	
			phdata tmp = argin[k];
			argin[k] = argin[k+c];
			argin[k+c] = tmp;
			// swap(&argin[k], &argin[k+c]);
		}
	}
	else
	{
		phdata tmp = argin[0];
		for (k=0; k<c;k++)
		{
			argin[k] = argin[c+k+1];
			argin[c+k+1] = argin[k+1];
		}
		argin[c]=tmp;
	}
};

void ifftshift(std::vector<phdata> &argin, std::size_t count){
	int k = 0;
	int c = (int)floor((float)count/2);
	
	//std::cout << c << std::endl;
	//std::cout << argin[0] << std::endl;
	//std::cout << argin->k << std::endl;
	
	if (count % 2 == 0)
	{
		std::cout << "IFFT Case Even." << std::endl;
		for (k=0; k<c;k++)
		{
			//std::cout << &argin << std::endl;
			//std::cout << "Swapping " << k << " Index." << std::endl;
			//std::cout << argin[k] << std::endl;
			//std::cout << argin[k+c] << std::endl;
			
			phdata tmp = argin[k];
			//std::cout << tmp << std::endl;
			
			argin[k] = argin[k+c];
			argin[k+c] = tmp;
			//swap(&argin[k], &argin[k+c]);
		}
	}
	else
	{
		std::cout << "IFFT Case Odd." << std::endl;
		phdata tmp = argin[count-1];
		for (k=c-1; k>=0; k--)
		{
			//std::cout << "Swapping " << k << " Index." << std::endl;
			argin[c+k+1] = argin[k];
			argin[k] = argin[c+k];
		}
		argin[c] = tmp;
	}
};

double dR(double txazim, double txgraz, double rxazim, double rxgraz,
	  double x, double y, double z){
	
	//std::cout << "Calling dR Function" << std::endl;
	
	return x*(cos(txgraz)*cos(txazim)+cos(rxgraz)*cos(rxazim)) + 
	y*(cos(txgraz)*sin(txazim)+cos(rxgraz)*sin(rxazim)) + 
	z*(sin(txgraz)+sin(rxgraz));
	
};

void backproject(rxdata& rx, rxdata& tx, imgdata& image){
	
	std::cout << "Running BackProjection" << std::endl;
	
	std::vector<phdata> filteredData;
	
	fftw_complex *in, *out;
	fftw_plan plan;
	
	double deltaFreq, diffFreq, sampleRange;
	phdata phaseCorr, interp;
	
	// these must match Nx, Ny, Nz...or find a better way to define this part...
	constexpr size_t xdim = 128, ydim=128, zdim = 1;
	
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image.Nfft );
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image.Nfft );
	plan = fftw_plan_dft_1d(image.Nfft, in, out, FFTW_BACKWARD, FFTW_MEASURE);
	
	std::cout << "FFTW Plan Created" << std::endl;
	
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
		image.r_vec.push_back((((-image.Nfft+(image.Nfft%2))/2)+i)*sampleRange);
	};
	
	// variables for min and max value of range vector for use in later comparison
	double rvecMin = *std::min_element(image.r_vec.begin(), image.r_vec.end());
	double rvecMax = *std::max_element(image.r_vec.begin(), image.r_vec.end());
	
	std::cout << "Image Range Vector Created" << std::endl;
	
	// Define the inverse Radon Transform frequency ramp filter
	for ( unsigned int i=0; i< tx.freq.size(); i++ )
	{
		radonFilter.push_back(fabs(tx.freq[i])); 
		// Assumes the frequencies are the same for every pulse
	};
	
	std::cout << "Radon Filter Created" << std::endl;
	
	// Determine size of desired zero-padding
	int L = image.Nfft - rx.freqbins;
	
	int zeros = (L+L%2)/2;
	std::vector<phdata> zeroPad(zeros);
	//vector<phdata>::iterator it;
	std::vector<phdata> outData;
	
	// Determine the Amplitude scaling factor based on zero padding
	double ampScaleFactor = image.Nfft/rx.freqbins/fcenter;
	
	// iterate through pulses to backproject across the aperture
	for ( int i=0; i<=pulses; i++ )
	{
		// Reset vectors for each pulse
		fx.clear();
		fy.clear();
		fz.clear();
		filteredData.clear();
		outData.clear();
		
		// Define Bistatic collection manifold/geometry
		/*
		x = (cos(tx.az[i]*M_PI/180)*cos(tx.grazing*M_PI/180) 
		     + cos(rx.az[i]*M_PI/180)*cos(rx.grazing*M_PI/180) )/2;
		y = (sin(tx.az[i]*M_PI/180)*cos(tx.grazing*M_PI/180) 
		     + sin(rx.az[i]*M_PI/180)*cos(rx.grazing*M_PI/180) )/2;
		z = (sin(tx.grazing*M_PI/180) + sin(rx.grazing*M_PI/180) )/2;
		*/
		x = (cos(tx.az[i])*cos(tx.grazing) 
		     + cos(rx.az[i])*cos(rx.grazing) )/2;
		y = (sin(tx.az[i])*cos(tx.grazing) 
		     + sin(rx.az[i])*cos(rx.grazing) )/2;
		z = (sin(tx.grazing) + sin(rx.grazing) )/2;
		// Need to figure out how to handle the Bistatic look angle
		for ( unsigned int jj=0; jj< tx.freq.size(); jj++ )
		{
			fx.push_back(tx.freq[jj] * x);
			fy.push_back(tx.freq[jj] * y);
			fz.push_back(tx.freq[jj] * z);
			//elBiLook.pushback(atan(fz(j),sqrt(fx(j)^2,fy(j)^2)));
		}
		
		for ( int ii=0; ii<radonFilter.size(); ii++ )
		{
			filteredData.push_back(radonFilter[ii]*rx.rxphase[ii][i]);
		}
		
		//zero pad the data for the FFT
		//it=filteredData.begin();
		filteredData.insert( filteredData.begin(), zeroPad.begin(), zeroPad.end()); 
		filteredData.insert( filteredData.end(), zeroPad.begin(), zeroPad.end());
		
		std::cout << "Pulse " << i << " Data Zero Padded." << std::endl;
		
		/*
		for (int ii=0; ii<filteredData.size(); ii++)
		{
			std::cout << filteredData[ii] ;
		}
		std::cout << std::endl;
		*/
		
		//std::cout << &filteredData << std::endl;
		//std::cout << &filteredData.front() << std::endl;
		
		// Shift the data to put the center frequency at DC position (zero index)
		ifftshift(filteredData, filteredData.size());
		
		std::cout << "Pulse " << i << " IFFT Shift Complete." << std::endl;
		
		/*
		for (int ii=0; ii<filteredData.size(); ii++)
		{
			std::cout << filteredData[ii] ;
		}
		std::cout << std::endl;
		*/
		
		// Recast the vector to the type needed for FFTW and compute FFT
		in = reinterpret_cast<fftw_complex*>(&filteredData[0]);
		
		
		for (int ii=0; ii<filteredData.size(); ii++)
		{
			std::cout << "(" << (in[ii])[0] << "," << (in[ii])[1] << ");";
		}
		std::cout << std::endl;
		
		
		std::cout << "Pulse " << i << " Data Prepared for FFT." << std::endl;
		
		fftw_execute(plan);
		
		std::cout << "Pulse " << i << " FFT Executed" << std::endl;
		
		std::memcpy(&outData, &out, sizeof(filteredData) );
		
		// Shift the output to put zero range at middle of range profile
		fftshift(outData, outData.size());
		
		
		for (int ii=0; ii<filteredData.size(); ii++)
		{
			std::cout << "(" << (out[ii])[0] << "," << (out[ii])[1] << ");";
			//std::cout << outData[ii];
		}
		std::cout << std::endl;
		
		
		// Calculate differential range for each pixel in the image
		// for x, for y, for z...
		
		double diffRange[xdim][ydim][zdim];
		int cnt=0;
		
		for (size_t ii=0; ii != xdim; ii++)
		{
			for (size_t jj=0; jj != ydim; jj++)
			{
				for (size_t kk=0; kk != zdim; kk++)
				{
					diffRange[ii][jj][kk] = dR(tx.az[i], tx.grazing, 
						       rx.az[i], rx.grazing,
						       image.x[ii], image.y[jj],
						       image.z[kk]);
					// std::cout << diffRange[ii][jj][kk] << ", ";
					// if diffRange[ii][jj][kk] < image.r_vec[0] (extrapolate)
					// elseif diffRange[ii][jj][kk] <= image.r_vec.back() (interpolate)
					// Do i need to check for which direction the interpolation/extrapolation goes?
					// Case where rvec has finer spacing or case where image pixels have finer spacing...
					phaseCorr = std::exp((-I)*twoPI*(fcenter/c)*diffRange[ii][jj][kk]);
					if ( (diffRange[ii][jj][kk] >= rvecMin) && (diffRange[ii][jj][kk] <= rvecMax) )
					{
						// interpolate image value here
						// linear interpolation is used for now...
						// take indices of "out" near current diffRange and interpolate a value
						// find out which value of "out" are near current diffRange...
						while ( image.r_vec[cnt] < diffRange[ii][jj][kk] )
						{
							cnt++;
						}
						// Lagrange approximation of degree <=1
						interp = outData[cnt]*((diffRange[ii][jj][kk]-image.r_vec[cnt+1])/(image.r_vec[cnt]-image.r_vec[cnt+1])) + outData[cnt+1]*((diffRange[ii][jj][kk]-image.r_vec[cnt])/(image.r_vec[cnt+1]-image.r_vec[cnt]));
						//interp = out[cnt]*((diffRange[ii][jj][kk]-image.r_vec[cnt+1])/(image.r_vec[cnt]-image.r_vec[cnt+1])) + out[cnt+1]*((diffRange[ii][jj][kk]-image.r_vec[cnt])/(image.r_vec[cnt+1]-image.r_vec[cnt]));
						image.img_final[ii][jj][kk] = image.img_final[ii][jj][kk] + interp;
					}
				}
				
			}
			//std::cout << std::endl;
		}		
	}
	
	//fftw_cleanup();
	fftw_destroy_plan(plan);
	//fftw_free(in);
	//fftw_free(out);
};

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
	std::ofstream imgfile ("image.csv", std::ios::trunc | std::ios::out);
	if (imgfile.is_open())
	{
	// Write a header to the columns of the csv
	// imgdata << "'Index', 'Time', 'Frequency', 'Pulse'" << endl;
		for (int i = 0; i < image.Nx; i++)
		{
			for (int j = 0; j < image.Ny; j++)
			{
				for (int k=0; k<image.Nz; k++)
				{
					imgfile << image.img_final[i][j][k] << ";";
				}
				imgfile << std::endl;
			}
			imgfile << ":";
		}
	}
	else std::cout << "Unable to open file";

	return 0;
};

