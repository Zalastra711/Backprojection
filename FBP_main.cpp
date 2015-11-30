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
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <numeric>
#include <algorithm>
#include "fac_params.h"


// Compile with:
// g++ -I .. -fopenmp <filename.cc> <path to library> -lfftw2 -lfftw3_omp

/*
 KEEP AZIMUH POSITIONS (PULSES) AND FREQUENCY BINS IN POWERS OF 2! THE RECURSION
 FUNCTION IS NOT DESIGNED TO HANDLE ODD NUMBERS AT THE MOMENT.
 */

// Create FFTSHIFT and IFFTSHIFT functions
void fftshift(std::vector<phdata> &argin, std::size_t count){
	int k = 0;
	int plc = (int) floor((float)count/2);
	
	if (count % 2 == 0) //define case for even length vector
	{
		for (k=0; k<plc; k++)
		{	
			phdata tmp = argin[k];
			argin[k] = argin[k+plc];
			argin[k+plc] = tmp;
		}
	}
	else
	{
		phdata tmp = argin[0];
		for (k=0; k<plc;k++)
		{
			argin[k] = argin[plc+k+1];
			argin[plc+k+1] = argin[k+1];
		}
		argin[plc]=tmp;
	}
};

void ifftshift(std::vector<phdata> &argin, std::size_t count){
	int k = 0;
	int plc = (int)floor((float)count/2);
	
	if (count % 2 == 0)
	{
		//std::cout << "IFFT Case Even." << std::endl;
		for (k=0; k<plc;k++)
		{
			phdata tmp = argin[k];
			argin[k] = argin[k+plc];
			argin[k+plc] = tmp;
		}
	}
	else
	{
		//std::cout << "IFFT Case Odd." << std::endl;
		phdata tmp = argin[count-1];
		for (k=plc-1; k>=0; k--)
		{
			argin[plc+k+1] = argin[k];
			argin[k] = argin[plc+k];
		}
		argin[plc] = tmp;
	}
};

double dR(double txazim, double txgraz, double rxazim, double rxgraz,
	  double x, double y, double z){
	
	//std::cout << "Calling dR Function" << std::endl;
	
	return (x*((cos(txgraz)*cos(txazim))+(cos(rxgraz)*cos(rxazim)))) + 
	(y*((cos(txgraz)*sin(txazim))+(cos(rxgraz)*sin(rxazim)))) + 
	(z*(sin(txgraz)+sin(rxgraz)));
	
};

void rvec_Create(rxdata& rx, imgdata& image)
{
	double deltaFreq, diffFreq, sampleRange;
	for ( int i=0; i<rx.freqbins-1; i++ )
	{
		diffFreq = rx.freq[i+1] - rx.freq[i];
		deltaFreq += diffFreq;
	};
	
	// Assumes uniformly sampled frequency vector
	deltaFreq /= (rx.freqbins-1); 
	
	// Range sample spacing corresponding to collected frequency samples
	sampleRange = C/(image.Nfft*deltaFreq); 
	
	// Define the range at each sample
	for ( int i = 0; i <image.Nfft; i++ )
	{
		image.r_vec.push_back((((-image.Nfft+(image.Nfft%2))/2)+i)*sampleRange);
	}
	//std::cout << image.r_vec.size() << std::endl;
	
};


void factorize(rxdata& rx, rxdata& tx, imgdata& image, int numSteps, int currStep)
{
	if ( currStep == numSteps || image.Nx==1 || image.Ny==1 )
	{
		// Create x and y distance vectors
		image.makeX();
		image.makeY();
		
		//image.printX();
		//image.printY();
		
		// Call backproject function for each quadrant
		backproject(rx, tx, image);
	}
	else
	{
		// Flag to keep track of how deep into recursion we've gone
		currStep++;
		
		std::cout << "Factorization Level " << currStep << std::endl;
		
		imgdata imSub;
		
		rxdata rxtmp, txtmp;
		
		// Populate new rxdata members
		rxtmp.azpos = rx.azpos/2;
		rxtmp.freqbins = rx.freqbins;
		rxtmp.minaz = rx.minaz;
		rxtmp.maxaz = rx.maxaz;
		rxtmp.cFreq = rx.cFreq;
		rxtmp.bWidth = rx.bWidth;
		rxtmp.grazing = rx.grazing;
		rxtmp.range = rx.range;
		
		txtmp.azpos = tx.azpos/2;
		txtmp.minaz = tx.minaz;
		txtmp.maxaz = tx.maxaz;
		txtmp.grazing = tx.grazing;
		txtmp.range = tx.range;
		
		
		for (int i = 0; i<rxtmp.azpos; i++)
		{
			txtmp.az.push_back(txtmp.minaz + ((txtmp.maxaz-txtmp.minaz)/(txtmp.azpos-1) *i));
			rxtmp.az.push_back(rxtmp.minaz + ((rxtmp.maxaz-rxtmp.minaz)/(rxtmp.azpos-1) *i));
		}
		
		for  (int i = 0; i<rxtmp.freqbins; i++)
		{
			rxtmp.freq.push_back((rxtmp.cFreq-rxtmp.bWidth/2) + (rxtmp.bWidth/(rxtmp.freqbins-1))*i);
		}
		
		rxtmp.rxphase.assign(rxtmp.freqbins, std::vector<phdata> (rxtmp.azpos, {0,0}) );
		 
		// Break image/subimage into quadrants
		imSub.x_extent = image.x_extent/2;
		imSub.y_extent = image.y_extent/2;
		imSub.Nx = image.Nx/2;
		imSub.Ny = image.Ny/2;
		
		// Add pairs of phase history pulses, i.e. 0+1, 2+3, etc.
		for ( int ii=0; ii<rxtmp.freqbins; ii++ )
		{
			for ( int jj=0; jj<rxtmp.azpos; jj++ )
			{
				int idx = 2*jj;
				rxtmp.rxphase[ii][jj] = rx.rxphase[ii][idx] + rx.rxphase[ii][idx+1];
				//std::cout << rxtmp.rxphase[ii][jj] << ", ";
			}
			//std::cout << std::endl;
		}
		
		std::cout << "Phase History Decimated." << std::endl;
		
		// Create Iterators for beginning, midpoint, and end of image vectors
		auto xvecBegin = image.img_final.begin(); 
		auto xvecEnd = image.img_final.end(); 
		auto xvecMid = std::distance(xvecBegin, xvecEnd)/2; 

		// Assumed factorization parameter of 2
		// Iterate over quadrants - for loop
		for ( int ii=0 ; ii<4; ii++ )
		{
			// Assign quadrant center based on quadrant being evaluated
			switch (ii)
			{
				case 0: // Quadrant 1
					
					std::cout << "Quadrant 1 being processed." << std::endl;
					imSub.x0 = image.x0 + imSub.x_extent/2;
					imSub.y0 = image.y0 + imSub.y_extent/2;
					imSub.x.clear();
					imSub.y.clear();
					imSub.r_vec.clear();
					// Create subimage grid
					imSub.img_final.assign(imSub.Nx, std::vector<phdata> (imSub.Ny, {0,0}) );
					
					factorize(rxtmp, txtmp, imSub, numSteps, currStep);
					
					std::cout << "Sub-image for quadrant 1 created.\n";
					// Place imSub from factorize into image.imFinal at top right
					// First half of rows, second half of columns
					for ( int jj= 0; jj< xvecMid; jj++ )
					{
						auto yvecBegin = image.img_final[jj+xvecMid].begin(); 
						auto yvecEnd = image.img_final[jj+xvecMid].end(); 
						auto yvecMid = std::distance(yvecBegin, yvecEnd)/2; 
						std::copy(imSub.img_final[jj].begin(), imSub.img_final[jj].end(), yvecBegin + yvecMid);
					}
					break;
				case 1: // Quadrant 2
					std::cout << "Quadrant 2 being processed." << std::endl;
					imSub.x0 = image.x0 - imSub.x_extent/2;
					imSub.y0 = image.y0 + imSub.y_extent/2;
					imSub.x.clear();
					imSub.y.clear();
					imSub.r_vec.clear();
					// Create subimage grid
					imSub.img_final.assign(imSub.Nx, std::vector<phdata> (imSub.Ny, {0,0}) );
					
					factorize(rxtmp, txtmp, imSub, numSteps, currStep);
					
					std::cout << "Sub-image for quadrant 2 created.\n";
					// Place imSub from factorize into image.imFinal at top left
					for ( int jj= 0; jj< xvecMid; jj++ )
					{
						auto yvecBegin = image.img_final[jj].begin(); 
						auto yvecEnd = image.img_final[jj].end(); 
						auto yvecMid = std::distance(yvecBegin, yvecEnd)/2; 
						std::copy(imSub.img_final[jj].begin(), imSub.img_final[jj].end(), yvecBegin + yvecMid);
					}
					break;
				case 2: // Quadrant 3
					std::cout << "Quadrant 3 being processed." << std::endl;
					imSub.x0 = image.x0 + imSub.x_extent/2;
					imSub.y0 = image.y0 - imSub.y_extent/2;
					
					imSub.x.clear();
					imSub.y.clear();
					imSub.r_vec.clear();
					// Create subimage grid
					imSub.img_final.assign(imSub.Nx, std::vector<phdata> (imSub.Ny, {0,0}) );
					
					factorize(rxtmp, txtmp, imSub, numSteps, currStep);
					
					std::cout << "Sub-image for quadrant 3 created.\n";
					
					// Place imSub from factorize into image.imFinal at bottom right
					for ( int jj= 0; jj<xvecMid; jj++ )
					{
						auto yvecBegin = image.img_final[jj+xvecMid].begin();
						std::copy(imSub.img_final[jj].begin(), imSub.img_final[jj].end(), yvecBegin);
					}
					break;
				case 3: // Quadrant 4
					std::cout << "Quadrant 4 being processed." << std::endl;
					imSub.x0 = image.x0 - imSub.x_extent/2;
					imSub.y0 = image.y0 - imSub.y_extent/2;
					imSub.x.clear();
					imSub.y.clear();
					imSub.r_vec.clear();
					// Create subimage grid
					imSub.img_final.assign(imSub.Nx, std::vector<phdata> (imSub.Ny, {0,0}) );
					
					factorize(rxtmp, txtmp, imSub, numSteps, currStep);
					
					std::cout << "Sub-image for quadrant 4 created.\n";
					// Place imSub from factorize into image.imFinal at bottom left
					for ( int jj= 0; jj<xvecMid; jj++ )
					{
						auto yvecBegin = image.img_final[jj].begin(); 
						std::copy(imSub.img_final[jj].begin(), imSub.img_final[jj].end(), yvecBegin);
					}
					break;
				default:
					break;
			}
		}
	}
}

//Ingest and parse data from system parameter file
void sysParam(rxdata& rx, rxdata& tx){
	// Open parameter file
	std::ifstream paramFile ("./parameters.csv", std::ios::in);
	if (paramFile.is_open())
	{
		int i=0;
		while (!paramFile.eof())
		{
			std::string line;
			std::getline( paramFile, line );
			
			std::stringstream ss( line);
			std::string parse;
			int j=0;
			switch (i)
			{
				case 0:
					while ( std::getline( ss, parse, ',') )
					{
						std::stringstream val(parse);
						bool bistatic;
						switch (j)
						{
							case 0:
								val >> rx.cFreq;
								break;
							case 1:
								val >> rx.bWidth;
							case 2:
								val >> rx.freqbins;
								break;
							case 3:
								val >> bistatic;
								break;
							default:
								break;
						}
						j++;
					}
					break;
				case 1:
					while ( std::getline( ss, parse, ',') )
					{
						std::stringstream val(parse);
						switch (j)
						{
							case 0:
								val >> tx.minaz;
								break;
							case 1:
								val >> tx.maxaz;
							case 2:
								val >> tx.azpos;
								break;
							case 3:
								val >> tx.grazing;
								break;
							default:
								break;
						}
						j++;
					}
					
				case 2:
					while ( std::getline( ss, parse, ',') )
					{
						std::stringstream val(parse);
						switch (j)
						{
							case 0:
								val >> rx.minaz;
								break;
							case 1:
								val >> rx.maxaz;
							case 2:
								val >> rx.azpos;
								break;
							case 3:
								val >> rx.grazing;
								break;
							default:
								break;
						}
						j++;
					}
					
				default:
					break;
			}
			std::cout << "Parameter file line " << i+1 << " parsed." << std::endl;
			i++;
		}
	}
	else std::cout << "Unable to open file";
	paramFile.close();
}

void phsDataInput(rxdata& rx)
{
	// Open parameter file
	std::ifstream dataFile("phase_history.csv", std::ios::in);
	if (dataFile.is_open())
	{
		int ii=0;
		while (!dataFile.eof())
		{
			std::string line;
			std::getline( dataFile, line );
			
			std::stringstream ss( line);
			std::string parse;
			
			int jj=0;
			bool toggle = 1;
			double temp;
			while ( std::getline(ss, parse, ',') )
			{
				if (toggle)
				{
					std::stringstream val(parse);
					val >> temp;
					rx.rxphase[ii][jj].real(temp);
					//std::cout << real(rx.rxphase[ii][jj]) << ",";
					toggle = !toggle;
				}
				else
				{
					std::stringstream val(parse);
					val >> temp;
					rx.rxphase[ii][jj].imag(temp);
					//std::cout << imag(rx.rxphase[ii][jj]) << ",";
					toggle = !toggle;
					jj++;
				}
			}
			//std::cout << std::endl;
			//std::cout << "Line " << ii << " of phase history parsed." << std::endl;
			ii++;
		}
		
	}
	else std::cout << "Unable to open file";
	dataFile.close();
}


void backproject(rxdata& rx, rxdata& tx, imgdata& image){
	
	std::cout << "Running BackProjection" << std::endl;
	
	std::vector<phdata> filteredData;
	std::vector<phdata> outData;
	
	fftw_complex *in, *rawOut;
	fftw_plan plan;
	
	phdata phaseCorr, interp;
	
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image.Nfft );
	rawOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image.Nfft );
	plan = fftw_plan_dft_1d(image.Nfft, in, rawOut, FFTW_BACKWARD, FFTW_MEASURE);
	
	std::cout << "FFTW Plan Created" << std::endl;
	
	// Create image range vector
	rvec_Create(rx,image);
	
	std::vector<double> radonFilter;
	
	int pulses = rx.azpos;
	//std::cout << pulses << std::endl;
	
	// variables for min and max value of range vector for use in later comparison
	double rvecMin = *std::min_element(image.r_vec.begin(), image.r_vec.end());
	double rvecMax = *std::max_element(image.r_vec.begin(), image.r_vec.end());
	
	std::cout << "Image Range Vector Created" << std::endl;
	
	// Define the inverse Radon Transform frequency ramp filter
	for ( int i=0; i<rx.freqbins; i++ )
	{
		// Assumes the frequencies are the same for every pulse
		radonFilter.push_back(fabs(rx.freq[i])); 
		//std::cout << radonFilter[i] << ", ";
	}
	//std::cout << std::endl;
	
	std::cout << "Radon Filter Created" << std::endl;
	
	// Determine size of desired zero-padding
	int L = image.Nfft - rx.freqbins;
	
	int zeros = (L+L%2)/2;
	std::vector<phdata> zeroPad(zeros);
	
	// Determine the Amplitude scaling factor based on zero padding
	//std::cout << rx.freqbins << ", " << rx.cFreq << std::endl;
	
	double ampScaleFactor = (1.0/rx.freqbins/rx.cFreq);
	//std::cout.precision(std::numeric_limits< double >::max_digits10);
	//std::cout << ampScaleFactor << std::endl;
	
	// Calculate distance to subimage center
	double subImgCenter = sqrt((image.x0*image.x0)+(image.y0*image.y0));
	
	// Find the maximum distance from center for image
	//double distanceX = image.x0-image.x[0];
	//double distanceY = image.y0-image.y[0];
	//double maxDis = sqrt((distanceX*distanceX)+(distanceY*distanceY));
	//double extension = maxDis/2;
	
	// Shift the range profile to subimage center
	phdata rpshift = std::exp(J*twoPI*rx.cFreq/C*subImgCenter);
	
	//std::cout << rpshift << std::endl;
	
	//auto minBound = std::lower_bound(image.r_vec.begin(), image.r_vec.end(), -(maxDis*3));
	//auto maxBound = std::lower_bound(image.r_vec.begin(), image.r_vec.end(), maxDis*3);
	
	//std::cout << *minBound << ", " << *maxBound << std::endl;
	// iterate through pulses to backproject across the aperture
	for ( int i=0; i<pulses; i++ )
	{
		// Reset vectors for each pulse
		//fx.clear();
		//fy.clear();
		//fz.clear();
		filteredData.clear();
		outData.clear();
		
		// Define Bistatic collection manifold/geometry
		/*
		 x = (cos(tx.az[i]*M_PI/180)*cos(tx.grazing*M_PI/180) 
		 + cos(rx.az[i]*M_PI/180)*cos(rx.grazing*M_PI/180) )/2;
		 y = (sin(tx.az[i]*M_PI/180)*cos(tx.grazing*M_PI/180) 
		 + sin(rx.az[i]*M_PI/180)*cos(rx.grazing*M_PI/180) )/2;
		 z = (sin(tx.grazing*M_PI/180) + sin(rx.grazing*M_PI/180) )/2;
		 
		 x = (cos(tx.az[i])*cos(tx.grazing) 
		 + cos(rx.az[i])*cos(rx.grazing) )/2;
		 y = (sin(tx.az[i])*cos(tx.grazing) 
		 + sin(rx.az[i])*cos(rx.grazing) )/2;
		 z = (sin(tx.grazing) + sin(rx.grazing) )/2;
		 // Need to figure out how to handle the Bistatic look angle
		 for ( int jj=0; jj<=tx.freqbins; jj++ )
		 {
		 fx.push_back(tx.freq[jj] * x);
		 fy.push_back(tx.freq[jj] * y);
		 fz.push_back(tx.freq[jj] * z);
		 //elBiLook.pushback(atan(fz(j),sqrt(fx(j)^2,fy(j)^2)));
		 }
		 */
		
		for ( int ii=0; ii<rx.freqbins; ii++ )
		{
			filteredData.push_back(radonFilter[ii]*rx.rxphase[ii][i]);
			//std::cout << filteredData[ii] << ", ";
		}
		//std::cout << std::endl;
		
		//zero pad the data for the FFT
		//it=filteredData.begin();
		filteredData.insert( filteredData.begin(), zeroPad.begin(), zeroPad.end()); 
		filteredData.insert( filteredData.end(), zeroPad.begin(), zeroPad.end());
		
		// Shift the data to put the center frequency at DC position (zero index)
		ifftshift(filteredData, filteredData.size());
		
		// Recast the vector to the type needed for FFTW and compute FFT
		
		for (size_t ii=0; ii<filteredData.size(); ii++)
		{
			in[ii][0] = real(filteredData[ii]);
			in[ii][1] = imag(filteredData[ii]);
			//std::cout << "(" << (in[ii])[0] << "," << (in[ii])[1] << ");";
		}
		//std::cout << std::endl;
		
		//std::cout << "Pulse " << i << " Data Prepared for FFT." << std::endl;
		
		fftw_execute(plan);
		
		if ( i%100 == 0 )
		{
			std::cout << "Pulse " << i << " FFT Executed" << std::endl;
		}
		
		// Scale for ifft implementation of FFTW
		for (int ii=0; ii<image.Nfft; ii++)
		{
			outData.push_back(phdata(rawOut[ii][0],rawOut[ii][1]));
			outData[ii] *= ampScaleFactor * rpshift;
			//std::cout << outData[ii] << ";";
		}
		//std::cout << std::endl;
		
		// Shift the output to put zero range at middle of range profile
		fftshift(outData, outData.size());
		/*
		for (int ii=0; ii<image.Nfft; ii++)
		{
			outData[ii] *= rpshift;
			//std::cout << outData[ii] << ";";
		}
		//std::cout << std::endl;
		 */
		
		// Calculate differential range for each pixel in the image
		// for x, for y, for z...
		
		double diffRange = 0;
		int cnt;
		
		for (int ii=0; ii < image.Nx; ii++)
		{
			for (int jj=0; jj < image.Ny; jj++)
			{
				
				diffRange = dR(tx.az[i], tx.grazing, 
					       rx.az[i], rx.grazing,
					       image.x[ii], image.y[jj],
					       image.z0);
				
				//std::cout << diffRange << ", ";
				
				phaseCorr = std::exp(J*twoPI*(rx.cFreq/C)*diffRange);
				//std::cout << phaseCorr << std::endl;
				
				if ( (diffRange >= rvecMin) && (diffRange <= rvecMax) )
				{
					// interpolate image value here
					// linear interpolation is used for now...
					// take indices of "out" near current diffRange and interpolate a value
					// find out which value of "out" are near current diffRange...
					auto idx = std::lower_bound(image.r_vec.begin(), image.r_vec.end(), diffRange);
					//std::cout << *idx << std::endl;
					cnt = std::distance(image.r_vec.begin(), idx);
					//std::cout << cnt << ", ";
					// Lagrange approximation of degree <=1
					// interp = outData[cnt-1]*((diffRange-image.r_vec[cnt])/(image.r_vec[cnt-1]-image.r_vec[cnt])) 
					//       + outData[cnt]*((diffRange-image.r_vec[cnt-1])/(image.r_vec[cnt]-image.r_vec[cnt-1]));
					interp = outData[cnt-1]+ ((outData[cnt]-outData[cnt-1])*((diffRange-image.r_vec[cnt-1])/(image.r_vec[cnt]-image.r_vec[cnt-1])));
					//std::cout << interp << ", ";
					image.img_final[ii][jj] = image.img_final[ii][jj] + (interp * (phaseCorr/(double) pulses));
				}
			}
			//std::cout << std::endl;
			//std::cout << "Image column " << ii << " completed.\n";
		}
		//std::cout << "Image for pulse " << i << " completed.\n";
	}
	
	//fftw_destroy_plan(plan);
	//fftw_free(in);
	//fftw_free(rawOut);
	//fftw_cleanup();
};

int main()
{
	int imgdim = 2;
	
	rxdata rx,tx;
	
	sysParam(rx,tx);
	
	rx.rxphase.assign(rx.freqbins, std::vector<phdata> (rx.azpos, {0,0}) );
	
	phsDataInput(rx);
	std::cout << "Phase History Loaded.\n";
	
	// Define the azimuth vectors
	for (int i = 0; i<rx.azpos; i++)
	{
		rx.az.push_back(rx.minaz + ((rx.maxaz-rx.minaz)/(rx.azpos-1) *i));
		tx.az.push_back(tx.minaz + ((tx.maxaz-tx.minaz)/(tx.azpos-1) *i));
	}
	
	//rx.printAz();
	//tx.printAz();
	
	// Define the frequency vector
	for  (int i = 0; i<rx.freqbins; i++)
	{
		//tx.freq.push_back((rx.cFreq-rx.bWidth/2) + (rx.bWidth/(rx.freqbins-1))*i);
		rx.freq.push_back((rx.cFreq-rx.bWidth/2) + (rx.bWidth/(rx.freqbins-1))*i);
		//std::cout << tx.freq[i] << ",";
	}
	//std::cout << std::endl;
	
	imgdata imFinal;
	
	imFinal.img_final.assign(imFinal.Nx, std::vector<phdata> (imFinal.Ny, {0,0}) );
	//imFinal.makeX();
	//imFinal.makeY();
	
	int fac_steps = 3;
	int current_step = 0;
	
	// Call factorize recursive function	
	factorize(rx, tx, imFinal, fac_steps, current_step);

	fftw_destroy_plan(plan);
	fftw_cleanup();
	// Define the file to write the image data to
	// Note that the file will be overwritten through this function
	std::ofstream imgfile ("ffbp_image.csv", std::ios::trunc | std::ios::out);
	if (imgfile.is_open())
	{
		for (int i = 0; i < imFinal.Ny; i++)
		{
			for (int j = 0; j < imFinal.Nx; j++)
			{
				imgfile << real(imFinal.img_final[i][j]) << "," << imag(imFinal.img_final[i][j]) << ",";
			}
			imgfile << std::endl;
		}
	}
	else std::cout << "Unable to open file";
	imgfile.close();
	
	return 0;
};

