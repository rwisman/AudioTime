package edu.ius.fft;

public class FFT {

	public double[] fft(double data[]) {
		
		int samples=(int) Math.pow(2.0, Math.round(Math.log(data.length)/Math.log(2)));
		
		int halfSpectrum=samples/2;

		double[] re = new double[samples];				// in place FFT, input N samples, output N complex results i.e. (re, imaginary)
		double[] im = new double[samples];				// in place FFT, input N samples, output N complex results i.e. (re, imaginary)

		double mag[]=new double[halfSpectrum];
		
		for (int i = 0; i < samples; i++) {
			im[i] = 0;
			re[i]=data[i];
		}
		
		fft(re,im);
		
		for (int i = 0; i < halfSpectrum; i++) 
			mag[i] = Math.sqrt(re[i] * re[i] + im[i] * im[i]);

         /* Frequency = Fs * max / N
             Fs = sample rate (Hz)
             max = index of peak
             N = number of points in FFT
          */
		return mag;
	}
	/*
	 *  Copyright 2006-2007 Columbia University.
	 *
	 *  This file is part of MEAPsoft.
	 *
	 *  MEAPsoft is free software; you can redistribute it and/or modify
	 *  it under the terms of the GNU General Public License version 2 as
	 *  published by the Free Software Foundation.
	 *
	 *  MEAPsoft is distributed in the hope that it will be useful, but
	 *  WITHOUT ANY WARRANTY; without even the implied warranty of
	 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	 *  General Public License for more details.
	 *
	 *  You should have received a copy of the GNU General Public License
	 *  along with MEAPsoft; if not, write to the Free Software
	 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
	 *  02110-1301 USA
	 *
	 *  See the file "COPYING" for the text of the license.
	 */

	/***************************************************************
	 * fft.c
	 * Douglas L. Jones
	 * University of Illinois at Urbana-Champaign
	 * January 19, 1992
	 * http://cnx.rice.edu/content/m12016/latest/
	 *
	 *   fft: in-place radix-2 DIT DFT of a complex input
	 *
	 *   input:
	 * n: length of FFT: must be a power of two
	 * m: n = 2**m
	 *   input/output
	 * x: double array of length n with real part of data
	 * y: double array of length n with imag part of data
	 *
	 *   Permission to copy and use this program is granted
	 *   as long as this header is included.
	 ****************************************************************/


	  public void fft(double[] x, double[] y)
	  {
		  int n, m;
	  
		  double[] cos;
		  double[] sin;
		  
		  int i,j,k,n1,n2,a;
		  double c,s,t1,t2;
		  
	      n = x.length;
	      m = (int)(Math.log(n) / Math.log(2));
	       
	      // precompute tables
	      cos = new double[n/2];
	      sin = new double[n/2];
	      
	      for(i=0; i<n/2; i++) {
	          cos[i] = Math.cos(-2*Math.PI*i/n);
	          sin[i] = Math.sin(-2*Math.PI*i/n);
	      }
	  
	    // Bit-reverse
	    j = 0;
	    n2 = n/2;
	    for (i=1; i < n - 1; i++) {
	      n1 = n2;
	      while ( j >= n1 ) {
	        j = j - n1;
	        n1 = n1/2;
	      }
	      j = j + n1;
	    
	      if (i < j) {
	        t1 = x[i];
	        x[i] = x[j];
	        x[j] = t1;
	        t1 = y[i];
	        y[i] = y[j];
	        y[j] = t1;
	      }
	    }

	    // FFT
	    n1 = 0;
	    n2 = 1;
	  
	    for (i=0; i < m; i++) {
	      n1 = n2;
	      n2 = n2 + n2;
	      a = 0;
	    
	      for (j=0; j < n1; j++) {
	        c = cos[a];
	        s = sin[a];
	        a +=  1 << (m-i-1);

	        for (k=j; k < n; k=k+n2) {
	          t1 = c*x[k+n1] - s*y[k+n1];
	          t2 = s*x[k+n1] + c*y[k+n1];
	          x[k+n1] = x[k] - t1;
	          y[k+n1] = y[k] - t2;
	          x[k] = x[k] + t1;
	          y[k] = y[k] + t2;
	        }
	      }
	    }
	  }                          
}
