package edu.ius.fftresult;

import java.io.Serializable;

public class FFTresult implements Serializable {
	private static final long serialVersionUID = 1L;
	public double magnitudeFFT[];
	public int frequencyFFT[];
	public int sampleRate;
	public int samples;
	
	public FFTresult(double []mag, int newSampleRate, int newSamples) {

		sampleRate = newSampleRate;
		samples = newSamples;
		
		int halfSpectrum = sampleRate / 2;
		int spectrumSize = Math.min(mag.length, halfSpectrum);
				
		frequencyFFT = new int[spectrumSize];
		magnitudeFFT = new double[spectrumSize];

		int freq = 0;
		double max=0;
		
		if (mag.length < halfSpectrum)
			for (int i = 0; i < mag.length; i++) {
				frequencyFFT[i] = (int) (sampleRate * (i/(double)samples));
				magnitudeFFT[i] = mag[i];
				if(magnitudeFFT[i] > max)
					max = magnitudeFFT[i];
			}
		else
			for (int i = 0; i < mag.length; i++) {
				freq = (int) (sampleRate * (i/(double)samples));
				magnitudeFFT[freq] = Math.max(magnitudeFFT[freq], mag[i]);
				frequencyFFT[freq] = freq;
				if(magnitudeFFT[freq] > max)
					max = magnitudeFFT[freq];
			}
		
		for(int i=0; i< magnitudeFFT.length; i++)							// Normalize to 0..1
			magnitudeFFT[i]=magnitudeFFT[i]/max;
		
	}
	
	public void setMagnitudeFFT(double newMagnitudeFFT[]) {
		magnitudeFFT = newMagnitudeFFT;
	}
	public void setFrequencyFFT(int newFrequencyFFT[]) {
		frequencyFFT = newFrequencyFFT;
	}
	public void setSampleRate(int newSampleRate) {
		sampleRate = newSampleRate;
	}
	public void setSamples(int newSamples) {
		samples = newSamples;
	}
	public double[] getMagnitudeFFT() {
		return magnitudeFFT;
	}
	public int []  getFrequencyFFT() {
		return frequencyFFT;
	}
	public int getSampleRate() {
		return sampleRate;
	}
	public int getSamples() {
		return samples;
	}	
}
