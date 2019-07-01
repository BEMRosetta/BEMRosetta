#include "BEMRosetta.h"

/*
	double samplingFrecuency = 1/0.03;	
	String csvSep = ";";

	Vector<double> data;
	FileIn file(AppendFileName(GetDesktopFolder(), "data.txt"));
	
	file.GetLine();
	while (!file.IsEof()) {
		String line = file.GetLine();
		int pos = line.Find('\t');
		double freq = ScanDouble(line.Left(pos));
		double val = ScanDouble(line.Mid(pos+1));
		data << val;
	}
	
	int numData = data.GetCount();
	
    // Filling the data series
    VectorXd timebuf(numData);
    {
	    double t = 0;
	    for (int i = 0; i < numData; ++i) 
	       	timebuf[i] = data[i];
    }
    
    VectorXcd freqbuf;
    FFT<double> fft;
    fft.SetFlag(fft.HalfSpectrum);
    fft.fwd(freqbuf, timebuf);
	
	VectorXcd freqbuf2(freqbuf.size());
	{
	    for (int i = 0; i < freqbuf.size(); ++i) {
	        double freq = i*samplingFrecuency/numData;
	        double T = 1/freq;
	        if (T < 1.5)
	            freqbuf2[i] = 0;
	        else
	            freqbuf2[i] = freqbuf[i];
	    }
	}
	
	VectorXd timebuf2(numData);
	fft.inv(timebuf2, freqbuf2);
	
	// Saving original and filtered FFT
	{
	    String str;
	    str << "Frec" << csvSep << "T" << csvSep << "fft" << csvSep << "Filtered fft";
	    for (int i = 0; i < freqbuf.size(); ++i) {
	        double freq = i*samplingFrecuency/numData;
	        double T = 1/freq;
	        str << "\n" << freq << csvSep << (freq > 0 ? FormatDouble(T) : "") << csvSep 
	        			<< 2*std::abs(freqbuf[i])/numData << csvSep 
	        			<< 2*std::abs(freqbuf2[i])/numData;;
	    }
	    String fftFileName = AppendFileName(GetDesktopFolder(), "fft.csv");
	    Cout() << "\nFFT saved in '" << fftFileName << "'";
	    SaveFile(fftFileName, str);
	}
	
	// Saving original and filtered series
	{
	    String str;
	    str << "Time" << csvSep << "Data" << csvSep << "Filtered data";
	    double t = 0;
	    for (int i = 0; i < numData; ++i, t = i*1/samplingFrecuency) 
	       	str << "\n" << t << csvSep << timebuf[i] << csvSep << timebuf2[i];;
	    String dataFileName = AppendFileName(GetDesktopFolder(), "data.csv");
	    Cout() << "\nSource data saved in '" << dataFileName << "'";
	    SaveFile(dataFileName, str);
    }
}
*/