CC = g++
CFLAGS = -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 
LIBS        = -lsz -lz -lm
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


SIMULSRC = TETHYS_1D_Main_v131.cpp Tethys1DLib.cpp 
SIMULOBJ = $(SIMULSRC:.cpp = .o)

ANALYSISSRC = TETHYS_1D_ElectronicAnalysis.cpp Tethys1DLib.cpp
ANALYSISOBJ = $(ANALYSISSRC:.cpp = .o)

TIMESERIESSRC = TETHYS_1D_TimeSeries.cpp Tethys1DLib.cpp
TIMESERIESOBJ = $(TIMESERIESSRC:.cpp = .o)

#all: tethys

all: tethys \
     analysis \
     timeseries \



tethys: $(SIMULOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o RichtmyerHDF5 $(SIMULOBJ) $(H5LIBS)

analysis: $(ANALYSISOBJ)
	$(CC) $(CFLAGS) -o ElectronicAnalysis $(ANALYSISOBJ) $(FFTLIBS) $(CLIBS)

timeseries: $(TIMESERIESOBJ)
	$(CC) $(CFLAGS) -o TimeSeries $(TIMESERIESOBJ) $(CLIBS)



clean:
	rm -f *.o

