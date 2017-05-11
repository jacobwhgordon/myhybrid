CC=g++

IDIR1 =/home/gordon.661/root/root/include
IDIR2 =/home/gordon.661/anita/AnitaTools/include 

ROOT_FLAGS = `root-config --cflags`
CFLAGS=-I$(IDIR1) -I$(IDIR2) -g -O0 
ROOT_LIBS = `root-config --libs`
ANITA_LIBS = /home/gordon.661/anita/AnitaTools/lib/libAnitaAnalysis.so /home/gordon.661/anita/AnitaTools/lib/libAnitaAnalysisTools.so /home/gordon.661/anita/AnitaTools/lib/libAnitaCorrelator.so /home/gordon.661/anita/AnitaTools/lib/libAnitaEvent.so /home/gordon.661/anita/AnitaTools/lib/libAnitaMagicDisplay.so /home/gordon.661/anita/AnitaTools/lib/libRootFftwWrapper.so /home/gordon.661/anita/AnitaTools/lib/libUCorrelator.so


DEPS = myHybridFunctions.h FFTtools.h
OBJECT_FILES = myHybrid.o myHybridFunctions.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

myHybrid: ${OBJECT_FILES}
	$(CC) -o $@ $^ $(CFLAGS) ${ROOT_LIBS} ${ANITA_LIBS}
  

clean: 
	rm myHybrid myHybridFunctions.o myHybrid.o











