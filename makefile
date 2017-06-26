CC=g++

IDIR1 =${ROOTSYS}/include
IDIR2 =${ANITA_UTIL_INSTALL_DIR}/include 

ROOT_FLAGS = `root-config --cflags`
CFLAGS=-I$(IDIR1) -I$(IDIR2) -g -O0 
ROOT_LIBS = `root-config --libs`
ANITA_LIBS = ${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaAnalysis.so ${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaAnalysisTools.so ${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaCorrelator.so ${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaEvent.so ${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaMagicDisplay.so ${ANITA_UTIL_INSTALL_DIR}/lib/libRootFftwWrapper.so ${ANITA_UTIL_INSTALL_DIR}/lib/libUCorrelator.so


DEPS = myHybridFunctions.h FFTtools.h
OBJECT_FILES1 = myHybrid.o myHybridFunctions.o 
OBJECT_FILES2 = test.o myHybridFunctions.o 

  
all: myHybrid test

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

myHybrid: ${OBJECT_FILES1}
	$(CC) -o $@ $^ $(CFLAGS) ${ROOT_LIBS} ${ANITA_LIBS}
  
test: ${OBJECT_FILES2}
	$(CC) -o $@ $^ $(CFLAGS) ${ROOT_LIBS} ${ANITA_LIBS}

clean: 
	rm myHybrid myHybridFunctions.o myHybrid.o test test.o















