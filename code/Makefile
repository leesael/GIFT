CXX=g++

LIB_FLAGS = -larmadillo -llapack -lblas -DARMA_DONT_USE_WRAPPER

OPT = -O2 -mcmodel=medium  -fopenmp

CXXFLAGS = $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT)

all: ptf stf gift


ptf: PTucker.cpp 
	$(CXX) $(CXXFLAGS)  -o $@  $< $(LIB_FLAGS)

stf: Silenced-TF.cpp 
	$(CXX) $(CXXFLAGS)  -o $@  $< $(LIB_FLAGS)

gift: GIFT.cpp 
	$(CXX) $(CXXFLAGS)  -o $@  $< $(LIB_FLAGS) 



.PHONY: clean

clean:
	rm -f ptf stf gift

