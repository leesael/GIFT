CXX=g++

LIB_FLAGS = -larmadillo -llapack -lblas -DARMA_DONT_USE_WRAPPER

OPT = -O2 -mcmodel=medium  -fopenmp

CXXFLAGS = $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT)

all: ptf stf gift demo

demo: code/GIFT.cpp
	$(CXX) $(CXXFLAGS)  -o $@  $< $(LIB_FLAGS) 
	./demo sample/sample.data sample/sample.mask sample/ 3 10 10 50 2
	rm demo

ptf: code/PTucker.cpp 
	$(CXX) $(CXXFLAGS)  -o $@  $< $(LIB_FLAGS)

stf: code/Silenced-TF.cpp 
	$(CXX) $(CXXFLAGS)  -o $@  $< $(LIB_FLAGS)

gift: code/GIFT.cpp 
	$(CXX) $(CXXFLAGS)  -o $@  $< $(LIB_FLAGS) 



.PHONY: clean

clean:
	rm -f ptf stf gift

