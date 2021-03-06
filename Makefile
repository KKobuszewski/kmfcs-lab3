
CXX       = g++ -pipe
CXX_FLAGS = -mtune=native -march=native -m64 -O3 -fPIC -fopenmp
CXX_INC   = -I/usr/local/include -I.


all: bands

bands:
	$(CXX) -o bands.exe lab3.cpp $(CXX_FLAGS) $(CXX_INC)



clean:
	@rm *.exe
