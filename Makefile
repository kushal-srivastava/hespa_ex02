CXXFLAGS = -std=c++11 -O3
CUD = nvcc

CUDA:	
	mkdir tmp
	$(CUD) $(CXXFLAGS) cuda_m.cu -o hespa

run1:
	./hespa stable.par
run2:
	./hespa repell.par
run3:
	./hespa attract.par
run4:
	./hespa blocks.par

clean:
	 -rm hespa tmp/*.vtk
	 -rm -rf tmp
