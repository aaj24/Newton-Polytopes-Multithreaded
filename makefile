G=g++ -std=c++14 -O3 -w

main: main.cpp ArrayX.o PolyDecomp.o convex_hull.o priority_update.o
	$(G) -fopenmp $^ -o $@
PolyDecomp.o: PolyDecomp.cpp PolyDecomp.hpp
	$(G) -fopenmp $^ -c
ArrayX.o: ArrayX.cpp ArrayX.hpp
	$(G) $^ -c
convex_hull.o: convex_hull.cpp convex_hull.hpp
	$(G) convex_hull.cpp convex_hull.hpp -c
priority_update.o: priority_update.cpp priority_update.hpp
	$(G) priority_update.cpp priority_update.hpp -c
generator: generator.cpp convex_hull.cpp convex_hull.hpp
	$(G) $^ -o $@
clean:
	rm *.o *.gch main
