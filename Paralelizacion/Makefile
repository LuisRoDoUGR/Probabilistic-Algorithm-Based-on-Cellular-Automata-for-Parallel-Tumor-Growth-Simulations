ct: 
	g++ -std=c++11 crecimientoTumoral.cpp -fopenmp -o BIN/ct
	
rm:
	rm -f TABLES/*	
	rm -f IMAGES/*	
	
run: ct rm
	./BIN/ct
	python3 grafiquear.py
	
prueba:
	g++ -std=c++11 prueba.cpp -fopenmp -o BIN/prueba
	./BIN/prueba
	
clean:
	rm -f BIN/*
