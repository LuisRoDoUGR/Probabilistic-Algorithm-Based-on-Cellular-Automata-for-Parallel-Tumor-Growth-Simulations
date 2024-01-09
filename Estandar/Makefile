ct: 
	g++ -std=c++11 SOURCES/crecimientoTumoral.cpp -o BIN/ct
	
rm:
	rm -f TABLES/*	
	rm -f IMAGES/*	
	
run: ct rm
	./BIN/ct
	python3 SOURCES/grafiquear.py
	
clean:
	rm -f BIN/*
