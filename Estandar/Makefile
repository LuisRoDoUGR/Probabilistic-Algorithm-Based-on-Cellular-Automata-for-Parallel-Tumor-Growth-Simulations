ct: 
	g++ -std=c++11 *.cpp -o BIN/ct
	
rm:
	rm -f TABLES/*	
	rm -f IMAGES/*	
	
run: ct rm
	./BIN/ct
	python3 grafiquear.py
	
clean:
	rm -f BIN/*
