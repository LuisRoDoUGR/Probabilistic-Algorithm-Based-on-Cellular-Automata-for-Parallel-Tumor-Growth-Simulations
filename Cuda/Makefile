ct: 
	nvcc crecimientoTumoral.cu -o BIN/ct
	
rm:
	rm -f TABLES/*	
	rm -f IMAGES/*	
	
run: ct rm
	./BIN/ct
	python3 grafiquear.py
	
test:
	nvcc test1.cu -lcurand -o BIN/test
	./BIN/test
	
wat:
	nvcc wat.cu -lcurand -o BIN/wat
	./BIN/wat
		
dev:
	nvcc deviceRand.cu -lcurand -o BIN/dev
	./BIN/dev

clean:
	rm -f BIN/*
