import matplotlib.pyplot as plt
import numpy as np
from glob import glob

directorio = glob("TABLES/*")

for fichero in directorio: 
	with open(fichero, "r+") as file1:
		long = file1.readline()
		nueva_rejilla = np.zeros(( int(long), int(long) ))
		for line in file1:
			nueva_rejilla[ int(line.split()[0]) ][ int(line.split()[1]) ] = 1

	plt.imshow(nueva_rejilla, vmin = 0, vmax = 1, cmap="Greys")
	plt.xlabel('Tiempo')
	plt.ylabel('Celulas')
	plt.savefig('IMAGES'+fichero[6:-3]+'png')
