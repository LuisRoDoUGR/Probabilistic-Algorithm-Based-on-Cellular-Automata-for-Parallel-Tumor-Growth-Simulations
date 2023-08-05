#include <iostream>
#include "random.hpp"
using namespace std;

struct Celula {
	float cct;
	float ro;
	float mu;
	float alpha;
};
// Probablemente funciones
struct Opciones {
	float morir;
	float migrar;
	float proliferar;
};

struct OpcionesHijas {
	float Iguales;
	float Distintas;
};

static float PS = 0.1;
static float T = 1/24;

static int ROMAX = 10;
static float ALPHAMAX = 0.01;

static int LONG = 201;

void introducir_en_casilla( int casilla, Celula celula){
	if (casilla == 0)
        	rejilla[i - 1][j - 1] = celula;
    	else if (casilla == 1)
        	rejilla[i][j - 1] = celula;
    	else if (casilla == 2)
        	rejilla[i + 1][j - 1] = celula;
    	else if (casilla == 3)
        	rejilla[i + 1][j] = celula;
    	else if (casilla == 4)
        	rejilla[i + 1][j + 1] = celula;
    	else if (casilla == 5)
        	rejilla[i][j + 1] = celula;
    	else if (casilla == 6)
        	rejilla[i - 1][j + 1] = celula;
    	else if (casilla == 7)
        	rejilla[i - 1][j] = celula;
}

int* casillas_libres( int i, int j){
	int libres[8];
	for( int h = 0; h < 8; h++){
		libres[h] = h;
	}
	
	if ( i-1 >= 0 && j-1 >= 0 && i-1 < LONG && j-1 < LONG && !( type(rejilla[i-1][j-1]) == Celula) )
       	libres.append(0);
	if ( i >= 0 && j-1 >= 0 && i < LONG && j-1 < LONG && !( type(rejilla[i][j-1]) == Celula) )
       	libres.append(1);
	if ( i+1 >= 0 && j-1 >= 0 && i+1 < LONG && j-1 < LONG && !( type(rejilla[i+1][j-1]) == Celula) )
       	libres.append(2);
	if ( i+1 >= 0 && j >= 0 && i+1 < LONG && j < LONG && !( type(rejilla[i+1][j]) == Celula) )
       	libres.append(3);
	if ( i+1 >= 0 && j+1 >= 0 && i+1 < LONG && j+1 < LONG && !( type(rejilla[i+1][j+1]) == Celula) )
       	libres.append(4);
	if ( i >= 0 && j+1 >= 0 && i < LONG && j+1 < LONG && !( type(rejilla[i][j+1]) == Celula) )
       	libres.append(5);
	if ( i-1 >= 0 && j+1 >= 0 && i-1 < LONG && j+1 < LONG && !( type(rejilla[i-1][j+1]) == Celula) )
       	libres.append(6);
	if ( i-1 >= 0 && j >= 0 && i-1 < LONG && j < LONG && !( type(rejilla[i-1][j]) == Celula) )
       	libres.append(7);
	
	return libres;
}


int main(){
	cout << "Hello World" << endl;
	
	return 0;
}
/**
    "def seleccionar_casilla(tabla):\n",
    "\n",
    "    casilla_encontrada = False\n",
    "\n",
    "    while not casilla_encontrada:\n",
    "        casilla_i = np.random.randint(long)\n",
    "        casilla_j = np.random.randint(long)\n",
    "\n",
    "        if not tabla[casilla_i][casilla_j] == 1:\n",
    "            tabla[casilla_i][casilla_j] = 1\n",
    "            return [casilla_i, casilla_j]\n",
    "\n",
    "rejilla = np.empty((long, long), dtype=Celula)\n",
    "rejilla[int((long - 1) / 2)][int((long - 1) / 2)] = Celula(24, 100000000000, 100, 0)\n",
    "\n",
    "celulas = [0]\n",
    "tiempos = [0]\n",
    "\n",
    "tiempo_extra_total = 0\n",
    "start_time = time()\n",
    "for t in range(1, 51):\n",
    "    #t = 1 # dias (24*180) pasos\n",
    "    pasos = 24\n",
    "\n",
    "    for paso in range(pasos):\n",
    "        tabla_hash = np.zeros((long, long))\n",
    "        for indice in range(long*long):\n",
    "\n",
    "            casilla_elegida = seleccionar_casilla(tabla_hash)\n",
    "            i = casilla_elegida[0]\n",
    "            j = casilla_elegida[1]\n",
    "            celula = rejilla[i][j]\n",
    "\n",
    "            if type(celula) == Celula: # si hay una celula\n",
    "                libres = casillas_libres(i, j)\n",
    "\n",
    "                if len(libres) > 0:\n",
    "                    alpha = celula.alpha\n",
    "                    pd = (24/celula.cct)*T\n",
    "                    migrar = (1-pd)*(celula.mu*T)\n",
    "\n",
    "                    pesos = [alpha, migrar, pd]\n",
    "                    if celula.ro == 0:\n",
    "                        pesos[2] = 0\n",
    "\n",
    "                    opcion = np.random.choice(opciones, size = 1, p = (pesos/(np.sum(pesos))))[0]\n",
    "                 \n",
    "                    if opcion == 'morir':\n",
    "                        # Vaciamos casilla\n",
    "                        rejilla[i][j] = None\n",
    "                    elif opcion == 'migrar':\n",
    "                        # Buscamos aleatorio una casilla libre para migrar\n",
    "                        casilla = np.random.choice(libres)\n",
    "\n",
    "                        copia = copy.deepcopy(celula)\n",
    "                        # Vaciamos casilla\n",
    "                        rejilla[i][j] = None\n",
    "\n",
    "                        introducir_en_casilla(casilla, copia)\n",
    "\n",
    "                    elif opcion == 'proliferar':\n",
    "                        # Buscamos aleatorio una casilla libre para la célula hija\n",
    "                        casilla = np.random.choice(libres)\n",
    "\n",
    "                        if celula.alpha == 0: # Si célula madre\n",
    "                            # Vemos si serán dos celulas distintas o iguales\n",
    "                            pesos_hija = [ps, 1-ps]\n",
    "                            celula_hija = np.random.choice(opciones_hijas, size = 1, p = (pesos_hija/(np.sum(pesos_hija))))[0]\n",
    "\n",
    "                            if celula_hija == 'iguales':\n",
    "                                introducir_en_casilla(casilla, copy.deepcopy(celula))\n",
    "\n",
    "                            else:\n",
    "                                copia = copy.deepcopy(celula)\n",
    "\n",
    "                                copia.alpha = alpha_max\n",
    "                                copia.ro = romax\n",
    "\n",
    "                                introducir_en_casilla(casilla, copia)\n",
    "                        else:\n",
    "                            celula.ro -= 1\n",
    "                            introducir_en_casilla(casilla, copy.deepcopy(celula))\n",
    "\n",
    "    start_time_extra = time()\n",
    "    numero_celulas = 0\n",
    "\n",
    "    nueva_rejilla = np.zeros((long,long))\n",
    "\n",
    "    for i in range(0, len(rejilla)):  # para cada columna\n",
    "        for j in range(0, len(rejilla[i])):  # para cada celda de la columna\n",
    "            # Si hay una célula\n",
    "            if not rejilla[i][j] is None:\n",
    "                # Añadimos un 1 para representación gráfica\n",
    "                nueva_rejilla[i][j] = 1.0\n",
    "                # Sumamos 1 al número de células\n",
    "                numero_celulas += 1\n",
    "\n",
    "    celulas.append(numero_celulas)\n",
    "\n",
    "    # Imprimimos datos\n",
    "    print('\\nDía: ' + str(t))\n",
    "    print(' Número de células: ' + str(numero_celulas))\n",
    "\n",
    "    if t == 1 or t%5 == 0:\n",
    "        plt.imshow(nueva_rejilla, vmin = 0, vmax = 1, cmap=\"Greys\")\n",
    "        plt.title('Día ' + str(t))\n",
    "        plt.show()\n",
    "    end_time_extra = time()\n",
    "    tiempo_extra = end_time_extra - start_time_extra\n",
    "    tiempo_extra_total += tiempo_extra\n",
    "\n",
    "end_time = time()\n",
    "tiempo = end_time - start_time\n",
    "print('\\nTiempo: ' + str(tiempo - tiempo_extra_total))\n",
    "\n",
    "tiempos.append(tiempo)\n",
    "\n",
    "plt.figure()\n",
    "plt.plot(celulas)\n",
    "plt.xlabel(\"Días\")\n",
    "plt.ylabel(\"Células\")\n",
    "plt.show()\n",
    "plt.plot(tiempos)\n",
    "plt.xlabel(\"Días\")\n",
    "plt.ylabel(\"Segundos\")\n",
    "plt.show()"
   ]

**/
