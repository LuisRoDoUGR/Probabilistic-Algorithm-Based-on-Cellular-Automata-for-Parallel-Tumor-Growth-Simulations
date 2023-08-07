#include <iostream>
#include <vector>
#include "random.hpp"

using namespace std;

struct Celula {
	bool cancer = false;
	float cct;
	float ro;
	float mu;
	float alpha;
};

struct Coordenadas{
	int x;
	int y;
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
/*
void introducir_en_casilla( int i, int j, int casilla, Celula celula){
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

vector<int> casillas_libres( int i, int j){
	vector<int> libres;
	
	if ( i-1 >= 0 && j-1 >= 0 && i-1 < LONG && j-1 < LONG && !( type(rejilla[i-1][j-1]) == Celula) )
       	libres.push_back(0);
	if ( i >= 0 && j-1 >= 0 && i < LONG && j-1 < LONG && !( type(rejilla[i][j-1]) == Celula) )
       	libres.push_back(1);
	if ( i+1 >= 0 && j-1 >= 0 && i+1 < LONG && j-1 < LONG && !( type(rejilla[i+1][j-1]) == Celula) )
       	libres.push_back(2);
	if ( i+1 >= 0 && j >= 0 && i+1 < LONG && j < LONG && !( type(rejilla[i+1][j]) == Celula) )
       	libres.push_back(3);
	if ( i+1 >= 0 && j+1 >= 0 && i+1 < LONG && j+1 < LONG && !( type(rejilla[i+1][j+1]) == Celula) )
       	libres.push_back(4);
	if ( i >= 0 && j+1 >= 0 && i < LONG && j+1 < LONG && !( type(rejilla[i][j+1]) == Celula) )
       	libres.push_back(5);
	if ( i-1 >= 0 && j+1 >= 0 && i-1 < LONG && j+1 < LONG && !( type(rejilla[i-1][j+1]) == Celula) )
       	libres.append(6);
	if ( i-1 >= 0 && j >= 0 && i-1 < LONG && j < LONG && !( type(rejilla[i-1][j]) == Celula) )
       	libres.push_back(7);
	
	return libres;
}
*/
int main(){

	cout << "Hello World" << endl;
	
	vector <Coordenadas> matriz(LONG*LONG);
	Coordenadas coor;
		
	for( int i = 0; i < LONG; i++ ){
		coor.x = i;
		coor.y = i;
		matriz[i*LONG+i] = coor;
		for( int j = i+1; j < LONG; j++){
			coor.x = i;
			coor.y = j;
			matriz[i*LONG + j] = coor;
			coor.x = j;
			coor.y = i;
			matriz[j*LONG + i] = coor;
		}
	}	
	//for ( int i = 0; i < matriz.size(); i++)
	//	cout <<"("<< matriz[i].x << "," << matriz[i].y << ")-\t";
	
	vector < vector <Celula> > rejilla(LONG);
	for( int i = 0; i < LONG; i++)
		rejilla[i] = vector <Celula> (LONG);
	
	/*if( rejilla[20][3].cancer )	
		cout << "Cancer" << endl;
	else
		cout << "Sana" << endl;
		
	return 0;*/
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
