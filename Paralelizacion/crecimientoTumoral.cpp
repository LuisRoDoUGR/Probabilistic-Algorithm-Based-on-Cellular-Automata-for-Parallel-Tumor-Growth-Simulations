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
	*/
		
	Celula cell;
	cell.cancer = true;
	cell.cct = 24;
	cell.ro = 100000000000;
	cell.mu = 100;
	cell.alpha = 0;
	
	rejilla[(long - 1)/2][(long - 1)/2] = cell;
	
	return 0;
}

