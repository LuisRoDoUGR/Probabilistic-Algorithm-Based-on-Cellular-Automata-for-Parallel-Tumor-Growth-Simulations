#include <iostream>
#include <vector>
#include "random.hpp"

using namespace std;
using Random = effolkronium::random_static;

struct Celula {
	bool cancer = false;
	int cct;
	int ro;
	int mu;
	float alpha;
	Celula (bool c, int cc, int r, int m, float a){
		cancer = c;
		cct = cc;
		ro = r;
		mu = m;
		alpha = a;
	}
	
	Celula (const Celula &cell){
		cancer = cell.cancer;
		cct = cell.cct;
		ro = cell.ro;
		mu = cell.mu;
		alpha = cell.alpha;
	}
	
	Celula (){
		cancer = false;
	}
	
	Celula& operator=( const Celula* cell){
		cancer = cell->cancer;
		cct = cell->cct;
		ro = cell->ro;
		mu = cell->mu;
		alpha = cell->alpha;
		return *this;
	}
	
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

void introducir_en_casilla( int i, int j, int casilla, Celula celula, vector< vector <Celula> > &rejilla){
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

vector<int> casillas_libres( int i, int j, vector< vector <Celula> > &rejilla){
	vector<int> libres;
	
	if ( i-1 >= 0 && j-1 >= 0 && i-1 < LONG && j-1 < LONG && !rejilla[i-1][j-1].cancer )
       	libres.push_back(0);
	if ( i >= 0 && j-1 >= 0 && i < LONG && j-1 < LONG && !rejilla[i-1][j-1].cancer )
       	libres.push_back(1);
	if ( i+1 >= 0 && j-1 >= 0 && i+1 < LONG && j-1 < LONG && !rejilla[i-1][j-1].cancer )
       	libres.push_back(2);
	if ( i+1 >= 0 && j >= 0 && i+1 < LONG && j < LONG && !rejilla[i-1][j-1].cancer )
       	libres.push_back(3);
	if ( i+1 >= 0 && j+1 >= 0 && i+1 < LONG && j+1 < LONG && !rejilla[i-1][j-1].cancer )
       	libres.push_back(4);
	if ( i >= 0 && j+1 >= 0 && i < LONG && j+1 < LONG && !rejilla[i-1][j-1].cancer )
       	libres.push_back(5);
	if ( i-1 >= 0 && j+1 >= 0 && i-1 < LONG && j+1 < LONG && !rejilla[i-1][j-1].cancer )
       	libres.push_back(6);
	if ( i-1 >= 0 && j >= 0 && i-1 < LONG && j < LONG && !rejilla[i-1][j-1].cancer )
       	libres.push_back(7);
	
	return libres;
}

int accion_cell(float alpha, float migrar, float pd){
	float p_total = alpha + migrar + pd;
	float p_obtenida = Random::get(0.0f, p_total);
	
	if( p_obtenida < alpha ){
		return 1;
	} else if(p_obtenida < (alpha + migrar) ){
		return 2;
	} else 
		return 3;
	
}


void simulacion_cancer(vector <Coordenadas> &matriz, vector< vector <Celula> > &rejilla){
	Coordenadas casilla_elegida;
	int i, j, action;
	Celula cell, nueva_cell;
	vector<int> cas_libres;
	float alpha, pd, migrar, p_obtenida;

	for( int a = 0; a < 51; a++){
		Random::shuffle(matriz);
		for( int indice = 0; indice < LONG*LONG; indice++){
			casilla_elegida = matriz[indice];
			i = casilla_elegida.x;
			j = casilla_elegida.y;
			cell = rejilla[i][j];
			
			if( cell.cancer ){
				cas_libres = casillas_libres( i, j, rejilla );
				
				if( cas_libres.size() > 0 ){
					alpha = cell.alpha;
					pd = (24.0/cell.cct)*T; //Operaciones enteras
					migrar = (1-pd)*(cell.mu*T);
					
					if( cell.ro == 0)
						pd = 0;
						
					action = accion_cell(alpha, migrar, pd);
					
					if( action == 1 ){
						rejilla[i][j] = new Celula();
					}else if( action == 2 ){
						Random::shuffle(cas_libres);
						
						nueva_cell = new Celula(cell);
						
						rejilla[i][j] = new Celula();
						
						introducir_en_casilla( i, j, cas_libres[0], nueva_cell, rejilla);
					
					}else if( action == 3 ){
						Random::shuffle(cas_libres);
						
						if( cell.alpha == 0 ){
							p_obtenida = Random::get(0.0, 1.0);
							if( p_obtenida < PS ){
								introducir_en_casilla( i, j, cas_libres[0], cell, rejilla);
							}else{
								nueva_cell = new Celula(true, cell.cct, ROMAX, cell.mu, ALPHAMAX);
								introducir_en_casilla( i, j, cas_libres[0], nueva_cell, rejilla);
							}
						} else {
							cell.ro = cell.ro-1;
							introducir_en_casilla( i, j, cas_libres[0], cell, rejilla);
						}
					}
				}
			}
		}
	}
}


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
		
	Celula cell(true, 24 , 1000000000, 100, 0);
	
	rejilla[(LONG - 1)/2][(LONG - 1)/2] = cell;
	
	simulacion_cancer(matriz, rejilla);
	
	return 0;
}

