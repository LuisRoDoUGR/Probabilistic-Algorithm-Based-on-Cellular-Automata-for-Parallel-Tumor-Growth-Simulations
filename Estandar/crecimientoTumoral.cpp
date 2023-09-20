#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include "random.hpp"

using namespace std;
using namespace std::chrono;
using Random = effolkronium::random_static;

// Clase Celula
struct Celula {
	//Indicar si es cancerigena o no
	bool cancer = false;
	//Cálculo de migración
	double cct;
	//Nº de reproducciones
	double ro;
	//Cálculo de Reproducción
	double mu;
	//Prababilidad de muerte
	double alpha;
	
	//Constructor con parametros, para cada uno de ellos
	Celula (bool c, double cc, double r, double m, double a){
		cancer = c;
		cct = cc;
		ro = r;
		mu = m;
		alpha = a;
	}
	
	//Constructor copia, genera una instancia nueva con los mismos parametros
	Celula (const Celula &cell){
		cancer = cell.cancer;
		cct = cell.cct;
		ro = cell.ro;
		mu = cell.mu;
		alpha = cell.alpha;
	}
	
	//Constructor por defecto, genera celula sin cancer, el resto de parametros sin inicializar
	Celula (){
		cancer = false;
	}
	
	//Definición del operador de asignación, copia los valores de los parametros de la celula a la derecha del operador
	Celula& operator=( const Celula* cell){
		cancer = cell->cancer;
		cct = cell->cct;
		ro = cell->ro;
		mu = cell->mu;
		alpha = cell->alpha;
		return *this;
	}
	
};

//Clase de Coordenadas, par de interos, para representar coordenadas de una matriz
struct Coordenadas{
	int x;
	int y;
	
	Coordenadas& operator=( const Coordenadas* coor){
		x = coor->x;
		y = coor->y;
		return *this;
	}
};

//Probabilidad de reproducción identica o normal
static double PS = 0.1;
//Tiempo al día, cálculo de migración
static double T = (1/24.0);

//Nº máximo de reproducciones de una célula
static double ROMAX = 10.0;
//Probabilidad máxima de muerte
static double ALPHAMAX = 0.01;

//Tamaño del grid
static int LONG = 601;
//Nº de días
static int NUM_DIAS = 121;

/* Introduce el valor de la celula, en el sitio correspondiente, alrededor de las coordenadas que se pasan
	como parametro. Al final de la función, hay una nueva celula en una casilla adyacente a las coordenadas
	que se han pasado
	@param i: Valor x de las coordenadas
	@param j: Valor y de las coordenadas
	@param casilla: Indicativo de en cuál de las 8 casillas adyacentes se ha de introducir la nueva célula, este valor debe ser válido
	@param celula: Nueva celula que se va a introducir en el grid
	@param rejilla: Grid de células que se está simulando
*/
Coordenadas introducir_en_casilla( int i, int j, int casilla, Celula celula, vector< vector <Celula> > &rejilla){
	Coordenadas coor;
	//Se comprueba cúal de los posibles 8 casillas es la que se ha pasado, y se introduce la célula en las coordenadas correspondientes
	if (casilla == 0){
        	rejilla[i - 1][j - 1] = celula;
        	coor.x = i-1;
        	coor.y = j-1;
    	}else if (casilla == 1){
        	rejilla[i][j - 1] = celula;
        	coor.x = i;
        	coor.y = j-1;
    	}
    	else if (casilla == 2){
        	rejilla[i + 1][j - 1] = celula;
        	coor.x = i+1;
        	coor.y = j-1;
    	}
    	else if (casilla == 3){
        	rejilla[i + 1][j] = celula;
        	coor.x = i+1;
        	coor.y = j;
    	}
    	else if (casilla == 4){
        	rejilla[i + 1][j + 1] = celula;
        	coor.x = i+1;
        	coor.y = j+1;
    	}
    	else if (casilla == 5){
        	rejilla[i][j + 1] = celula;
        	coor.x = i;
        	coor.y = j+1;
    	}
    	else if (casilla == 6){
        	rejilla[i - 1][j + 1] = celula;
        	coor.x = i-1;
        	coor.y = j+1;
    	}
    	else if (casilla == 7){
        	rejilla[i - 1][j] = celula;
        	coor.x = i-1;
        	coor.y = j;
    	}
    	return coor;
}

/* Cálcula cuales de las 8 casillas adyacentes a una dada están libres, es decir, sin células cancerigenas.
	@param i: Valor x de las coordenadas de la casilla
	@param j: Valor y de las coordenadas de la casilla
	@param rejilla: Grid de células que se está simulando
	@return libres: Vector con las posiciones alrededor de la casilla indicada, donde no hay células cancerígenas, indicadas con números del 0 al 7.
*/
vector<int> casillas_libres( int i, int j, vector< vector <Celula> > &rejilla){
	vector<int> libres(0); //Vector donde se guardan las posiciones
	//Desde la esquina superior izquierda, recorriendo las 8 posiciones circundantes en sentido horario, se comprueba que no hay células cancerígenas.
	if ( i-1 >= 0 && j-1 >= 0 && i-1 < LONG && j-1 < LONG && !rejilla[i-1][j-1].cancer )
       		libres.push_back(0);
	if ( i >= 0 && j-1 >= 0 && i < LONG && j-1 < LONG && !rejilla[i][j-1].cancer )
       		libres.push_back(1);
	if ( i+1 >= 0 && j-1 >= 0 && i+1 < LONG && j-1 < LONG && !rejilla[i+1][j-1].cancer )
       		libres.push_back(2);
	if ( i+1 >= 0 && j >= 0 && i+1 < LONG && j < LONG && !rejilla[i+1][j].cancer )
       		libres.push_back(3);
	if ( i+1 >= 0 && j+1 >= 0 && i+1 < LONG && j+1 < LONG && !rejilla[i+1][j+1].cancer )
       		libres.push_back(4);
	if ( i >= 0 && j+1 >= 0 && i < LONG && j+1 < LONG && !rejilla[i][j+1].cancer )
       		libres.push_back(5);
	if ( i-1 >= 0 && j+1 >= 0 && i-1 < LONG && j+1 < LONG && !rejilla[i-1][j+1].cancer )
       		libres.push_back(6);
	if ( i-1 >= 0 && j >= 0 && i-1 < LONG && j < LONG && !rejilla[i-1][j].cancer )
       		libres.push_back(7);
	
	return libres;
}

/* Cálculo de la acción que toma una célula en una fase de tiempo, según las probabilidades que se le pasan
	@param alpha: Probabilidad de morir
	@param migrar: Probabilidad de migrar
	@param pd: Probabilidad de reproducción
	@return: Valor entero del 1 al 3, indicando una de las 3 posibles acciones que toma una célula
*/
int accion_cell(double alpha, double migrar, double pd){
	//Se suman las probabilidades de cada acción
	double p_total = alpha + migrar + pd;
	//Se busca un número aleatorio entre 0 y las probabilidades sumadas
	double p_obtenida = Random::get(0.0, p_total);
	
	//Dependiendo de el número obtenido se elige una de las opciones
	if( p_obtenida < alpha ){
		return 1;
	} else if(p_obtenida < (alpha + migrar) ){
		return 2;
	} else if(p_obtenida < p_total ){
		return 3;
	}
	return 4;
}

/* Simulación durante 50 días, con 24 pasos por día del crecimiento de un tumor en un grid.
	@param matriz: Vector con todas las posibles coordenadas que tiene el grid
	@param rejilla: Grid de células donde se va a simular el crecimiento del tumor
*/
int simulacion_cancer(vector <Coordenadas> &matriz, vector <Coordenadas> &futuras, vector< vector <Celula> > &rejilla){
	// Declaración de variables auxiliares
	Coordenadas casilla_elegida;
	int i, j, action, contador_cells = 1, pasos = 24;
	Celula cell, nueva_cell;
	vector<int> cas_libres;
	double alpha, pd, migrar, p_obtenida;
	ofstream myfile;
	string file;
	vector <Coordenadas>::iterator it;
	
	auto new_extra = high_resolution_clock::now(); //Declaramos los valores que vamos a usar para el calculo de tiempo
	auto fin_extra = high_resolution_clock::now();
	auto extra = duration_cast<chrono::milliseconds>(fin_extra - fin_extra);
	int time = 0;

	//Durante 50 días
	for( int dia = 1; dia < NUM_DIAS; dia++){
		//En cada día 24 pasos
		for (int paso = 0; paso < pasos; paso++){
			//Aleatorizamos el acceso al grid
			Random::shuffle(matriz);
			it = futuras.begin();
			//Para el número de elementos del grid
			for( int indice = 0; indice < matriz.size(); indice++){
				//Vemos que casilla vamos a acceder
				casilla_elegida = matriz[indice];
				i = casilla_elegida.x; //Guardamos las coordenadas
				j = casilla_elegida.y;
				cell = rejilla[i][j]; //Y guardamos la célula en esa coordenada
				
				//Innecesario
				if( cell.cancer ){ //Si la célula es cancerígena
					cas_libres = casillas_libres( i, j, rejilla ); //Cálculamos las casillas libres alrededor
					
					if( cas_libres.size() > 0 ){ //Si tiene casillas libres alrededor
						//Se cálcula las probabilidades de la célula de cada acción según sus parametros
						alpha = cell.alpha;	
						pd = (24.0/cell.cct)*T; 
						migrar = (1-pd)*(cell.mu*T);
						if( cell.ro <= 0){
							pd = 0;
						}
						//Se escoge qué acción va a hacer la célula
						action = accion_cell(alpha, migrar, pd);
						
						if( action == 1 ){ //Si muere
							rejilla[i][j] = new Celula(); //Se reemplaza en el grid por una celula no cancerígena
							contador_cells --; //Disminuye el contador de células
						}else if( action == 2 ){ //Si migra
							Random::shuffle(cas_libres); //Aleatorizamos el acceso a las casilla adyacentes
							
							nueva_cell = new Celula(cell); //Se genera una copia de la célula
							
							rejilla[i][j] = new Celula(); //Eliminamos la célula del lugar actual en el que está
							
							(*it) = introducir_en_casilla( i, j, cas_libres[0], nueva_cell, rejilla); //Introducimos la copia en la casilla adyacente
							it ++;
							
								
						}else if( action == 3 ){ // Si se reproduce
							Random::shuffle(cas_libres); //Aleatorizamos el acceso a las casilla adyacentes
							
							if( cell.alpha == 0 ){ //Si no puede morir (Célula madre)
								p_obtenida = Random::get(0.0, 1.0); //Cálculamos si sera copia exacta o hija
								if( p_obtenida < PS ){ //Si es copia exacta
									nueva_cell = new Celula(cell); //Hacemos la copia
									(*it) = introducir_en_casilla( i, j, cas_libres[0], nueva_cell, rejilla); //Introducimos la copia al lado
									it ++;
									
								}else{ //Si es hija
									nueva_cell = new Celula(true, cell.cct, ROMAX, cell.mu, ALPHAMAX); //Creamos una hija con nuevos parametros
									(*it) = introducir_en_casilla( i, j, cas_libres[0], nueva_cell, rejilla); //Introducimos la hija al lado
									it ++;
								}
							} else { //Si no es célula madre
								// ESTO ESTÁ ASÍ EN EL ORIGINAL Y NO SE SI DEBERÍA ESTAR AL REVES
								cell.ro = cell.ro-1; //Reducimos el número de reproducciones de la célula
								nueva_cell = new Celula(cell); //Hacemos una copia
								(*it) = introducir_en_casilla( i, j, cas_libres[0], nueva_cell, rejilla); //Introducimos la hija al lado
								it ++;
								
							}
							(*it) = casilla_elegida;
							it ++;
							
							contador_cells ++; //Aumentamos el número de células
						}
					}
					
				}
			}
			matriz.assign(futuras.begin(), it);
		}
		new_extra = high_resolution_clock::now();		
		//Cada día indicamos cuantas células hay
		cout << "Día: " << dia << endl;
		cout << "Numero de celulas: "<< matriz.size() << endl;
		
		if (dia % 5 == 0){
			
			file = "TABLES/" + to_string(dia) + "dia.txt";
						
			myfile.open( file );
	  		if (myfile.is_open()){
	  			myfile << LONG << "\n";
	  			for( int filas = 0; filas < rejilla.size(); filas++ ){
	  				for ( int columnas = 0; columnas < rejilla[filas].size(); columnas++){
		  				if( rejilla[filas][columnas].cancer )
			  				myfile << filas << " " << columnas << "\n";
		  			}
	  			}
	    			myfile.close();
	  		} else cout << "Unable to open file";
	  	}
	  	
	  	fin_extra = high_resolution_clock::now();
	  	extra = duration_cast<milliseconds>(fin_extra - new_extra);
		time += extra.count();
	}
	return time;
}


int main(){
	//Inicializamos el vector con las coordenadas
	vector <Coordenadas> matriz; // (LONG*LONG);
	vector <Coordenadas> futuras(LONG*LONG);
	Coordenadas coor;
	/*	
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
	}*/
	//for ( int i = 0; i < matriz.size(); i++)
	//	cout <<"("<< matriz[i].x << "," << matriz[i].y << ")-\t";
	
	// Inicializamos el grid de Células
	vector < vector <Celula> > rejilla(LONG);
	for( int i = 0; i < LONG; i++)
		rejilla[i] = vector <Celula> (LONG);
		
	//Introducimos en la mitad del grid una célula madre	
	Celula cell(true, 24.0 , 1000000000.0 , 100.0, 0.0);
	
	rejilla[(LONG - 1)/2][(LONG - 1)/2] = cell;
	
	coor.x = (LONG - 1)/2;
	coor.y = (LONG - 1)/2;
	matriz.push_back(coor);
	
	//Simulamos el grid cálculando el tiempo que tarda
	auto start = high_resolution_clock::now(); //Declaramos los valores que vamos a usar para el calculo de tiempo
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	int tiempo = 0;
	
	start = high_resolution_clock::now();
	tiempo -= simulacion_cancer(matriz, futuras, rejilla);
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	tiempo += duration.count();
	
	cout << "Tiempo: " << tiempo << endl;
	
	return 0;
}

