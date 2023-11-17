#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include "random.hpp"

using namespace std;
using namespace std::chrono;
using Random = effolkronium::random_static;


//Probabilidad de reproducción identica o normal
static double PS = 0.1;
//Tiempo al día, cálculo de migración
static double T = (1/24.0);

//Nº máximo de reproducciones de una célula
static double ROMAX = 10.0;
//Probabilidad máxima de muerte
static double ALPHAMAX = 0.01;

//Tamaño del grid
static int LONG = 751;
//Nº de días
static int NUM_DIAS = 11;

static vector< int > FACTORIAL(9);


class Celula {
	
	public:
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
	//Numero de vecinos
	int n_vecinos;
	//Numero de vecinos
	Celula * vecinos[8];
	
	double P_morir;
	double P_migra;
	double P_repro;
	// HAY QUE INICIALIZAR ESTOS DOS EN ALGUN MOMENTO
	
	
	//Constructor con parametros, para cada uno de ellos
	Celula (bool c, double cc, double r, double m, double a){
		cancer = c;
		cct = cc;
		ro = r;
		mu = m;
		alpha = a;
		n_vecinos = 0;
		
		vecinos[0] = NULL;
		vecinos[1] = NULL;
		vecinos[2] = NULL;
		vecinos[3] = NULL;
		vecinos[4] = NULL;
		vecinos[5] = NULL;
		vecinos[6] = NULL;
		vecinos[7] = NULL;
		
		calculoProbabilidades();
	}
	
	//Constructor copia, genera una instancia nueva con los mismos parametros
	Celula (const Celula &cell){
		cancer = cell.cancer;
		cct = cell.cct;
		ro = cell.ro;
		mu = cell.mu;
		alpha = cell.alpha;
		n_vecinos = cell.n_vecinos;
		P_morir = cell.P_morir;
		P_migra = cell.P_migra;
		P_repro = cell.P_repro;
		
		vecinos[0] = NULL;
		vecinos[1] = NULL;
		vecinos[2] = NULL;
		vecinos[3] = NULL;
		vecinos[4] = NULL;
		vecinos[5] = NULL;
		vecinos[6] = NULL;
		vecinos[7] = NULL;
	}
	
	//Constructor por defecto, genera celula sin cancer, el resto de parametros sin inicializar
	Celula (){
		cancer = false;
		ro = 0.0;
		n_vecinos = 0;
		P_morir = 0;
		P_migra = 0;
		P_repro = 0;
		vecinos[0] = NULL;
		vecinos[1] = NULL;
		vecinos[2] = NULL;
		vecinos[3] = NULL;
		vecinos[4] = NULL;
		vecinos[5] = NULL;
		vecinos[6] = NULL;
		vecinos[7] = NULL;
	}
	
	//Definición del operador de asignación, copia los valores de los parametros de la celula a la derecha del operador
	Celula& operator=( const Celula* cell){
		cancer = cell->cancer;
		cct = cell->cct;
		ro = cell->ro;
		mu = cell->mu;
		alpha = cell->alpha;
		n_vecinos = 0;
		P_morir = cell->P_morir;
		P_migra = cell->P_migra;
		P_repro = cell->P_repro;
		vecinos[0] = NULL;
		vecinos[1] = NULL;
		vecinos[2] = NULL;
		vecinos[3] = NULL;
		vecinos[4] = NULL;
		vecinos[5] = NULL;
		vecinos[6] = NULL;
		vecinos[7] = NULL;
		return *this;
	}
		
	void calculoProbabilidades(){
		double pd = (24.0/cct)*T; 
		double migrar = (1-pd)*(mu*T);
		
		if( ro <= 0){
			pd = 0.0;
		}
		//Se escoge qué acción va a hacer la célula
		double p_total = alpha + migrar + pd;
		P_morir = alpha/p_total;
		P_migra = migrar/p_total;
		P_repro = pd/p_total;
	}
	
	void decreaseRo(){
		ro --;
		if(ro <= 0)
			calculoProbabilidades();
	}

	
	/* Cálcula cuales de las 8 casillas adyacentes a una dada están libres, es decir, sin células cancerigenas.
	@param i: Valor x de las coordenadas de la casilla
	@param j: Valor y de las coordenadas de la casilla
	@param rejilla: Grid de células que se está simulando
	@return libres: Vector con las posiciones alrededor de la casilla indicada, donde no hay células cancerígenas, indicadas con números del 0 al 7.
	*/	
	void actualizarVecinos( int i, int j, vector< vector <Celula> > &rejilla){
		n_vecinos = 0;
		vecinos[0] = NULL;
		vecinos[1] = NULL;
		vecinos[2] = NULL;
		vecinos[3] = NULL;
		vecinos[4] = NULL;
		vecinos[5] = NULL;
		vecinos[6] = NULL;
		vecinos[7] = NULL;
		 //Vector donde se guardan las posiciones
		//Desde la esquina superior izquierda, recorriendo las 8 posiciones circundantes en sentido horario, se comprueba que no hay células cancerígenas.
		if ( i-1 >= 0  && j-1 >= 0 && i-1 < LONG && j-1 < LONG && rejilla[i-1][j-1].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[0] = &rejilla[i-1][j-1];
	       		n_vecinos++;
		}if ( i >= 0   && j-1 >= 0 && i < LONG   && j-1 < LONG && rejilla[i][j-1].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[1] = &rejilla[i][j-1];
	       		n_vecinos++;
		}if ( i+1 >= 0 && j-1 >= 0 && i+1 < LONG && j-1 < LONG && rejilla[i+1][j-1].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[2] = &rejilla[i+1][j-1];
	       		n_vecinos++;
		}if ( i+1 >= 0 && j >= 0   && i+1 < LONG && j < LONG   && rejilla[i+1][j].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[3] = &rejilla[i+1][j];
	       		n_vecinos++;
		}if ( i+1 >= 0 && j+1 >= 0 && i+1 < LONG && j+1 < LONG && rejilla[i+1][j+1].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[4] = &rejilla[i+1][j+1];
	       		n_vecinos++;
		}if ( i >= 0   && j+1 >= 0 && i < LONG   && j+1 < LONG && rejilla[i][j+1].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[5] = &rejilla[i][j+1];
	       		n_vecinos++;
		}if ( i-1 >= 0 && j+1 >= 0 && i-1 < LONG && j+1 < LONG && rejilla[i-1][j+1].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[6] = &rejilla[i-1][j+1];
	       		n_vecinos++;
		}if ( i-1 >= 0 && j >= 0   && i-1 < LONG && j < LONG   && rejilla[i-1][j].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[7] = &rejilla[i-1][j];
	       		n_vecinos++;
		}
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
    	else if (casilla == 8){
        	rejilla[i][j] = celula;
        	coor.x = i;
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



double combinaciones(int elegido, int inicio, int vecinos, int nivel, vector< double > &probabilidad_v){
	double p = 0;
	for( int j = inicio; j < vecinos; j++){
		if( j != elegido)
			p += (1-probabilidad_v[j])*( FACTORIAL[vecinos - 1 - nivel]*FACTORIAL[nivel] + combinaciones(elegido, j+1, vecinos, nivel+1, probabilidad_v) );
	}
	return p;
}

/* Cálculo de la acción que toma una célula en una fase de tiempo, según las probabilidades que se le pasan
	@param alpha: Probabilidad de morir
	@param migrar: Probabilidad de migrar
	@param pd: Probabilidad de reproducción
	@return: Valor entero del 1 al 3, indicando una de las 3 posibles acciones que toma una célula
*/
vector <double> calculoProbabilidadVecinos(Celula &cell){
	
	vector< double > probabilidad_v;
	vector< double > probabilidad_f;
	for( int i = 0 ; i < 8; i++){
		if( cell.vecinos[i] != NULL ){
			probabilidad_v.push_back( (cell.vecinos[i]->P_migra + cell.vecinos[i]->P_repro) / (8.0 - cell.vecinos[i]->n_vecinos) );
		}
	}
		
	for( int i = 0; i < cell.n_vecinos; i++){
		probabilidad_f.push_back( probabilidad_v[i]*(FACTORIAL[cell.n_vecinos-1] + combinaciones( i, 0, cell.n_vecinos, 1, probabilidad_v) ) / (double)cell.n_vecinos );
		 
	}
	
	return probabilidad_f;
}


/* Simulación durante 50 días, con 24 pasos por día del crecimiento de un tumor en un grid.
	@param matriz: Vector con todas las posibles coordenadas que tiene el grid
	@param rejilla: Grid de células donde se va a simular el crecimiento del tumor
*/
int simulacion_cancer(vector < vector< vector <Celula> > > &rejillas){
	int min_x = ((LONG - 1)/2) - 1, max_x = ((LONG - 1)/2) + 1, min_y = ((LONG - 1)/2) - 1, max_y = ((LONG - 1)/2) + 1;
	int indice_rejilla = 0;
	
	// Declaración de variables auxiliares
	vector <Coordenadas> celulas_madre;
	vector <Coordenadas> futuras;
	Coordenadas coor;
	coor.x = (LONG - 1)/2;
	coor.y = (LONG - 1)/2;
	celulas_madre.push_back( coor );
	
	int i, j, indice_libre, vecino, contador = 1, pasos = 24;
	Celula cell, nueva_cell;
	vector<int> cas_libres;
	vector<double> p_vecinos;
	double alpha, pd, migrar, p_obtenida, p_acumulada, p_auxiliar;
	double p_migrar, p_morir, p_repro;
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
			indice_rejilla = paso%2;
			for (int madre = 0; madre < celulas_madre.size(); madre++){
				cell = rejillas[indice_rejilla][celulas_madre[madre].x][celulas_madre[madre].y];
				cas_libres = casillas_libres( celulas_madre[madre].x, celulas_madre[madre].y, rejillas[indice_rejilla]);
				p_obtenida = Random::get(0.0, 1.0);
				if( cas_libres.size() != 0){
					vecino = Random::get(0, int(cas_libres.size()-1));
					
					if( p_obtenida < cell.P_migra ){
						coor = introducir_en_casilla( celulas_madre[madre].x, celulas_madre[madre].y, cas_libres[vecino], cell, rejillas[indice_rejilla]);
						futuras.push_back(coor);
						rejillas[indice_rejilla][celulas_madre[madre].x][celulas_madre[madre].y] = new Celula();
						
					}else{
						futuras.push_back(celulas_madre[madre]);	
						if( Random::get(0.0, 1.0) < PS){
							coor = introducir_en_casilla( celulas_madre[madre].x, celulas_madre[madre].y, cas_libres[vecino], cell, rejillas[indice_rejilla]);
							futuras.push_back( coor );
							contador ++;
						}else{
							nueva_cell = new Celula(true, cell.cct, ROMAX, cell.mu, ALPHAMAX);
							introducir_en_casilla( celulas_madre[madre].x, celulas_madre[madre].y, cas_libres[vecino], nueva_cell, rejillas[indice_rejilla]);
						}
					}
					
					if( coor.x == min_x && min_x != 0){
						min_x--;
					}else if( coor.x == max_x && max_x != LONG-1){
						max_x++;
					}else if( coor.y == min_y && min_y != 0){
						min_y--;
					}else if( coor.y == max_y && max_y != LONG-1)
						max_y++;
				} else {
					futuras.push_back(celulas_madre[madre]);
				}
			}
			celulas_madre = futuras;
			futuras.clear();
			
			for(int x = min_x; x <= max_x ; x++){ // Cálculo de vecinos
				for(int y = min_y; y <= max_y; y++){
					rejillas[indice_rejilla][x][y].actualizarVecinos(x, y, rejillas[indice_rejilla]);
				}
			}
			
			for(int x = min_x; x <= max_x ; x++){ // Cálculo de función de transición
				for(int y = min_y; y <= max_y; y++){	
					//Vemos que casilla vamos a acceder
					
					cell = rejillas[indice_rejilla][x][y]; //Y guardamos la célula en esa coordenada
					
					if( cell.cancer ){ //Si la célula es cancerígena
						if( cell.ro < 11 ){ //Células madre se tratan al final, en serie
							if( cell.n_vecinos > 0 ){ //Si tiene casillas libres alrededor
								if( cell.n_vecinos == 8){ // P=1 Quiescencia, se mantiene igual
									rejillas[(indice_rejilla+1)%2][x][y] = cell;
								} else {
									//Se cálcula las probabilidades de la célula de cada acción según sus parametros
									p_auxiliar = 0.0;
									p_vecinos = calculoProbabilidadVecinos(cell);	
									for( int vec = 0; vec < p_vecinos.size(); vec ++){
										p_auxiliar += p_vecinos[vec]; }
									
									p_acumulada = cell.P_repro + (1-cell.P_repro)*p_auxiliar;			
																
									p_obtenida = Random::get(0.0, 1.0);
									if( p_obtenida < cell.P_repro ){
										cell.decreaseRo();
										rejillas[(indice_rejilla+1)%2][x][y] = cell;
									} else if( p_obtenida < p_acumulada ){
										
										p_auxiliar = cell.P_repro;
										for(int vec = 0; vec < p_vecinos.size(); vec++){
											p_auxiliar += (1-cell.P_repro)*p_vecinos[vec];
											if( p_obtenida < p_auxiliar ){
												vecino = vec;
												break;
											}
											}
										for(int indice_vecinos = 0; indice_vecinos < 8; indice_vecinos ++){
											if( cell.vecinos[indice_vecinos] != NULL && vecino == 0){
												rejillas[(indice_rejilla+1)%2][x][y] = *cell.vecinos[indice_vecinos];
												
												if( cell.vecinos[indice_vecinos]->ro > 10 ){
													rejillas[(indice_rejilla+1)%2][x][y] = new Celula(true, cell.vecinos[indice_vecinos]->cct, ROMAX, cell.vecinos[indice_vecinos]->mu, ALPHAMAX);
												}
												indice_vecinos = 8;
																							
												if( Random::get(0.0, 1.0) < nueva_cell.P_repro )
													rejillas[(indice_rejilla+1)%2][x][y].decreaseRo();
												
											}else if(cell.vecinos[indice_vecinos] != NULL)
												vecino --;
										}
									}else{
										rejillas[(indice_rejilla+1)%2][x][y] = new Celula();
										contador --;
									}
								}
							} else { // Tiene cancer y no tiene vecinas cancerígenas
								//Se cálcula las probabilidades de la célula de cada acción según sus parametros
								p_obtenida = Random::get(0.0, 1.0);
								
								if( p_obtenida < cell.P_repro ){
									cell.decreaseRo();
									rejillas[(indice_rejilla+1)%2][x][y] = cell;
								}else{
									rejillas[(indice_rejilla+1)%2][x][y] = new Celula();
									contador --;
								}
									
							}
						} else {
							rejillas[(indice_rejilla+1)%2][x][y] = cell;
						}
					} else { // Si es célula no cancerígena		
						if( cell.n_vecinos > 0 ){ //Si tiene celulas cancerígenas alrededor
							
							p_auxiliar = 0.0;
							p_vecinos = calculoProbabilidadVecinos(cell);	
							for( int vec = 0; vec < p_vecinos.size(); vec ++){
								p_auxiliar += p_vecinos[vec];  }
							
							p_acumulada = p_auxiliar;
							p_obtenida = Random::get(0.0, 1.0);
							if( p_obtenida < p_acumulada ){	
								p_auxiliar = 0.0;
								for(int vec = 0; vec < p_vecinos.size(); vec++){
									p_auxiliar += p_vecinos[vec];
									if( p_obtenida < p_auxiliar ){
										vecino = vec;
										vec = p_vecinos.size();
									}
								}
								for(int indice_vecinos = 0; indice_vecinos < 8; indice_vecinos ++){
									if( cell.vecinos[indice_vecinos] != NULL && vecino == 0){
										nueva_cell = new Celula(*cell.vecinos[indice_vecinos]);
										if( cell.vecinos[indice_vecinos]->ro > 10 ){
											nueva_cell = new Celula(true, cell.vecinos[indice_vecinos]->cct, ROMAX, cell.vecinos[indice_vecinos]->mu, ALPHAMAX);
										}
										rejillas[(indice_rejilla+1)%2][x][y] = nueva_cell;
										contador ++;
										indice_vecinos = 8;									
										if( Random::get(0.0, 1.0) < nueva_cell.P_repro )
											rejillas[(indice_rejilla+1)%2][x][y].decreaseRo();
										
									}else if(cell.vecinos[indice_vecinos] != NULL)
										vecino --;
								}
								if( x == min_x && min_x != 0){
									min_x--;
								}else if( x == max_x && max_x != LONG-1){
									max_x++;
								}else if( y == min_y && min_y != 0){
									min_y--;
								}else if( y == max_y && max_y != LONG-1)
									max_y++;
							}else{
								rejillas[(indice_rejilla+1)%2][x][y] = cell;
							}
						} else { // P = 0; Sin cancer cercano, se mantiene igual
							rejillas[(indice_rejilla+1)%2][x][y] = cell;	
						}
					}
				}
			}
		}
		
		new_extra = high_resolution_clock::now();
		//Cada día indicamos cuantas células hay
		cout << "Día: " << dia << endl;
		cout << "Numero de celulas: "<< contador << ", de las cuales madre: " << celulas_madre.size() << endl;
		
		if (dia % 5 == 0 || dia == 1){
			file = "TABLES/" + to_string(dia) + "dia.txt";
						
			myfile.open( file );
	  		if (myfile.is_open()){
	  			myfile << LONG << "\n";
	  			for( int filas = 0; filas < LONG; filas++ ){
	  				for ( int columnas = 0; columnas < LONG; columnas++){
		  				if( rejillas[indice_rejilla][filas][columnas].cancer ){
			  				myfile << filas << " " << columnas << " " << rejillas[indice_rejilla][filas][columnas].ro << "\n";
			  			}
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

int calc_factorial( int n ){
	return (n == 1 || n == 0) ? 1 : calc_factorial(n - 1) * n;
}

int main(){
	for(int i = 0; i < FACTORIAL.size(); i++ )
		FACTORIAL[i] = calc_factorial(i);

	//Inicializamos el vector con las coordenadas
	Coordenadas coor;
	Random::seed(1);
	// Inicializamos el grid de Células
	vector < vector < vector < Celula > > > rejillas;
	rejillas.push_back( vector < vector <Celula> > (LONG) );
	rejillas.push_back( vector < vector <Celula> > (LONG) );
	for( int i = 0; i < LONG; i++){
		rejillas[0][i] = vector <Celula> (LONG);
		rejillas[1][i] = vector <Celula> (LONG);
	}
		
	//Introducimos en la mitad del grid una célula madre	
	Celula cell(true, 24.0 , 1000.0 , 100.0, 0.0);
	
	coor.x = (LONG - 1)/2;
	coor.y = (LONG - 1)/2;
	
	rejillas[0][coor.x][coor.y] = cell;

	
	//Simulamos el grid cálculando el tiempo que tarda
	auto start = high_resolution_clock::now(); //Declaramos los valores que vamos a usar para el calculo de tiempo
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	int tiempo = 0;
	
	start = high_resolution_clock::now();
	tiempo -= simulacion_cancer(rejillas);
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	tiempo += duration.count();
	
	cout << "Tiempo: " << tiempo << endl;
	
	return 0;
}

