#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <curand.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "random.hpp"

using namespace std;
using namespace std::chrono;
using Random = effolkronium::random_static;


//Probabilidad de reproducción identica o normal
static const double PS = 0.1;
//Tiempo al día, cálculo de migración
static const double T = (1/24.0);

//Nº máximo de reproducciones de una célula
static const double ROMAX = 10.0;
//Probabilidad máxima de muerte
static const double ALPHAMAX = 0.01;

//Tamaño del grid
static const int LONG = 512;
//Nº de días
static const int NUM_DIAS = 25;

//Vector con los valores de n!, para los numeros del 0 al 8
__device__ __managed__ int FACTORIAL[9]; // = [1, 1, 2, 6, 24, 120, 720, 5040, 40320]

//Clase que representan los elementos de la rejilla que usamos para simular el crecimiento tumoral
class Celula {
	
	public:
	//Indica si es cancerigena o no
	bool cancer = false;
	//Cálculo de migración
	double cct;
	//Nº de reproducciones
	double ro;
	//Cálculo de Reproducción
	double mu;
	//Prababilidad de muerte
	double alpha;
	//Numero de vecinos cancerigenos
	int n_vecinos;
	//Referencias a los vecinos cancerigenos o NULL en otro caso
	Celula * vecinos[8];
	
	// Probabilidades de esta célula de morir, migrar o reproducirse respectivamente
	double P_morir;
	double P_migra;
	double P_repro;	
	
	/* Inicializa la célula con los valores que se le pasa por referencia, se establecen los vecinos a 0 y se calcula las probabilidades de morir, migrar o reproducirse
	@param c: Booleano que indica si la célula es cancerígena
	@param cc: Parametro usado para el calculo de la probabilidad de migración
	@param r: Numero de reproducciones de la celula
	@param m: Parametro usado para el cálculo de la probabilidad de reproducción
	@param a: Parametro usado para el cálculo de la probabilidad de morir
	*/
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
	
	/* Inicializa la célula copiando los valores de la célula pasada como parámetro
	@param cell: Célula de la que se va a hacer la copia
	*/
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
	
	/* Inicializa la célula como no cancerígena y con todos los valores a 0
	*/
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
	
	
	/* Usa los valores de cct, mu y alpha para calcular las probabilidades de morir, migrar y reproducirse
	*/	
	__host__ __device__ void calculoProbabilidades(){
		double pd = (24.0/cct)*T; 
		double migrar = (1-pd)*(mu*T);
		
		if( ro <= 0){
			pd = 0.0;
		}
		//Con la suma de los valores de pd, migrar y alpha se calculan las probabilidades de reproducirse, migrar y morir respectivamente
		double p_total = alpha + migrar + pd;
		P_morir = alpha/p_total;
		P_migra = migrar/p_total;
		P_repro = pd/p_total;
	}
	
	/* Reduce el valor de ro en uno y, si ha llegado a 0, se recalculan las probabilidades
	*/
	__device__ void decreaseRo(){
		ro --;
		if(ro == 0)
			calculoProbabilidades();
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
	@return: Coordenadas donde se ha insertado la célula
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

/*_________________________________________________________
	CUDA Declarations
  _________________________________________________________	
*/
	//Grid donde se simula el tumor
	__device__ __managed__ Celula *grid[2][LONG][LONG];
	
/* Kernel para inicializar el generador de números aleatorios
	@param *state: Vector de direcciones donde guardar el estado del generador para cada hilo
*/
__global__ void setup_kernel(curandState *state){
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	/* Each thread gets same seed, a different sequence
	number, no offset */
	curand_init(1, id, 0, &state[id]);
}

/* Kernel que actualiza los vectores de vecinos y el número de vecinos de cada célula
	@param estado: Entero 0 o 1 que nos indica en cual de los dos grid trabajamos
*/
__global__ void actualizarVecinos( int estado ){
	int x = threadIdx.x; //Para cada célula
	int y = blockIdx.x;
	
	//Reiniciamos los valores a 0 y NULL
	grid[estado][x][y]->n_vecinos = 0;
	grid[estado][x][y]->vecinos[0] = NULL;
	grid[estado][x][y]->vecinos[1] = NULL;
	grid[estado][x][y]->vecinos[2] = NULL;
	grid[estado][x][y]->vecinos[3] = NULL;
	grid[estado][x][y]->vecinos[4] = NULL;
	grid[estado][x][y]->vecinos[5] = NULL;
	grid[estado][x][y]->vecinos[6] = NULL;
	grid[estado][x][y]->vecinos[7] = NULL;
	
	//Desde la esquina superior izquierda, recorriendo las 8 posiciones circundantes en sentido horario, se comprueba que no hay células cancerígenas.
	if ( x-1 >= 0  && y-1 >= 0 && x-1 < LONG && y-1 < LONG && grid[estado][x-1][y-1]->cancer && grid[estado][x-1][y-1]->ro < 11 ){
       		grid[estado][x][y]->vecinos[0] = grid[estado][x-1][y-1];
       		grid[estado][x][y]->n_vecinos++;
	}if ( x >= 0  && y-1 >= 0 && x < LONG && y-1 < LONG && grid[estado][x][y-1]->cancer && grid[estado][x][y-1]->ro < 11 ){
       		grid[estado][x][y]->vecinos[1] = grid[estado][x][y-1];
       		grid[estado][x][y]->n_vecinos++;
	}if ( x+1 >= 0  && y-1 >= 0 && x+1 < LONG && y-1 < LONG && grid[estado][x+1][y-1]->cancer && grid[estado][x+1][y-1]->ro < 11 ){
       		grid[estado][x][y]->vecinos[2] = grid[estado][x+1][y-1];
       		grid[estado][x][y]->n_vecinos++;
	}if ( x+1 >= 0  && y >= 0 && x+1 < LONG && y < LONG && grid[estado][x+1][y]->cancer && grid[estado][x+1][y]->ro < 11 ){
       		grid[estado][x][y]->vecinos[3] = grid[estado][x+1][y];
       		grid[estado][x][y]->n_vecinos++;
	}if ( x+1 >= 0  && y+1 >= 0 && x+1 < LONG && y+1 < LONG && grid[estado][x+1][y+1]->cancer && grid[estado][x+1][y+1]->ro < 11 ){
       		grid[estado][x][y]->vecinos[4] = grid[estado][x+1][y+1];
       		grid[estado][x][y]->n_vecinos++;
	}if ( x >= 0  && y+1 >= 0 && x < LONG && y+1 < LONG && grid[estado][x][y+1]->cancer && grid[estado][x][y+1]->ro < 11 ){
       		grid[estado][x][y]->vecinos[5] = grid[estado][x][y+1];
       		grid[estado][x][y]->n_vecinos++;
	}if ( x-1 >= 0  && y+1 >= 0 && x-1 < LONG && y+1 < LONG && grid[estado][x-1][y+1]->cancer && grid[estado][x-1][y+1]->ro < 11 ){
       		grid[estado][x][y]->vecinos[6] = grid[estado][x-1][y+1];
       		grid[estado][x][y]->n_vecinos++;
	}if ( x-1 >= 0  && y >= 0 && x-1 < LONG && y < LONG && grid[estado][x-1][y]->cancer && grid[estado][x-1][y]->ro < 11 ){
       		grid[estado][x][y]->vecinos[7] = grid[estado][x-1][y];
       		grid[estado][x][y]->n_vecinos++;
	}
}


/* Función que añade la probabilidad de que uno de los vecinos no ocurra, y se llama a si misma para añadir el resto de vecinos
	@param elegido: Número del vecino que tiene probabilidad positiva
	@param inicio: Número del vecino desde el que empezar a iterar
	@param vecinos: Número total de vecinos
	@param nivel: Número de elementos múltiplicados en este nivel
	@probabilidad_v: Vector con las probabilidades de migrar y reproducirse de los vecinos
	@return: Probabilidad de que los vecinos a partir del inicio, y sin contar el elegido, no migren ni se reproduzcan en una célula vecina
*/
__device__ double combinaciones(int elegido, int inicio, int vecinos, int nivel, double * probabilidad_v){
	double p = 0;
	for( int j = inicio; j < vecinos; j++){
		if( j != elegido)
			p += (1-probabilidad_v[j])*( FACTORIAL[vecinos - 1 - nivel]*FACTORIAL[nivel] + combinaciones(elegido, j+1, vecinos, nivel+1, probabilidad_v) );
	}
	return p;
}

/* Cálculo de la probabilidad de que los vecinos de la celula migren o se reproduzcan donde la célula se encuentra
	@param cell: Célula en la posición que estamos calculando
	@return: Vector con los valores de las probabilidades de los vecinos de migrar o reproducirse en el sitio de la célula
*/
__device__ void calculoProbabilidadVecinos(Celula &cell, double * probabilidad_f){
	
	double probabilidad_v[8]; //Vector donde se van a guardar las probabilidades de los vecinos
	 //Vector final con las probabilidades múltiplicadas
	int vec = 0;
	
	for( int i = 0 ; i < 8; i++){
		if( cell.vecinos[i] != NULL ){ 
			probabilidad_v[vec] = (cell.vecinos[i]->P_migra + cell.vecinos[i]->P_repro) / (8.0 - cell.vecinos[i]->n_vecinos);
			vec++;
		}
	}
	
	for( int i = 0; i < cell.n_vecinos; i++){
		probabilidad_f[i] = probabilidad_v[i]*(FACTORIAL[cell.n_vecinos-1] + combinaciones( i, 0, cell.n_vecinos, 1, probabilidad_v) ) / (double)cell.n_vecinos;
		 
	}
}

/* Función de transición de cada célula, utilizando las células adyacentes se cálcula el valor en la siguiente iteración
	@param estado: Entero 0 o 1 que nos indica en cual de los dos grid trabajamos
	@param *state: Vector de direcciones donde guardar el estado del generador para cada hilo  
*/
__global__ void funcionTransicion(int estado, curandState *state){
	//Variables auxiliares
	double p_auxiliar, p_acumulada, p_obtenida;
	double p_vecinos[8];
	int vecino;
	
	int x = threadIdx.x; //Para cada célula
	int y = blockIdx.x;
	int rand = threadIdx.x + blockIdx.x * blockDim.x;
	
	if( grid[estado][x][y]->cancer ){ //Si la célula es cancerígena
		if( grid[estado][x][y]->ro < 11 ){ //Si no es célula madre
			if( grid[estado][x][y]->n_vecinos > 0 ){ //Si tiene vecinos cancerígenos
				if( grid[estado][x][y]->n_vecinos == 8){ // P=1 Quiescencia, se mantiene igual
					grid[(estado+1)%2][x][y]->cancer = true;
					grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->cct;
					grid[(estado+1)%2][x][y]->ro = grid[estado][x][y]->ro;
					grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->mu;
					grid[(estado+1)%2][x][y]->alpha = grid[estado][x][y]->alpha;
					grid[(estado+1)%2][x][y]->calculoProbabilidades();
				} else {
					//Se cálcula las probabilidades de cada vecino
					p_auxiliar = 0.0;
					calculoProbabilidadVecinos(*grid[estado][x][y], p_vecinos);	
					for( int vec = 0; vec < grid[estado][x][y]->n_vecinos; vec ++){
						p_auxiliar += p_vecinos[vec]; }
					
					//Se calcula la probabilidad total de la posición de que haya una célula en la iteración siguiente
					p_acumulada = grid[estado][x][y]->P_repro + (1-grid[estado][x][y]->P_repro)*p_auxiliar;
					p_obtenida = curand_uniform(&state[rand]);
					
					if( p_obtenida < grid[estado][x][y]->P_repro ){ //Si la célula actual se reproduce
						grid[estado][x][y]->decreaseRo();
						grid[(estado+1)%2][x][y]->cancer = true;
						grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->cct;
						grid[(estado+1)%2][x][y]->ro = grid[estado][x][y]->ro;
						grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->mu;
						grid[(estado+1)%2][x][y]->alpha = grid[estado][x][y]->alpha;
						grid[(estado+1)%2][x][y]->calculoProbabilidades();
					} else if( p_obtenida < p_acumulada ){ //Si otra célula ocupa este sitio
						
						p_auxiliar = grid[estado][x][y]->P_repro;
						for(int vec = 0; vec < grid[estado][x][y]->n_vecinos; vec++){ //Vemos que vecino ocupa este sitio
							p_auxiliar += (1-grid[estado][x][y]->P_repro)*p_vecinos[vec];
							if( p_obtenida < p_auxiliar ){
								vecino = vec;
								vec = grid[estado][x][y]->n_vecinos;
							}
						} //Copiamos la célula 
						for(int indice_vecinos = 0; indice_vecinos < 8; indice_vecinos ++){
							if( grid[estado][x][y]->vecinos[indice_vecinos] != NULL && vecino == 0){
								grid[(estado+1)%2][x][y]->cancer = grid[estado][x][y]->vecinos[indice_vecinos]->cancer;
								grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->vecinos[indice_vecinos]->cct;
								grid[(estado+1)%2][x][y]->ro = grid[estado][x][y]->vecinos[indice_vecinos]->ro;
								grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->vecinos[indice_vecinos]->mu;
								grid[(estado+1)%2][x][y]->alpha = grid[estado][x][y]->vecinos[indice_vecinos]->alpha;
								grid[(estado+1)%2][x][y]->calculoProbabilidades();
								
								if( grid[estado][x][y]->vecinos[indice_vecinos]->ro > 10 ){
									grid[(estado+1)%2][x][y]->cancer = true;
									grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->vecinos[indice_vecinos]->cct;
									grid[(estado+1)%2][x][y]->ro = ROMAX;
									grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->vecinos[indice_vecinos]->mu;
									grid[(estado+1)%2][x][y]->alpha = ALPHAMAX;
									grid[(estado+1)%2][x][y]->calculoProbabilidades();
								}
								indice_vecinos = 8;
																			
								if( curand_uniform(&state[rand]) < grid[(estado+1)%2][x][y]->P_repro ) //Posibilidad de quitarlo
									grid[(estado+1)%2][x][y]->decreaseRo();
								
							}else if(grid[estado][x][y]->vecinos[indice_vecinos] != NULL)
								vecino --;
						}
					}else{
						grid[(estado+1)%2][x][y]->cancer = false;
						grid[(estado+1)%2][x][y]->cct = 0.0;
						grid[(estado+1)%2][x][y]->ro = 0.0;
						grid[(estado+1)%2][x][y]->mu = 0.0;
						grid[(estado+1)%2][x][y]->alpha = 0.0;
						grid[(estado+1)%2][x][y]->calculoProbabilidades();
					}
				}
			} else { // Tiene cancer y no tiene vecinas cancerígenas
				//Se cálcula las probabilidades de la célula de cada acción según sus parametros
				p_obtenida = curand_uniform(&state[rand]);
				
				if( p_obtenida < grid[estado][x][y]->P_repro ){
					grid[estado][x][y]->decreaseRo();
					grid[(estado+1)%2][x][y]->cancer = true;
					grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->cct;
					grid[(estado+1)%2][x][y]->ro = grid[estado][x][y]->ro;
					grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->mu;
					grid[(estado+1)%2][x][y]->alpha = grid[estado][x][y]->alpha;
					grid[(estado+1)%2][x][y]->calculoProbabilidades();
				}else{
					grid[(estado+1)%2][x][y]->cancer = false;
					grid[(estado+1)%2][x][y]->cct = 0.0;
					grid[(estado+1)%2][x][y]->ro = 0.0;
					grid[(estado+1)%2][x][y]->mu = 0.0;
					grid[(estado+1)%2][x][y]->alpha = 0.0;
					grid[(estado+1)%2][x][y]->calculoProbabilidades();
				}
					
			}
		} else {
			// Comprobar si alrededor hay celulas madre, y si no -> hacer copia forzada. 
			// Problema: Dos celulas madre juntas se convierten en una
			
			// 
			
			grid[(estado+1)%2][x][y]->cancer = true;
			grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->cct;
			grid[(estado+1)%2][x][y]->ro = grid[estado][x][y]->ro;
			grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->mu;
			grid[(estado+1)%2][x][y]->alpha = grid[estado][x][y]->alpha;
			grid[(estado+1)%2][x][y]->calculoProbabilidades();
		}
	} else { // Si es célula no cancerígena
		if( grid[estado][x][y]->n_vecinos > 0 ){ //Si tiene celulas cancerígenas alrededor
			//Se cálcula las probabilidades de cada vecino
			p_auxiliar = 0.0;
			calculoProbabilidadVecinos(*grid[estado][x][y], p_vecinos);	
			for( int vec = 0; vec < grid[estado][x][y]->n_vecinos; vec ++){
				p_auxiliar += p_vecinos[vec];  }
			//Se calcula la probabilidad total de la posición de que haya una célula en la iteración siguiente
			p_acumulada = p_auxiliar;
			p_obtenida = curand_uniform(&state[rand]);
			if( p_obtenida < p_acumulada ){//Si otra célula ocupa este sitio
				p_auxiliar = 0.0;
				for(int vec = 0; vec < grid[estado][x][y]->n_vecinos; vec++){//Vemos que vecino ocupa este sitio
					p_auxiliar += p_vecinos[vec];
					if( p_obtenida < p_auxiliar ){
						vecino = vec;
						vec = grid[estado][x][y]->n_vecinos;
					}
				} //Copiamos la célula 
				for(int indice_vecinos = 0; indice_vecinos < 8; indice_vecinos ++){
					if( grid[estado][x][y]->vecinos[indice_vecinos] != NULL && vecino == 0){
						grid[(estado+1)%2][x][y]->cancer = grid[estado][x][y]->vecinos[indice_vecinos]->cancer;
						grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->vecinos[indice_vecinos]->cct;
						grid[(estado+1)%2][x][y]->ro = grid[estado][x][y]->vecinos[indice_vecinos]->ro;
						grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->vecinos[indice_vecinos]->mu;
						grid[(estado+1)%2][x][y]->alpha = grid[estado][x][y]->vecinos[indice_vecinos]->alpha;
						grid[(estado+1)%2][x][y]->calculoProbabilidades();
						
						if( grid[estado][x][y]->vecinos[indice_vecinos]->ro > 10 ){
							grid[(estado+1)%2][x][y]->cancer = true;
							grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->vecinos[indice_vecinos]->cct;
							grid[(estado+1)%2][x][y]->ro = ROMAX;
							grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->vecinos[indice_vecinos]->mu;
							grid[(estado+1)%2][x][y]->alpha = ALPHAMAX;
							grid[(estado+1)%2][x][y]->calculoProbabilidades();
						}
						indice_vecinos = 8;									
						if( curand_uniform(&state[rand]) < grid[(estado+1)%2][x][y]->P_repro )
							grid[(estado+1)%2][x][y]->decreaseRo();
						
					}else if(grid[estado][x][y]->vecinos[indice_vecinos] != NULL)
						vecino --;
				}
			}else{
				grid[(estado+1)%2][x][y]->cancer = false;
				grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->cct;
				grid[(estado+1)%2][x][y]->ro = grid[estado][x][y]->ro;
				grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->mu;
				grid[(estado+1)%2][x][y]->alpha = grid[estado][x][y]->alpha;
				grid[(estado+1)%2][x][y]->calculoProbabilidades();
			}
		} else { // P = 0; Sin cancer cercano, se mantiene igual
			grid[(estado+1)%2][x][y]->cancer = false;
			grid[(estado+1)%2][x][y]->cct = grid[estado][x][y]->cct;
			grid[(estado+1)%2][x][y]->ro = grid[estado][x][y]->ro;
			grid[(estado+1)%2][x][y]->mu = grid[estado][x][y]->mu;
			grid[(estado+1)%2][x][y]->alpha = grid[estado][x][y]->alpha;
			grid[(estado+1)%2][x][y]->calculoProbabilidades();	
		}
	}	
}




/* Simulación durante 50 días, con 24 pasos por día del crecimiento de un tumor en un grid.
	@param matriz: Vector con todas las posibles coordenadas que tiene el grid
	@param rejilla: Grid de células donde se va a simular el crecimiento del tumor
*/
int simulacion_cancer(vector < vector< vector <Celula> > > &rejillas){
	int indice_rejilla = 0;
	
	// Declaración de variables auxiliares
	vector <Coordenadas> celulas_madre;
	vector <Coordenadas> futuras;
	Coordenadas coor;
	coor.x = (LONG - 1)/2;
	coor.y = (LONG - 1)/2;
	celulas_madre.push_back( coor );
	
	int vecino, contador = 1, pasos = 24;
	Celula cell, nueva_cell;
	vector<int> cas_libres;
	vector<double> p_vecinos;
	double  p_obtenida;
	ofstream myfile;
	string file;
	
	auto new_extra = high_resolution_clock::now(); //Declaramos los valores que vamos a usar para el calculo de tiempo
	auto fin_extra = high_resolution_clock::now();
	auto extra = duration_cast<chrono::milliseconds>(fin_extra - fin_extra);
	int time = 0;
	
	const unsigned int threadsPerBlock = LONG;
	const unsigned int blockCount = LONG;
	const unsigned int totalThreads = threadsPerBlock * blockCount;
	curandState *devStates;
	cudaMalloc((void **)&devStates, totalThreads * sizeof(curandState));
	
	setup_kernel<<<blockCount, threadsPerBlock>>>(devStates);
	
	//Durante 50 días
	for( int dia = 1; dia < NUM_DIAS; dia++){
		//En cada día 24 pasos
		for (int paso = 0; paso < pasos; paso++){
			
			//Actualizamos el indice de rejillas
			indice_rejilla = paso % 2;
			
			cudaDeviceSynchronize();
			//Primero se actualizan todas las células madres
			for (int madre = 0; madre < celulas_madre.size(); madre++){
				//Guardamos la célula actual y vemos si tiene vecinos
				cell = rejillas[indice_rejilla][celulas_madre[madre].x][celulas_madre[madre].y];
				cas_libres = casillas_libres( celulas_madre[madre].x, celulas_madre[madre].y, rejillas[indice_rejilla]);
				
				if( cas_libres.size() != 0){
					//Si hay huecos alrededor, se calcula la probabilidad y a que vecino se va a realizar la acción
					p_obtenida = Random::get(0.0, 1.0);
					vecino = Random::get(0, int(cas_libres.size()-1));
					
					if( p_obtenida < cell.P_migra ){ //Si migra
						//Se introduce la célula en el vecino libre
						coor = introducir_en_casilla( celulas_madre[madre].x, celulas_madre[madre].y, cas_libres[vecino], cell, rejillas[indice_rejilla]);
						futuras.push_back(coor);
						//Se limpia la casilla actual
						rejillas[indice_rejilla][celulas_madre[madre].x][celulas_madre[madre].y] = new Celula();
						
					}else{ //Si se reproduce
						futuras.push_back(celulas_madre[madre]);	
						if( Random::get(0.0, 1.0) < PS){//Se comprueba si va a ser copia o hijo normal
							coor = introducir_en_casilla( celulas_madre[madre].x, celulas_madre[madre].y, cas_libres[vecino], cell, rejillas[indice_rejilla]);
							futuras.push_back( coor );
							contador ++;
						}else{
							nueva_cell = new Celula(true, cell.cct, ROMAX, cell.mu, ALPHAMAX);
							introducir_en_casilla( celulas_madre[madre].x, celulas_madre[madre].y, cas_libres[vecino], nueva_cell, rejillas[indice_rejilla]);
						}
					}
					
				} else {
					futuras.push_back(celulas_madre[madre]);
				}
			} //Se actualiza el vector de celulas madre
			
			celulas_madre = futuras;
			futuras.clear();
						
			actualizarVecinos <<< blockCount, threadsPerBlock>>> (indice_rejilla);			
					
			funcionTransicion <<< blockCount, threadsPerBlock>>> (indice_rejilla, devStates);
			
		}
		contador = 0;
		for(int x = 0; x < LONG ; x++){ // Se actualizan los vecinos de cada célula
			for(int y = 0; y < LONG; y++){
				if(rejillas[indice_rejilla][x][y].cancer)
					contador ++;
			}
		}
		
		new_extra = high_resolution_clock::now();
		//Cada día indicamos cuantas células hay
		cout << "Día: " << dia << endl;
		cout << "Numero de celulas: "<< contador << endl;
		
		
		if (dia % 4 == 0 || dia == 1){
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
	cudaFree(devStates);
	return time;
}





int calc_factorial( int n ){
	return (n == 1 || n == 0) ? 1 : calc_factorial(n - 1) * n;
}

int main(){
	for(int i = 0; i < 9; i++ )
		FACTORIAL[i] = calc_factorial(i);

	//Inicializamos 
	Random::seed(1);
	// Inicializamos el grid de Células
	vector < vector < vector < Celula > > > rejillas;
	rejillas.push_back( vector < vector <Celula> > (LONG) );
	rejillas.push_back( vector < vector <Celula> > (LONG) );
	for( int i = 0; i < LONG; i++){
		rejillas[0][i] = vector <Celula> (LONG);
		rejillas[1][i] = vector <Celula> (LONG);
	}
	
	for(int i = 0; i < LONG; i++){
		for(int j = 0; j < LONG; j++){
			grid[0][i][j] = &rejillas[0][i][j];
			grid[1][i][j] = &rejillas[1][i][j];
		}
	}
		
	//Introducimos en la mitad del grid una célula madre	
	Celula cell(true, 24.0 , 1000.0 , 100.0, 0.0);
		
	rejillas[0][(LONG - 1)/2][(LONG - 1)/2] = cell;

	
	//Simulamos el grid cálculando el tiempo que tarda
	auto start = high_resolution_clock::now(); //Declaramos los valores que vamos a usar para el calculo de tiempo
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	int tiempo = 0;
	
	start = high_resolution_clock::now();
	//
	tiempo -= simulacion_cancer(rejillas);
	//
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	tiempo += duration.count();
	
	cout << "Tiempo: " << tiempo << endl;
	
	return 0;
}
