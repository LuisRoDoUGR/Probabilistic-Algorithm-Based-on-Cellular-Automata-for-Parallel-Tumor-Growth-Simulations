//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
#include <iostream>
#include <math.h>
#include <vector>

#include <cuda_runtime.h>
#include <curand.h>

#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
  printf("Error at %s:%d\n",__FILE__,__LINE__);\
  return EXIT_FAILURE; }} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
  printf("Error at %s:%d\n",__FILE__,__LINE__);\
  return EXIT_FAILURE; }} while(0)  

using namespace std;

static double T = (1/24.0);
//Tamaño del grid
static const int LONG = 32;

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
	void calculoProbabilidades(){
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
	void decreaseRo(){
		ro --;
		if(ro == 0)
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
		//Desde la esquina superior izquierda, recorriendo las 8 posiciones circundantes en sentido horario, se comprueba que no hay células cancerígenas.
		if ( i-1 >= 0  && j-1 >= 0 && i-1 < LONG && j-1 < LONG && rejilla[i-1][j-1].cancer && rejilla[i-1][j-1].ro < 11 ){
	       		vecinos[0] = &rejilla[i-1][j-1];
	       		n_vecinos++;
		}if ( i >= 0   && j-1 >= 0 && i < LONG   && j-1 < LONG && rejilla[i][j-1].cancer && rejilla[i][j-1].ro < 11 ){
	       		vecinos[1] = &rejilla[i][j-1];
	       		n_vecinos++;
		}if ( i+1 >= 0 && j-1 >= 0 && i+1 < LONG && j-1 < LONG && rejilla[i+1][j-1].cancer && rejilla[i+1][j-1].ro < 11 ){
	       		vecinos[2] = &rejilla[i+1][j-1];
	       		n_vecinos++;
		}if ( i+1 >= 0 && j >= 0   && i+1 < LONG && j < LONG   && rejilla[i+1][j].cancer && rejilla[i+1][j].ro < 11 ){
	       		vecinos[3] = &rejilla[i+1][j];
	       		n_vecinos++;
		}if ( i+1 >= 0 && j+1 >= 0 && i+1 < LONG && j+1 < LONG && rejilla[i+1][j+1].cancer && rejilla[i+1][j+1].ro < 11 ){
	       		vecinos[4] = &rejilla[i+1][j+1];
	       		n_vecinos++;
		}if ( i >= 0   && j+1 >= 0 && i < LONG   && j+1 < LONG && rejilla[i][j+1].cancer && rejilla[i][j+1].ro < 11 ){
	       		vecinos[5] = &rejilla[i][j+1];
	       		n_vecinos++;
		}if ( i-1 >= 0 && j+1 >= 0 && i-1 < LONG && j+1 < LONG && rejilla[i-1][j+1].cancer && rejilla[i-1][j+1].ro < 11 ){
	       		vecinos[6] = &rejilla[i-1][j+1];
	       		n_vecinos++;
		}if ( i-1 >= 0 && j >= 0   && i-1 < LONG && j < LONG   && rejilla[i-1][j].cancer && rejilla[i-1][j].ro < 11 ){
	       		vecinos[7] = &rejilla[i-1][j];
	       		n_vecinos++;
		}
	}
};

	__device__ __managed__ int N = 1 << 20;
	__device__ __managed__ float x[12], y[12];
	__device__ __managed__ Celula *a[2][LONG][LONG];
	
// Kernel function to add the elements of two 
// big arrays
__global__ void add(int n, float *x, float *y)
{
	for (int i = 0; i < n; i++)
		y[i] = x[i] + y[i];
}

__global__ void roUp(int n){
	int tid = threadIdx.x;
	
	if(tid < n)
		a[0][blockIdx.x][tid]->ro = blockDim.y;
}

int main()
{
	size_t n = 100;
	size_t i;
	curandGenerator_t gen;
	float *devData, *hostData;
	
	hostData = (float *)calloc(n, sizeof(float));
	
	CUDA_CALL(cudaMalloc((void **)&devData, n*sizeof(float)));
	
	CURAND_CALL(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));
	
	CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, 1234ULL));
	
	CURAND_CALL(curandGenerateUniform(gen, devData, n));
	
	CUDA_CALL(cudaMemcpy(hostData, devData, n*sizeof(float), cudaMemcpyDeviceToHost));
	
	for(i=0; i<n; i++){
		printf("%1.4f ", hostData[i]);
	}
	printf("\n\n");
	
	CURAND_CALL(curandGenerateUniform(gen, devData, n));
	
	CUDA_CALL(cudaMemcpy(hostData, devData, n*sizeof(float), cudaMemcpyDeviceToHost));
	
	for(i=0; i<n; i++){
		printf("%1.4f ", hostData[i]);
	}
	printf("\n\n");
	
	CURAND_CALL(curandDestroyGenerator(gen));
	CUDA_CALL(cudaFree(devData));
	free(hostData);

	//vector <int> a;
	//int N = 1 << 20;
	//float *x, *y;
		
	Celula cell(true, 24.0 , 123.0 , 100.0, 0.0);
	
	vector < vector < vector < Celula > > > rejillas;
	rejillas.push_back( vector < vector <Celula> > (LONG) );
	rejillas.push_back( vector < vector <Celula> > (LONG) );
	for( int i = 0; i < LONG; i++){
		rejillas[0][i] = vector <Celula> (LONG);
		rejillas[1][i] = vector <Celula> (LONG);
	}
	
	for(int i = 0; i < LONG; i++){
		for(int j = 0; j < LONG; j++){
			a[0][i][j] = &rejillas[0][i][j];
			a[1][i][j] = &rejillas[1][i][j];
			//cout << a[0][i][j] << " "; 
		}
		//cout << endl;
	}
	
// Allocate Unified Memory accessible from CPU
// or GPU

	//cudaMallocManaged(&x, N * sizeof(float));
	//cudaMallocManaged(&y, N * sizeof(float));
	
// initialize x and y arrays on the host
	for (int i = 0; i < 12; i++) {
		x[i] = 1.0f;
		y[i] = 2.0f;
	}

// Run kernel on 1M elements on the GPU
// this is almost sequential code: 1 block and
//only 1 thread

	roUp << <16, 32 >> > ( LONG);

// Wait for GPU to finish before accessing 
// on host
	cudaDeviceSynchronize();
	
	for(int i = 0; i < LONG; i++){
		for(int j = 0; j < LONG; j++){
				cout << a[0][i][j]->ro << " ";
		}
		cout << endl;
	}

// Check for errors (all values should be 3.0f)
	float maxError = 0.0f;
	for (int i = 0; i < 12; i++)
		maxError = fmax(maxError, fabs(y[i] - 3.0f));
	cout << "Max error: " << maxError << endl;

// Free memory
	cudaFree(x);
	cudaFree(y);	
    	return 0;
}
/*
CUDA_HOME=/usr/local/cuda
PATH=${CUDA_HOME}/bin:${PATH}
LD_LIBRARY_PATH=${CUDA_HOME}/lib64:$LD_LIBRARY_PATH
*/
