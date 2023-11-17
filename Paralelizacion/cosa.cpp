#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "random.hpp"

using namespace std;
using namespace std::chrono;
using Random = effolkronium::random_static;

int main(){
	vector < int > Vector_aleatorio(7);
	int Indice_aleatorio;
	
	for (int i = 0; i < 7; i++ ){
		Vector_aleatorio[i] = 0;
	}
	
	for (int i = 0; i < 10000; i++ ){
		Vector_aleatorio[Random::get(0, 6)] += 1;
	}
	
	for (int i = 0; i < 7; i++ ){
		cout << Vector_aleatorio[i] << endl;
	}
}
