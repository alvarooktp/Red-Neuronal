// g++ mired.cpp -o test && ./test mnist.txt
#include <iostream> //cout
#include <cstdlib> //exit
#include<sstream> //stringstream
#include <string>
#include <fstream> // ifstream
#include <ctime>
#include <cmath>
using namespace std;
#define N 784
#define DATA 900
void getFromFile(char * archivo,int x[][N],double *clasif ){
  string line;
  ifstream myfile (archivo);
	int j=0,i=0;
  if (myfile.is_open()){
		while ( getline (myfile,line) && i < DATA){
      std::stringstream   data(line);
			for (j = 0;j<N;j++)
				data >> x[i][j];
			data>>clasif[i];
      i++;
    }
    myfile.close();
		cout<<"Total de reglones leidos = "<<i<<endl;
  }
  else cout << "Unable to open file";
}
void printPosX(int x[][N],double *clasif ,int pos){
	const int miPos = pos -1;
	for (int j = 0; j<N;j++)
		cout<<x[miPos][j]<<", ";
	cout<<clasif[miPos]<<endl;
}

void getRands(int *posicion){
	int*arreglo = NULL;
	arreglo = new int[DATA];
	for(int i=0; i<DATA; i++){
		arreglo[i]=i;
		posicion[i]=-1;
	}

	for(int i=0; i<DATA; i++){
		int num=rand()%(DATA);
		if(posicion[num]==-1){
			posicion[num] = arreglo[i];
		}else{
			int k = i, potencia=1, asignado = 0;

			if(k%2==0){
				potencia=2;}
				do{
						k = k + pow(-1,potencia);
						if(k==-1)
							k=DATA-1;
						if(k==DATA+1)
							k=0;
						if(posicion[k]==-1){
							asignado = 1;
							posicion[k] = arreglo[i];
						}
				}while(asignado==0);
		}
	}
	// for(int i = 0; i<DATA;i++)
	// 	cout<<posicion[i]<<endl;
	delete[] arreglo;
}
int main(int argc, char *argv[]) {
	if (argc < 2){
    cerr << "Faltan Parámetros:\n\n\tUso: examen3v2 <Archivo a leer>\n\n";
    exit (EXIT_FAILURE);
  }
	srand(time(NULL));
  int x[DATA][N];

  double clasif[DATA];
  getFromFile(argv[1],x,clasif);//proceder a leer del archivo

	// printPosX(x,clasif,3);

	int *posicionRands = NULL;
	posicionRands = new int[DATA];
	getRands(posicionRands);

	delete[] posicionRands;
	return 0;
}
