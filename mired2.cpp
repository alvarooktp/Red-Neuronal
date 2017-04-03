// g++ mired.cpp -o test && ./test mnist.txt
#include <iostream> //cout
#include <cstdlib> //exit
#include <sstream> //stringstream
#include <string>
#include <fstream> // ifstream
#include <ctime>
#include <cmath>
using namespace std;
#define N 784 //input
#define DATA 900
#define M 10 //output
static double eps .000000001//se puede asi?

//Funciones prototipos
void getFromFile(char * archivo,int x[][N],double *clasif );
void printPosX(int x[][N],double *clasif ,int pos);
void getRands(int *posicion);
void randTheta(double *thetas, int long1, int long2);


int main(int argc, char *argv[]) {
	if (argc < 2){
    cerr << "Faltan Parámetros:\n\n\tUso: test <Archivo a leer>\n\n";
    exit (EXIT_FAILURE);
  }

	srand(time(NULL));
	int x[DATA][N];

	double clasif[DATA];
	getFromFile(argv[1],x,clasif);//proceder a leer del archivo

  //Random thetas
	double new thetas1 = NULL;
	double new thetas2 = NULL;
	double new thetas3 = NULL;
  
  //theta[numero de la neurona][elegir la theta]
	int r = pow(N/M, 1/3);
	thetas1 = new thetas1[N][M*r*r];//datos->capa oculta 1
	thetas2 = new thetas2[M*r*r][M*r]//capa oculta 1 -> capa oculta 2
	thetas3 = new thetas3[M*r][M];//capa oculta 2 -> salida

  //inicializa thetas random
	randTheta(thetas1, N, M*r*r);
	randTheta(thetas2, M*r*r, M*r);
	randTheta(thetas3, M*r, M);


	int *posicionRands = NULL;
	posicionRands = new int[DATA];
	getRands(posicionRands);

	delete[] posicionRands;
	delete[] thetas1;
	delete[] thetas2;
	delete[] thetas3;
	return 0;
}


//Fucniones usadas
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
		//cout<<x[miPos][j]<<", ";
	cout<<clasif[miPos]<<endl;
}

void randTheta(double *thetas, int long1, int long2){
	for(int i=0; i<long1; i++){
		for(int j=0; j<long2; j++){
			thetas[i][j] = rand()%(101-1)*2*eps-eps;
		}
	}
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

	delete[] arreglo;
}


