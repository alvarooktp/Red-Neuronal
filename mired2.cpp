// g++ mired2.cpp -o test && ./test mnist.txt
#include <iostream> //cout
#include <cstdlib> //exit
#include <sstream> //stringstream
#include <string>
#include <fstream> // ifstream
#include <ctime>
#include <cmath>
using namespace std;
const int N = 784;//input
const int DATA = 600;
const int M = 10;//output
const double eps = 0.000000001;//se puede asi?

//Funciones prototipos
void getFromFile(char * archivo,double x[][N],double *clasif );
void printPosX(double x[][N],double *clasif ,int pos);
void getRands(int *posicion);
void randTheta(double **thetas, int long1, int long2);
double getH(int *posicionRands,double *theta, double *entrada, int tamEntrada);
void ForwardPropagation(double **theta, double *entrada, int tamEntrada, double *neuronasSigCapa,int numNeuronas);


int main(int argc, char *argv[]) {
	if (argc < 2){
    cerr << "Faltan Parámetros:\n\n\tUso: test <Archivo a leer>\n\n";
    exit (EXIT_FAILURE);
  }

	srand(time(NULL));
	double x[DATA][N];

	double clasif[DATA];
	getFromFile(argv[1],x,clasif);//proceder a leer del archivo

  //Random thetas
	double **thetas1 = NULL;
	double **thetas2 = NULL;
	double **thetas3 = NULL;

  //theta[numero de la neurona][elegir la theta]
	// const int N = 784;//input
	// const int DATA = 900;
	// const int M = 10;//out
	static int r = pow(N/M, 1/3.0);
	// cout<<"r = "<<r<<" r2 = "<<pow(double(N/M), 1/3.0)<<endl;
	thetas1 = new double *[N];//datos->capa oculta 1
	thetas2 = new double *[M*r*r+1];//capa oculta 1 -> capa oculta 2
	thetas3 = new double *[M*r+1];//capa oculta 2 -> salida

	for (int pos = 0; pos < N;pos++)
		thetas1[pos] = new double [M*r*r +1 ];

	for (int pos = 0; pos < M*r*r +1 ;pos++)
		thetas2[pos] = new double [M*r +1 ];

	for (int pos = 0; pos < M*r +1 ;pos++)
		thetas3[pos] = new double [M];

	// Neuronas
	double *neuronas1 = NULL;
	double *neuronas2 = NULL;
	double *salida = NULL;

	neuronas1 = new double [M*r*r]; //de la capa 1
	neuronas2 = new double [N*r]; //de la capa 2
	salida = new double [N];

  //inicializa thetas random
	randTheta(thetas1, N, M*r*r +1);
	randTheta(thetas2, M*r*r + 1, M*r +1);
	randTheta(thetas3, M*r + 1, M);
	ForwardPropagation(thetas1, (x[DATA]),DATA, neuronas1, M*r*r);

	for (int pos = 0; pos < N;pos++)
		 delete[] thetas1[pos];
	delete[] thetas1;

	for (int pos = 0; pos < M*r*r+1;pos++)
		delete[] thetas2[pos];
	delete[] thetas2;

	for (int pos = 0; pos < M*r +1;pos++)
		delete[] thetas3[pos];
	delete[] thetas3;
	return 0;
}


//Fucniones usadas
void getFromFile(char * archivo,double x[][N],double *clasif ){
  string line;
  ifstream myfile (archivo);
	int j=0,i=0;
  if (myfile.is_open()){
		while ( getline (myfile,line) && i < DATA){
      std::stringstream   data(line);
			for (j = 0;j<N;j++){
				data >> x[i][j];
				x[i][j]/=255;
			}
			data>>clasif[i];
      i++;
    }
    myfile.close();
		cout<<"Total de reglones leidos = "<<i<<endl;
  }
  else cout << "Unable to open file";
}

void printPosX(double x[][N],double *clasif ,int pos){
	const int miPos = pos -1;
	for (int j = 0; j<N;j++)
		//cout<<x[miPos][j]<<", ";
	cout<<clasif[miPos]<<endl;
}

void randTheta(double **thetas, int long1, int long2){
	int i,j;
	for(i=0; i<long1; i++){
		for(j=0; j<long2; j++){
			thetas[i][j] = rand()%(101-1)*2*eps-eps;
		}
	}
	cout<<i*j<<endl;
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
double getH(int *posicionRands,double *theta, double *entrada, int tamEntrada){
  double H = 0.0;
	double Z = theta[0];

	for (int pos = 0; pos < tamEntrada ;pos++){
		Z += entrada[posicionRands[pos]] * theta [posicionRands[pos+1]];
		cout<<"entreeeeeeeeeeee "<< endl;
	}
  H = 1.0/(1.0+exp(-Z));
  cout<<"H = "<<H<<" tamEntrada "<<tamEntrada<< endl;
  return H;
}
// ForwardPropagation() funciona para una unica capa indicado por neuronasSigCapa
void ForwardPropagation(double **theta, double *entrada, int tamEntrada, double *neuronasSigCapa,int numNeuronas){
	int *posicionRands = NULL;
	posicionRands = new int[tamEntrada];
	getRands(posicionRands);

	double *ptrTheta = NULL;//*ptrEntrada = NULL;
	int capaActualTheta = 0;

	for (int pos = 0; pos < numNeuronas;pos++){
		ptrTheta = theta[capaActualTheta];
		cout<<"entre a forward "<<pos<<endl;
		// ptrEntrada = entrada[capaActualTheta];
		neuronasSigCapa[pos] = getH(posicionRands,ptrTheta,entrada,tamEntrada);
		capaActualTheta++;
	}
	// Liberación de memoria
	delete[] posicionRands;
}
// ***** Cost() calcula el costo según la definición derivada
// Se apolla de getH() al final no fue utilizada pero se deja por si en un
// futuro sea requerida
double Cost(double *theta,double x[][DATA],int i,double clasif){
  // double costo = 0.0;
  // double h = getH(theta,x,i);
  // if(clasif == 1){
  //   costo = - log(h + eps);
  // }else{
  //   costo = - log(1 - h + eps);
  // }
  // // cout<<"Costo = "<<costo<<endl;
  // return costo;
}

// ***** getGradient() utilizada a Cost()
// Al final ésta función no fue utilizada pero se deja por si en un futuro sea
// requerida
double getGradient(double *theta,double x[][DATA],double *clasif){

  double J = 0.0, suma =0.0;
  // size_t m = sizeof(x)/sizeof(x[0][0])/2;
  for (int i = 0; i < DATA; i++){
    // cout<<"Entre "<<endl;
    int clasifACtual = clasif[i];
    suma += Cost(theta,x,i,clasifACtual);
  }
  J = suma/double(DATA);
  // cout<<"Sali "<<endl;
  return J;
}
