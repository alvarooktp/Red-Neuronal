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
const int DATA = 6;
const int M = 10;//output
const double eps = 0.000000001;//se puede asi?
const int EPOCAS = 200;
const float ALFA = 0.5;

//Funciones prototipos
void getFromFile(char * archivo,double **x,double *clasif );
void printPosX(double **x,double *clasif ,int pos);
void getRands(int *posicion,int tam);
void randTheta(double **thetas, int long1, int long2);
double getH(double **theta, double *entrada,int capaActual);
void ForwardCapa(double **theta, double *entrada, double *neuronasSigCapa,int numNeuronas);
// void ForwardPropagation(double **theta1,double **theta2,double **theta3,int *cantidadNeuronasCapa,double **entrada,int *posicionRands);
double Cost(double **theta,double *entrada,int capaActualTheta,double clasif);
void ParameterVector(double **theta, double *entrada, double *clasif, int numNeuAct, int numNeuPrev);
void BackwardCapa(void);


int main(int argc, char *argv[]) {
	if (argc < 2){
    cerr << "Faltan Parámetros:\n\n\tUso: test <Archivo a leer>\n\n";
    exit (EXIT_FAILURE);
  }

	srand(time(NULL));

	double **x = NULL;//[DATA][N];
	x = new double *[DATA];
	for (int pos = 0; pos < DATA;pos++)
		x[pos] = new double [N];

	double clasif[DATA];

	getFromFile(argv[1],x,clasif);//proceder a leer del archivo

  //Random thetas
	double **thetas1 = NULL;
	double **thetas2 = NULL;
	double **thetas3 = NULL;
	double **granDelta1 = NULL;
	double **granDelta2 = NULL;
	double **granDelta3 = NULL;

  //theta[numero de la neurona][elegir la theta]
	// const int N = 784;//input
	// const int DATA = 900;
	// const int M = 10;//out
	static int r = pow(N/M, 1/3.0);
	// cout<<"r = "<<r<<" r2 = "<<pow(double(N/M), 1/3.0)<<endl;
	thetas1 = new double *[N];//datos->capa oculta 1
	thetas2 = new double *[M*r*r+1];//capa oculta 1 -> capa oculta 2
	thetas3 = new double *[M*r+1];//capa oculta 2 -> salida
	granDelta1 = new double *[N];//datos->capa oculta 1
	granDelta2 = new double *[M*r*r+1];//capa oculta 1 -> capa oculta 2
	granDelta3 = new double *[M*r+1];//capa oculta 2 -> salida


	for (int pos = 0; pos < N;pos++){
		thetas1[pos] = new double [M*r*r +1 ];
		granDelta1[pos] = new double [M*r*r +1 ];
	}

	for (int pos = 0; pos < M*r*r +1 ;pos++){
		thetas2[pos] = new double [M*r +1 ];
		granDelta2[pos] = new double [M*r +1 ];
	}

	for (int pos = 0; pos < M*r +1 ;pos++){
		thetas3[pos] = new double [M];
		granDelta3[pos] = new double [M];}

	// Neuronas
	double *neuronas1 = NULL;
	double *neuronas2 = NULL;
	double *salida = NULL;
	// for (int r=1;r<DATA+1;r++)
	// 	printPosX(x,clasif ,r);

	neuronas1 = new double [M*r*r]; //de la capa 1
	neuronas2 = new double [M*r]; //de la capa 2
	salida = new double [M];


  //inicializa thetas random
	randTheta(thetas1, N, M*r*r +1);


	//Inicializa las grandes deltas
	for(int k=0; k<N; k++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M*r*r; j++){
			granDelta1[j][k] = 0;
		}
	}

	for(int k=0; k<M*r*r; k++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M*r; j++){
			granDelta2[j][k] = 0;
		}
	}

	for(int k=0; k<M*r; k++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M; j++){
			granDelta3[j][k] = 0;	
		}
	}






	// for(int pos = 0; pos<N;pos++)
	// 	cout<<"thetas1[pos][0] "<<thetas1[pos][0]<<endl;
	cout<<"Cantidad de elementos en thetas1[0] "<<N<<endl;
	randTheta(thetas2, M*r*r + 1, M*r +1);
	randTheta(thetas3, M*r + 1, M);

	int *posicionRands = NULL;
	posicionRands = new int[DATA];
	getRands(posicionRands,DATA);
	cout<<"posicionRands[0] "<<posicionRands[0]<<endl;

//For de las epocas
for(int i = 0; i<EPOCAS; i++){
	for(int j = 0; j < DATA; j++){
		ForwardCapa(thetas1, x[posicionRands[j]], neuronas1, M*r*r);
		ForwardCapa(thetas2, x[posicionRands[j]], neuronas2,M*r);
		ForwardCapa(thetas3, x[posicionRands[j]], salida, M);
		backPropagation(theta1, theta2, theta3,  entrada, neuronas1, neuronas2, salida, posicionRands[j], clasif, granDelta2, granDelta3, granDelta4);
	}

	//Actualizacion de valores de theta
	for(int k=0; k<N; k++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M*r*r; j++){
			thetas1[j][k] = granDelta1[j][k]/DATA;
			granDelta1[j][k] = 0;
		}
	}

	for(int k=0; k<M*r*r; k++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M*r; j++){
			granDelta2[j][k] = granDelta2[j][k]/DATA;
			granDelta2[j][k] = 0;
		}
	}

	for(int k=0; k<M*r; k++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M; j++){
			granDelta3[j][k] = granDelta3[j][k]/DATA;	
			granDelta3[j][k] = 0;	
		}
	}
}





	for (int pos = 0; pos < M; pos++){
		cout<<"Salida "<<salida[pos]<<endl;
	}

// ForwardCapa(double **theta, double **entrada, int tamEntrada, double *neuronasSigCapa,int numNeuronas);
	for (int pos = 0; pos < DATA;pos++)
		 delete[] x[pos];
	delete[] x;

	for (int pos = 0; pos < N;pos++){
		delete[] thetas1[pos];
		delete[] granDelta1;
	}
	delete[] thetas1;

	for (int pos = 0; pos < M*r*r+1;pos++){
		delete[] thetas2[pos];
		delete[] granDelta2;
	}
	delete[] thetas2;

	for (int pos = 0; pos < M*r +1;pos++){
		delete[] thetas3[pos];
		delete[] granDelta3;
	}
	delete[] thetas3;

	delete[] posicionRands;
	delete[] neuronas1;
	delete[] neuronas2;
	delete[] salida;
	
	
	

	
	return 0;
}
//Funciones usadas

void getFromFile(char * archivo,double **x,double *clasif ){
  string line;
  ifstream myfile (archivo);
	int j=0,i=0;
  if (myfile.is_open()){
		while ( getline (myfile,line) && i < DATA){
      std::stringstream   data(line);
			for (j = 0;j<N;j++){
				data >> x[i][j];
				// x[i][j]/=255;
			}
			data>>clasif[i];
      i++;
    }
    myfile.close();
		cout<<"Total de reglones leidos = "<<i<<endl;
  }
  else cout << "Unable to open file";
}

void printPosX(double **x,double *clasif ,int pos){
	const int miPos = pos -1;
	for (int j = 0; j<N;j++)
		cout<<x[miPos][j]<<", ";
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

void getRands(int *posicion,int tam){
	int*arreglo = NULL;
	// tam = tam +1;
	arreglo = new int[tam];
	for(int i=0; i<tam; i++){
		arreglo[i]=i;
		posicion[i]=-1;
	}

	for(int i=0; i<tam; i++){
		int num=rand()%(tam);
		if(posicion[num]==-1){
			posicion[num] = arreglo[i];
		}else{
			int k = i, potencia=1, asignado = 0;

			if(k%2==0){
				potencia=2;}
				do{
						k = k + pow(-1,potencia);
						if(k==-1)
							k=tam-1;
						if(k==tam+1)
							k=0;
						if(posicion[k]==-1){
							asignado = 1;
							posicion[k] = arreglo[i];
						}
				}while(asignado==0);
		}
	}
	for(int i = 0; i<tam;i++){
	cout<<posicion[i]<<endl<<endl;}
	delete[] arreglo;
}
double getH(double **theta, double *entrada, int neuronaActual){
  double H = 0.0;
	double Z = theta[0][neuronaActual];

	for (int position = 0; position < N;position++){
		// cout<<"  "<<endl;
		// Z += entrada[capaActual][0] ;
		// cout<<"despues de entrada "<<position<<endl;
		Z+= theta [position][neuronaActual];//x[mantener][iterar]
		// cout<<"despues DE THETA =  "<<pos<<endl;
	}
  H = 1.0/(1.0+exp(-Z));
  // cout<<"H = "<<H<<" tamEntrada "<<tamEntrada<< endl;
  return H;
}
// void ForwardPropagation(double **theta1,double **theta2,double **theta3,int *cantidadNeuronasCapa,double **entrada,int *posicionRands, double **salida){
// 	int thetaActual = 0;
// 	int capaActual, ;
// 	for(thetaActual = 0; thetaActual<3;thetaActual++){
// 		// ForwardCapa(double **theta, double *entrada, double *neuronasSigCapa,int numNeuronas)
// 		// ForwardCapa(double **theta, double *entrada, double *neuronasSigCapa,int numNeuronas)
// 		// ForwardCapa(double **theta, double *entrada, double *neuronasSigCapa,int numNeuronas)
// 	}
// }
// ForwardCapa() funciona para una unica capa indicado por neuronasSigCapa
void ForwardCapa(double **theta, double *entrada, double *neuronasSigCapa,int numNeuronas){
	int capaActual = 0;
	double *ptrTheta = NULL;//*ptrEntrada = NULL;
	int capaActualTheta = 0;
	int pos;
	// for(int pos = 0; pos<N;pos++)
	// 	cout<<"thetas1[pos][0] "<<thetas1[pos][0]<<endl;

	for (pos = 0; pos < numNeuronas;pos++){
		ptrTheta = theta[capaActualTheta];
		// ptrEntrada = entrada[capaActualTheta];
		// cout<<"Invocando a H "<<pos<<endl;
		neuronasSigCapa[pos] = getH(theta,entrada,capaActualTheta);
		capaActualTheta++;
	}
	 cout<<"Neuronas en la capa "<<pos<<" Cantidad de caracteristicas "<<N<<endl;
}
// ***** Cost() calcula el costo según la definición derivada
// Se apolla de getH() al final no fue utilizada pero se deja por si en un
// futuro sea requerida
double Cost(double **theta,double *entrada,int capaActualTheta,double clasif){
  double costo = 0.0;
  double h = getH(theta,entrada,capaActualTheta);
  if(clasif == 1){
    costo = - log(h + eps);
  }else{
    costo = - log(1 - h + eps);
  }
  // cout<<"Costo = "<<costo<<endl;
  return costo;
}
void ParameterVector(double **theta,double *entrada,double *clasif, int numNeuAct, int numNeuPrev){
	double *gradiente = NULL;
	gradiente = new double[numNeuAct];

	for(int j=0; j<numNeuAct; j++){//neurona en la que estamos
			gradiente[j]=0;
		for(int i= 0; i<numNeuPrev; i++){
			//Theta en al que estamos
			//for(int k=0; k<DATA; k++){//?????
				theta[j][i] = theta[j][i]+eps;
	 			gradiente[j] = (clasif[i]*log(getH(theta,entrada,j, clasif)+eps) + (1-clasif[i])*log(1-getH(theta,entrada,j, clasif)+eps))+gradiente[j];
				theta[j][i] = theta[j][i]-2*eps;
				gradiente[j] = (clasif[i]*log(getH(theta,entrada,j, clasif)+eps) + (1-clasif[i])*log(1-getH(theta,entrada,j, clasif)+eps))+gradiente[j];
			//}
		}

			gradiente[j] = (-(1.0/double(numNeuACt)))*gradiente[i]/(2*eps);
			theta[j][i] = theta[j][i]+eps;

		//ACTUALIZACION DE THETAS
		for(int i=0; i<numNeu; i++){
			theta[j][i] = gradiente[i];
		}
	}

	delete[] gradiente;
}

void backPropagation(double **theta1, double **theta2, double **theta3, double **entrada, double *neuronas1, double *neuronas2, double *salida, int posicionRand, double *clasif, double *granDelta2, double *granDelta3, double *granDelta4){
	
	double delta4 = NULL, delta3 = NULL, aux3A = NULL, delta2 = NULL, aux2A = NULL, aux3B = NULL, aux2B = NULL;
	delta4 = new double[M];
	delta3 = new double[M*r];
	delta2 = new double[M*r*r];

	aux3A = new double[M*r];
	aux2A = new double[M*r*r];

	aux3B = new double[M*r];
	aux2B = new double[M*r*r];

	for(int i= 0; i<M, i++){
		delta4[i] = salida[i]-clasif[posicionRand];}
	

	for(int i = 0; i<M*r; i++){
			aux3A[i]=0;
			aux3B[i]=neuronas2[i]*(1-neuronas2[i]);
		for(int j=0; j<M, j++){
			aux3A[i] = aux3A[i] + theta3[j][i]*neuronas2[i];
		}	
	}

	for(int i=0; i<M*r*r; i++){
		aux2A[i]=0;
		aux2B[i]=neuronas1[i]*(1-neuronas2[i]);
		for(int j = 0; j<M*r; j++){
			aux2A[i] = aux2A[i] + theta2[j][i]*neuronas1[i];
		}	
	}
	
	for(int i=0; i<M*r; i++){
		delta3[i] = aux3A[i]*aux3B;}

	for(int i=0; i<M*r*r; i++){
		delta2[i] = aux2A[i]*aux2B;}

	for(int i=0; i<N; i++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M*r*r; j++){
			granDelta1[j][i] = granDelta1[j][i] + neuronas1[j]*delta2[j];
		}
	}

	for(int i=0; i<M*r*r; i++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M*r; j++){
			granDelta2[j][i] = granDelta2[j][i] + neuronas2[j]*delta3[j];
		}
	}

	for(int i=0; i<M*r; i++){ //N, M*r*r; M*r*r, M*r; M*r, M
		for(int j=0; j<M; j++){
			granDelta3[j][i] = granDelta3[j][i] + salida[j]*delta4[j];
		}
	}


	delete[] delta4;
	delete[] delta3;
	delete[] delta2;
	delete[] aux3A;
	delete[] aux2A;
	delete[] aux3B;
	delete[] aux2B;

}


// ***** getGradient() utilizada a Cost()
// Al final ésta función no fue utilizada pero se deja por si en un futuro sea
// requerida
// double getGradient(double *theta,double x[][DATA],double *clasif){
//
//   double J = 0.0, suma =0.0;
//   // size_t m = sizeof(x)/sizeof(x[0][0])/2;
//   for (int i = 0; i < DATA; i++){
//     // cout<<"Entre "<<endl;
//     int clasifACtual = clasif[i];
//     suma += Cost(theta,x,i,clasifACtual);
//   }
//   J = suma/double(DATA);
//   // cout<<"Sali "<<endl;
//   return J;
// }