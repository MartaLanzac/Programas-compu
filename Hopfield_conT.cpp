#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "gsl_rng.h"

using namespace std;

#define N 10 // N*N dimensón de la red
#define pasos_montecarlo 10000
#define T 0.005
#define P 17 //num de patrones que introducimos

//definimos las variables aleatorias
gsl_rng *tau;

int main(void)
{
	//introducimos los patrones a memorizar
	int I[P][N]; //matriz que almacena todos los patrones
	ifstream fpatr;
    string npatr;

	for(int mu=0; mu<P; mu++)
    {
        //Abro el patron correspondiente y lo vuelco en I
        npatr="PATRON"+to_string(mu)+".txt";
        fpatr.open(npatr);
        for(int j=0; j<N; j++)
        {
            fpatr >> I[mu][j];
        }
        fpatr.close();
    }

    //Para calcular la diferencis de energía necesitamos conocer varios parámetros que están 
    //relacionados con los patrones introducidos. Los calculamos:

    //Calculamos a

    double a;
    a=0.0;

    for (int mu=0; mu<P; mu++)
    {
    	for (int i = 0; i < N; ++i)
    	{
    		a=a+I[mu][i]*1.0;
    	}
    }

    a=a/(N*P);

    cout << "a vale " << a << endl;

    //calculamos las matrices J

    double J[N][N];

    for (int i=0; i<N; i++)
    {
    	for(int j=0; j<N; j++)
    	{
    		J[i][j]=0.0;

    		if (i==j) J[i][j]=0.0;
    		else
    		{
    			for (int mu=0; mu<P; mu++)
    			{
    				J[i][j]=J[i][j]+(I[mu][i]*1.0-a)*(I[mu][j]*1.0-a);
    			}	

    			J[i][j]=J[i][j]/(a*(1-a)*N*1.0);
    		}
    	}
    }


    //Calculamos las matrices theta

    double theta[N];

    for (int i=0; i<N; i++)
    {
    	theta[i]=0.0;

    	for(int j=0; j<N; j++)
    	{
    		theta[i]=theta[i]+J[i][j];
    	}

    	theta[i]=theta[i]*0.5;
    }

    //Ahora que tenemos todos los parámetros introducimos a nuestra red la conf inicial

    int s[N];
    ifstream c_inicial;

    c_inicial.open("Conf_inicial.txt");

    for (int i=0; i<N; i++)
    {
    	c_inicial >> s[i];
    }

    //Ahora evolucionamos el sistema;

    extern gsl_rng *tau;
	int semilla1=98474389;

	tau=gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(tau,semilla1);

	int n; //nodo que elegimos al azar
	double dH, suma, p, xi;
	int ds;

    for (int pMC=0; pMC<pasos_montecarlo; pMC++)
    {
    	for (int k=0; k<N; k++)
    	{
    		//elegimos el primer punto al azar
    		n=gsl_rng_uniform_int(tau,N);

    		//cout << "sn inicial " << s[n] <<endl;

    		//calculamos el delta de energía teniendo en cuenta las condiciones de contorno
    		suma=0.0;

    		for (int i=0; i<N; i++)
    		{
    			if(i != n)
    			{
    				suma=suma+J[n][i]*s[i];
    			} 		
    		}

    		if (s[n]==1) ds=-1;
    		if (s[n]==0) ds=1;

    		//cout << "su ds " << ds << "	";

    		dH=ds*(theta[n]-suma);

    		//cout << "su dH " << dH << "	";

    		//evaluamos p=min(1,exp(-de/T));

			if(1.0 < exp(-dH/T)) p=1.0;
			else p=exp(-dH/T);

			//cout << "su p " << p << "	";

			//generamos un número aleatorio uniforme entre 0 y 1

			xi = gsl_rng_uniform(tau);

			//cout << "la xi " << xi << endl;

			if (xi < p)
			{
				if(s[n]==1) s[n]=0;
				else s[n]=1;
			} 

			//cout << "sn final " << s[n] << endl << endl;
    	}
    }

    //Calculamos el overlap

    double overlap[P];
    ofstream salida;

    salida.open("salida.txt");

    for (int mu=0; mu<P; mu++)
    {
    	overlap[mu]=0.0;
    	for (int i=0; i<N; i++)
    	{
    		overlap[mu]=overlap[mu]+(I[mu][i]-a)*(s[i]-0.5);
    	}

    	overlap[mu]=overlap[mu]/(a*(1-a)*N*1.0);

     	salida << overlap[mu] << "	";
    }

    salida.close();

    return 0;
}