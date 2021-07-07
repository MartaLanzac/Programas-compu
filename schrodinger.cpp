#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <complex>
#include "gsl_rng.h"

using namespace std;

#define N 1000
#define nciclos 100
#define lambda 1.0
#define pi 3.1415926535
#define mfinal 1 // número de medidas totales que realizamos
#define nD 580

//definicmos el puntero para numeros aleatorios como variable externa
gsl_rng*tau;
gsl_rng*tau2;

int main (void)
{
	//variables que definimos para los numeros aleatorios
	extern gsl_rng*tau;
	extern gsl_rng*tau2;

	int semilla=98635888;
	tau=gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(tau,semilla);

	int semilla2=83472364;
	tau2=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(tau2,semilla2);


	float xD, xI; //números aleatorios para la prob derecha e izq respectivamente

	cout << "N = " << N << ", nciclos = " << nciclos << ", lambda = " << lambda << "y nD = " << nD <<endl;

	//generamos en primer lugar todos los parámetros que no dependen del tiempo

	float k0, s, V[N+1];

	
	k0=2*pi*nciclos/(N*1.0);

	s=1.0/(4*k0*k0);

	//calculamos el potencial

	for (int j=0; j<N+1; j++)
	{
		if ( (j > 2.0*N/(5.0)) && (j < 3.0*N/(5.0)))
		{
			V[j]=lambda*k0*k0;

		}

		else V[j]=0.0;
	}

	//calculamos la función de onda inicial

	complex<float> phi[N+1], phi_0[N+1];

	double real, im;

	for (int j = 0; j < N+1; j++)
	{
		if (j == 0 || j == N) phi[j]= 0.0;

		else
		{	
			real = cos(k0*j) * exp(-8*(4*j-N)*(4*j-N)/(1.0*N*N));

			if ((2*nciclos) % N == 0)
			{
				im = 0.0;
			}
			else im=sin(k0*j) * exp(-8*(4*j-N)*(4*j-N)/(1.0*N*N));

			phi[j] = complex<float> (real,im);
			phi_0[j] = phi[j];
		}
	}

	//Calculamos A j^0

	complex<float> Aj0[N+1];

	im=2/s;

	for (int j=0; j<N+1; j++)
	{
		real=-2-V[j];

		Aj0[j]=complex<float>(real,im);
	}

	//calculamos alfa
	complex<float> cte_1=complex<float>(1.,0.);
	complex<float> alfa[N];
	complex<float> denom;
	
	alfa[N-1]=0.0;
	alfa[0]=0.0;

	for (int j=N-2; j > 0; j--)
	{
		denom=Aj0[j+1]+alfa[j+1];

		alfa[j]=-cte_1/denom;
	}

	//para el cálculo de las b
	complex<float> b[N+1];
	complex<float> i=complex<float>(0.,1.);
	complex<float> cte_4=complex<float>(4.,0.);

	//para el cálculo de las bets
	complex<float> beta[N];
	complex<float> num;

	//definimos las chi
	complex<float> chi[N+1];


	complex<float> newphi[N+1];
	float norma, amplitud[N+1], prob;

	norma=0.0;
	prob=0.0;

	float PD, PI;
	PD = PI = 0.0; 

	int mT, no;
	mT=0;
	no=0;

	ofstream fich;

	fich.open("salida.txt");

	for (int m=0; m<mfinal; m++)
	{
		//en primer lugar hacemos al sistema evolucionar de formal normal nD veces

		for (int k=0; k < nD; k++) // k es en realidad el parámetro que marca el tiempo en estas iteraciones, empezamos en k=0 y calculamos las ctes a k=0, para calcula la phi a k=1
		{
			//calculamos las bs y betas

			for (int j=0; j < N+1; j++)
			{
				b[j] = cte_4*i*phi[j]/s;
			}

			beta[N-1]=0.0;
			beta[0]=0.0;

			for (int j=N-2; j > 0; j--)
			{
				num=b[j+1]-beta[j+1];
				denom=Aj0[j+1]+alfa[j+1];

				beta[j]=num/denom;
			}

			//calculamos las chi

			chi[0]=0.0;
			chi[N]=0.0;

			for (int j=1; j < N; j++)
			{
				chi[j] = alfa[j-1] * chi[j-1] + beta[j-1];
			}

			//una vez que tenemos todos los parámetros podemos calcular las nuevas funciones de onda

			for (int j=0; j < N+1; j++)
			{
				newphi[j] = chi[j] - phi[j];
			}

			//estas phis deben estar normalizadas por lo que calculamos su norma

			norma=0.0;

			for (int j=0; j < N+1; j++)
			{
				amplitud[j]=norm(newphi[j]);

				norma=norma+amplitud[j];
			}

			//normalizamos las newphi, que las renombramos phi
			for (int j=0; j < N+1; j++)
			{
				phi[j]=newphi[j]/sqrt(norma);
			}

		} //aqui aacaba la evol temporal

		for (int j=0; j<N+1; j++)
		{
			fich << j*1.0/1000.0 << "	" << norm(phi[j]) << endl; 
		}

		//ahora empezamos a medir

			//1-> Calculamos PD

			PD=0.0;

			for (int j=4*N/5 ; j<N+1; j++)
			{
				PD=PD+norm(phi[j]);
			}

			//2-> generamos un numero aleatorio entere 0 y 1 

			xD=gsl_rng_uniform(tau);

			if(xD < PD) //3-> aumentamos mT y volvemos a 1
			{
				mT=mT+1;

				for (int j=0; j<N+1; j++)
				{
					phi[j] = phi_0[j];
				}
			} //si se ha medido mT

			else //4-> cambiamos la función de onda
			{
				for (int j=4*N/5; j<N+1; j++)
				{
					phi[j]=0.0;
				}

				norma=0.0;
				for (int j=0; j<N+1; j++)
				{
					norma = norma + norm(phi[j]);
				}

				for (int j=0; j<N+1; j++)
				{
					phi[j]=phi[j]/sqrt(norma);
				}

				//5-> Calculamos PI

				PI = 0.0;

				for (int j=0; j<N/5; j++)
				{
					PI = PI + norm(phi[j]);
				}

				//6-> Generamos un número alatorio entre 0 y 1

				xI=gsl_rng_uniform(tau2);

				if(xI < PI) // 7-> ya ha sido detectada por lo que volvemos al inicio del experimento
				{
					no=no+1;

					for (int j=0; j<N+1; j++)
					{
						phi[j] = phi_0[j];
					}
				}
				else // 8-> generamos una nueva función de onda que evoluciona y que volveremos a medir
				{
					m=m-1; 
					for (int j=0; j<N/5 +1; j++)
					{
						phi[j]=0.0;
					}

					norma=0.0;
					for (int j=0; j<N+1; j++)
					{
						norma = norma + norm(phi[j]);
					}

					prob=0.0;

					for (int j=0; j<N+1; j++)
					{
						phi[j]=phi[j]/sqrt(norma);

						prob = prob + norm(phi[j]);
					}
				}
			}
	}

	cout << "mT= " << mT << " no= " << no << endl;


	fich.close();
}
