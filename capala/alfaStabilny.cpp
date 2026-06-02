#include <iostream>
#include <random>
#include <cmath>
#include "alfaStabilny.h"
using namespace std;

#define PI 3.14159265
#define PI_DRUGICH 1.57079633
#define JEDEN_NAD_PI 0.318309886

alfaStabilny::alfaStabilny(double Alfa, double Beta, double Mi,double Sigma,  mt19937* Gen)
{	
	alfa=Alfa;
	beta=Beta;
	mi=Mi;
	sigma=Sigma;
	jedenNadAlfa = 1./alfa;
	jedenMinusAlfaNadAlfa = (1.-alfa) * jedenNadAlfa;
	gen=Gen;
	
	if(alfa>2.||alfa<0.)
	{
		cout<<"Alfa poza przedzialem (0,2]"<<endl;
		exit(0);
	}
	else if(alfa<0.99999 || alfa>1.00001)
	{
		C=atan(beta * tan(PI*alfa*0.5) ) * jedenNadAlfa;
		D=sigma*pow(cos(atan(beta* tan(0.5 * PI * alfa) ) ),-jedenNadAlfa);
	}
	else
	{
		C=0.;
		D=0.;	
	}
};

double alfaStabilny::losuj()
{
	exponential_distribution<double> exp_dyst(1.);
	double V=PI*(generate_canonical<double,10>(*gen)-0.5);
	double W=exp_dyst(*gen);

	double liczba_losowa(0.);
		
	if(alfa<0.99999 || alfa>1.00001)
	{
		liczba_losowa= D * sin( alfa * (V+C) )/pow(cos(V),jedenNadAlfa)*pow(cos(V - alfa * (V+C))/W,jedenMinusAlfaNadAlfa);
	}
	else
	{
		liczba_losowa= 2. * sigma * JEDEN_NAD_PI * ( (PI_DRUGICH + beta * V ) * tan(V) - beta  * log(PI_DRUGICH * W * cos(V) / (PI_DRUGICH + beta * V)) ) + mi;
	}
	
	return liczba_losowa;
		
}
