#ifndef alfaStabilny_h
#define alfaStabilny_h
#include <random>
#include <cmath>
using namespace std;

class alfaStabilny
{
	public:
        alfaStabilny(double Alfa, double Beta, double Mi,double Sigma, mt19937* Gen);
        double losuj();
        double dL(double dt);
    private:
    	double C;
    	double D;
    	double jedenNadAlfa;
    	double jedenMinusAlfaNadAlfa;
    	double alfa;
    	double beta;
    	double mi;
    	double sigma;
		mt19937* gen;
};

#endif
