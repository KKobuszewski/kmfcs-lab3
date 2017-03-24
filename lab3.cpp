#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


#include <Eigen/Dense>


double V0 = 1.0;
double a = 10;
double b= 1;
double Ecut = 15;

double kmin = -M_PI/a;
double kmax =  M_PI/a;

int npw = (int) round(sqrt(Ecut)/(2.0*M_PI/a)+0.5);


// inline double G(const unsigned n)
// {
//     return 2.0*M_PI*n/a;
// }


void construct_hamiltonian(Eigen::MatrixXd &H, const double k, double* G)
{
    for (unsigned i=0; i<npw; i++)
    for (unsigned j=0; j<npw; j++)
    {
        
        if (i == j)  H(i,j) = pow( k + G[i],2) - V0/a/b;
        else         H(i,j) = -2.0 * V0 * sin( 0.5*b*(G[j] - G[i]) ) / ( a*(G[j] - G[i]) );
    }
}


int main(int argc, char* argv[])
{
    std::cout << "Ecut: " << Ecut << std::endl;
    std::cout << "npw:  " << npw  << std::endl;
    
    double* G = (double*) malloc( npw*sizeof(double) );
    G[0] = 0.;
    for (unsigned i=1; i < npw; i+= 2)
    {
        G[i]   =  M_PI*((double)i+1.0)/a;
        G[i+1] = -M_PI*((double)i+1.0)/a;
    }
    
    std::cout << "Reciprocal lattice: " << std::endl;
    for (unsigned i=0; i<npw; i++)
        std::cout << G[i] << "  ";
    std::cout << std::endl;
    std::cout << std::endl;
    
    Eigen::MatrixXd H(npw,npw);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
    
    std::ofstream file;
    file.open("bands.out");
    file << std::setprecision(15) << std::fixed;
    std::cout << std::setprecision(5) << std::fixed;
    
    double N = 50.0;
    double dk = (kmax - kmin)/(N-1);
    for (double k = kmin; k <= kmax; k+= dk)
    {
        construct_hamiltonian(H, k, G);
        if (k == kmin) std::cout << H << std::endl << std::endl;
        eigensolver.compute(H);
        
        std::cout << k << "\t";
        std::cout << eigensolver.eigenvalues().transpose() << std::endl;
        file << k << "  " << eigensolver.eigenvalues().transpose() << std::endl;
    }
    
    file.close();
    free(G);
    
    return EXIT_SUCCESS;
}