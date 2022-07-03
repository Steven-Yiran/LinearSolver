#include <iomanip> 
#include "spmat.h"
#include "vec.h"
#include "vec.hpp"
#include "poisson.h"

namespace flux{

void
get_model_problem( int n , spmat<double>& A , vecd<double>& f , int IC_TYPE ) {
   
    flux_assert( f.m() == (n+1)*(n+1) );
    double h = 1.0/n;
    

    // lower boundary condition coefficients
    for (int i = 0; i < n+1 ; i++ ) {
        A(i,i) = 1;
    }
    
    for (int i = n+1; i < (n+1)*n; i++) {
        // left/right side boundary condition coefficients
        if( i % (n+1) == 0 || (i+1) % (n+1) == 0) {
            A(i,i) = 1;
            continue;
        }

        A(i,i-1) = -1 / (h*h);
        A(i,i)   =  4 / (h*h);
        if( i+1 < (n+1)*(n+1) ){
            A(i,i+1) = -1 / (h*h);
        }
        if( n+i < (n+1)*(n+1) ){
            A(i,n+i) = -1 / (h*h);
        }
        if( i >= (n+1) ){
            A(i, i-n) = -1 / (h*h);
        }

        //f(i) = 1;
    }

    // upper boundary condition coefficients
    for (int i = (n+1)*n; i < (n+1)*(n+1); i++) {
        A(i,i) = 1;
    }
    
    #if IC_TYPE
    for (int i = 0; i < (n+1)*(n+1); i++) {
        int r = i / (n+1);
        int c = i % (n+1);
        if(r == 0 || r == n || c == 0 || c == n) continue;
        double xi = double(r) / n;
        double yi = double(c) / n;
        f(i) = 2 * M_PI * M_PI * sin(M_PI * xi) * cos(M_PI * yi); 
    }
    #else
    for (int i = 0; i < (n+1)*(n+1); i++) {
        int r = i / (n+1);
        int c = i % (n+1);
        double xi = double(r) / n;
        double yi = double(c) / n;
        if(r == 0 || r == n) {
            f(i) = 0.5*yi*(1 - yi);
            continue;
        }
        if(c == 0 || c == n) {
            f(i) = 0.5*xi*(1 - xi);
            continue;
        }
        f(i) = 2; 
    }
    #endif
}

vecd<double>
solve_model_problem(int n, double tol, int max_iter , int IC_TYPE ){

    int l = (n+1)*(n+1);
    spmat<double> A(l,l);
    vecd<double> f(l);

    get_model_problem(n,A,f,IC_TYPE);
    if (n < 30) A.print_full();

    vecd<double> x(l);
    vecd<double> r(l);
    int iter = 0;
    double e = norm(A*x - f);
    // implement jacobi or gauss-siedel algorithm here
    while (e > tol && iter++ < max_iter) {

        r = f - A*x;
        for (int i = 0; i < f.m(); i++) {
            // skip boundary
            x(i) = (r(i) + A(i,i)*x(i)) / A(i,i);
        }
        e = norm(r);
        //printf("[iter %d], error = %g\n", iter, e);
    }
    printf("--> converged to %g in %d iterations\n",e,iter);

    return x;
}

vecd<double>
solve_model_problem_cg( int n , int IC_TYPE ) {

    spmat<double> A((n+1)*(n+1),(n+1)*(n+1));
    vecd<double> f((n+1)*(n+1));
    vecd<double> x((n+1)*(n+1));

    get_model_problem(n,A,f,IC_TYPE);

    A.solve_nl(f,x,false);

    return x;
}

}