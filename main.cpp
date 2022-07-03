/* main.cpp */

#include <iostream>
#include "element.h"
#include "grid.h"
#include "vec.h"
#include "vec.hpp"
#include "webgl.h"
#include "spmat.h"
#include "poisson.h"
#include "convection.h"

#include <iomanip> 

using namespace flux;

// define which test case to use
// setting RHS_SINE to 1 uses the u(x) = sin(pi*x)*sin(pi*y), f(x) = 2 pi^2 sin(pi*x) sin(pi*y) problem
// setting RHS_SINE to 0 uses the u(x) = 0.5*x*(1 - x) + 0.5*y*(1 - y), f(x) = 2 problem 
#define RHS_SINE 1

// define which types of equation to use
// setting EQ_TYPE to 1 uses Poisson equation solver
// setting EQ_TYPE to 0 uses Convection-Diffusion equation solver
#define EQ_TYPE 0



int 
main() {
    clock_t tStart = clock();

    Viewer viewer;

    int n = 100;
    Grid<Quad> mesh( {n,n} , 2 );

    #if EQ_TYPE
    vecd<double> u = solve_model_problem(n,1e-10,1e2,RHS_SINE);
    vecd<double> u_nl = solve_model_problem_cg(n,RHS_SINE);
    #else
    vecd<double> u = solve_model_problem_con(n,1e-10,1e2,RHS_SINE);
    vecd<double> u_nl = solve_model_problem_cg_con(n,RHS_SINE);
    #endif

    std::vector<double> res;
    for(int k = 0;k < (n+1)*(n+1); k++) {
        res.push_back(u(k));
    }
    
    std::vector<double> res_nl;
    for(int k = 0;k < (n+1)*(n+1); k++) {
        res_nl.push_back(u_nl(k));
    }

    // calculate error
    double error = 0.0;

    for (int i = 0; i < (n+1)*(n+1); i++) {
        int r = i / (n+1);
        int c = i % (n+1);
        double xi = double(c) / n;
        double yi = double(r) / n;

        #if RHS_SINE
        double ua = sin(M_PI * xi) * cos(M_PI * yi);
        #else
        double ua = 0.5*xi*(1 - xi) + 0.5*yi*(1 - yi);
        #endif
        error += std::pow( ua - u(i) , 2.0 );
        printf("u[%d,%d] = %g, ua = %g, error = %g\n",r,c,u(i),ua,std::fabs(ua - u(i)));
    }

    error = std::sqrt(error/((n+1)*(n+1)));
    printf("error = %g\n",error);

    mesh.create_vertex_field("solution", res);
    mesh.create_vertex_field("nl_solution", res_nl);  

    viewer.add(mesh);
    printf("runtime: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    //https://stackoverflow.com/questions/31905157/how-to-print-the-execution-time-of-my-code
    viewer.run();
    
}

