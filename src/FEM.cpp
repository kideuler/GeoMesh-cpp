#include "GeoMesh.hpp"
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef  Eigen::Triplet<double> T;

vector<double> Poisson_2d(Mesh* msh, vector<double> &frhs, vector<int> &bndnodes, vector<double> &dvals){
    int ndofs = msh->coords.size();

    // setting up matrix and rhs vector
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    
    SpMat A(ndofs,ndofs);
    A.reserve(Eigen::VectorXi::Constant(ndofs,10));
    Eigen::VectorXd b(ndofs);
    vector<double> u(ndofs);

    vector<vector<double>> elem_Stiff = Zeros<double>(3,3);
    vector<double> elem_load(3);
    vector<vector<double>> ps = Zeros<double>(3,2);
    double area, rhs_coeff;
    auto start = chrono::high_resolution_clock::now();
    for (int ii = 0; ii<msh->nelems; ii++){
        ps[0] = msh->coords[msh->elems[ii][0]];
        ps[1] = msh->coords[msh->elems[ii][1]];
        ps[2] = msh->coords[msh->elems[ii][2]];

        area = ((ps[1][0]-ps[0][0])*(ps[2][1]-ps[0][1]) - \
        (ps[1][1]-ps[0][1])*(ps[2][0]-ps[0][0]))/2;
        rhs_coeff = area/3;

        for (int jj = 0; jj<3; jj++){
            // inserting into stiffness matrix
            for (int kk = 0; kk<3; kk++){

            }

            // inserting right hand side
            b(msh->elems[ii][jj]) = frhs[msh->elems[ii][jj]]*(rhs_coeff);
        }
    }
    A.makeCompressed();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Assembled stiffness matrix in: " << duration.count()/1e6 << " seconds" << endl;


}
