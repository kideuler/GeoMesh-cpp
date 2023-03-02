#include "GeoMesh.hpp"
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef  Eigen::Triplet<double> T;

void apply_dbc(SpMat &A, Eigen::VectorXd &b, vector<int> &bndnodes, vector<double> &dvals, vector<bool> &is_dir);
void eval_shapefunc(int degree, vector<double> &cs, vector<double> &phi, vector<vector<double>> &dphi);
void Tri_3(double xi, double eta, vector<double> &phi, vector<vector<double>> &dphi);
void Tri_6(double xi, double eta, vector<double> &phi, vector<vector<double>> &dphi);
void Quadrature_rule(int degree, vector<vector<double>> &cs, vector<double> &ws);

vector<double> Poisson_2d(Mesh* msh, vector<double> &frhs, double kappa, vector<int> &bndnodes, vector<double> &dvals){
    int ndofs = msh->coords.size();
    int nnodes = (msh->degree+1)*(msh->degree+2)/2;
    assert(frhs.size() == ndofs);
    vector<bool> is_dir_node(ndofs,false);
    for (int n = 0; n<bndnodes.size(); n++){
        is_dir_node[bndnodes[n]] = true;
    }

    // setting up matrix and rhs vector
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    
    SpMat A(ndofs,ndofs);
    A.reserve(Eigen::VectorXi::Constant(ndofs,msh->degree*10));
    Eigen::VectorXd b(ndofs),sol(ndofs-bndnodes.size());
    vector<double> u(ndofs);

    auto start = chrono::high_resolution_clock::now();

    // P1 triangular elements
    if (false){ // fast assembly for P1 elements
        vector<vector<double>> elem_Stiff = Zeros<double>(3,3);
        vector<vector<double>> D = Zeros<double>(2,3);
        vector<double> elem_load(3);
        vector<vector<double>> ps = Zeros<double>(3,2);
        double area, rhs_coeff;
        for (int ii = 0; ii<msh->nelems; ii++){
            ps[0] = msh->coords[msh->elems[ii][0]];
            ps[1] = msh->coords[msh->elems[ii][1]];
            ps[2] = msh->coords[msh->elems[ii][2]];

            area = ((ps[1][0]-ps[0][0])*(ps[2][1]-ps[0][1]) - \
            (ps[1][1]-ps[0][1])*(ps[2][0]-ps[0][0]))/2;
            rhs_coeff = area/3;

            D = {{ps[2][0]-ps[1][0], ps[0][0]-ps[2][0], ps[1][0]-ps[0][0]}, \
            {ps[2][1]-ps[1][1], ps[0][1]-ps[2][1], ps[1][1]-ps[0][1]}};

            elem_Stiff = (Transpose(D)*D) / (4*area);

            for (int jj = 0; jj<3; jj++){
                // inserting into stiffness matrix
                for (int kk = 0; kk<3; kk++){
                    A.coeffRef(msh->elems[ii][jj], msh->elems[ii][kk]) += kappa*elem_Stiff[jj][kk];
                }

                // inserting right hand side
                b(msh->elems[ii][jj]) = frhs[msh->elems[ii][jj]]*(rhs_coeff);
            }
        }
    } else{
        vector<vector<double>> Jt = Zeros<double>(2,2);
        vector<vector<double>> Jtinv = Zeros<double>(2,2);

        // quadrature rules
        vector<vector<double>> cs;
        vector<double> ws;
        Quadrature_rule(msh->degree, cs,ws);
        int nq = cs.size();
        
        double detJ;

        vector<vector<double>> phi = Zeros<double>(nq,nnodes);
        vector<vector<vector<double>>> dphi(nq);
        vector<vector<double>> grad_sol = Zeros<double>(nnodes,2);
        for(int n = 0; n<nq; n++){
            dphi[n].assign(nnodes, vector<double>(2,0));
            eval_shapefunc(msh->degree, cs[n], phi[n], dphi[n]);
        }

        vector<vector<double>> ps = Zeros<double>(nnodes,2);
        vector<double> fs(nnodes);
        double v;
        for (int ii = 0; ii<msh->nelems; ii++){

            for (int nn = 0; nn<nnodes; nn++){
                ps[nn] = msh->coords[msh->elems[ii][nn]];
                fs[nn] = frhs[msh->elems[ii][nn]];
            }

            for (int q = 0; q<nq; q++){

                // compute Jacobian
                Jt[0][0] = 0.0; Jt[0][1] = 0.0; Jt[1][0] = 0.0; Jt[1][1] = 0.0;
                for (int n = 0; n<nnodes; n++){
                    for (int jj = 0; jj<2; jj++){
                        for (int kk = 0; kk<2; kk++){
                            Jt[jj][kk] += ps[n][kk]*dphi[q][n][jj];
                        }
                    }
                }
                detJ = (Jt[0][0]*Jt[1][1] - Jt[1][0]*Jt[0][1]);

                Jtinv[0][0] = Jt[1][1]/detJ;
                Jtinv[0][1] = -Jt[0][1]/detJ;
                Jtinv[1][0] = -Jt[1][0]/detJ;
                Jtinv[1][1] = Jt[0][0]/detJ;
                for (int n = 0; n<nnodes; n++){
                    grad_sol[n] = Jtinv*dphi[q][n];
                }


                for (int jj = 0; jj<nnodes; jj++){
                    // inserting into stiffness matrix
                    for (int kk = 0; kk<nnodes; kk++){
                        for (int dim = 0; dim<2; dim++){
                            A.coeffRef(msh->elems[ii][jj], msh->elems[ii][kk]) += kappa*ws[q]*detJ*grad_sol[kk][dim]*grad_sol[jj][dim];
                        }
                    }
                    // inserting right hand side
                    v = inner(fs,phi[q]);
                    b(msh->elems[ii][jj]) += v*ws[q]*detJ*phi[q][jj];
                }
            }
        }
        


    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Assembled stiffness matrix in: " << duration.count()/1e6 << " seconds" << endl;

    // apply dirichlet boundary conditions
    start = chrono::high_resolution_clock::now();
    apply_dbc(A,b,bndnodes,dvals,is_dir_node);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Applied boundary conditions in: " << duration.count()/1e6 << " seconds" << endl;
    
    A.makeCompressed();

    start = chrono::high_resolution_clock::now();
    solver.analyzePattern(A);
    solver.factorize(A); 
    if (solver.info() == Eigen::Success){
        sol = solver.solve(b); 
    } else {
        cout << "Eigen failed to solve" << endl;
    }
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Solved matrix in: " << duration.count()/1e6 << " seconds" << endl;

    for (int ii=0; ii<ndofs; ii++){
        u[ii] = sol(ii);
    }

    return u;
}

void apply_dbc(SpMat &A, Eigen::VectorXd &b, vector<int> &bndnodes, vector<double> &dvals, vector<bool> &is_dir){
    int ndbc = bndnodes.size();
    assert(ndbc == dvals.size());
    int nv = b.size();
    double val;
    for (int n = 0; n<nv; n++){
        if (!is_dir[n]){
            val = 0.0;
            for (int jj = 0; jj<ndbc; jj++){
                val += A.coeff(n,bndnodes[jj])*dvals[jj];
                if (abs(A.coeff(n,bndnodes[jj])) > 1e-12){
                    A.coeffRef(n,bndnodes[jj]) = 0.0;
                }
            }
            b(n) = b(n) - val;
        }
    }
    for (int n = 0; n<ndbc; n++){
        b(bndnodes[n]) = dvals[n];
    }
    for (int ii = 0; ii<nv; ii++){
        if (is_dir[ii]){
            for (int jj = 0; jj<nv; jj++){
                if (abs(A.coeff(ii,jj)) > 1e-12){
                    A.coeffRef(ii,jj) = 0.0;
                }
            }
            A.coeffRef(ii,ii) = 1.0;
        }
    }
    A.prune(1e-6);
    return;
}

void eval_shapefunc(int degree, vector<double> &cs, vector<double> &phi, vector<vector<double>> &dphi){
    switch (degree){
        case 1:
            Tri_3(cs[0],cs[1],phi,dphi);
            break;
        case 2:
            Tri_6(cs[0],cs[1],phi,dphi);
            break;
    }
}

void Tri_3(double xi, double eta, vector<double> &phi, vector<vector<double>> &dphi){
    phi[0] = 1-xi-eta;
    phi[1] = xi;
    phi[2] = eta;

    dphi[0] = {-1,-1};
    dphi[1] = {1,0};
    dphi[2] = {0,1};
    return;
}

void Tri_6(double xi, double eta, vector<double> &phi, vector<vector<double>> &dphi){
    double t2,t3,t4,t5,t8,t9,t10,t11;

    t2 = eta*4.0;
    t3 = pow(eta,2);
    t4 = xi*4.0;
    t5 = pow(xi,2);
    t8 = t3*2.0;
    t9 = t5*2.0;
    t10 = t2-3.0;
    t11 = -t2+4.0;

    phi = {eta*-3.0+t8+t9+t10*xi+1.0, \
    t9-1.0*xi, \
    -eta+t8, \
    t5*-4.0-xi*(t10-1.0), \
    t2*xi, \
    t2-t3*4.0-t2*xi};

    dphi = {{t4+t10, t4+t10},
        {t4-1.0, 0.0},
        {0.0, t2-1.0},
        {t11-8.0*xi, -t4},
        {eta*4.0, xi*4.0},
        {-t2, eta*-8.0-t4+4.0}};
}

void Quadrature_rule(int degree, vector<vector<double>> &cs, vector<double> &ws){
    switch (degree){
        case 1:
            cs = {{1.0/3.0,1.0/3.0}};
            ws = {0.5};
            break;
        case 2:
            cs = {{2.0/3.0,1.0/6.0},{1.0/6.0,2.0/3.0},{1.0/6.0,1.0/6.0}};
            ws = {1.0/6.0,1.0/6.0,1.0/6.0};
            break;
    }
}