#include "GeoMesh.hpp"
#include <unistd.h>
#include <Eigen/Sparse>
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat;
typedef  Eigen::Triplet<double> T;

static Spline ndegree_spline(const vector<vector<double>> &xs, int degree);
static Spline Cubic_spline(const vector<vector<double>> &xs);
static Spline Quadratic_spline(const vector<vector<double>> &xs);

function<double(vector<double>)> create_curvature_hfunction(Spline* spl, int npoints, double theta, double h, double h_min, double hgrad){
    vector<double> K(npoints), H(npoints);
    vector<vector<double>> ps = Zeros<double>(npoints,2);
    for (int i = 0; i<npoints; i++){
        ps[i] = spline_var(spl, double(i)/double(npoints));
        K[i] = spline_curvature(spl, double(i)/double(npoints));
        H[i] = (theta*M_PI/180.0) / K[i];
        H[i] = min(max(h*h_min, H[i]), h);
    }

    function<double(vector<double>)> hF = [ps,H,h,hgrad](vector<double> xs){
        int nv = ps.size();
        if (xs.size() == 0){
            return mean(H);
        }
        double r,xi,alpha;
        alpha = 1000;
        for (int i = 0; i<nv; i++){
            r = (pow(xs[0]-ps[i][0],2)+pow(xs[1]-ps[i][1],2));
            xi = r/(hgrad*h);
            alpha = min((1-min(xi,1.0))*H[i] + h*min(xi,1.0),alpha);
        }
        return alpha;
    };
    return hF;
    
}

vector<double> spline_point_segment(Spline* spl, double a, double b, double ratio){
    double arclength=0.0;
    double q[5] = {-(1/3)*sqrt(5+2*sqrt(10/7)), -(1/3)*sqrt(5-2*sqrt(10/7)),0,(1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))};
    double w[5] = {(322-13*sqrt(70))/900.0, (322+13*sqrt(70))/900.0, 125/225, (322+13*sqrt(70))/900.0,(322-13*sqrt(70))/900.0};

    // find total arclength
    for (int k=0; k<5; k++){
        arclength += w[k]*norm(spline_var(spl,q[k]*(b-a)/2 + (a+b)/2, 1));
    }

    // implement full bisection search method later if needed
    if (abs(b-a) > 0.5){
        if (b<a){
            b +=1;
        } else {
            a +=1;
        }
    }
    vector<double> xs = spline_var(spl, (a+b)/2,0);
    return xs;
}

vector<double> spline_var(Spline* spl, double t, int order){
    assert(order < 3);
    if (t<0){
        t = 1-t;
    }
    if (t>1){
        t = fmod(t,1);
    }
    int n = (int) floor(double(spl->nv)*t);
    
    if (n==spl->nv){
        n--;
    }
    bool stop = false;
    while (!stop){
        if (n == 0){
            if (spl->params[n+1] >= t){
                stop = true;
            } else {
                n++;
            }
        } else if (n == spl->nv-1){
            if (spl->params[n] <= t){
                stop = true;
            } else {
                n--;
            }
        } else {
            if (spl->params[n] <= t && spl->params[n+1] >= t){
                stop = true;
            } else if (spl->params[n] > t){
                n--;
            } else {
                n++;
            }
        }
    }
    
    // computing function values or derivatives
    double x_j = (t-spl->params[n])/(spl->params[n+1]-spl->params[n]);
    vector<double> xs = {0.0,0.0}; double coeff;
    for (int i = order; i<spl->degree+1; i++){
        coeff = 1.0;
        for (int j = 0; j<order; j++){
            coeff = coeff*double(i-j);
        }
        xs[0] += coeff*pow(x_j,double(i-order))*spl->xweights[n][i];
        xs[1] += coeff*pow(x_j,double(i-order))*spl->yweights[n][i];
    }

    return xs;
}

double spline_curvature(Spline* spl, double t){
    vector<double> D1 = spline_var(spl,t,1);
    vector<double> D2 = spline_var(spl,t,2);
    double K = abs(D1[0]*D2[1] - D1[1]*D2[0])/pow(D1[0]*D1[0] + D1[1]*D1[1],1.5);
    return K;
}

Spline spline_init(const vector<vector<double>> &xs, int degree){

    Spline spl;
    
    // special efficient algorithm for degrees 2 and 3 
    if (degree == 2){
        spl = Quadratic_spline(xs);
    } else if (degree == 3){
        spl = Cubic_spline(xs);
    } else {
        spl = ndegree_spline(xs,degree);
    }


    double arclength=0.0;
    double q[5] = {-(1/3)*sqrt(5+2*sqrt(10/7)), -(1/3)*sqrt(5-2*sqrt(10/7)),0,(1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))};
    double w[5] = {(322-13*sqrt(70))/900.0, (322+13*sqrt(70))/900.0, 125/225, (322+13*sqrt(70))/900.0,(322-13*sqrt(70))/900.0};

    double a,b,I;
    double* temp = new double[spl.nv+1];
    temp[0] = 0.0;
    for (int i = 0; i<spl.nv; i++){
        a = spl.params[i];
        b = spl.params[i+1];
        I = 0.0;
        for (int j = 0; j<5; j++){
            I += w[j]*norm(spline_var(&spl,q[j]*(b-a)/2 + (a+b)/2, 1));
        }
        I = I*(b-a)/2;
        arclength += I;
        temp[i+1] = arclength;
    }

    for (int i = 1; i<spl.nv+1; i++){
        spl.params[i] = temp[i]/arclength;
    }

    delete temp;
    return spl;
}

vector<vector<double>> construct_CVM(double t[],int nv, int degree){
    vector<vector<double>> V = Zeros<double>(nv*degree, degree+1);

    int k = 0;
    for (int i = 0; i<degree; i++){
        for (int j = i; j<degree+1; j++){
            for (int n = 0; n<nv; n++){
                V[k+n][j] = (tgamma(j+1)/tgamma(j-i+1))*pow(t[n],j-i);
            }
        }
        k+=nv;
    }

    return V;
}

static Spline ndegree_spline(const vector<vector<double>> &xs, int degree){

    // Boolean array for which points are corners
    int nv = xs.size();
    int ii,jj;

    // Setting up spline data structure
    Spline spl;
    spl.degree = degree;
    spl.nv = nv;
    spl.coords = Zeros<double>(nv,2);
    spl.xweights = Zeros<double>(nv, degree+1);
    spl.yweights = Zeros<double>(nv, degree+1);
    spl.params.resize(nv+1);

    int ndofs = (degree+1)*nv;
    // set up eigen parameters
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    
    SpMat A(ndofs,ndofs);
    A.reserve(Eigen::VectorXi::Constant(ndofs,3*(degree+1)));
    Eigen::VectorXd bx(ndofs), by(ndofs), Dxs(ndofs), Dys(ndofs);

    for (int i = 0; i<nv; i++){
        spl.coords[i][0] = xs[i][0];
        spl.coords[i][1] = xs[i][1];
        spl.params[i] = (double) i / ((double) nv+1);
    }
    spl.params[nv] = 1;

    double t[2] = {0.0,1.0};
    vector<vector<double>> V = construct_CVM(t,2,degree);

    int i,j,k;
    k = 1;
    for (i = 0; i<nv; i++){
        for (j = 0; j<(degree+1); j++){
            A.insert(2*(k-1), (degree+1)*(k-1)+j) = V[0][j];
            A.insert(2*(k-1)+1, (degree+1)*(k-1)+j) = V[1][j];
        }

        bx(2*(k-1)) = spl.coords[i][0];
        bx(2*(k-1)+1) = spl.coords[(i+1)%nv][0];
        by(2*(k-1)) = spl.coords[i][1];
        by(2*(k-1)+1) = spl.coords[(i+1)%nv][1];
        k++;
    }

    int p, v2;
    for (int deg = 1; deg<degree; deg++){
        p=1;
        k = nv + (deg)*nv;
        for (i=0; i<nv; i++){
            for (j=0; j<degree+1; j++){
                v2 = (degree+1)*(p-1)+j;
                A.insert(k, v2) = V[2*(deg)+1][j];
                v2 = (degree+1)*p+j;
                A.insert(k, v2%(nv*(degree+1))) = -V[2*(deg)][j];
            }
            k++;
            p++;
        }
    }
    A.makeCompressed();

    solver.analyzePattern(A);
    solver.factorize(A); 
    Dxs = solver.solve(bx); 
    Dys = solver.solve(by);

    for (i=0; i<nv; i++){
        for (j=0; j<degree+1; j++){
            spl.xweights[i][j] = Dxs((degree+1)*(i)+j);
            spl.yweights[i][j] = Dys((degree+1)*(i)+j);
        }
    }

    return spl;
}


static Spline Cubic_spline(const vector<vector<double>> &xs){

    // Boolean array for which points are corners
    int nv = xs.size();
    int ii,jj;

    // Setting up spline data structure
    Spline spl;
    spl.degree = 3;
    spl.nv = nv;
    spl.coords = Zeros<double>(nv,2);
    spl.xweights = Zeros<double>(nv, 4);
    spl.yweights = Zeros<double>(nv, 4);
    spl.params.resize(nv+1);

    // set up eigen parameters
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    SpMat A(nv,nv);
    A.reserve(Eigen::VectorXi::Constant(nv,3));
    Eigen::VectorXd bx(nv), by(nv), Dxs(nv), Dys(nv);

    for (int i = 0; i<nv; i++){
        spl.coords[i][0] = xs[i][0];
        spl.coords[i][1] = xs[i][1];
        spl.params[i] = (double) i / ((double) nv+1);
    }
    spl.params[nv] = 1;

    ii = 0;
    while (ii < nv){
        A.insert(ii,ii) = 4.0;
        A.insert(ii,(ii+nv-1)%nv) = 1.0;
        A.insert(ii,(ii+1)%nv) = 4.0;
        bx(ii) = 3*(xs[(ii+1)%nv][0] - xs[(ii+nv-1)%nv][0]);
        by(ii) = 3*(xs[(ii+1)%nv][1] - xs[(ii+nv-1)%nv][1]);
        ii++;
    }
    A.makeCompressed();

    solver.analyzePattern(A);
    solver.factorize(A); 
    Dxs = solver.solve(bx); 
    Dys = solver.solve(by); 

    for (ii=0; ii<nv; ii++){
        spl.xweights[ii][0] = xs[ii][0];
        spl.yweights[ii][0] = xs[ii][1];
        spl.xweights[ii][1] = Dxs(ii);
        spl.yweights[ii][1] = Dys(ii);
        spl.xweights[ii][2] = 3*(xs[(ii+1)%nv][0] - xs[ii][0]) - 2*Dxs(ii) - Dxs((ii+1)%nv);
        spl.yweights[ii][2] = 3*(xs[(ii+1)%nv][1] - xs[ii][1]) - 2*Dys(ii) - Dys((ii+1)%nv);
        spl.xweights[ii][3] = 2*(xs[ii][0] - xs[(ii+1)%nv][0]) + Dxs(ii) + Dxs((ii+1)%nv);
        spl.yweights[ii][3] = 2*(xs[ii][1] - xs[(ii+1)%nv][1]) + Dys(ii) + Dys((ii+1)%nv);
    }

    return spl;
}

static Spline Quadratic_spline(const vector<vector<double>> &xs){

    // Boolean array for which points are corners
    int nv = xs.size();
    int ii,jj;

    // Setting up spline data structure
    Spline spl;
    spl.degree = 2;
    spl.nv = nv;
    spl.coords = Zeros<double>(nv,2);
    spl.xweights = Zeros<double>(nv, 3);
    spl.yweights = Zeros<double>(nv, 3);
    spl.params.resize(nv+1);

    // set up eigen parameters
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    SpMat A(nv,nv);
    A.reserve(Eigen::VectorXi::Constant(nv,3));
    Eigen::VectorXd bx(nv), by(nv), Dxs(nv), Dys(nv);

    for (int i = 0; i<nv; i++){
        spl.coords[i][0] = xs[i][0];
        spl.coords[i][1] = xs[i][1];
        spl.params[i] = (double) i / ((double) nv+1);
    }
    spl.params[nv] = 1;

    ii = 0;
    while (ii < nv){
        A.insert(ii,ii) = 1.0;
        A.insert(ii,(ii+1)%nv) = 1.0;
        bx(ii) = 2*(xs[(ii+1)%nv][0] - xs[ii][0]);
        by(ii) = 2*(xs[(ii+1)%nv][1] - xs[ii][1]);
        ii++;
    }
    A.makeCompressed();

    solver.analyzePattern(A);
    solver.factorize(A); 
    Dxs = solver.solve(bx); 
    Dys = solver.solve(by); 

    for (ii=0; ii<nv; ii++){
        spl.xweights[ii][0] = xs[ii][0];
        spl.yweights[ii][0] = xs[ii][1];
        spl.xweights[ii][1] = Dxs(ii);
        spl.yweights[ii][1] = Dys(ii);
        spl.xweights[ii][2] = (Dxs((ii+1)%nv) - Dxs(ii))/2;
        spl.yweights[ii][2] = (Dys((ii+1)%nv) - Dys(ii))/2;
    }

    return spl;
}