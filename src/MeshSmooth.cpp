#include "GeoComp.hpp"
#include <unistd.h>
using namespace std;

double isometry_energy_tri(vector<vector<double>>& xs, vector<vector<double>> &Grad, vector<vector<vector<double>>> &Hess, double refarea, double mu);
double mesh_smoothing_tri_2d_iter(Mesh* mesh, vector<bool> no_move, vector<double> refareas, double mu, double h_max);

void mesh_smoothing_2d(Mesh* mesh, vector<bool> no_move, function<double(vector<double>)> r_ref, double mu){

    // calculate reference areas
    vector<double> refareas(mesh->nelems);
    int i,j;
    vector<double> centroid(2);
    double h,h_max;
    h=0;
    if (mu > 0){
    for (i=0; i<mesh->nelems; i++){
        centroid = (mesh->coords[mesh->elems[i][0]]+mesh->coords[mesh->elems[i][1]]+mesh->coords[mesh->elems[i][2]])/3;
        h = (3/sqrt(3))*r_ref(centroid);
        h_max = max(h,h_max);
        refareas[i] = h*h/2;
    }
    } else {
        for (i=0; i<mesh->nelems; i++){
            for (j=0; j<3; j++){
                h = norm(mesh->coords[mesh->elems[i][(j+1)%3]]-mesh->coords[mesh->elems[i][j]]);
                h_max = max(h,h_max);
            }
        }
    }

    // main smoothing algorithm
    double tol = 0.3*h_max;
    double errs = 1e6;
    int iter = 0;
    while (errs > tol && iter < 10){
        // iteration of Energy based smoothing
        errs = mesh_smoothing_tri_2d_iter(mesh, no_move, refareas, mu, h_max);
        iter++;
    }
}


double mesh_smoothing_tri_2d_iter(Mesh* mesh, vector<bool> no_move, vector<double> refareas, double mu, double h_max){

    int nv = mesh->coords.size();
    Mat* Grads = new Mat;
    *Grads = Zeros(2,nv);
    vector<vector<vector<double>>>* Hess = new vector<vector<vector<double>>>;
    (*Hess).resize(2);
    (*Hess)[0].resize(2);
    (*Hess)[1].resize(2);
    (*Hess)[0][0].resize(nv);
    (*Hess)[1][0].resize(nv);
    (*Hess)[0][1].resize(nv);
    (*Hess)[1][1].resize(nv);
    Mat ps = Zeros(3,2);
    Mat Grad_elem = Zeros(2,3);
    vector<vector<vector<double>>> Hess_elem;
    Hess_elem.resize(2);
    Hess_elem[0].resize(2);
    Hess_elem[1].resize(2);
    Hess_elem[0][0].resize(3);
    Hess_elem[1][0].resize(3);
    Hess_elem[0][1].resize(3);
    Hess_elem[1][1].resize(3);

    // Accurmulating Energy
    double Energy;
    int v;
    for (int n = 0; n<mesh->nelems; n++){
        ps[0] = mesh->coords[mesh->elems[n][0]];
        ps[1] = mesh->coords[mesh->elems[n][1]];
        ps[2] = mesh->coords[mesh->elems[n][2]];

        Energy = isometry_energy_tri(ps, Grad_elem, Hess_elem, refareas[n], 0.0);
        for (int ii = 0; ii<3; ii++){
            v = mesh->elems[n][ii];
            (*Grads)[0][v] = (*Grads)[0][v] + Grad_elem[0][ii];
            (*Grads)[1][v] = (*Grads)[1][v] + Grad_elem[1][ii];

            (*Hess)[0][0][v] = (*Hess)[0][0][v] + Hess_elem[0][0][ii];
            (*Hess)[0][1][v] = (*Hess)[0][1][v] + Hess_elem[0][1][ii];
            (*Hess)[1][0][v] = (*Hess)[1][0][v] + Hess_elem[1][0][ii];
            (*Hess)[1][1][v] = (*Hess)[1][1][v] + Hess_elem[1][1][ii];
        }
    }


    // Moving nodes
    vector<vector<double>> H = Zeros(2,2);
    vector<vector<double>> H_inv = Zeros(2,2);
    vector<double> G(2);
    vector<double> xs_smooth(2);
    double det;
    double delta = 0.0;
    for (int n = 0; n<nv; n++){
        if (!no_move[n]){
        H[0][0] = (*Hess)[0][0][n];
        H[0][1] = (*Hess)[0][1][n];
        H[1][0] = (*Hess)[1][0][n];
        H[1][1] = (*Hess)[1][1][n];
        G[0] = (*Grads)[0][n];
        G[1] = (*Grads)[1][n];

        det = H[0][0]*H[1][1] - H[1][0]*H[0][1];
        if (abs(det) > 1e-3){
            H_inv[0][0] = H[1][1]/det;
            H_inv[1][0] = -H[1][0]/det;
            H_inv[0][1] = -H[0][1]/det;
            H_inv[1][1] = H[0][0]/det;
            xs_smooth = H_inv*(-G);
            if (norm(xs_smooth) < h_max){
            mesh->coords[n] = mesh->coords[n] + xs_smooth;
            delta += xs_smooth[0]*xs_smooth[0] + xs_smooth[1]*xs_smooth[1];
            }
        }
        }
    }
    delete Grads, Hess;

    delta = sqrt(delta);
    return delta;
}




double isometry_energy_tri(vector<vector<double>>& xs, vector<vector<double>> &Grad, vector<vector<vector<double>>> &Hess, double refarea, double mu){

// computing angle based energy
double e12[2],e12_orth[2],e23[2],e23_orth[2],e31[2],e31_orth[2], sql12, sql23, sql31;
e12[0] = xs[1][0]-xs[0][0]; e12[1] = xs[1][1]-xs[0][1];
sql12 = e12[0]*e12[0] + e12[1]*e12[1];
e23[0] = xs[2][0]-xs[1][0]; e23[1] = xs[2][1]-xs[1][1];
sql23 = e23[0]*e23[0] + e23[1]*e23[1];
e31[0] = xs[0][0]-xs[2][0]; e31[1] = xs[0][1]-xs[2][1];
sql31 = e31[0]*e31[0] + e31[1]*e31[1];
double area2 = (e12[0]*e23[1] - e12[1]*e23[0]);
double area = 0.5*area2;
e12_orth[0] = -e12[1]; e12_orth[1] = e12[0];
e23_orth[0] = -e23[1]; e23_orth[1] = e23[0];
e31_orth[0] = -e31[1]; e31_orth[1] = e31[0];
double cts[3] = {1/(sqrt(3)*area),1/(sqrt(3)*area),1/(sqrt(3)*area)};

double Energy = (cts[2]*sql12 + cts[0]*sql23 + cts[1]*sql31);


// computing angle based gradient
int i,j;
for (i=0; i<2; i++){
    for(j=0; j<3; j++){
        Grad[i][j] = 0.0;
    }
}

Grad[0][0] = -2*cts[2]*e12[0];
Grad[1][0] = -2*cts[2]*e12[1];
Grad[0][1] = 2*cts[2]*e12[0];
Grad[1][1] = 2*cts[2]*e12[0];

Grad[0][1] += -2*cts[0]*e23[0];
Grad[1][1] += -2*cts[0]*e23[1];
Grad[0][2] = 2*cts[0]*e23[0];
Grad[1][2] = 2*cts[0]*e23[1];

Grad[0][2] += -2*cts[1]*e31[0];
Grad[1][2] += -2*cts[1]*e31[1];
Grad[0][0] += 2*cts[1]*e31[0];
Grad[1][0] += 2*cts[1]*e31[1];
double energy_a = Energy/area2;

Grad[0][0] -= energy_a*e23_orth[0];
Grad[1][0] -= energy_a*e23_orth[1];
Grad[0][1] -= energy_a*e31_orth[0];
Grad[1][1] -= energy_a*e31_orth[1];
Grad[0][2] -= energy_a*e12_orth[0];
Grad[1][2] -= energy_a*e12_orth[1];

// Computing Hessian
double c = 2*(cts[2]+cts[1]);
double C[2][2];
C[0][0] = c;
C[0][1] = 0.0;
C[1][0] = 0.0;
C[1][1] = c;

for (i = 0; i<2; i++){
    for (j = 0; j<2; j++){
        Hess[i][j][0] = C[i][j] - (Grad[i][0]*e23_orth[j]/area2 + Grad[j][0]*e23_orth[i]/area2);
    }
}
for (i = 0; i<2; i++){
    for (j = 0; j<2; j++){
        Hess[i][j][1] = C[i][j] - (Grad[i][1]*e31_orth[j]/area2 + Grad[j][1]*e31_orth[i]/area2);
    }
}
for (i = 0; i<2; i++){
    for (j = 0; j<2; j++){
        Hess[i][j][2] = C[i][j] - (Grad[i][2]*e12_orth[j]/area2 + Grad[j][2]*e12_orth[i]/area2);
    }
}


// compute area energy
if(mu > 0){
    c = 2*refarea/area2;
    Energy = (1-mu)*Energy + mu*(area2/(2*refarea) + c);

    c = c / area2;
    double coeff = mu* (1/(2*refarea) - c);
    Grad[0][0] = (1-mu)*Grad[0][0] + coeff*e23_orth[0];
    Grad[1][0] = (1-mu)*Grad[1][0] + coeff*e23_orth[1];

    Grad[0][1] = (1-mu)*Grad[0][1] + coeff*e31_orth[0];
    Grad[1][1] = (1-mu)*Grad[1][1] + coeff*e31_orth[1];

    Grad[0][2] = (1-mu)*Grad[0][2] + coeff*e12_orth[0];
    Grad[1][2] = (1-mu)*Grad[1][2] + coeff*e12_orth[1];

    c = c / area2; 
    coeff = mu*2*c;
    for (i = 0; i<2; i++){
        for (j = 0; j<2; j++){
            Hess[i][j][0] = (1-mu)*Hess[i][j][0] + coeff*(e23_orth[i]*e23_orth[j]);
        }
    }
    for (i = 0; i<2; i++){
        for (j = 0; j<2; j++){
            Hess[i][j][1] = (1-mu)*Hess[i][j][1] + coeff*(e31_orth[i]*e31_orth[j]);
        }
    }
    for (i = 0; i<2; i++){
        for (j = 0; j<2; j++){
            Hess[i][j][2] = (1-mu)*Hess[i][j][2] + coeff*(e12_orth[i]*e12_orth[j]);
        }
    }
}

return Energy;
}