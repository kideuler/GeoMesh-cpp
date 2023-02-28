#include "GeoMesh.hpp"
#include <unistd.h>
using namespace std;

double isometry_energy_tri(vector<vector<double>>& xs, double refarea, double mu);
double isometry_energy_tri(vector<vector<double>>& xs, vector<vector<double>> &Grad, vector<vector<vector<double>>> &Hess, double refarea, double mu);

void Mesh::mesh_smoothing_2d(vector<bool> no_move, int niters, double mu,vector<double> refareas){

    // calculate CRS one-ring elems
    if (stl.index.size() == 0){
        compute_Onering();
    }

    // main smoothing algorithm
    int iter = 0;
    double Energy = 1e10;
    double Energy_old = 1e12;
    
    while (iter < niters && Energy_old > Energy){
        // iteration of Energy based smoothing
        Energy_old = Energy;
        Energy = mesh_smoothing_tri_2d_iter(no_move, mu, refareas, Energy_old);

        iter++;
    }
}


double Mesh::mesh_smoothing_tri_2d_iter(vector<bool> no_move, double mu, vector<double> refareas, double Energy_old){

    int nv = coords.size();
    vector<vector<double>>* Grads = new vector<vector<double>>;
    *Grads = Zeros<double>(2,nv);
    vector<vector<vector<double>>>* Hess = new vector<vector<vector<double>>>;
    (*Hess).resize(2);
    (*Hess)[0].resize(2);
    (*Hess)[1].resize(2);
    (*Hess)[0][0].resize(nv);
    (*Hess)[1][0].resize(nv);
    (*Hess)[0][1].resize(nv);
    (*Hess)[1][1].resize(nv);
    vector<vector<double>> ps = Zeros<double>(3,2);
    vector<vector<double>> Grad_elem = Zeros<double>(2,3);
    vector<vector<vector<double>>> Hess_elem;
    Hess_elem.resize(2);
    Hess_elem[0].resize(2);
    Hess_elem[1].resize(2);
    Hess_elem[0][0].resize(3);
    Hess_elem[1][0].resize(3);
    Hess_elem[0][1].resize(3);
    Hess_elem[1][1].resize(3);

    // Accurmulating Energy
    double Energy_total = 0.0;
    double Energy;
    int v;
    bool use_area = refareas.size() != 0;
    for (int n = 0; n<nelems; n++){
        ps[0] = coords[elems[n][0]];
        ps[1] = coords[elems[n][1]];
        ps[2] = coords[elems[n][2]];

        if (use_area){
            Energy = isometry_energy_tri(ps, Grad_elem, Hess_elem, refareas[n], mu);
        } else {
            Energy = isometry_energy_tri(ps, Grad_elem, Hess_elem, 0.0, mu);
        }
        for (int ii = 0; ii<3; ii++){
            v = elems[n][ii];
            (*Grads)[0][v] = (*Grads)[0][v] + Grad_elem[0][ii];
            (*Grads)[1][v] = (*Grads)[1][v] + Grad_elem[1][ii];

            (*Hess)[0][0][v] = (*Hess)[0][0][v] + Hess_elem[0][0][ii];
            (*Hess)[0][1][v] = (*Hess)[0][1][v] + Hess_elem[0][1][ii];
            (*Hess)[1][0][v] = (*Hess)[1][0][v] + Hess_elem[1][0][ii];
            (*Hess)[1][1][v] = (*Hess)[1][1][v] + Hess_elem[1][1][ii];
        }
        Energy_total += Energy;
    }


    // Moving nodes
    // Energy_total < Energy_old
    vector<vector<double>> H = Zeros<double>(2,2);
    vector<vector<double>> H_inv = Zeros<double>(2,2);
    vector<double> G(2);
    vector<double> xs_smooth(2);
    vector<vector<double>> xs_diff = Zeros<double>(nv,2);
    double det;
    double alpha;
    int iter;
    for (int n = 0; n<nv; n++){
        if (!no_move[n]){
        H[0][0] = (*Hess)[0][0][n];
        H[0][1] = (*Hess)[0][1][n];
        H[1][0] = (*Hess)[1][0][n];
        H[1][1] = (*Hess)[1][1][n];
        G[0] = -(*Grads)[0][n];
        G[1] = -(*Grads)[1][n];

        det = H[0][0]*H[1][1] - H[1][0]*H[0][1];
        if (abs(det) > 1e-3){
            H_inv[0][0] = H[1][1]/det;
            H_inv[1][0] = -H[1][0]/det;
            H_inv[0][1] = -H[0][1]/det;
            H_inv[1][1] = H[0][0]/det;
            xs_smooth = H_inv*(G);
            alpha = 0.5;
            iter = 0;
            while (!check_jacobians_node(n,xs_smooth*alpha) && iter<8){
                alpha = alpha/2;
                //cout << "naegative jacobian found while trying to move node: " << n << " with alpha: " << alpha << endl;
                iter++;
            }
            if (iter < 8 && norm(xs_smooth*alpha) < 1){
                xs_diff[n] = xs_smooth*alpha;
            }
        }
        }
    }
    delete Grads, Hess;


    // evaulating total energy for new nodal placement
    double Energy_new = 0;
    for (int n = 0; n<nelems; n++){
        ps[0] = coords[elems[n][0]]+xs_diff[elems[n][0]];
        ps[1] = coords[elems[n][1]]+xs_diff[elems[n][1]];
        ps[2] = coords[elems[n][2]]+xs_diff[elems[n][2]];
        if (use_area){
            Energy = isometry_energy_tri(ps, refareas[n], mu);
        } else {
            Energy = isometry_energy_tri(ps, 1.0, mu);
        }
        Energy_new += Energy;
    }
    Energy_new = Energy_new / (double) nelems;
    Energy_old = Energy_total / (double) nelems;
    if (Energy_new < Energy_old){
        coords += xs_diff;
    }
    cout << "Energy old: " << Energy_old << " Energy new: " << Energy_new << endl;
    return Energy_new;
}

bool Mesh::check_jacobians_node(int vid, vector<double> dir){
    const vector<vector<double>> dphi = {{-1,-1},{1,0},{0,1}};
    vector<vector<double>> ps =Zeros<double>(3,2);
    vector<vector<double>> J = Zeros<double>(2,2);
    int hvid,eid,lvid;
    double detJ;
    for (int v = stl.index[vid]; v<stl.index[vid+1]-1; v++){
        hvid = stl.hvids[v];
        eid = hfid2eid(hvid)-1;
        lvid = hfid2lid(hvid)-1;
        ps[0] = coords[vid]+dir;
        ps[1] = coords[elems[eid][(lvid+1)%3]];
        ps[2] = coords[elems[eid][(lvid+2)%3]];
        J = Transpose(ps)*dphi;
        detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (detJ < 0){
            return false;
        }
    }

    return true;
}

double isometry_energy_tri(vector<vector<double>>& xs, double refarea, double mu){

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

if (mu > 0){
    double c = 2*refarea/area2;
    Energy = (1-mu)*Energy + mu*(area2/(2*refarea) + c);
}
return Energy;
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