#include "GeoMesh.hpp"
#include <Eigen/Sparse>
#include <queue>
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat;
typedef  Eigen::Triplet<double> T;

static void Circle_param(double t, double &x, double &y);
static void Square_param(double t, double &x, double &y);
static void Harmonic_map(vector<vector<double>> &elemat, vector<vector<double>> &ps);
static void Convex_Combo_Map(vector<vector<double>> &elemat);
static bool find_enclosing_tri(Mesh* DT,vector<vector<double>> &params, int* tri, vector<double> ps);

void Parametric2Surface(Mesh *Surf, vector<vector<double>> &params, Mesh *msh){
    int nv = msh->coords.size();
    bool use_cubic = Surf->normals.size() > 0;
    vector<vector<double>> coords_new = Zeros<double>(nv,3);

    // main loop
    int tri = 0;
    bool passed;
    vector<vector<double>> xi_u = {{0,0},{0,0}};
    vector<double> xi_eta(2);
    vector<vector<double>> ps = {{0,0},{0,0},{0,0}};
    vector<vector<double>> xs = {{0,0,0},{0,0,0},{0,0,0}};
    double det_xi_u,xi,eta,zeta;
    double w;
    vector<double> a300(3), a030(3), a003(3), a210(3), a120(3), a021(3), a012(3), a102(3), a201(3), E(3), V(3), a111(3), n1(3), n2(3), n3(3);
    for (int i = 0; i<nv; i++){
        passed = find_enclosing_tri(Surf ,params, &tri, msh->coords[i]);
        if (passed){
            ps[0] = params[Surf->elems[tri][0]];
            ps[1] = params[Surf->elems[tri][1]];
            ps[2] = params[Surf->elems[tri][2]];

            xs[0] = Surf->coords[Surf->elems[tri][0]];
            xs[1] = Surf->coords[Surf->elems[tri][1]];
            xs[2] = Surf->coords[Surf->elems[tri][2]];
            det_xi_u = (ps[1][0]-ps[0][0])*(ps[2][1]-ps[0][1]) - (ps[1][1]-ps[0][1])*(ps[2][0]-ps[0][0]);
            xi_u[0][0] = (ps[2][1]-ps[0][1])/det_xi_u;
            xi_u[0][1] = -(ps[2][0]-ps[0][0])/det_xi_u;
            xi_u[1][0] = -(ps[1][1]-ps[0][1])/det_xi_u;
            xi_u[1][1] = (ps[1][0]-ps[0][0])/det_xi_u;
            xi_eta = xi_u*(msh->coords[i] - ps[0]);

            // Linear interpolation for now
            if (use_cubic){
                n1 = Surf->normals[Surf->elems[tri][0]];
                n2 = Surf->normals[Surf->elems[tri][1]];
                n3 = Surf->normals[Surf->elems[tri][2]];
                a300 = xs[0]; a030 = xs[1]; a003 = xs[2];

                w = inner(xs[1]-xs[0],n1);
                a210 = (2.0*xs[0] + xs[1] - w*n1)/3.0;
                w = inner(xs[0]-xs[1],n2);
                a120 = (2.0*xs[1] + xs[0] - w*n2)/3.0;

                w = inner(xs[2]-xs[1],n2);
                a021 = (2.0*xs[1] + xs[2] - w*n2)/3.0;
                w = inner(xs[1]-xs[2],n3);
                a012 = (2.0*xs[2] + xs[1] - w*n3)/3.0;

                w = inner(xs[0]-xs[2],n3);
                a102 = (2.0*xs[2] + xs[0] - w*n3)/3.0;
                w = inner(xs[2]-xs[0],n1);
                a201 = (2.0*xs[0] + xs[2] - w*n1)/3.0;

                E = (a210 + a120 + a021 + a012 + a102 + a201)/6.0;
                V = (xs[0] + xs[1] + xs[2])/3.0;
                a111[0] = E[0] + (E[0]-V[0])/2;
                a111[1] = E[1] + (E[1]-V[1])/2;
                a111[2] = E[2] + (E[2]-V[2])/2;
                xi = xi_eta[0]; eta = xi_eta[1];
                zeta = 1-xi-eta;

                coords_new[i] = zeta*a300 + xi*a030 + eta*a300 + \
                3.0*zeta*zeta*xi*a210 +  3.0*zeta*xi*xi*a120 + 3.0*zeta*zeta*eta*a201 + \
                3.0*xi*xi*eta*a021 + 3.0*zeta*eta*eta*a102 + 3.0*xi*eta*eta*a201 + 6.0*xi*eta*zeta*a111;
                cout << coords_new[i][0] << " " << coords_new[i][1] << " " << coords_new[i][2] << endl;
            } else{
                coords_new[i] = (1-xi_eta[0]-xi_eta[1])*xs[0] + xi_eta[0]*xs[1] + xi_eta[1]*xs[2];
            }
        }
    }

    msh->coords = coords_new;
    return;
}

vector<vector<double>> Parametric_Mapping(Mesh* Surf, int domain, int algo){

    int nv = Surf->coords.size();
    vector<vector<double>> params = Zeros<double>(nv,2);
    vector<vector<int>> bdy = find_boundary(Surf, true);
    int nb = bdy.size();
    bool* bdymask = new bool[nv];
    for(int i = 0; i<nv; i++){ bdymask[i]=false; }
    for(int i = 0; i<nb; i++){ bdymask[bdy[i][0]] = true; bdymask[bdy[i][1]] = true; }

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Upper|Eigen::Lower, Eigen::COLAMDOrdering<int> >   solver;
    SpMat A(nv,nv);
    A.reserve(Eigen::VectorXi::Constant(nv,10));
    Eigen::VectorXd bx(nv), by(nv), u(nv), v(nv);

    int vn;
    double t,x,y;
    for(int i = 0; i<nb; i++){
        vn = bdy[i][0];
        t = (((double) i) / ((double) nb));
        if (domain == 0){
            Circle_param(t,x,y);
        } else if(domain == 1){
            Square_param(t,x,y);
        }
        
        bx(vn) = x;
        by(vn) = y;
        A.coeffRef(vn,vn) = 1.0;
    }
    cout << "setup rhs vectors" << endl;
    vector<vector<double>> ps = {{0,0,0},{0,0,0},{0,0,0}};
    vector<vector<double>> elemat = {{0,0,0},{0,0,0},{0,0,0}};
    for(int ii = 0; ii<Surf->nelems; ii++){
        ps[0] = Surf->coords[Surf->elems[ii][0]];
        ps[1] = Surf->coords[Surf->elems[ii][1]];
        ps[2] = Surf->coords[Surf->elems[ii][2]];

        if (algo == 0){
            Harmonic_map(elemat,ps);
        } else if (algo == 1){
            Convex_Combo_Map(elemat);
        }
        for(int i = 0; i<3; i++){
            if(!bdymask[Surf->elems[ii][i]]){
                for(int j = 0; j<3; j++){
                    A.coeffRef(Surf->elems[ii][i], Surf->elems[ii][j]) += elemat[i][j];
                }
            }
        }
    }
    A.makeCompressed();
    cout << "assembled matrix" << endl;

    solver.compute(A);

    u = solver.solve(bx); 
    v = solver.solve(by);

    for(int i = 0; i<nv; i++){
        params[i][0] = u(i);
        params[i][1] = v(i);
    }

    delete bdymask;
    return params;
}

vector<double> Average_nodal_edgelength(Mesh* Surf, vector<vector<double>> &params){
    int nv = params.size();

    vector<double> hnode(nv);
    double* nnode = new double[nv];
    for(int i = 0; i<nv; i++) {hnode[i]=0.0; nnode[i]=0.0; }
    double l;
    int v1,v2;
    for (int i = 0; i<Surf->nelems; i++){
        for (int j = 0; j<3; j++){
            v1 = Surf->elems[i][j];
            v2 = Surf->elems[i][(j+1)%3];
            l = norm(params[v2]-params[v1]);
            hnode[v1] += l;
            hnode[v2] += l;   
            nnode[v1] += 1;
            nnode[v2] += 1;         
        }
    }   

    for(int i = 0; i<nv; i++){
        hnode[i] = hnode[i]/nnode[i];
    }

    delete nnode;
    return hnode;
}

void Compute_normals(Mesh* Surf){
    int nv = Surf->coords.size();
    Surf->normals = Zeros<double>(nv,3);
    vector<double> u(3);
    vector<double> v(3);
    vector<double> n(3);
    double dotuv, ul, vl, w, nl;

    for(int i = 0; i<Surf->nelems; i++){
        for(int j = 0; j<3; j++){
            u = Surf->coords[Surf->elems[i][(j+1)%3]]-Surf->coords[Surf->elems[i][j]];
            v = Surf->coords[Surf->elems[i][(j+2)%3]]-Surf->coords[Surf->elems[i][j]];
            dotuv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
            ul = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
            vl = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            w = acos(dotuv/(ul*vl));
            n = {u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2] ,u[0]*v[1] - u[1]*v[0]};
            nl = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

            Surf->normals[Surf->elems[i][j]][0] += (w/(2*M_PI))*n[0]/nl;
            Surf->normals[Surf->elems[i][j]][1] += (w/(2*M_PI))*n[1]/nl;
            Surf->normals[Surf->elems[i][j]][2] += (w/(2*M_PI))*n[2]/nl;
        }
    }
    return;
}

static void Harmonic_map(vector<vector<double>> &elemat, vector<vector<double>> &ps){
    // can be much more efficient but works
    double x1 = ps[0][0]; double y1 = ps[0][1]; double z1 = ps[0][2];
    double x2 = ps[1][0]; double y2 = ps[1][1]; double z2 = ps[1][2];
    double x3 = ps[2][0]; double y3 = ps[2][1]; double z3 = ps[2][2];
    double area = area_tri(ps);
    vector<vector<double>> dphi = {{-1,-1},{1,0},{0,1}};
    vector<vector<double>> M_inv = Zeros<double>(2,2);
    M_inv[0][0] = pow(x1-x3,2.0)+pow(y1-y3,2.0)+pow(z1-z3,2.0);
    M_inv[0][1] = -(x1-x2)*(x1-x3)-(y1-y2)*(y1-y3)-(z1-z2)*(z1-z3);
    M_inv[1][0] = -(x1-x2)*(x1-x3)-(y1-y2)*(y1-y3)-(z1-z2)*(z1-z3);
    M_inv[1][1] = pow(x1-x2,2.0)+pow(y1-y2,2.0)+pow(z1-z2,2.0);
    elemat = (dphi*M_inv*Transpose(dphi)) * (1/(4*area));
    return;
}

static void Convex_Combo_Map(vector<vector<double>> &elemat){
    elemat = {{2,-1,-1},{-1,2,-1},{-1,-1,2}};
    return;
}

static void Circle_param(double t, double &x, double &y){
    x = 0.5*cos(2*M_PI*t)+0.5;
    y = 0.5*sin(2*M_PI*t)+0.5;
    return;
}
static void Square_param(double t, double &x, double &y){
    if (t<=0.25){
        x = 4*t;
        y = 0;
    } else if (t<=0.5){
        x = 1.0;
        y = 4*(t-0.25);
    } else if (t<=0.75){
        x = 1-4*(t-0.5);
        y = 1.0;
    } else {
        x = 0.0;
        y = 1-4*(t-0.75);
    }
    return;
}

bool find_enclosing_tri(Mesh* DT, vector<vector<double>> &params, int* tri, vector<double> ps){
    int v1,v2,v3,i,hfid;
    bool stop;
    int iters = 0;
    i = 0;
    stop = false;
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    while (!stop){
        v1 = DT->elems[*tri][0];
        v2 = DT->elems[*tri][1];
        v3 = DT->elems[*tri][2];
        if (DT->delete_elem[*tri]){
            cout << "deleted tri passed ran through in find_enclosing_tri, should not be" << endl;

        }
        xs[0] = params[v1];
        xs[1] = params[v2];
        xs[2] = params[v3];
        if (inside_tri(xs,ps)){
            stop = true;
            return true;
        } else {
            double AB[2] = {params[v2][0]-params[v1][0], params[v2][1]-params[v1][1]};
            double BC[2] = {params[v3][0]-params[v2][0], params[v3][1]-params[v2][1]};
            double CA[2] = {params[v1][0]-params[v3][0], params[v1][1]-params[v3][1]};
            double AP[2] = {ps[0]-params[v1][0], ps[1]-params[v1][1]};
            double BP[2] = {ps[0]-params[v2][0], ps[1]-params[v2][1]};
            double CP[2] = {ps[0]-params[v3][0], ps[1]-params[v3][1]};
            double N1[2] = {AB[1],-AB[0]};
            double N2[2] = {BC[1],-BC[0]};
            double N3[2] = {CA[1],-CA[0]};
            double S1 = AP[0]*N1[0]+AP[1]*N1[1];
            double S2 = BP[0]*N2[0]+BP[1]*N2[1];
            double S3 = CP[0]*N3[0]+CP[1]*N3[1];
            if ((S1>0)&&(S1>=S2)&&(S1>=S3)){
                hfid = DT->sibhfs[*tri][0];
                if (hfid != 0){
                    *tri = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tri = elids2hfid(*tri+1,1);
                    return false;
                }
            } else if ((S2>0)&&(S2>=S1)&&(S2>=S3)) {
                hfid = DT->sibhfs[*tri][1];
                if (hfid != 0){
                    *tri = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tri =  elids2hfid(*tri+1,2);
                    return false;
                }
            } else if ((S3>0)&&(S3>=S1)&&(S3>=S2)){
                hfid = DT->sibhfs[*tri][2];
                if (hfid != 0){
                    *tri = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tri = elids2hfid(*tri+1,3);
                    return false;
                }
            } 
            
        }
        iters++;
    }
    *tri = -1;
    return false;
    
}