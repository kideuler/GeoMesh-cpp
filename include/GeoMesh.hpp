#include "Matrix.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <functional>
using namespace std;

struct Spline {
    int nv;
    int degree;
    vector<vector<double>> coords;
    vector<vector<double>> xweights;
    vector<vector<double>> yweights;
    vector<double> params;
};

struct Stencil {
    vector<int> index;
    vector<int> hvids;
};

struct Mesh {
    vector<vector<double>> coords;
    vector<vector<double>> normals;
    vector<double> param;
    vector<vector<int>> elems;
    vector<vector<int>> sibhfs;
    int nelems;
    vector<vector<int>> facets;
    vector<bool> delete_elem;
    vector<int> badtris;
    vector<int> vedge;
    vector<bool> bwork;
    vector<bool> on_boundary;
    Spline spl;
    Stencil stl;

    // major member functions
    void compute_AHF();
    void compute_Onering(int maxne = 10);
    void Delaunay_refine(function<double(vector<double>)> r_ref);
    void Delaunay_refine(double r_ref);
    void mesh_smoothing_2d(vector<bool> no_move, int niters = 5, double mu = 0.0, vector<double> refareas={});
    double mesh_smoothing_tri_2d_iter(vector<bool> no_move,  double mu, vector<double> refareas, double Energy_old);

    // minor utility functions
    void delete_tris();
    vector<bool> find_boundary_nodes();
    bool find_enclosing_tri(int *starting_tri, vector<double> &ps);
    void Flip_Insertion(int* vid, int starting_tri);
    void Flip_Insertion_segment(int vid, int hfid);
    void flip_edge(int eid, int lid);
    bool inside_diametral(int hfid, vector<double> &ps);
    bool check_jacobians_node(int vid, vector<double> dir = {0.0,0.0});
};

// Surface remeshing
void Parametric2Surface(Mesh *Surf, vector<vector<double>> &params, Mesh *msh);
vector<vector<double>> Parametric_Mapping(Mesh* Surf, int domain=0, int algo=0);
void Compute_normals(Mesh* Surf);
vector<double> Average_nodal_edgelength(Mesh* Surf, vector<vector<double>> &params);

// spline functions
Spline spline_init(const vector<vector<double>> &xs, int degree = 3);
vector<double> spline_var(Spline* spl, double t, int order=0);
double spline_curvature(Spline* spl, double t);
vector<double> spline_point_segment(Spline* spl, double a, double b, double ratio);

// delaunay Mesh functions 2d
Mesh GeoMesh_Delaunay_Mesh(vector<vector<double>> &xs, vector<double> &params);
Mesh GeoMesh_Delaunay_Mesh(const vector<vector<int>> &segs, vector<vector<double>> &xs);
Mesh GeoMesh_Delaunay_Mesh(vector<vector<double>> &xs);
void Bowyer_watson2d(Mesh* DT, int vid, int tri_s,bool refine);
double check_minangle(Mesh* DT);
double area_tri(const vector<vector<double>> &xs);
vector<bool> find_boundary_nodes(Mesh* DT);
bool inside_tri(const vector<vector<double>> &xs, const vector<double> &ps);

// Meshutils functions
vector<vector<int>> find_boundary(Mesh* msh, bool findloop);
void WrtieVtk_tri(const Mesh &msh);
void WrtieVtk_tri(const Mesh &msh, const vector<double> &data);
Mesh ReadObj_tri(string filename);
bool check_sibhfs(Mesh* DT);
bool check_jacobians(Mesh* DT);

// small functions to be used in multiple files
vector<double> min_array(const vector<vector<double>> &xs);
vector<double> max_array(const vector<vector<double>> &xs);
int hfid2eid(int hfid);
int hfid2lid(int hfid);
int elids2hfid(int eid, int lid);

struct stack
{
    int hfid;
    struct stack *next;
};
void push_stack(stack** head, int hfid);
void pop_stack(stack** head);


// Functions necessary for Blossom algorithm of triangle stitching
struct Blossom {
    int n, m;
    vector<int> mate;
    vector<vector<int>> b;
    vector<int> p, d, bl;
    vector<vector<int>> G;
    vector<vector<int>> E;
    Blossom(int n) : n(n) {
        m = n + n / 2;
        mate.assign(n, -1);
        b.resize(m);
        p.resize(m);
        d.resize(m);
        bl.resize(m);
        G.assign(m, vector<int>(m, -1));
        E.assign(m, vector<int>(m, -1));
    }
    void add_edge(int u, int v);
    void add_edge(int u, int v, int lid);
    void match(int u, int v);
    vector<int> trace(int x);
    void contract(int c, int x, int y, vector<int> &vx, vector<int> &vy);
    vector<int> lift(vector<int> &vx);
    int solve();
};

Blossom Mesh2Graph(Mesh* msh);
void Tris2quads_blossom(Mesh *msh, Blossom* B);