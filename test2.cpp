#include "GeoMesh.hpp"
using namespace std;

vector<vector<double>> Circle(int npoints){
    vector<vector<double>> xs = Zeros<double>(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.49*cos(t)+0.5;
        xs[i][1] = 0.49*sin(t)+0.5;
    }
    return xs;
}

vector<vector<double>> Ellipse(int npoints){
    vector<vector<double>> xs = Zeros<double>(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.2*cos(t)+0.5;
        xs[i][1] = 0.1*sin(t)+0.5;
    }
    return xs;
}

vector<vector<double>> Flower(int npoints){
    vector<vector<double>> xs = Zeros<double>(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = (0.25 + 0.1*sin(5*t))*cos(t)+0.5;
        xs[i][1] = (0.25 + 0.1*sin(5*t))*sin(t)+0.5;
    }
    return xs;
}
/* TODO:
* convert delaunay refinement to struct function

*/

int main(){
    int n = 80;
    Spline spl;
    spl.nv = 0; 

    vector<vector<double>> xs = Circle(n);
    vector<vector<int>> segs = Zeros<int>(n,2);
    double h = 0.0;
    for (int i = 0; i<n; i++){
        segs[i][0] = i;
        segs[i][1] = (i+1)%n;
        h = h + norm(xs[(i+1)%n]-xs[i]);
    }
    h = 1*(sqrt(3)/3)*(h/(double(n)-1));

    Mesh msh = GeoMesh_Delaunay_Mesh(xs);
    GeoMesh_refine(&msh, h, &spl);
    WrtieVtk_tri(msh);
}