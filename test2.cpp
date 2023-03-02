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
    /*
    // loading bunny mesh
    string filename = "data/bunny.obj";
    Mesh bunny = ReadObj_tri(filename);
    bunny.compute_AHF();
    vector<double> refareas(bunny.nelems);
    vector<vector<double>> ps = Zeros<double>(3,2);
    double total=0.0;
    for (int i = 0; i<bunny.nelems; i++){
        ps[0] = bunny.coords[bunny.elems[i][0]];
        ps[1] = bunny.coords[bunny.elems[i][1]];
        ps[2] = bunny.coords[bunny.elems[i][2]];
        refareas[i] = area_tri(ps);
        total += refareas[i];
    }
    refareas = refareas / total;
    vector<vector<double>> P =  Parametric_Mapping(&bunny);
    refareas *= M_PI*0.25;
    Mesh S = bunny;
    S.coords = P;
    cout << "finished writing to file" << endl;
    vector<bool> bnd;
    bnd.assign(S.coords.size(), false);
    cout << "minimum angle: " << check_minangle(&S) << endl;
    check_jacobians(&S);
    S.mesh_smoothing_2d(S.find_boundary_nodes(), 1500, 0.8, refareas);
    cout << "minimum angle: " << check_minangle(&S) << endl;
    WrtieVtk_tri(S);
    cout << "finished writing to file" << endl;
    */

    int n = 300;
    vector<vector<double>> xs = Circle(n);
    Spline spl;
    spl.nv = 0; 
    double h = 0.0;
    vector<vector<int>> segs = Zeros<int>(n,2);
    for (int i = 0; i<n; i++){
        segs[i][0] = i;
        segs[i][1] = (i+1)%n;
        h = h + norm(xs[(i+1)%n]-xs[i]);
    }
    h = 1*(sqrt(3)/3)*(h/(double(n)-1));

    Mesh msh = GeoMesh_Delaunay_Mesh(xs);
    msh.Delaunay_refine(h);
    cout << "minimum angle: " << check_minangle(&msh) << endl;
    vector<bool> free_bnd(msh.coords.size());
    msh.mesh_smoothing_2d(free_bnd, 50, 0.0);
    cout << "minimum angle: " << check_minangle(&msh) << endl;
    check_jacobians(&msh);
    
    
    // writing
    WrtieVtk_tri(msh, "test.vtk");
    cout << "finished writing to file" << endl;
}