#include <GeoComp.hpp>
#include <cmath>
#include <chrono>


#define FAIL "\033[0;31mFAIL\033[0m"
#define PASS "\033[0;32mPASS\033[0m"

using namespace std;

vec uexpr(Mat xs){
    int nv = xs.size();
    vec u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = xs[i][0]*xs[i][1];
    }
    return u;
}

Mat Circle(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.49*cos(t)+0.5;
        xs[i][1] = 0.49*sin(t)+0.5;
    }
    return xs;
}

Mat Ellipse(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.2*cos(t)+0.5;
        xs[i][1] = 0.1*sin(t)+0.5;
    }
    return xs;
}

Mat Flower(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = (0.25 + 0.1*sin(5*t))*cos(t)+0.5;
        xs[i][1] = (0.25 + 0.1*sin(5*t))*sin(t)+0.5;
    }
    return xs;
}

Mat Star(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.25*pow(cos(t),3)+0.5;
        xs[i][1] = 0.25*pow(sin(t),3)+0.5;
    }
    return xs;
}

Mat Cardiod(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.25*(2*cos(t) - cos(2*t))+0.5;
        xs[i][1] = 0.25*(2*sin(t) - sin(2*t))+0.5;
    }
    return xs;
}


Mat Box(int npoints){
    Mat xs = Zeros(npoints*npoints,2);
    int k = 0;
    for (int i =0;i<npoints;i++){
        for (int j=0; j<npoints;j++){
            xs[k][0] = ((double) j)/((double) npoints-1);
            xs[k][1] = ((double) i)/((double) npoints-1);
            k++;
        }
    }
    return xs;
}

Mat rBox(int npoints, double bmin=0.0, double bmax=1.0){
    double h = (bmax-bmin)/sqrt((double)npoints);
    Mat xs = rMat(npoints-4,2,bmin+2*h,bmax-2*h);
    xs.push_back({bmin,bmin});
    xs.push_back({bmax,bmin});
    xs.push_back({bmax,bmax});
    xs.push_back({bmin,bmax});
    return xs;
}

double Gradiate(vec xs){
    if ((pow(xs[0]-0.5,2)+pow(xs[1]-0.5,2))< .1*.1 || abs(xs[0]-xs[1]) < 0.05){
        return 0.7*(sqrt(3)/3)*M_PI/(199);
    } else{
        return (sqrt(3)/3)*M_PI/(199);
    }
}

function<double(vec)> create_grad(const Mat &ps, double h, double ratio, double hgrad){
    function<double(vec)> hF = [ps,h,ratio,hgrad](vec xs){
        int nv = ps.size();
        double r,xi,alpha;
        alpha = 1000;
        for (int i = 0; i<nv; i++){
            r = (pow(xs[0]-ps[i][0],2)+pow(xs[1]-ps[i][1],2));
            xi = r/(hgrad*h);
            alpha = min((1-min(xi,1.0))*ratio*h + h*min(xi,1.0),alpha);
        }
        return alpha;
    };
    return hF;
}
function<double(vec)> create_grad(const Mat &ps, double h, vector<double> ratio, double hgrad){
    function<double(vec)> hF = [ps,h,ratio,hgrad](vec xs){
        int nv = ps.size();
        double r,xi,alpha;
        alpha = 1000;
        for (int i = 0; i<nv; i++){
            r = (pow(xs[0]-ps[i][0],2)+pow(xs[1]-ps[i][1],2));
            xi = r/(hgrad*h);
            alpha = min((1-min(xi,1.0))*ratio[i]*h + h*min(xi,1.0),alpha);
        }
        return alpha;
    };
    return hF;
}
function<double(vec)> create_curvature_grad(const Mat &ps, vector<double> K, double theta, double h, double h_min, double hgrad){
    int nv = ps.size();
    vector<double> H(nv);
    for (int i=0; i<nv; i++){
        H[i] = (theta*M_PI/180.0) / K[i];
        H[i] = min(max(h*h_min, H[i]), h);
    }

    function<double(vec)> hF = [ps,H,h,hgrad](vec xs){
        int nv = ps.size();
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

/**
 * STILL TO DO:
 *  - Take out as many std::vectors as possible to minimize allocations
 *  - Add Test cases from CAD from scratch and get them working
 *  - Add constriained delaunay Mesh support
 *  - Add spline geometry processing (general degree, NEED TO USE EIGEN FOR MATRIX STUFF)
 *  - Add BLossom implementation for graphs
 *  - Add recombination algorithm to generate quad meshes and add quad smoothing
 *  - 3D tetgen (STALLED FOR NOW)
 */
Mat picture();
int main(){
    // loading bunny mesh
    string filename = "bunny.obj";
    Mesh bunny = ReadObj_tri(filename);
    bunny.compute_AHF();
    Mat P =  Parametric_Mapping(&bunny);

    Mesh S = bunny;
    S.coords = P;
    WrtieVtk_tri(S);
    cout << "finished writing to file" << endl;
    vector<bool> bnd = find_boundary_nodes(&S);
    cout << "minimum angle: " << check_minangle(&S) << endl;
    vector<double> hnode = Average_nodal_edgelength(&bunny, P);
    function<double(vec)> H;
    mesh_smoothing_2d(&S, bnd, H, 0);
    cout << "minimum angle: " << check_minangle(&S) << endl;
    WrtieVtk_tri(S);
    cout << "finished writing to file" << endl;

    /*
    // compute average mesh density for a uniform surface mesh distribution
    vector<double> hnode = Average_nodal_edgelength(&bunny, P);
    // spline setup
    Mat sps = Circle(100);
    auto start = chrono::high_resolution_clock::now();
    Spline spl = spline_init(sps,3);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "created Spline in " << duration.count()/1e6 << " seconds" << endl;
    int n = 500;
    Mat xs  = Zeros(n,2);
    vector<double> param(n);
    for (int i=0; i<n; i++){
        param[i] = double(i)/double(n);
        xs[i] = spline_var(&spl, param[i]);
    }
    
    int kk = 251;
    Mat ps = Star(kk);
    vector<double> K(kk);
    for (int i=0; i<ps.size(); i++){
        double p = double(i)/double(kk);
        //K[i] = spline_curvature(&spl, p);
    }
    
    vector<vector<int>> segs = Zerosi(n,2);
    double h = 0.0;
    for (int i = 0; i<n; i++){
        segs[i][0] = i;
        segs[i][1] = (i+1)%n;
        h = h + norm(xs[(i+1)%n]-xs[i]);
    }
    h = 1*(sqrt(3)/3)*(h/(double(n)-1));
    hnode = hnode/h;
    for (int i = 0; i<hnode.size(); i++) { hnode[i] =  max(0.02, hnode[i]);}
    function<double(vec)> H = create_grad(P, h, hnode, 0.001);
    // initial Mesh
    start = chrono::high_resolution_clock::now();
    Mesh DT = GeoComp_Delaunay_Mesh(xs);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "finished initial Mesh in " << duration.count()/1e6 << " seconds" << endl;
    // refinement
    Spline S;
    start = chrono::high_resolution_clock::now();
    GeoComp_refine(&DT, h, &spl);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "finished delaunay refinement in " << duration.count()/1e6 << " seconds" << endl;
    cout << "minimum angle: " << check_minangle(&DT) << endl;
    vector<bool> bnd = find_boundary_nodes(&DT);
    // smoothing
    start = chrono::high_resolution_clock::now();
    mesh_smoothing_2d(&DT, bnd, H, 0);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "finished mesh smoothing in " << duration.count()/1e6 << " seconds" << endl;
    cout << "minimum angle: " << check_minangle(&DT) << endl;
    // surface remeshing
    Parametric2Surface(&bunny, P, &DT);
    // writing
    WrtieVtk_tri(DT);
    cout << "finished writing to file" << endl;
    */
}