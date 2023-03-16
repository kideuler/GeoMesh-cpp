#include "GeoMesh.hpp"
#include <cmath>

using namespace std;

vector<double> uexpr(vector<vector<double>> xs){
    int nv = xs.size();
    vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = 10*sin(2*xs[i][0])*cos(2*xs[i][1]);
    }
    return u;
}

vector<double> ulap(vector<vector<double>> xs, double kappa){
    int nv = xs.size();
    vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = -80*sin(2*xs[i][0])*cos(2*xs[i][1])*(-kappa);
    }
    return u;
}

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

vector<vector<double>> picture();

vector<vector<double>> Box(int npoints){
    vector<vector<double>> xs = Zeros<double>(npoints*npoints,2);
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

double Gradiate(vector<double> xs){
    if ((pow(xs[0]-0.5,2)+pow(xs[1]-0.5,2))< .1*.1 || abs(xs[0]-xs[1]) < 0.05){
        return 0.7*(sqrt(3)/3)*M_PI/(199);
    } else{
        return (sqrt(3)/3)*M_PI/(199);
    }
}

function<double(vector<double>)> create_grad(const vector<vector<double>> &ps, double h, double ratio, double hgrad){
    function<double(vector<double>)> hF = [ps,h,ratio,hgrad](vector<double> xs){
        int nv = ps.size();
        if (xs.size() == 0){
            return h*ratio;
        }
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

function<double(vector<double>)> create_grad(const vector<vector<double>> &ps, double h, vector<double> ratio, double hgrad){
    function<double(vector<double>)> hF = [ps,h,ratio,hgrad](vector<double> xs){
        int nv = ps.size();
        if (xs.size() == 0){
            return mean(ratio)*h;
        }
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

Mesh Ellipse_mesh(int npoints);

void test1(); // mesh generation and fem on ellipse
void test2(); // mesh generation and fem on curvature adaptive flower hole
void test3(); // mesh generation and fem on circle with mona-lisa embedded 
void test4(); // surface reconstruction of stanford bunny
void test5(); // testing branch for edge coloring and quad mesh gen
void test6(); // testing branch for constraine deltri

/**
 * STILL TO DO: create tests print timings and angles
 * - offcenters
 * - cubic hermite interpolation of surface points
 */
int main(int argc, char *argv[]){
    istringstream iss(argv[1]);
    int arg;
    iss >> arg;
    if (arg == 1){
        test1();
    } else if (arg == 2){
        test2();
    } else if (arg == 3){
        test3();
    } else if (arg == 4){
        test4();
    } else if (arg == 5){
        test5();
    } else if (arg == 6){
        test6();
    }
}

void test1(){
    cout << "performing test 1: meshing an ellipse" << endl << endl;
    Mesh msh = Ellipse_mesh(100);

        // FEM goes here
        cout << "Performing FEM with P1 elements" << endl;
        double kappa = 2.0;
        int nv = msh.coords.size();
        vector<double> u = uexpr(msh.coords); // true solution
        vector<double> frhs = ulap(msh.coords, kappa); // rhs expression
        vector<int> bndnodes = msh.boundary_nodes(); // finding boundary nodes
        int nbnd = bndnodes.size();
        vector<double> dvals(nbnd); // dirichlet values
        for (int n = 0; n<nbnd; n++) {dvals[n] = u[bndnodes[n]];}
        vector<double> u_h = Poisson_2d(&msh, frhs, kappa, bndnodes, dvals);
        vector<double> error(nv);
        double maxE = 0.0;
        for (int n = 0; n<nv; n++) {error[n] = abs(u[n]-u_h[n]); if(error[n]>maxE){maxE=error[n];}}
        double P1_error = norm(error)/norm(u);
        cout << "l2 relative error from Linear mesh: " << P1_error << endl;
        cout << "l_inf relative error from Linear mesh: " << maxE << endl << endl;
        
        cout << "Performing FEM with P2 elements" << endl;
        msh.make_quadratic();
        nv = msh.coords.size();
        u = uexpr(msh.coords); // true solution
        frhs = ulap(msh.coords, kappa); // rhs expression
        bndnodes = msh.boundary_nodes(); // finding boundary nodes
        nbnd = bndnodes.size();
        dvals.resize(nbnd);
        for (int n = 0; n<nbnd; n++) {dvals[n] = u[bndnodes[n]];}
        u_h = Poisson_2d(&msh, frhs, kappa, bndnodes, dvals);
        error.resize(nv);
        maxE = 0.0;
        for (int n = 0; n<nv; n++) {error[n] = abs(u[n]-u_h[n]); if(error[n]>maxE){maxE=error[n];}}
        double P2_error = norm(error)/norm(u);
        cout << "l2 relative error from Quadratic mesh: " << P2_error << endl;
        cout << "l_inf relative error from Linear mesh: " << maxE << endl << endl;

        // writing
        WrtieVtk_tri(msh, error, "test1.vtk");
        cout << "finished writing to file" << endl;
}

void test2(){
    cout << "performing test 2: meshing flower hole with curvature adaptivity" << endl;
        int npoints_spline = 250;

        // generating points on ellipse for the spline
        vector<vector<double>> Spline_points = Flower(npoints_spline);

        // creating spline degree 5
        auto start = chrono::high_resolution_clock::now();
        Spline spl = spline_init(Spline_points,5);
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "created Spline in " << duration.count()/1e6 << " seconds" << endl << endl;

        // interpolating points using spline
        int npoints = 80;
        vector<vector<double>> points = Zeros<double>(npoints,2);
        vector<vector<int>> segments = Zeros<int>(npoints,2);
        vector<double> param(npoints);
        for (int i = 0; i<npoints; i++){
            param[i] = double(i)/double(npoints);
            points[i] = spline_var(&spl, param[i]);
            segments[i][1] = i;
            segments[i][0] = (i+1)%npoints;
        }
        points.push_back({0,0});
        points.push_back({1,0});
        points.push_back({1,1});
        points.push_back({0,1});

        // calculating average edge length
        double h = 0.0;
        for (int i = 0; i<npoints; i++){
            h = h + norm(points[(i+1)%npoints]-points[i]);
        }
        h = 1*(sqrt(3)/3)*(h/(double(npoints)-1));

        // calculating local edgelength function
        function<double(vector<double>)> H = create_curvature_hfunction(&spl, 200, 10.0, h, 0.2, 0.2);
        /*
        Note: this is not an efficient way to do local refinement as the time complexity is O(n^2)
        In the future with more time I would like to implement a Marching squares or quadtree structure to
        interpolate local edge lengths
        */

        // initial Mesh
        start = chrono::high_resolution_clock::now();
        Mesh msh = GeoMesh_Delaunay_Mesh(segments,points);
        stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "finished initial Mesh in " << duration.count()/1e6 << " seconds" << endl << endl;
        msh.spl = spl;
        msh.param = param;
        
        // mesh refinement
        start = chrono::high_resolution_clock::now();
        msh.Delaunay_refine(H);
        stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "finished delaunay refinement in " << duration.count()/1e6 << " seconds" << endl;
        cout << "minimum angle: " << check_minangle(&msh) << endl << endl;

        // mesh smoothing
        start = chrono::high_resolution_clock::now();
        msh.mesh_smoothing_2d(msh.find_boundary_nodes(), 50, 0.0);;
        stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "finished mesh smoothing in " << duration.count()/1e6 << " seconds" << endl;
        cout << "minimum angle: " << check_minangle(&msh) << endl << endl;

        // FEM goes here
        cout << "Performing FEM with P1 elements" << endl;
        double kappa = 2.0;
        int nv = msh.coords.size();
        vector<double> u = uexpr(msh.coords); // true solution
        vector<double> frhs = ulap(msh.coords, kappa); // rhs expression
        vector<int> bndnodes = msh.boundary_nodes(); // finding boundary nodes
        int nbnd = bndnodes.size();
        vector<double> dvals(nbnd); // dirichlet values
        for (int n = 0; n<nbnd; n++) {dvals[n] = u[bndnodes[n]];}
        vector<double> u_h = Poisson_2d(&msh, frhs, kappa, bndnodes, dvals);
        vector<double> error(nv);
        double maxE = 0.0;
        for (int n = 0; n<nv; n++) {error[n] = abs(u[n]-u_h[n]); if(error[n]>maxE){maxE=error[n];}}
        double P1_error = norm(error)/norm(u);
        cout << "l2 relative error from Linear mesh: " << P1_error << endl;
        cout << "l_inf absolute error from Linear mesh: " << maxE << endl << endl;

        cout << "Performing FEM with P2 elements" << endl;
        msh.make_quadratic();
        nv = msh.coords.size();
        u = uexpr(msh.coords); // true solution
        frhs = ulap(msh.coords, kappa); // rhs expression
        bndnodes = msh.boundary_nodes(); // finding boundary nodes
        nbnd = bndnodes.size();
        dvals.resize(nbnd);
        for (int n = 0; n<nbnd; n++) {dvals[n] = u[bndnodes[n]];}
        u_h = Poisson_2d(&msh, frhs, kappa, bndnodes, dvals);
        error.resize(nv);
        maxE = 0.0;
        for (int n = 0; n<nv; n++) {error[n] = abs(u[n]-u_h[n]); if(error[n]>maxE){maxE=error[n];}}
        double P2_error = norm(error)/norm(u);
        cout << "l2 relative error from Quadratic mesh: " << P2_error << endl;
        cout << "l_inf absolute error from Quadratic mesh: " << maxE << endl << endl;

        // writing
        WrtieVtk_tri(msh, error, "test2.vtk");
        cout << "finished writing to file" << endl;
}

void test3(){
    cout << "performing test 3: meshing circle with mona lisa gradiation" << endl;
        int npoints_spline = 250;

        // generating points on ellipse for the spline
        vector<vector<double>> Spline_points = Circle(npoints_spline);

        // creating spline degree 3
        auto start = chrono::high_resolution_clock::now();
        Spline spl = spline_init(Spline_points,3);
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "created Spline in " << duration.count()/1e6 << " seconds" << endl << endl;

        // interpolating points using spline
        int npoints = 100;
        vector<vector<double>> points = Zeros<double>(npoints,2);
        vector<vector<int>> segments = Zeros<int>(npoints,2);
        vector<double> param(npoints);
        for (int i = 0; i<npoints; i++){
            param[i] = double(i)/double(npoints);
            points[i] = spline_var(&spl, param[i]);
            segments[i][0] = i;
            segments[i][1] = (i+1)%npoints;
        }

        // calculating average edge length
        double h = 0.0;
        for (int i = 0; i<npoints; i++){
            h = h + norm(points[(i+1)%npoints]-points[i]);
        }
        h = 1*(sqrt(3)/3)*(h/(double(npoints)-1));

        // calculating local edgelength function
        vector<vector<double>> ps = read_array("data/picture.dat");
        function<double(vector<double>)> H = create_grad(ps, h, 0.15, 0.01);
        /*
        Note: this is not an efficient way to do local refinement as the time complexity is O(n^2)
        In the future with more time I would like to implement a Marching squares or quadtree structure to
        interpolate local edge lengths
        */

        // initial Mesh
        start = chrono::high_resolution_clock::now();
        Mesh msh = GeoMesh_Delaunay_Mesh(points);
        stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "finished initial Mesh in " << duration.count()/1e6 << " seconds" << endl << endl;
        msh.spl = spl;
        msh.param = param;
        
        // mesh refinement
        start = chrono::high_resolution_clock::now();
        msh.Delaunay_refine(H);
        stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "finished delaunay refinement in " << duration.count()/1e6 << " seconds" << endl;
        cout << "minimum angle: " << check_minangle(&msh) << endl << endl;

        // mesh smoothing
        start = chrono::high_resolution_clock::now();
        msh.mesh_smoothing_2d(msh.find_boundary_nodes(), 50, 0.0);;
        stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "finished mesh smoothing in " << duration.count()/1e6 << " seconds" << endl;
        cout << "minimum angle: " << check_minangle(&msh) << endl << endl;

        // FEM goes here
        cout << "Performing FEM with P1 elements" << endl;
        double kappa = 2.0;
        int nv = msh.coords.size();
        vector<double> u = uexpr(msh.coords); // true solution
        vector<double> frhs = ulap(msh.coords, kappa); // rhs expression
        vector<int> bndnodes = msh.boundary_nodes(); // finding boundary nodes
        int nbnd = bndnodes.size();
        vector<double> dvals(nbnd); // dirichlet values
        for (int n = 0; n<nbnd; n++) {dvals[n] = u[bndnodes[n]];}
        vector<double> u_h = Poisson_2d(&msh, frhs, kappa, bndnodes, dvals);
        vector<double> error(nv);
        double maxE = 0.0;
        for (int n = 0; n<nv; n++) {error[n] = abs(u[n]-u_h[n]); if(error[n]>maxE){maxE=error[n];}}
        double P1_error = norm(error)/norm(u);
        cout << "l2 relative error from Linear mesh: " << P1_error << endl;
        cout << "l_inf absolute error from Linear mesh: " << maxE << endl << endl;

        cout << "Performing FEM with P2 elements" << endl;
        msh.make_quadratic();
        nv = msh.coords.size();
        u = uexpr(msh.coords); // true solution
        frhs = ulap(msh.coords, kappa); // rhs expression
        bndnodes = msh.boundary_nodes(); // finding boundary nodes
        nbnd = bndnodes.size();
        dvals.resize(nbnd);
        for (int n = 0; n<nbnd; n++) {dvals[n] = u[bndnodes[n]];}
        u_h = Poisson_2d(&msh, frhs, kappa, bndnodes, dvals);
        error.resize(nv);
        maxE = 0.0;
        for (int n = 0; n<nv; n++) {error[n] = abs(u[n]-u_h[n]); if(error[n]>maxE){maxE=error[n];}}
        double P2_error = norm(error)/norm(u);
        cout << "l2 relative error from Quadratic mesh: " << P2_error << endl;
        cout << "l_inf absolute error from Quadratic mesh: " << maxE << endl << endl;

        // writing
        WrtieVtk_tri(msh, error, "test3.vtk");
        cout << "finished writing to file" << endl;
}

void test4(){
    // loading bunny mesh
        string filename = "data/bunny.obj";
        Mesh Surface = ReadObj_tri(filename);
        Surface.compute_AHF();
        Surface.compute_normals();

        auto start = chrono::high_resolution_clock::now();
        vector<vector<double>> Params =  Parametric_Mapping(&Surface);
        auto stop = chrono::high_resolution_clock::now();
        //Surface.coords = Params;

        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "finished parametric mapping of surface in " << duration.count()/1e6 << " seconds" << endl << endl;

        start = chrono::high_resolution_clock::now();
        Surface.MIPS_minimize(Params,10);
        stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << "finished MIPS minimization in " << duration.count()/1e6 << " seconds" << endl << endl;

        Surface.coords = Params;
        WrtieVtk_tri(Surface, "test4.vtk");
        cout << "finished writing to file" << endl;
}

void test5(){
    cout << "performing test 5: edge coloring a mesh" << endl << endl;
    Mesh msh = Ellipse_mesh(50);

    // edge coloring
    msh.Graph_init();
    queue<int> Q = msh.Graph_Color_greedy();
    check_graph(msh);
    WrtieVtk_Graph(msh, "test5.vtk");

}

void test6(){
    cout << "performing test 6: using constrained delaunay to mesh lake superior" << endl << endl;
    vector<vector<double>> points;
    vector<vector<int>> segments;
    load_lake(points, segments);
    double h = 0.0;
    for (int i = 0; i<segments.size(); i++){
        h = h + norm(points[segments[i][1]]-points[segments[i][0]]);
    }
    h = (sqrt(3)/3)*(h/(double(segments.size())-1));

    auto start = chrono::high_resolution_clock::now();
    Mesh msh = GeoMesh_Delaunay_Mesh(segments,points,true); 
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "created initial constrained mesh in " << duration.count()/1e6 << " seconds" << endl << endl;

    // creating mesh refinement function
    int nspace = 20;
    for (int i = 0; i<segments.size(); i++){
        for (int j = 0; j<nspace; j++){
            points.push_back((points[segments[i][1]]-points[segments[i][0]])*(double(j+1)/double(nspace+2)) + points[segments[i][0]]);
        }
    }
    function<double(vector<double>)> H = create_grad(points, h, 0.2, 0.5);

    // mesh refinement
    start = chrono::high_resolution_clock::now();
    msh.Delaunay_refine(H);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "finished delaunay refinement in " << duration.count()/1e6 << " seconds" << endl;
    cout << "minimum angle: " << check_minangle(&msh) << endl << endl;

    // mesh smoothing
    start = chrono::high_resolution_clock::now();
    msh.mesh_smoothing_2d(msh.find_boundary_nodes(), 50, 0.0);;
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "finished mesh smoothing in " << duration.count()/1e6 << " seconds" << endl;
    cout << "minimum angle: " << check_minangle(&msh) << endl << endl;

    // FEM goes here
    cout << "Performing FEM with P1 elements" << endl;
    double kappa = 2.0;
    int nv = msh.coords.size();
    vector<double> u = uexpr(msh.coords); // true solution
    vector<double> frhs = ulap(msh.coords, kappa); // rhs expression
    vector<int> bndnodes = msh.boundary_nodes(); // finding boundary nodes
    int nbnd = bndnodes.size();
    vector<double> dvals(nbnd); // dirichlet values
    for (int n = 0; n<nbnd; n++) {dvals[n] = u[bndnodes[n]];}
    vector<double> u_h = Poisson_2d(&msh, frhs, kappa, bndnodes, dvals);
    vector<double> error(nv);
    double maxE = 0.0;
    for (int n = 0; n<nv; n++) {error[n] = abs(u[n]-u_h[n]); if(error[n]>maxE){maxE=error[n];}}
    double P1_error = norm(error)/norm(u);
    cout << "l2 relative error from Linear mesh: " << P1_error << endl;
    cout << "l_inf absolute error from Linear mesh: " << maxE << endl << endl;

    cout << "Performing FEM with P2 elements" << endl;
    msh.make_quadratic();
    nv = msh.coords.size();
    u = uexpr(msh.coords); // true solution
    frhs = ulap(msh.coords, kappa); // rhs expression
    bndnodes = msh.boundary_nodes(); // finding boundary nodes
    nbnd = bndnodes.size();
    dvals.resize(nbnd);
    for (int n = 0; n<nbnd; n++) {dvals[n] = u[bndnodes[n]];}
    u_h = Poisson_2d(&msh, frhs, kappa, bndnodes, dvals);
    error.resize(nv);
    maxE = 0.0;
    for (int n = 0; n<nv; n++) {error[n] = abs(u[n]-u_h[n]); if(error[n]>maxE){maxE=error[n];}}
    double P2_error = norm(error)/norm(u);
    cout << "l2 relative error from Quadratic mesh: " << P2_error << endl;
    cout << "l_inf absolute error from Quadratic mesh: " << maxE << endl << endl;
    
    WrtieVtk_tri(msh, error, "test6.vtk");
}


Mesh Ellipse_mesh(int npoints){
    int npoints_spline = 250;

    // generating points on ellipse for the spline
    vector<vector<double>> Spline_points = Ellipse(npoints_spline);

    // creating spline
    auto start = chrono::high_resolution_clock::now();
    Spline spl = spline_init(Spline_points,3);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "created Spline in " << duration.count()/1e6 << " seconds" << endl << endl;

    // interpolating points using spline
    vector<vector<double>> points = Zeros<double>(npoints,2);
    vector<double> param(npoints);
    for (int i = 0; i<npoints; i++){
        param[i] = double(i)/double(npoints);
        points[i] = spline_var(&spl, param[i]);
    } 

    // calculating average edge length
    double h = 0.0;
    for (int i = 0; i<npoints; i++){
        h = h + norm(points[(i+1)%npoints]-points[i]);
    }
    h = 1*(sqrt(3)/3)*(h/(double(npoints)-1));

    // initial Mesh
    start = chrono::high_resolution_clock::now();
    Mesh msh = GeoMesh_Delaunay_Mesh(points);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "finished initial Mesh in " << duration.count()/1e6 << " seconds" << endl << endl;;
    msh.spl = spl;
    msh.param = param;
        
    // mesh refinement
    start = chrono::high_resolution_clock::now();
    msh.Delaunay_refine(h);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "finished delaunay refinement in " << duration.count()/1e6 << " seconds" << endl;
    cout << "minimum angle: " << check_minangle(&msh) << endl << endl;

    // mesh smoothing
    start = chrono::high_resolution_clock::now();
    msh.mesh_smoothing_2d(msh.find_boundary_nodes(), 50, 0.0);;
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "finished mesh smoothing in " << duration.count()/1e6 << " seconds" << endl;
    cout << "minimum angle: " << check_minangle(&msh) << endl << endl;

    return msh;
}