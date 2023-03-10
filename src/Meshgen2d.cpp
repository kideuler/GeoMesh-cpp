#include "GeoMesh.hpp"
#include <unistd.h>
using namespace std;

/// subfunctions for delaunay triangulation
double eval_alpha(const vector<vector<double>> xs,double r_ref);
static vector<double> circumcenter(const vector<vector<double>> xs);
static vector<double> off_circumcenter(const vector<vector<double>> xs, double beta);
static bool inside_circumtri(const vector<vector<double>> xs, const vector<double> ps);
void Recursive_tri_delete(Mesh* DT, int hfid);
static bool Line_cross(const vector<double> &p1, const vector<double> &p2, const vector<double> &p3, const vector<double> &p4);
static bool Ray_in_triangle(Mesh* DT, int eid, int nid, int vid);
static double min_angle(const vector<vector<double>> &xs);
static vector<double> find_center(const vector<vector<double>> &xs);

/// Structure Member functions

/**
 * @brief Computes AHF array for triangular mesh
 * 
 */
void Mesh::compute_AHF(){
    bool oriented = true;
    bool manifold = true;
    int nelems = elems.size();
    int nv = coords.size();
    int nf = elems[0].size();
    sibhfs = Zeros<int>(nelems, nf);
    int* is_index = new int[nv+1];
    int v,c,vn;

    for(int ii=0; ii<nv+1; ii++){ is_index[ii] = 0; }
    for(int i = 0; i<nelems; i++){
        for(int j = 0; j<nf; j++){
            v = elems[i][j]+1;
            is_index[v]++;
        }
    }
    is_index[0] = 0;

    for(int ii=0; ii<nv; ii++){ is_index[ii+1] = is_index[ii]+is_index[ii+1]; }
    int ne = nelems*nf;

    int* v2nv = new int[ne];
    int* v2he_fid = new int[ne];
    int* v2he_leid = new int[ne];

    for(int i = 0; i<nelems; i++){
        for(int j = 0; j<nf; j++){
            v = elems[i][j];
            c = is_index[v];
            v2nv[c] = elems[i][(j+1)%3];
            v2he_fid[c] = i;
            v2he_leid[c] = j;
            is_index[v]++;
        }
    }

    for(int ii=nv-1; ii>=0; ii--){ is_index[ii+1] = is_index[ii]; }
    is_index[0] = 0;

    int first_heid_fid, first_heid_leid, prev_heid_fid, prev_heid_leid, nhes;
    for(int i = 0; i<nelems; i++){
        for(int j = 0; j<nf; j++){
            if(sibhfs[i][j] > 0){ continue; }
            v = elems[i][j];
            vn = elems[i][(j+1)%3];
            first_heid_fid = i;
            first_heid_leid = j;
            prev_heid_fid = i;
            prev_heid_leid = j;
            nhes = 0;
            
            // locate index in v2nv
            for(int index = is_index[vn]; index<is_index[vn+1]; index++){
                if(v2nv[index] == v){
                    sibhfs[prev_heid_fid][prev_heid_leid] = elids2hfid(v2he_fid[index]+1,v2he_leid[index]+1);
                    prev_heid_fid = v2he_fid[index];
                    prev_heid_leid = v2he_leid[index];
                    nhes++;
                }
            }
            // check for halfedges in the same orientation
            for(int index = is_index[v]; index<is_index[v+1]; index++){
                if(v2nv[index] == vn && v2he_fid[index] != i){
                    sibhfs[prev_heid_fid][prev_heid_leid] = elids2hfid(v2he_fid[index]+1,v2he_leid[index]+1);
                    prev_heid_fid = v2he_fid[index];
                    prev_heid_leid = v2he_leid[index];
                    nhes++;
                    oriented = true;
                }
            }
            if(prev_heid_fid != first_heid_fid){
                sibhfs[prev_heid_fid][prev_heid_leid] = elids2hfid(first_heid_fid+1,first_heid_leid+1);
                nhes++;
            }

            if(manifold && nhes>2){ manifold=false; oriented=false; }
        }
    }

    if (!manifold){
        cout << "Mesh is not a manifold" << endl;
    }
    if (!oriented){
        cout << "Mesh is not oriented" << endl;
    }
    delete is_index, v2nv, v2he_fid, v2he_leid;
}

/**
 * @brief raises linear mesh to quadratic mesh
 * 
 */
void Mesh::make_quadratic(){
    assert(degree == 1);
    vector<vector<int>> facets = obtain_trifacets(2);
    int nv = coords.size();
    if (spl.nv > 0){
        assert(param.size() == nv);
    }
    coords.resize(10*nv);
    param.resize(10*nv);

    // making space in elems
    int ii,jj;
    for (ii = 0; ii<nelems; ii++){
        elems[ii].resize(6);
    }
    
    int eid, lid;
    double a,b;
    for (ii = 0; ii<nelems; ii++){
        for (jj = 0; jj<3; jj++){
            eid = hfid2eid(sibhfs[ii][jj])-1;
            lid = hfid2lid(sibhfs[ii][jj])-1;
            if (eid < 0){ // boundary facet
                elems[ii][facets[jj][2]] = nv;
                a = param[elems[ii][facets[jj][0]]];
                b = param[elems[ii][facets[jj][1]]];
                if (spl.nv > 0 && abs(a-b)>1e-6){
                    param[nv] = (a+b)/2;
                    coords[nv] = spline_point_segment(&spl, a, b, 0.5);
                } else {
                    coords[nv] = (coords[elems[ii][facets[jj][0]]] + coords[elems[ii][facets[jj][1]]])/2.0;
                }
                nv++;
            } else { // interior facet
                if (elems[ii][facets[jj][2]] == 0){ // if not used
                    elems[ii][facets[jj][2]] = nv;
                    elems[eid][facets[lid][2]] = nv;
                    coords[nv] = (coords[elems[ii][facets[jj][0]]] + coords[elems[ii][facets[jj][1]]])/2.0;
                    nv++;
                }
            }
        }
    }
    degree = 2;
    coords.resize(nv);
    param.resize(nv);
}

/**
 * @brief decompose quadratic mesh to linear
 * 
 */
void Mesh::decompose_to_linear(){
    if (degree == 1){
        return;
    }
    vector<vector<int>> elems_new = Zeros<int>(degree*degree*nelems,3);
    vector<vector<int>> ho2lin = {{0,3,5},{3,1,4},{5,4,2},{3,4,5}};

    int kk = 0;
    for (int ii = 0; ii<nelems; ii++){
        for (int jj = 0; jj<degree*degree; jj++){
            elems_new[kk][0] = elems[ii][ho2lin[jj][0]];
            elems_new[kk][1] = elems[ii][ho2lin[jj][1]];
            elems_new[kk][2] = elems[ii][ho2lin[jj][2]];
            kk++;
        }
    }

    nelems = elems_new.size();
    elems = elems_new;
    degree = 1;
    compute_AHF();
}

/**
 * @brief compute the one ring stencil for the mesh
 * 
 * @param maxne max number of elements that can be around each point
 */
void Mesh::compute_Onering(int maxne, bool compress){
    int nv = coords.size();
    int* hvids = new int[maxne*nv];
    vector<int> nv2e;
    nv2e.assign(nv,0);
    int size = 0;

    int v;
    for (int i = 0; i<nelems; i++){
        for (int j = 0; j<3; j++){
            v = elems[i][j];
            nv2e[v]++;
            if (nv2e[v] >= maxne){
                delete hvids;
                cout << "Mesh::compute_Onering buffers too small, enlarging ring" << endl;
                compute_Onering(2*maxne, compress);
                return;
            }
            hvids[maxne*v+nv2e[v]-1] = elids2hfid(i+1,j+1);
            size++;
        }
    }

    // putting in structure
    if (compress){
        stl.index.assign(nv+1,0);
        stl.hvids.assign(size,0);
        int sz = 0;
        for (int vid = 0; vid<nv; vid++){
            for (int i = 0; i<nv2e[vid]; i++){
                stl.hvids[sz] = hvids[maxne*vid+i];
                sz++;
            }
            stl.index[vid+1] = sz;
        }
    } else {
        stl.index.assign(nv+1,0);
        stl.hvids.assign(maxne*nv,-1);
        int sz = 0;
        for (int vid = 0; vid<nv; vid++){
            for (int i = 0; i<nv2e[vid]; i++){
                stl.hvids[maxne*vid+i] = hvids[maxne*vid+i];
                sz++;
            }
            stl.index[vid+1] = stl.index[vid]+maxne;
        }
    }
    delete hvids;

    return;
}

/**
 * @brief delete triangles in dele
 * 
 */
void Mesh::delete_tris(){
    int i,j;
    int nelems2 = 0;
    int sz = sibhfs[0].size();
    vector<int> idx(nelems);
    vector<int> idx_rev(nelems);
    fill(idx_rev.begin(), idx_rev.end(),-1);

    // delete triangles to be deleted
    for (i = 0; i<nelems; i++){
        if (!delete_elem[i]){
            elems[nelems2] = elems[i];
            sibhfs[nelems2] = sibhfs[i];
            delete_elem[nelems2] = false;
            idx[nelems2] = i;
            nelems2++;
        } else {
            delete_elem[i] = false;
        }
    }
    nelems = nelems2;

    for (i = 0; i<nelems; i++){
        idx_rev[idx[i]] = i;
    }

    int nside;
    int hfid, eid, lid;
    for (i = 0; i<nelems2; i++){
        nside = 0;
        for (j = 0; j < sz; j++){
            hfid = sibhfs[i][j];
            if (!hfid == 0){
                eid = hfid2eid(hfid);
                lid = hfid2lid(hfid);
                if (!delete_elem[eid-1]){
                    if (idx_rev[eid-1] == -1){
                        sibhfs[i][j] = 0;
                    } else {
                        sibhfs[i][j] = elids2hfid(idx_rev[eid-1]+1,lid);
                    }
                } else {
                    sibhfs[i][j] = 0;
                }
                nside++;
            } else {
                sibhfs[i][j] = 0;
            }
        }
        if (nside == sz){
            on_boundary[i] = false;
        } else {
            on_boundary[i] = true;
        }
        delete_elem[i] = false;
    }
    elems.resize(nelems);
    sibhfs.resize(nelems);
}

void Mesh::Delaunay_refine(double r_ref){
    auto f = [r_ref](vector<double> xs){ return r_ref; };
    Delaunay_refine(f);
}

int find_hfid(vector<vector<int>> sibhfs, int eid){
    int hfid = 0;
    for (int i = 0; i<3; i++){
        if(sibhfs[eid][i] == 0){
            hfid = elids2hfid(eid+1,i+1);
            return hfid;
        }
    }
    cout << "not found" << endl;
    return hfid;
}
void Mesh::Delaunay_refine(function<double(vector<double>)> r_ref){
    // random number generator
    default_random_engine re;
    uniform_real_distribution<double> unif(-1, 1);

    int n,i,nelems2;
    bool exitl;
    double alpha,theta;
    int nv = coords.size();
    vector<vector<double>> ps = Zeros<double>(3,2);
    int tri,eid,lid,ub;

    // finding total area and estimate to define maximum size bounds
    double total_area = 0.0;
    vector<double> mid;
    double minr = r_ref(mid);
    for (n=0; n<nelems; n++){
        ps[0] = coords[elems[n][0]];
        ps[1] = coords[elems[n][1]];
        ps[2] = coords[elems[n][2]];
        total_area += area_tri(ps);
    }
    double area_single = max(minr*minr/2,1e-6);
    ub = elems.size(); 
    if (1.2*total_area/area_single > ub){
        ub = (int) 1.2*total_area/area_single;
    }
    ub = min(ub,1000000);
    coords.resize(ub);
    param.resize(ub);
    elems.resize(ub);
    sibhfs.resize(ub);
    delete_elem.resize(ub);
    on_boundary.resize(ub);
    vector<double> C;
    nelems2 = nelems;
    vector<int> order(ub);
    for (i = 0; i<ub; i++){
        order[i] = i;
    }

    // main loop 
    n=0;
    bool inside_domain;
    int e;
    bool freed = false;
    while (n<nelems){
        e = order[n];
        if (!delete_elem[e]){
            ps[0] = coords[elems[e][0]];
            ps[1] = coords[elems[e][1]];
            ps[2] = coords[elems[e][2]];
            C = off_circumcenter(ps,sqrt(2));
            C[0] += 1e-4*unif(re);
            C[1] += 1e-4*unif(re);
            // sanity check
            if (abs(C[0] > 1e6 || abs(C[1]) > 1e6)){
                cout << "abnormally large point added from points " << C[0] << "," << C[1] << endl;
                printMat(ps); 
            }
            alpha = eval_alpha(ps,r_ref(C));
            theta = min_angle(ps);
            if (alpha > 1+0.1*double(3*n/(ub))){

                // add circumcircle to Mesh
                coords[nv] = C;
                tri = e;
                inside_domain = find_enclosing_tri(&tri, C);
                //cout << alpha << endl;
                if (!inside_domain){
                    if (tri == -1){
                        cout << "find triangle location failed: deleting point" << endl;
                        nv--;
                    } else {
                        if (inside_diametral(tri, C)){
                            Flip_Insertion_segment(nv, tri);
                        } else {
                            if (on_boundary[e]){
                                Flip_Insertion_segment(nv, tri);
                            } else{
                            if (n == nelems-1 || e == nelems-1){
                                Flip_Insertion_segment(nv, tri);
                            } else {
                                nv--;
                                order[n] = nelems-1;
                                order[nelems-1] = e;
                                n--;
                            }
                            }
                        }
                    }
                } else {
                    bool stop = false;
                    if (on_boundary[tri]){
                        int hfid = find_hfid(sibhfs,tri);
                        if (inside_diametral(hfid, C)){
                            Flip_Insertion_segment(nv, hfid);
                            stop = true;
                        }
                    } else if(sibhfs[tri][0] > 0) {
                        if(on_boundary[hfid2eid(sibhfs[tri][0])-1]){
                        int hfid = find_hfid(sibhfs,hfid2eid(sibhfs[tri][0])-1);
                        if (inside_diametral(hfid, C)){
                            Flip_Insertion_segment(nv, hfid);
                            stop = true;
                        }
                        }
                    } else if(sibhfs[tri][1] > 0) {
                        if(on_boundary[hfid2eid(sibhfs[tri][1])-1]){
                        int hfid = find_hfid(sibhfs,hfid2eid(sibhfs[tri][1])-1);
                        if (inside_diametral(hfid, C)){
                            Flip_Insertion_segment(nv, hfid);
                            stop = true;
                        }
                        }
                    } else if(sibhfs[tri][2] > 0) {
                        if(on_boundary[hfid2eid(sibhfs[tri][2])-1]){
                        int hfid = find_hfid(sibhfs,hfid2eid(sibhfs[tri][2])-1);
                        if (inside_diametral(hfid, C)){
                            Flip_Insertion_segment(nv, hfid);
                            stop = true;
                        }
                        }
                    }
                    


                    if (!stop){
                        Flip_Insertion(&nv,tri);
                    }
                }
                nv++;

                if ((double) nelems >= (0.95)*((double) elems.size())){
                    cout << "approaching size bound, freeing up space" << endl;
                    ub = ub*2;
                    coords.resize(ub);
                    param.resize(ub);
                    elems.resize(ub);
                    sibhfs.resize(ub);
                    delete_elem.resize(ub);
                    on_boundary.resize(ub);
                    order.resize(ub);
                }
            }
        }
        n++;
    }
    coords.resize(nv);
    param.resize(nv);
    delete_tris();
    cout << "created " << nelems-nelems2 << " triangles from refining the mesh" << endl;
}


/**
 * @brief Converts half-facet id to element id
 * 
 * @param hfid Half facet id
 * @return element id
 */
int hfid2eid(int hfid){
    return (hfid >> 8);
}
/**
 * @brief Convert half-facet id to local edge id
 * 
 * @param hfid Half facet id
 * @return local edge id
 */
int hfid2lid(int hfid){
    return (hfid&255) + 1;
}
/**
 * @brief Converts element id and local edge id to half-facet id
 * 
 * @param eid element id
 * @param lid local edge id
 * @return half-facet id
 */
int elids2hfid(int eid, int lid){
    return ((eid << 8) + lid - 1);
}

// need stack and push and pop functions
void push_stack(stack** head, int hfid){
    stack* new_stack = new stack;
    new_stack->next = *head;
    new_stack->hfid = hfid;
    (*head) = new_stack;
    return;
}
void pop_stack(stack** head){
    if (*head == NULL){
        return;
    }
    stack* temp = *head;
    *head = (*head)->next;
    delete temp;
    return;
}

/**
 * @brief returns a vector of booleans where boundary nodes are true
 * 
 * @return vector<bool> 
 */
vector<bool> Mesh::find_boundary_nodes(){
    int nv = coords.size();
    vector<vector<int>> facets = obtain_trifacets(degree);
    int nnodesf = facets[0].size();
    vector<bool> bnd;
    bnd.assign(nv,false);
    for (int i = 0; i<nelems; i++){
        for (int j = 0; j<3; j++){
            if (sibhfs[i][j] == 0){
                for (int k = 0; k<nnodesf; k++){
                    bnd[elems[i][facets[j][k]]] = true;
                }
            }
        }
    }

    return bnd;
}

/**
 * @brief returns a vector of int of boundary nodes
 * 
 * @return vector<int> 
 */
vector<int> Mesh::boundary_nodes(){
    int nv = coords.size();
    vector<vector<int>> facets = obtain_trifacets(degree);
    int nnodesf = facets[0].size();
    vector<bool> bnd;
    bnd.assign(nv,false);
    vector<int> bdy;
    bdy.reserve(nv);
    int count=0;
    int v;
    for (int i = 0; i<nelems; i++){
        for (int j = 0; j<3; j++){
            if (sibhfs[i][j] == 0){
                for (int k = 0; k<nnodesf; k++){
                    v = elems[i][facets[j][k]];
                    if (!bnd[v]){
                        bnd[v] = true;
                        bdy.push_back(v);
                        count++;
                    }
                }
            }
        }
    }
    bdy.resize(count);
    return bdy;
}

/**
 * @brief Create constrained Delaunay Mesh from PSLG
 * 
 * @param segs Constrained segments (nsegs -by- 2)
 * @param xs Point data (nv -by- 2)
 * @return Constrained delaunay Mesh
 */
Mesh GeoMesh_Delaunay_Mesh(const vector<vector<int>> &segs, vector<vector<double>> &xs, bool delete_exterior){

    // segs define boundary segments for the mesh
    Mesh DT = GeoMesh_Delaunay_Mesh(xs);
    int maxne = 10;
    DT.compute_Onering(maxne,false);
    DT.facets.resize(DT.nelems);
    DT.bwork.assign(DT.nelems, false);

    int v1, v2, hvid, kk, eid, nid, nf, hfid, lid, oppeid, opplid;
    nf = 0;
    bool inray, connected, convex;
    for (int S = 0; S<segs.size(); S++){
        v1 = segs[S][0];
        v2 = segs[S][1];
        
        
        connected = false;
        while (!connected){ // segment not done
            inray = false;
            convex = false;
            hvid = 1;
            kk = 0;
            while (!inray && !connected && kk < maxne){ // finding segment connection or line that obstructs segment
                hvid = DT.stl.hvids[DT.stl.index[v1]+kk];
                if (hvid <= 0){
                    break;
                }

                eid = hfid2eid(hvid)-1;
                nid = hfid2lid(hvid)-1;
                if (DT.elems[eid][(nid+1)%3] == v2){
                    connected = true;
                    hfid = DT.sibhfs[eid][nid];
                    if (hfid > 0){
                        DT.facets[nf][0] = eid;
                        DT.facets[nf][1] = nid;
                        DT.bwork[eid] = true;
                        nf++;
                    }
                    break;
                }
                if (DT.elems[eid][(nid+2)%3] == v2){
                    connected = true;
                    hfid = DT.sibhfs[eid][(nid+2)%3];
                    if (hfid > 0){
                        DT.facets[nf][0] = hfid2eid(hfid)-1;
                        DT.facets[nf][1] = hfid2lid(hfid)-1;
                        DT.bwork[hfid2eid(hfid)-1] = true;
                        nf++;
                    }
                    break;
                }
                
                if (Line_cross(DT.coords[DT.elems[eid][(nid+1)%3]], DT.coords[DT.elems[eid][(nid+2)%3]], DT.coords[v1], DT.coords[v2])){
                    inray = true;
                    break;
                }
                kk++;
            }

            if (inray) { // flip segment and correct onering if convex quad
                lid = (nid+1)%3;
                hfid = DT.sibhfs[eid][lid];
                while (!DT.convex_quad(eid,lid)){ // keep finding segments till one is convex
                    oppeid = hfid2eid(hfid)-1;
                    opplid = hfid2lid(hfid)-1;
                    if (Line_cross(DT.coords[DT.elems[oppeid][(opplid+1)%3]], DT.coords[DT.elems[oppeid][(opplid+2)%3]], DT.coords[v1], DT.coords[v2])){
                        eid = oppeid;
                        lid = (opplid+1)%3;
                    } else if (Line_cross(DT.coords[DT.elems[oppeid][(opplid+2)%3]], DT.coords[DT.elems[oppeid][opplid]], DT.coords[v1], DT.coords[v2])){
                        eid = oppeid;
                        lid = (opplid+2)%3;
                    }
                }
                DT.flip_edge(eid, lid);
                // get working inefficient (highly inefficient for big mesh, need code that simply replaces nodes that are changed not all nodes)
                DT.compute_Onering(maxne,false);
            }

        }


    }


    // delete elements on opposite side of boundary;
    if (delete_exterior){
        for (int i=0;i<nf;i++){
            eid = DT.facets[i][0];
            lid = DT.facets[i][1];
            hfid = DT.sibhfs[eid][lid];
            Recursive_tri_delete(&DT, hfid);
        }
    }
    DT.delete_tris();
    DT.stl.index.resize(0);
    return DT;
}
void Recursive_tri_delete(Mesh* DT, int hfid){
    int eid = hfid2eid(hfid)-1;
    int lid = hfid2lid(hfid)-1;
    if (!DT->delete_elem[eid]){
        DT->delete_elem[eid] = true;
        int hfid1 = DT->sibhfs[eid][(lid+1)%3];
        if (hfid1 != 0){
            if (!DT->bwork[hfid2eid(hfid1)-1]){
                Recursive_tri_delete(DT, hfid1);
            }
        }
        int hfid2 = DT->sibhfs[eid][(lid+2)%3];
        if (hfid2 != 0){
            if (!DT->bwork[hfid2eid(hfid2)-1]){
                Recursive_tri_delete(DT, hfid2);
            }
        }
    } 
    return;
}
static bool Line_cross(const vector<double> &p1, const vector<double> &p2, const vector<double> &p3, const vector<double> &p4){
    bool a1 = (p4[1]-p1[1])*(p3[0]-p1[0]) > (p3[1]-p1[1])*(p4[0]-p1[0]);
    bool a2 = (p4[1]-p2[1])*(p3[0]-p2[0]) > (p3[1]-p2[1])*(p4[0]-p2[0]);
    bool b1 = (p3[1]-p1[1])*(p2[0]-p1[0]) > (p2[1]-p1[1])*(p3[0]-p1[0]);
    bool b2 = (p4[1]-p1[1])*(p2[0]-p1[0]) > (p2[1]-p1[1])*(p4[0]-p1[0]);

    return (a1 != a2 && b1 != b2);
}
static bool Ray_in_triangle(Mesh* DT, int eid, int nid, int vid){
    vector<double> v1 = DT->coords[DT->elems[eid][(nid+1)%3]] - DT->coords[DT->elems[eid][nid]];
    double n1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
    v1 = v1/n1;
    vector<double> v2 = DT->coords[DT->elems[eid][(nid+2)%3]] - DT->coords[DT->elems[eid][nid]];
    double n2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1]);
    v2 = v2/n2;
    vector<double> v3 = DT->coords[vid] - DT->coords[DT->elems[eid][nid]];
    double n3 = sqrt(v3[0]*v3[0] + v3[1]*v3[1]);
    v3 = v3/n3;

    double cross1 = v1[0]*v3[1] - v1[1]*v3[0];
    double cross2 = v3[0]*v2[1] - v3[1]*v2[0];

    return (cross1 > 0 && cross2 > 0);
}

/**
 * @brief Create constrained Delaunay Mesh and applies spline parameters to the mesh
 * 
 * @param xs Point data (nv -by- 2)
 * @param params parameters of the spline
 * @return Constrained delaunay Mesh
 */
Mesh GeoMesh_Delaunay_Mesh(vector<vector<double>> &xs, vector<double> &params){
    Mesh DT = GeoMesh_Delaunay_Mesh(xs);
    int nv = DT.coords.size();
    assert(nv == params.size());
    DT.param = params;
    return DT;
}

/**
 * @brief create Delaunay Mesh from point set
 * 
 * @param xs Point data (nv -by- 2)
 * @return Delaunay Mesh 
 */
Mesh GeoMesh_Delaunay_Mesh(vector<vector<double>> &xs){
    // size checking
    default_random_engine re;
    uniform_real_distribution<double> unif(-1, 1);
    int nv = xs.size();
    int n,i,j;

    Mesh DT;
    int ub = 2*nv*nv;
    DT.coords = Zeros<double>(nv+3,2);
    DT.elems = Zeros<int>(ub,3);
    DT.sibhfs = Zeros<int>(ub,3);
    DT.facets = Zeros<int>(ub,2);
    DT.delete_elem.resize(ub);
    DT.on_boundary.resize(ub);
    vector<double> a = min_array(xs);
    vector<double> b = max_array(xs);
    double dx = (b[0]-a[0])/10;
    double dy = (b[1]-a[1])/10;
    double d = max(dx,dy);

    // reorder points for optimal triangle search O(N*log(N))
    for (n=0; n<nv; n++){
        DT.coords[n][0] = (xs[n][0]-a[0])/d + 1e-4*unif(re);
        DT.coords[n][1] = (xs[n][1]-a[1])/d + 1e-4*unif(re);
    }
    
    // sort points by proximity
    int nbin = (int)ceil(pow(nv,0.5));
    int* bins = new int[nv];
    int* order = new int[nv];
    for (n = 0; n<nv; n++){
        int p = (int)(DT.coords[n][1]*nbin*0.999);
        int q = (int)(DT.coords[n][0]*nbin*0.999);
        if (p%2){
            bins[n] = (p+1)*nbin-q;
        } else {
            bins[n] = p*nbin+q+1;
        }
        order[n] = n;
    }

    int key;
    int temp;
    for (i = 1; i<nv; i++){
        key = bins[i];
        temp = order[i];
        j = i-1;
        while(j>=0 && bins[j]>key){
            bins[j+1] = bins[j];
            order[j+1] = order[j];
            j--;
        }
        bins[j+1] = key;
        order[j+1] = temp;
    }

    delete bins;



    DT.coords[nv][0] = -100;
    DT.coords[nv][1] = -100;
    DT.coords[nv+1][0] = 100;
    DT.coords[nv+1][1] = -100;
    DT.coords[nv+2][0] = 0;
    DT.coords[nv+2][1] = 100;

    DT.elems[0] = {nv,nv+1,nv+2};
    DT.sibhfs[0] = {0,0,0};

    // loop inserting each point into Mesh
    DT.nelems = 1;
    int tri = -1;
    bool exitl,inside;
    int vid;
    for (int n = 0; n < nv; n++){
        vid = order[n];
        tri = DT.nelems-1;
        inside = DT.find_enclosing_tri(&tri, DT.coords[vid]);
        if (!inside){
            cout << "no enclosing tri found" << endl;
        }
        // inserting node into the Mesh using Bowyer-Watson algorithm
        DT.Flip_Insertion(&vid,tri);
        if ((double) DT.nelems >= (0.9)*((double) ub)){
            cout << "approaching size bound, freeing up space" << endl;
            DT.delete_tris();
        }
    }
    delete order;

    for (int n = 0; n<DT.nelems; n++){
        for (int i = 0; i<3; i++){
            if (DT.elems[n][i] >= nv){
                DT.delete_elem[n] = true;
            }
        }
    }

    DT.delete_tris();
    DT.coords.resize(nv);
    for (n=0; n<nv; n++){
        DT.coords[n][0] = xs[n][0];
        DT.coords[n][1] = xs[n][1];
    }
    cout << "created " << DT.nelems << " triangles from initial points" << endl;
    DT.degree = 1;
    return DT;
}

/**
 * @brief Insert node into mesh using Lawson flipping algorithm
 * 
 * @param vid Node to be inserted
 * @param tri_s Triangle that encloses the node
 */
void Mesh::Flip_Insertion(int* vid, int tri_s){
    delete_elem[tri_s] = true;
    int hfid,eid,lid;

    vector<int> tri = {elems[tri_s][0], elems[tri_s][1], elems[tri_s][2]};
    vector<int> sib = {sibhfs[tri_s][0],sibhfs[tri_s][1],sibhfs[tri_s][2]};
    
    // splitting triangle and adding subtriangles to the stack if they are not on the boundary
    stack* head = NULL;
    vector<int> eids = {nelems, nelems+1, nelems+2};
    for (int i = 0; i<3; i++){
        elems[eids[i]] = {*vid, tri[i], tri[(i+1)%3]};
        hfid = sib[i];
        sibhfs[eids[i]] = {elids2hfid(eids[(i+2)%3]+1 ,3), hfid, elids2hfid(eids[(i+1)%3]+1,1)};
        if (hfid2eid(hfid) > 0){
            sibhfs[hfid2eid(hfid)-1][hfid2lid(hfid)-1] = elids2hfid(eids[i]+1, 2);
            push_stack(&head, elids2hfid(eids[i]+1, 2));
            on_boundary[eids[i]] = false;
        } else {
            on_boundary[eids[i]] = true;
        }
    }

    nelems+=3;
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    int oppeid,opplid;

    // loop through stack and flip edges
    while (head != NULL){
        hfid = head->hfid;
        pop_stack(&head);
        eid = hfid2eid(hfid)-1;
        lid = hfid2lid(hfid)-1;
        oppeid = hfid2eid(sibhfs[eid][lid])-1;
        opplid = hfid2lid(sibhfs[eid][lid])-1;
        xs[0] = coords[elems[oppeid][0]];
        xs[1] = coords[elems[oppeid][1]];
        xs[2] = coords[elems[oppeid][2]];
        
        if (inside_circumtri(xs, coords[elems[eid][0]])){
            flip_edge(eid,1);

            hfid = sibhfs[oppeid][1];
            if (hfid2eid(hfid) > 0){
                push_stack(&head, elids2hfid(oppeid+1,2));
            }
            hfid = sibhfs[eid][1];
            if (hfid2eid(hfid) > 0){
                push_stack(&head, elids2hfid(eid+1,2));
            }
        }

    }
    return;
}
/**
 * @brief Insert node into mesh on boundary segment by splitting segment using Lawson flipping algorithm
 * 
 * @param vid Node to be inserted
 * @param hfid Triangle and local edge id of the segment to be split
 */
void Mesh::Flip_Insertion_segment(int vid, int hfid){

    // adding two triangles instead of three
    int nvS = spl.nv;
    stack* head = NULL;
    int eid = hfid2eid(hfid)-1;
    int lid = hfid2lid(hfid)-1;
    double a = param[elems[eid][lid]];
    double b = param[elems[eid][(lid+1)%3]];
    if (abs(b-a) < 1e-6 || nvS == 0){
        coords[vid] = (coords[elems[eid][lid]] + coords[elems[eid][(lid+1)%3]])*0.5;
    } else {
        if (abs(b-a) > 0.5){
            if (b<a){
                b +=1;
            } else {
                a +=1;
            }
        }
        param[vid] = (a+b)/2;
        coords[vid] = spline_point_segment(&spl, param[elems[eid][lid]],  param[elems[eid][(lid+1)%3]], 0.5);
    }
    double radius = 0.5*(norm(coords[elems[eid][lid]] - coords[elems[eid][(lid+1)%3]]));
    delete_elem[eid] = true;
    elems[nelems] = {vid, elems[eid][(lid+2)%3], elems[eid][lid]};
    sibhfs[nelems] = {elids2hfid(nelems+2,3), sibhfs[eid][(lid+2)%3], 0};
    elems[nelems+1] = {vid, elems[eid][(lid+1)%3] ,elems[eid][(lid+2)%3]};
    sibhfs[nelems+1] = {0, sibhfs[eid][(lid+1)%3], elids2hfid(nelems+1,1)};
    if ( sibhfs[eid][(lid+1)%3] != 0){
        sibhfs[hfid2eid(sibhfs[eid][(lid+1)%3])-1][hfid2lid(sibhfs[eid][(lid+1)%3])-1] = elids2hfid(nelems+2, 2);
        push_stack(&head, elids2hfid(nelems+2, 2));
    }
    if ( sibhfs[eid][(lid+2)%3] != 0){
        sibhfs[hfid2eid(sibhfs[eid][(lid+2)%3])-1][hfid2lid(sibhfs[eid][(lid+2)%3])-1] = elids2hfid(nelems+1, 2);
        push_stack(&head, elids2hfid(nelems+1, 2));
    }
    on_boundary[nelems] = true;
    on_boundary[nelems+1] = true;

    nelems+=2;
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    int oppeid,opplid;

    // loop through stack and flip edges
    while (head != NULL){
        hfid = head->hfid;
        pop_stack(&head);
        eid = hfid2eid(hfid)-1;
        lid = hfid2lid(hfid)-1;
        oppeid = hfid2eid(sibhfs[eid][lid])-1;
        opplid = hfid2lid(sibhfs[eid][lid])-1;
        xs[0] = coords[elems[oppeid][0]];
        xs[1] = coords[elems[oppeid][1]];
        xs[2] = coords[elems[oppeid][2]];
        
        if (inside_circumtri(xs, coords[elems[eid][0]])){
            flip_edge(eid,lid);

            hfid = sibhfs[oppeid][1];
            if (hfid2eid(hfid) > 0){
                push_stack(&head, elids2hfid(oppeid+1,2));
            }
            hfid = sibhfs[eid][1];
            if (hfid2eid(hfid) > 0){
                push_stack(&head, elids2hfid(eid+1,2));
            }
        }

    }
    return;
}

/**
 * @brief Flip the edge in a Delaunay mesh
 * 
 * @param eid element to be flipped
 * @param lid edge to be flipped across
 */
void Mesh::flip_edge(int eid, int lid){
    vector<int> tri1(3);
    vector<int> tri2(3);
    vector<int> sib1(3);
    vector<int> sib2(3);

    int hfid = sibhfs[eid][lid];
    int oppeid = hfid2eid(hfid)-1;
    int opplid = hfid2lid(hfid)-1;

    int v1 = elems[eid][(lid+2)%3];
    int v2 = elems[eid][lid];
    int v3 = elems[eid][(lid+1)%3];
    int v4 = elems[oppeid][(opplid+2)%3];

    tri1 = {v1,v2,v4};
    tri2 = {v1,v4,v3};
    sib1 = {sibhfs[eid][(lid+2)%3], sibhfs[oppeid][(opplid+1)%3], elids2hfid(oppeid+1,1)};
    sib2 = {elids2hfid(eid+1,3), sibhfs[oppeid][(opplid+2)%3], sibhfs[eid][(lid+1)%3]};

    elems[eid] = tri1;
    elems[oppeid] = tri2;
    sibhfs[eid] = sib1;
    sibhfs[oppeid] = sib2;
    
    bool sib11 = false;
    bool sib12 = false;
    bool sib21 = false;
    bool sib22 = false;
    if (hfid2eid(sib1[0]) > 0){
        sibhfs[hfid2eid(sib1[0])-1][hfid2lid(sib1[0])-1] = elids2hfid(eid+1,1);
    } else {
        sib11 = true;
    }
    if (hfid2eid(sib1[1]) > 0){
        sibhfs[hfid2eid(sib1[1])-1][hfid2lid(sib1[1])-1] = elids2hfid(eid+1,2);
    } else {
        sib12 = true;
    }
    if (hfid2eid(sib2[1]) > 0){
        sibhfs[hfid2eid(sib2[1])-1][hfid2lid(sib2[1])-1] = elids2hfid(oppeid+1,2);
    } else {
        sib21 = true;
    }
    if (hfid2eid(sib2[2]) > 0){
        sibhfs[hfid2eid(sib2[2])-1][hfid2lid(sib2[2])-1] = elids2hfid(oppeid+1,3);
    } else {
        sib22 = true;
    }
    on_boundary[eid] = sib11 || sib12;
    on_boundary[oppeid] = sib21 || sib22;

    return;
}

/**
 * @brief find if triangles attatcched to eid,lid form convex quad
 * 
 * @param eid element id
 * @param lid edge id
 * @return true 
 * @return false 
 */
bool Mesh::convex_quad(int eid, int lid){
    int hfid = sibhfs[eid][lid];
    int oppeid = hfid2eid(hfid)-1;
    int opplid = hfid2lid(hfid)-1;

    int v1 = elems[eid][(lid+2)%3];
    int v2 = elems[eid][lid];
    int v4 = elems[eid][(lid+1)%3];
    int v3 = elems[oppeid][(opplid+2)%3];

    double sides[4][2] = {{coords[v2][0]-coords[v1][0],coords[v2][1]-coords[v1][1]},
    {coords[v3][0]-coords[v2][0],coords[v3][1]-coords[v2][1]},
    {coords[v4][0]-coords[v3][0],coords[v4][1]-coords[v3][1]},
    {coords[v1][0]-coords[v4][0],coords[v1][1]-coords[v4][1]}};

    double cp;
    bool sign = false;
    for (int i = 0; i<4; i++){
        cp = sides[i][0]*sides[(i+1)%4][1] - sides[i][1]*sides[(i+1)%4][0];
        if (i==0){
            sign = (cp>0);
        } else if (sign != (cp>0)){
            return false;
        }
    }

    return true;
}

/**
 * @brief Finding the triangle that encloses a node
 * 
 * @param tri Starting triangle
 * @param ps vector<double> containing location of the query point
 * @return true 
 * @return false 
 */
bool Mesh::find_enclosing_tri(int* tri, vector<double> &ps){
    int v1,v2,v3,i,hfid;
    bool stop;
    int iters = 0;
    i = 0;
    stop = false;
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    while (!stop){
        v1 = elems[*tri][0];
        v2 = elems[*tri][1];
        v3 = elems[*tri][2];
        if (delete_elem[*tri]){
            cout << "deleted tri passed ran through in find_enclosing_tri, should not be" << endl;

        }
        xs[0] = coords[v1];
        xs[1] = coords[v2];
        xs[2] = coords[v3];
        if (inside_tri(xs,ps)){
            stop = true;
            return true;
        } else {
            double AB[2] = {coords[v2][0]-coords[v1][0], coords[v2][1]-coords[v1][1]};
            double BC[2] = {coords[v3][0]-coords[v2][0], coords[v3][1]-coords[v2][1]};
            double CA[2] = {coords[v1][0]-coords[v3][0], coords[v1][1]-coords[v3][1]};
            double AP[2] = {ps[0]-coords[v1][0], ps[1]-coords[v1][1]};
            double BP[2] = {ps[0]-coords[v2][0], ps[1]-coords[v2][1]};
            double CP[2] = {ps[0]-coords[v3][0], ps[1]-coords[v3][1]};
            double N1[2] = {AB[1],-AB[0]};
            double N2[2] = {BC[1],-BC[0]};
            double N3[2] = {CA[1],-CA[0]};
            double S1 = AP[0]*N1[0]+AP[1]*N1[1];
            double S2 = BP[0]*N2[0]+BP[1]*N2[1];
            double S3 = CP[0]*N3[0]+CP[1]*N3[1];
            if ((S1>0)&&(S1>=S2)&&(S1>=S3)){
                hfid = sibhfs[*tri][0];
                if (hfid != 0){
                    *tri = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tri = elids2hfid(*tri+1,1);
                    return false;
                }
            } else if ((S2>0)&&(S2>=S1)&&(S2>=S3)) {
                hfid = sibhfs[*tri][1];
                if (hfid != 0){
                    *tri = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tri =  elids2hfid(*tri+1,2);
                    return false;
                }
            } else if ((S3>0)&&(S3>=S1)&&(S3>=S2)){
                hfid = sibhfs[*tri][2];
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


// sub-functions necessary for delaunay
/// evalute alpha for delaunay refinement
double eval_alpha(const vector<vector<double>> xs,double r_ref){
    double a = norm(xs[1]-xs[0]);
    double b = norm(xs[2]-xs[1]);
    double c = norm(xs[2]-xs[0]);
    double s = 0.5*(a+b+c);
    double A = sqrt(s*(s-a)*(s-b)*(s-c));
    double r = a*b*c/(4*A);
    return r/r_ref;
}
/// find circumcenter for a triangle
static vector<double> circumcenter(const vector<vector<double>> xs){
    double ax = xs[0][0];
    double ay = xs[0][1];
    double bx = xs[1][0];
    double by = xs[1][1];
    double cx = xs[2][0];
    double cy = xs[2][1];
    double D = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by));
    double ux = (ax*ax + ay*ay)*(by-cy) + \
        (bx*bx + by*by)*(cy-ay) + \
        (cx*cx + cy*cy)*(ay-by);
    double uy = (ax*ax + ay*ay)*(cx-bx) + \
        (bx*bx + by*by)*(ax-cx) + \
        (cx*cx + cy*cy)*(bx-ax);
    return {ux/D, uy/D};
}
/// find offcenter steiner point outlined in https://doi.org/10.1016/j.comgeo.2008.06.002
static vector<double> off_circumcenter(const vector<vector<double>> xs, double beta){
    vector<double> c1 = circumcenter(xs);
    vector<vector<double>> ps = Zeros<double>(3,2);
    vector<double> m(2);
    double distpq = 1e6;
    double temp;
    for (int i = 0; i<3; i++){
        temp = sqrt(pow(xs[(i+1)%3][0]-xs[i][0],2)+pow(xs[(i+1)%3][1]-xs[i][1],2));
        if (temp < distpq){
            distpq = temp;
            m = {(xs[(i+1)%3][0]+xs[i][0])/2, (xs[(i+1)%3][1]+xs[i][1])/2};
            ps[0] = xs[i];
            ps[1] = xs[(i+1)%3];
        }
    }
    ps[2] = c1;
    vector<double> c2 = circumcenter(ps);
    double distc1c2 = norm(c1-c2);
    vector<double> c;
    if (distc1c2 <= beta*distpq){
        c = c1;
    } else{
        c = c2 + 0.95*beta*distpq*(c2-m)/norm(c2-m);
    }
    return c;
}

/// find whether point is inside the circumcircle of a triangle
static bool inside_circumtri(const vector<vector<double>> xs, const vector<double> ps){
    vector<double> C = circumcenter(xs);
    double R = pow(xs[0][0]-C[0],2) + pow(xs[0][1] - C[1],2);
    bool D = pow(ps[0]-C[0],2) + pow(ps[1] - C[1],2) < R;
    return (D);
}
/// find whether or not point lies inside diametral circle of line segment
bool Mesh::inside_diametral(int hfid, vector<double> &ps){
    int eid = hfid2eid(hfid) - 1;
    int lid = hfid2lid(hfid) - 1;
    double v1[2] = {coords[elems[eid][lid]][0],coords[elems[eid][lid]][1]};
    double v2[2] = {coords[elems[eid][(lid+1)%3]][0],coords[elems[eid][(lid+1)%3]][1]};
    double r = sqrt(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2))/2;
    double p[2] = {(v1[0]+v2[0])/2, (v1[1]+v2[1])/2}; 
    double dist = sqrt(pow(ps[0]-p[0],2) + pow(ps[1]-p[1],2));
    return dist < r;
}
/// find whether point is inside a triangle
bool inside_tri(const vector<vector<double>> &xs, const vector<double> &ps){
    double val1 = (ps[0]-xs[1][0])*(xs[0][1]-xs[1][1]) - (xs[0][0]-xs[1][0])*(ps[1]-xs[1][1]);
    double val2 = (ps[0]-xs[2][0])*(xs[1][1]-xs[2][1]) - (xs[1][0]-xs[2][0])*(ps[1]-xs[2][1]);
    double val3 = (ps[0]-xs[0][0])*(xs[2][1]-xs[0][1]) - (xs[2][0]-xs[0][0])*(ps[1]-xs[0][1]);
    bool has_neg = (val1 <= 0) || (val2<=0) || (val3<=0);
    bool has_pos = (val1 >= 0) || (val2>=0) || (val3>=0);
    return !(has_neg && has_pos);
}

// tiny private functions for array operations
/// find smallest x and y values in point set
vector<double> min_array(const vector<vector<double>> &xs){
    int ndims = xs[0].size();
    vector<double> min_(ndims);

    for (int j = 0; j<ndims; j++){
        min_[j] = xs[0][j];
    }

    for (int i = 1; i<xs.size(); i++){
        for (int j = 0; j<ndims; j++){
            if (xs[i][j] < min_[j]){
                min_[j] = xs[i][j];
            }
        }
    }

    return min_;
}
/// find largest x and y values in point set
vector<double> max_array(const vector<vector<double>> &xs){
    int ndims = xs[0].size();
    vector<double> min_(ndims);

    for (int j = 0; j<ndims; j++){
        min_[j] = xs[0][j];
    }

    for (int i = 1; i<xs.size(); i++){
        for (int j = 0; j<ndims; j++){
            if (xs[i][j] > min_[j]){
                min_[j] = xs[i][j];
            }
        }
    }

    return min_;
}
/// find center in a point set
static vector<double> find_center(const vector<vector<double>> &xs){
    int nv = xs.size();
    int ndims = 2;
    vector<double> center(2);
    for (int i=0;i<nv; i++){
        for (int j=0;j<ndims;j++){
            center[j] += xs[i][j];
        }
    }
    for (int j = 0; j<ndims; j++){
        center[j] = center[j]/((double) nv);
    }
    return center;
}

/// find area of a triangle
double area_tri(const vector<vector<double>> &xs){
    int sz = xs[0].size();
    vector<double> u(sz);
    vector<double> v(sz);
    vector<double> N(3);
    u =  xs[1]-xs[0];
    v =  xs[2]-xs[0];
    if (sz == 2){
        N = {0,0,u[0]*v[1] - u[1]*v[0]};
    } else {
        N = {u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]};
    }
    
    return norm(N)/2;
}
/// find minimal angle in a triangle
static double min_angle(const vector<vector<double>> &xs){
    double e1[2] = {xs[1][0]-xs[0][0], xs[1][1]-xs[0][1]};
    double n1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1]);
    double e2[2] = {xs[2][0]-xs[1][0], xs[2][1]-xs[1][1]};
    double n2 = sqrt(e2[0]*e2[0] + e2[1]*e2[1]);
    double e3[2] = {xs[0][0]-xs[2][0], xs[0][1]-xs[2][1]};
    double n3 = sqrt(e3[0]*e3[0] + e3[1]*e3[1]);

    double theta1 = (180/M_PI)*acos((e1[0]*(-e2[0]) + e1[1]*(-e2[1]))/(n1*n2));
    double theta2 = (180/M_PI)*acos((e2[0]*(-e3[0]) + e2[1]*(-e3[1]))/(n2*n3));
    double theta3 = (180/M_PI)*acos((e3[0]*(-e1[0]) + e3[1]*(-e1[1]))/(n3*n1));
    return min(min(theta1,theta2),theta3);
}

double check_minangle(Mesh* DT){
    vector<vector<double>> ps = {{0,0},{0,0},{0,0}};
    double theta = 360.0;
    for (int i = 0; i<DT->nelems; i++){
        ps[0] = DT->coords[DT->elems[i][0]];
        ps[1] = DT->coords[DT->elems[i][1]];
        ps[2] = DT->coords[DT->elems[i][2]];
        theta = min(theta, min_angle(ps));
    }
    return theta;
}

/// facet nodes for high order triangular meshes
vector<vector<int>> obtain_trifacets(int degree){
    vector<vector<int>> facets = Zeros<int>(3,degree+1);
    facets[0][0] = 0;
    facets[0][1] = 1;
    facets[1][0] = 1;
    facets[1][1] = 2;
    facets[2][0] = 2;
    facets[2][1] = 0;
    if (degree > 1){
        int kk = 3;
        for (int i = 0; i<3; i++){
            for (int j = 2; j<2+degree-1; j++){
                facets[i][j] = kk;
                kk++;
            }
        }
    }
    return facets;
}