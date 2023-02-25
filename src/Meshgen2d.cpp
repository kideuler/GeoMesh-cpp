#include "GeoMesh.hpp"
#include <unistd.h>
using namespace std;

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

/// subfunctions for delaunay triangulation
double eval_alpha(const vector<vector<double>> xs,double r_ref);
static vector<double> circumcenter(const vector<vector<double>> xs);
static bool inside_circumtri(const vector<vector<double>> xs, const vector<double> ps);
static bool inside_diametral(Mesh* DT, int hfid, int vid);
void Recursive_tri_delete(Mesh* DT, int hfid);
static bool Line_cross(const vector<double> &p1, const vector<double> &p2, const vector<double> &p3, const vector<double> &p4);
static bool Ray_in_triangle(Mesh* DT, int eid, int nid, int vid);
static double area_tri(const vector<vector<double>> &xs);
static double min_angle(const vector<vector<double>> &xs);
void recursive_delauney_flip(Mesh* DT, int eid, int lid);
void find_bad_tri_recursive(Mesh* DT, int tri, int *nbad, int vid);
void find_bad_tri_recursive(Mesh* DT, int tri, int *nbad, int vid, bool* exitf);
void find_bad_tri(Mesh* DT, int tri, int *nbad, int vid);
void flip_edge(Mesh* DT, int eid, int lid);
vector<int> facet_reorder(Mesh *DT, int *nsegs);
static vector<double> find_center(const vector<vector<double>> &xs);
void reorder(vector<vector<double>> &xs);

vector<bool> find_boundary_nodes(Mesh* DT){
    int nv = DT->coords.size();
    vector<bool> bnd(nv);
    for (int i = 0; i<DT->nelems; i++){
        for (int j = 0; j<3; j++){
            if (DT->sibhfs[i][j] == 0){
                bnd[DT->elems[i][j]] = true;
                bnd[DT->elems[i][(j+1)%3]] = true;
            }
        }
    }

    return bnd;
}

int find_hfid(Mesh* DT, int eid){
    int hfid = 0;
    for (int i = 0; i<3; i++){
        if(DT->sibhfs[eid][i] == 0){
            hfid = elids2hfid(eid+1,i+1);
            return hfid;
        }
    }
    cout << "not found" << endl;
    return hfid;
}
/**
 * @brief Refine a delaunay Mesh using Chews second algorithm
 * 
 * @param DT Triangultion data structure
 * @param r_ref radius of circumcircle (double)
 */
void GeoMesh_refine(Mesh* DT, double r_ref, Spline* spl){
    auto f = [r_ref](vector<double> xs) {return r_ref; };
    GeoMesh_refine(DT, f, spl);
}

/**
 * @brief Refine a delaunay Mesh using Rupperts algorithm
 * 
 * @param DT Triangultion data structure
 * @param r_ref radius of circumcircle (function of position)
 */
void GeoMesh_refine(Mesh* DT, function<double(vector<double>)> r_ref, Spline* spl){
    // random number generator
    default_random_engine re;
    uniform_real_distribution<double> unif(-1, 1);

    int n,i,nelems;
    bool exitl;
    double alpha,theta;
    int nv = (*DT).coords.size();
    vector<vector<double>> ps = Zeros<double>(3,2);
    int tri,eid,lid,ub;

    // finding total area and estimate to define maximum size bounds
    double total_area = 0.0;
    double minr = 1e6;
    vector<double> mid;
    for (n=0; n<DT->nelems; n++){
        ps[0] = (*DT).coords[(*DT).elems[n][0]];
        ps[1] = (*DT).coords[(*DT).elems[n][1]];
        ps[2] = (*DT).coords[(*DT).elems[n][2]];
        mid = (ps[0]+ps[1]+ps[2])/3.0;
        total_area += area_tri(ps);
        minr = min(r_ref(mid),minr);
    }
    double area_single = minr*minr/2;
    ub = DT->elems.size(); 
    if (1.2*total_area/area_single > ub){
        ub = (int) 1.2*total_area/area_single;
    }
    cout << ub << " " << total_area << " " << area_single << endl;
    DT->coords.resize(ub);
    DT->param.resize(ub);
    DT->elems.resize(ub);
    DT->sibhfs.resize(ub);
    DT->delete_elem.resize(ub);
    DT->on_boundary.resize(ub);
    vector<double> C;
    nelems = DT->nelems;
    int* order = new int[ub];
    for (i = 0; i<ub; i++){
        order[i] = i;
    }

    // main loop 
    n=0;
    bool inside_domain;
    int e;
    bool freed = false;
    while (n<(*DT).nelems){
        e = order[n];
        if (!(*DT).delete_elem[e]){
            ps[0] = (*DT).coords[(*DT).elems[e][0]];
            ps[1] = (*DT).coords[(*DT).elems[e][1]];
            ps[2] = (*DT).coords[(*DT).elems[e][2]];
            C = circumcenter(ps);
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
                (*DT).coords[nv] = C;
                tri = e;
                inside_domain = find_enclosing_tri(DT, &tri, nv);
                //cout << alpha << endl;
                if (!inside_domain){
                    if (tri == -1){
                        cout << "find triangle location failed: deleting point" << endl;
                        nv--;
                    } else {
                        if (inside_diametral(DT,tri, nv)){
                            Flip_Insertion_segment(DT, nv, tri, spl);
                        } else {
                            if (DT->on_boundary[e]){
                                Flip_Insertion_segment(DT, nv, tri, spl);
                            } else{
                            if (n == DT->nelems-1 || e == DT->nelems-1){
                                Flip_Insertion_segment(DT, nv, tri, spl);
                            } else {
                                nv--;
                                order[n] = DT->nelems-1;
                                order[DT->nelems-1] = e;
                                n--;
                            }
                            }
                        }
                    }
                } else {
                    bool stop = false;
                    if (DT->on_boundary[tri]){
                        int hfid = find_hfid(DT,tri);
                        if (inside_diametral(DT, hfid, nv)){
                            Flip_Insertion_segment(DT, nv, hfid, spl);
                            stop = true;
                        }
                    } else if(DT->sibhfs[tri][0] > 0) {
                        if(DT->on_boundary[hfid2eid(DT->sibhfs[tri][0])-1]){
                        int hfid = find_hfid(DT,hfid2eid(DT->sibhfs[tri][0])-1);
                        if (inside_diametral(DT, hfid, nv)){
                            Flip_Insertion_segment(DT, nv, hfid, spl);
                            stop = true;
                        }
                        }
                    } else if(DT->sibhfs[tri][1] > 0) {
                        if(DT->on_boundary[hfid2eid(DT->sibhfs[tri][1])-1]){
                        int hfid = find_hfid(DT,hfid2eid(DT->sibhfs[tri][1])-1);
                        if (inside_diametral(DT, hfid, nv)){
                            Flip_Insertion_segment(DT, nv, hfid, spl);
                            stop = true;
                        }
                        }
                    } else if(DT->sibhfs[tri][2] > 0) {
                        if(DT->on_boundary[hfid2eid(DT->sibhfs[tri][2])-1]){
                        int hfid = find_hfid(DT,hfid2eid(DT->sibhfs[tri][2])-1);
                        if (inside_diametral(DT, hfid, nv)){
                            Flip_Insertion_segment(DT, nv, hfid, spl);
                            stop = true;
                        }
                        }
                    }
                    


                    if (!stop){
                        Flip_Insertion(DT,&nv,tri);
                    }
                }
                nv++;

                if ((double) DT->nelems >= (0.95)*((double) DT->elems.size())){
                    cout << "approaching size bound, freeing up space" << endl;
                    if (!freed) {
                        delete_tris(DT,&n);
                        freed = true;
                        ub = ub*1.5;
                        DT->coords.resize(ub);
                        DT->param.resize(ub);
                        DT->elems.resize(ub);
                        DT->sibhfs.resize(ub);
                        DT->delete_elem.resize(ub);
                        DT->on_boundary.resize(ub);
                    } else {
                        break;
                    }
                }
            }
        }
        n++;
    }
    delete order;
    (*DT).coords.resize(nv);
    DT->param.resize(nv);
    delete_tris(DT);
    cout << "created " << (*DT).nelems-nelems << " triangles from refining the mesh" << endl;
}

/**
 * @brief Create constrained Delaunay Mesh from PSLG
 * 
 * @param segs Constrained segments (nsegs -by- 2)
 * @param xs Point data (nv -by- 2)
 * @return Constrained delaunay Mesh
 */
Mesh GeoMesh_Delaunay_Mesh(const vector<vector<int>> &segs, vector<vector<double>> &xs, vector<double> &params){

    // segs define boundary segments for the mesh
    Mesh DT = GeoMesh_Delaunay_Mesh(xs,params);
    DT.bwork.resize(DT.nelems);
    DT.facets.resize(DT.nelems);
    int nv = xs.size();
    int nsegs = segs.size();
    vector<vector<int>> onering(nv);
    vector<int> numonering(nv);

    int vid,i,j,k,hfid;
    for (i = 0; i<DT.nelems; i++){
        for (j = 0; j<3; j++){
            vid = DT.elems[i][j];
            onering[vid].push_back(elids2hfid(i+1,j+1));
        }
    }

    int oppeid,eid,lnid,vid2,lid,opplid,nf;
    nf = 0;
    bool exitf,exiti;
    for (i = 0; i<nsegs; i++){
        stack* head = NULL;
        vid = segs[i][0];
        vid2 = segs[i][1];
        j=0;
        exitf = false;
        while(j<onering[vid].size() && !exitf){
            eid = hfid2eid(onering[vid][j])-1;
            lnid = hfid2lid(onering[vid][j])-1;
            if (DT.elems[eid][(lnid+1)%3] == vid2){
                hfid = DT.sibhfs[eid][(lnid)%3];
                if (hfid != 0){
                    DT.facets[nf][0] = eid;
                    DT.facets[nf][1] = lnid;
                    DT.bwork[eid] = true;
                    nf++;
                }
                exitf = true;
            } else if (Line_cross(DT.coords[DT.elems[eid][(lnid+1)%3]], DT.coords[DT.elems[eid][(lnid+2)%3]], DT.coords[vid], DT.coords[vid2])){
                cout << "found ray with triangle with node " << eid << " " << lnid << endl; 
                cout << "changing Mesh to satisfy constrained edge " << vid<< "," << vid2 << endl;
                // Do constrained Del alg here
                hfid = DT.sibhfs[eid][(lnid+1)%3];
                exiti = false;
                while (!exiti){
                    push_stack(&head, hfid);
                    eid = hfid2eid(hfid)-1;
                    lid = hfid2lid(hfid)-1;
                    cout << eid << " " << lid << endl;
                    k = 0;
                    cout << DT.elems[eid][0] << " " << DT.elems[eid][1] << " " << DT.elems[eid][2] << endl;
                    cout << vid2 << endl;
                    while (k<3 && !exiti){
                        if (DT.elems[eid][k] == vid2){
                            exiti = true;
                        }
                        k++;
                    }

                    if (Line_cross(DT.coords[DT.elems[eid][(lid+1)%3]], DT.coords[DT.elems[eid][(lid+2)%3]], DT.coords[vid], DT.coords[vid2])){
                        hfid = DT.sibhfs[eid][(lid+1)%3];
                    } else if (Line_cross(DT.coords[DT.elems[eid][(lid+2)%3]], DT.coords[DT.elems[eid][(lid)%3]], DT.coords[vid], DT.coords[vid2])){
                        hfid = DT.sibhfs[eid][(lid+2)%3];
                    } else {
                        cout << "no line cross found, error in geometry" << endl;
                    }
                }
                exitf = true;
            }
            j++;
        }
        if (!exitf){
            cout << "no constrained edge found for edge: " << segs[i][0] << "," << segs[i][1] << endl;
        }

        // Flipping triangles
        while (head != NULL){
            hfid = head->hfid;
            pop_stack(&head);
            eid = hfid2eid(hfid)-1;
            lid = hfid2lid(hfid)-1;
            oppeid = hfid2eid(DT.sibhfs[eid][lid])-1;
            opplid = hfid2lid(DT.sibhfs[eid][lid])-1;
            cout << "flipping edge " << eid << " " << lid << " and " << oppeid << " " << opplid << endl;
            cout << "eid " << DT.elems[eid][0] << " " << DT.elems[eid][1] << " " << DT.elems[eid][2] << endl;
            cout << "oppeid " << DT.elems[oppeid][0] << " " << DT.elems[oppeid][1] << " " << DT.elems[oppeid][2] << endl;
            flip_edge(&DT,eid,lid);
            cout << "after swapping " << endl;
            cout << "eid " << DT.elems[eid][0] << " " << DT.elems[eid][1] << " " << DT.elems[eid][2] << endl;
            cout << "oppeid " << DT.elems[oppeid][0] << " " << DT.elems[oppeid][1] << " " << DT.elems[oppeid][2] << endl;
        }

    }

    // delete elements on opposite side of boundary;
    for (i=0;i<nf;i++){
        eid = DT.facets[i][0];
        lid = DT.facets[i][1];
        if (DT.sibhfs[eid][lid] != 0){
            Recursive_tri_delete(&DT, DT.sibhfs[eid][lid]);
        }
    }

    delete_tris(&DT);

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
        inside = find_enclosing_tri(&DT, &tri, vid);
        if (!inside){
            cout << "no enclosing tri found" << endl;
        }
        // inserting node into the Mesh using Bowyer-Watson algorithm
        Flip_Insertion(&DT,&vid,tri);
        if ((double) DT.nelems >= (0.8)*((double) ub)){
            cout << "approaching size bound, freeing up space" << endl;
            delete_tris(&DT);
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

    delete_tris(&DT);
    DT.coords.resize(nv);
    for (n=0; n<nv; n++){
        DT.coords[n][0] = xs[n][0];
        DT.coords[n][1] = xs[n][1];
    }
    cout << "created " << DT.nelems << " triangles from initial points" << endl;
    return DT;
}

/**
 * @brief Insert node into mesh using Lawson flipping algorithm
 * 
 * @param DT Delaunay Mesh passed by reference
 * @param vid Node to be inserted
 * @param tri_s Triangle that encloses the node
 */
void Flip_Insertion(Mesh* DT, int* vid, int tri_s){
    DT->delete_elem[tri_s] = true;
    int hfid,eid,lid;

    vector<int> tri = {DT->elems[tri_s][0], DT->elems[tri_s][1], DT->elems[tri_s][2]};
    vector<int> sib = {DT->sibhfs[tri_s][0],DT->sibhfs[tri_s][1],DT->sibhfs[tri_s][2]};
    
    // splitting triangle and adding subtriangles to the stack if they are not on the boundary
    stack* head = NULL;
    vector<int> eids = {DT->nelems, DT->nelems+1, DT->nelems+2};
    for (int i = 0; i<3; i++){
        DT->elems[eids[i]] = {*vid, tri[i], tri[(i+1)%3]};
        hfid = sib[i];
        DT->sibhfs[eids[i]] = {elids2hfid(eids[(i+2)%3]+1 ,3), hfid, elids2hfid(eids[(i+1)%3]+1,1)};
        if (hfid2eid(hfid) > 0){
            DT->sibhfs[hfid2eid(hfid)-1][hfid2lid(hfid)-1] = elids2hfid(eids[i]+1, 2);
            push_stack(&head, elids2hfid(eids[i]+1, 2));
            DT->on_boundary[eids[i]] = false;
        } else {
            DT->on_boundary[eids[i]] = true;
        }
    }

    DT->nelems+=3;
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    int oppeid,opplid;

    // loop through stack and flip edges
    while (head != NULL){
        hfid = head->hfid;
        pop_stack(&head);
        eid = hfid2eid(hfid)-1;
        lid = hfid2lid(hfid)-1;
        oppeid = hfid2eid(DT->sibhfs[eid][lid])-1;
        opplid = hfid2lid(DT->sibhfs[eid][lid])-1;
        xs[0] = DT->coords[DT->elems[oppeid][0]];
        xs[1] = DT->coords[DT->elems[oppeid][1]];
        xs[2] = DT->coords[DT->elems[oppeid][2]];
        
        if (inside_circumtri(xs, DT->coords[DT->elems[eid][0]])){
            flip_edge(DT,eid,1);

            hfid = DT->sibhfs[oppeid][1];
            if (hfid2eid(hfid) > 0){
                push_stack(&head, elids2hfid(oppeid+1,2));
            }
            hfid = DT->sibhfs[eid][1];
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
 * @param DT Delaunay Mesh passed by reference
 * @param vid Node to be inserted
 * @param hfid Triangle and local edge id of the segment to be split
 */
void Flip_Insertion_segment(Mesh* DT, int vid, int hfid, Spline* spl){

    // adding two triangles instead of three
    int nvS = spl->nv;
    stack* head = NULL;
    int eid = hfid2eid(hfid)-1;
    int lid = hfid2lid(hfid)-1;
    double a = DT->param[DT->elems[eid][lid]];
    double b = DT->param[DT->elems[eid][(lid+1)%3]];
    if (abs(b-a) < 1e-6 || nvS == 0){
        DT->coords[vid] = (DT->coords[DT->elems[eid][lid]] + DT->coords[DT->elems[eid][(lid+1)%3]])*0.5;
    } else {
        if (abs(b-a) > 0.5){
            if (b<a){
                b +=1;
            } else {
                a +=1;
            }
        }
        DT->param[vid] = (a+b)/2;
        DT->coords[vid] = spline_point_segment(spl, DT->param[DT->elems[eid][lid]],  DT->param[DT->elems[eid][(lid+1)%3]], 0.5);
    }
    double radius = 0.5*(norm(DT->coords[DT->elems[eid][lid]] - DT->coords[DT->elems[eid][(lid+1)%3]]));
    DT->delete_elem[eid] = true;
    DT->elems[DT->nelems] = {vid, DT->elems[eid][(lid+2)%3], DT->elems[eid][lid]};
    DT->sibhfs[DT->nelems] = {elids2hfid(DT->nelems+2,3), DT->sibhfs[eid][(lid+2)%3], 0};
    DT->elems[DT->nelems+1] = {vid, DT->elems[eid][(lid+1)%3] ,DT->elems[eid][(lid+2)%3]};
    DT->sibhfs[DT->nelems+1] = {0, DT->sibhfs[eid][(lid+1)%3], elids2hfid(DT->nelems+1,1)};
    if ( DT->sibhfs[eid][(lid+1)%3] != 0){
        DT->sibhfs[hfid2eid(DT->sibhfs[eid][(lid+1)%3])-1][hfid2lid(DT->sibhfs[eid][(lid+1)%3])-1] = elids2hfid(DT->nelems+2, 2);
        push_stack(&head, elids2hfid(DT->nelems+2, 2));
    }
    if ( DT->sibhfs[eid][(lid+2)%3] != 0){
        DT->sibhfs[hfid2eid(DT->sibhfs[eid][(lid+2)%3])-1][hfid2lid(DT->sibhfs[eid][(lid+2)%3])-1] = elids2hfid(DT->nelems+1, 2);
        push_stack(&head, elids2hfid(DT->nelems+1, 2));
    }
    DT->on_boundary[DT->nelems] = true;
    DT->on_boundary[DT->nelems+1] = true;

    DT->nelems+=2;
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    int oppeid,opplid;

    // loop through stack and flip edges
    while (head != NULL){
        hfid = head->hfid;
        pop_stack(&head);
        eid = hfid2eid(hfid)-1;
        lid = hfid2lid(hfid)-1;
        oppeid = hfid2eid(DT->sibhfs[eid][lid])-1;
        opplid = hfid2lid(DT->sibhfs[eid][lid])-1;
        xs[0] = DT->coords[DT->elems[oppeid][0]];
        xs[1] = DT->coords[DT->elems[oppeid][1]];
        xs[2] = DT->coords[DT->elems[oppeid][2]];
        
        if (inside_circumtri(xs, DT->coords[DT->elems[eid][0]])){
            flip_edge(DT,eid,lid);

            hfid = DT->sibhfs[oppeid][1];
            if (hfid2eid(hfid) > 0){
                push_stack(&head, elids2hfid(oppeid+1,2));
            }
            hfid = DT->sibhfs[eid][1];
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
 * @param DT Delaunay Mesh passed by reference
 * @param eid element to be flipped
 * @param lid edge to be flipped across
 */
void flip_edge(Mesh* DT, int eid, int lid){
    vector<int> tri1(3);
    vector<int> tri2(3);
    vector<int> sib1(3);
    vector<int> sib2(3);

    int hfid = DT->sibhfs[eid][lid];
    int oppeid = hfid2eid(hfid)-1;
    int opplid = hfid2lid(hfid)-1;

    int v1 = DT->elems[eid][(lid+2)%3];
    int v2 = DT->elems[eid][lid];
    int v3 = DT->elems[eid][(lid+1)%3];
    int v4 = DT->elems[oppeid][(opplid+2)%3];

    tri1 = {v1,v2,v4};
    tri2 = {v1,v4,v3};
    sib1 = {DT->sibhfs[eid][(lid+2)%3], DT->sibhfs[oppeid][(opplid+1)%3], elids2hfid(oppeid+1,1)};
    sib2 = {elids2hfid(eid+1,3), DT->sibhfs[oppeid][(opplid+2)%3], DT->sibhfs[eid][(lid+1)%3]};

    DT->elems[eid] = tri1;
    DT->elems[oppeid] = tri2;
    DT->sibhfs[eid] = sib1;
    DT->sibhfs[oppeid] = sib2;
    
    bool sib11 = false;
    bool sib12 = false;
    bool sib21 = false;
    bool sib22 = false;
    if (hfid2eid(sib1[0]) > 0){
        DT->sibhfs[hfid2eid(sib1[0])-1][hfid2lid(sib1[0])-1] = elids2hfid(eid+1,1);
    } else {
        sib11 = true;
    }
    if (hfid2eid(sib1[1]) > 0){
        DT->sibhfs[hfid2eid(sib1[1])-1][hfid2lid(sib1[1])-1] = elids2hfid(eid+1,2);
    } else {
        sib12 = true;
    }
    if (hfid2eid(sib2[1]) > 0){
        DT->sibhfs[hfid2eid(sib2[1])-1][hfid2lid(sib2[1])-1] = elids2hfid(oppeid+1,2);
    } else {
        sib21 = true;
    }
    if (hfid2eid(sib2[2]) > 0){
        DT->sibhfs[hfid2eid(sib2[2])-1][hfid2lid(sib2[2])-1] = elids2hfid(oppeid+1,3);
    } else {
        sib22 = true;
    }
    DT->on_boundary[eid] = sib11 || sib12;
    DT->on_boundary[oppeid] = sib21 || sib22;

    return;
}

/**
 * @brief Finding the triangle that encloses a node
 * 
 * @param DT Delaunay Mesh passed by reference
 * @param tri Starting triangle
 * @param vid node query
 * @return true 
 * @return false 
 */
bool find_enclosing_tri(Mesh* DT, int* tri, int vid){
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
        xs[0] = DT->coords[v1];
        xs[1] = DT->coords[v2];
        xs[2] = DT->coords[v3];
        if (inside_tri(xs,DT->coords[vid])){
            stop = true;
            return true;
        } else {
            double AB[2] = {DT->coords[v2][0]-DT->coords[v1][0], DT->coords[v2][1]-DT->coords[v1][1]};
            double BC[2] = {DT->coords[v3][0]-DT->coords[v2][0], DT->coords[v3][1]-DT->coords[v2][1]};
            double CA[2] = {DT->coords[v1][0]-DT->coords[v3][0], DT->coords[v1][1]-DT->coords[v3][1]};
            double AP[2] = {DT->coords[vid][0]-DT->coords[v1][0], DT->coords[vid][1]-DT->coords[v1][1]};
            double BP[2] = {DT->coords[vid][0]-DT->coords[v2][0], DT->coords[vid][1]-DT->coords[v2][1]};
            double CP[2] = {DT->coords[vid][0]-DT->coords[v3][0], DT->coords[vid][1]-DT->coords[v3][1]};
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


/// mod function for Bowyer-Watson
static int modi(int a, int b){
    int z = a - b*(a/b);
    if (z<0){
        z=z+b;
    }
    return z;
}
/**
 * @brief Add point to Mesh using Bowyer-Watson Algorithm
 * 
 * @param DT Mesh (Pass by reference)
 * @param vid Point in Mesh to be added
 * @param tri_s Starting triangle
 * @param refine whether the mesh is being refined flag
 */
void Bowyer_watson2d(Mesh* DT, int vid, int tri_s,bool refine){
    int nbad = 0;
    int i,j;
    vector<vector<double>> ps = {{0,0},{0,0},{0,0}};
    if (refine) {
        bool exitf = false;
        int eid,lid;
        find_bad_tri_recursive(DT, tri_s, &nbad, vid, &exitf);
        if (exitf) {
            eid = hfid2eid(nbad)-1;
            lid = hfid2lid(nbad)-1;
            DT->coords[vid] = (DT->coords[DT->elems[eid][(lid+1)%3]] + DT->coords[DT->elems[eid][lid]])/2.0;
            nbad = 0;
            find_bad_tri_recursive(DT, eid, &nbad, vid);
        }
    } else {
        find_bad_tri_recursive(DT, tri_s, &nbad, vid);
    }

    
    for (i = 0; i<nbad; i++){
        DT->bwork[i] = true;
    }

    // checking for duplicate elements
    for (i=0; i<nbad; i++){
        if ((*DT).bwork[i]){
            for (j=i+1; j<nbad; j++){
                if (DT->elems[DT->badtris[i]][0] == DT->elems[DT->badtris[j]][0] && DT->elems[DT->badtris[i]][1] == DT->elems[DT->badtris[j]][1] && DT->elems[DT->badtris[i]][2] == DT->elems[DT->badtris[j]][2]){
                    cout << "found duplicate elements, sibhfs is: " << endl;
                    cout <<  DT->sibhfs[DT->badtris[i]][0] << " " << DT->sibhfs[DT->badtris[i]][1] << " " << DT->sibhfs[DT->badtris[i]][2] << endl;
                    cout <<  DT->sibhfs[DT->badtris[j]][0] << " " << DT->sibhfs[DT->badtris[j]][1] << " " << DT->sibhfs[DT->badtris[j]][2] << endl;
                    DT->bwork[j] = false;
                }   
            }
        }
    }

    // adding facets to buffer
    int nsegs = 0;
    for (i=0; i<nbad; i++){
        for (j=0; j<3; j++){
            (*DT).facets[nsegs][0] = (*DT).elems[(*DT).badtris[i]][j];
            (*DT).facets[nsegs][1] = (*DT).elems[(*DT).badtris[i]][(j+1)%3];
            (*DT).vedge[nsegs] = (*DT).sibhfs[(*DT).badtris[i]][j];
            nsegs++;
        }
    }
    // find unique facets and reordering
    vector<int> order = facet_reorder(DT, &nsegs);

    for (i = 0; i<nsegs; i++){
        for (j=i+1; j<nsegs; j++){
            if (DT->facets[i][0] == DT->facets[i][1] && DT->facets[i][1] == DT->facets[i][1]){
                cout << "duplicate segments found: " << DT->facets[i][0] << "," << DT->facets[i][1] << endl;
            }
        }
    }

    int hfid;
    int ntri_b = (*DT).nelems;
    // Adding new triangles from Polygonal hole
    for (i=0; i<nsegs; i++){
        ps[0] = DT->coords[(*DT).facets[order[i]][0]];
        ps[1] = DT->coords[(*DT).facets[order[i]][1]];
        ps[2] = DT->coords[vid];
        (*DT).elems[(*DT).nelems] = {(*DT).facets[order[i]][0],(*DT).facets[order[i]][1],vid};
        (*DT).sibhfs[(*DT).nelems][0] = (*DT).vedge[order[i]];
        if ((*DT).vedge[order[i]] != 0){
            hfid = (*DT).vedge[order[i]];
            (*DT).sibhfs[hfid2eid(hfid)-1][hfid2lid(hfid)-1] = elids2hfid((*DT).nelems+1,1);
        } else {
            (*DT).on_boundary[(*DT).nelems] = true;
        }
        int eid1 = ntri_b + modi(i+1,nsegs) + 1;
        int eid2 = ntri_b + modi(i-1,nsegs) + 1;
        (*DT).sibhfs[(*DT).nelems][1] = elids2hfid(ntri_b + modi(i+1,nsegs) + 1, 3);
        (*DT).sibhfs[(*DT).nelems][2] = elids2hfid(ntri_b + modi(i-1,nsegs) + 1, 2);
        (*DT).nelems++;
        
    }
}

/**
 * @brief Find triangles which have point in their circumcircle (optimal recursive algorithm)
 * 
 * @param DT Trangulation passed by reference
 * @param tri triangle to be focused on
 * @param nbad number of bad triangle passed by reference
 * @param vid point in question
 */
void find_bad_tri_recursive(Mesh* DT, int tri, int *nbad, int vid){
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    vector<double> ps = (*DT).coords[vid];
    xs[0] = (*DT).coords[(*DT).elems[tri][0]];
    xs[1] = (*DT).coords[(*DT).elems[tri][1]];
    xs[2] = (*DT).coords[(*DT).elems[tri][2]];
    int eid;
    if (inside_circumtri(xs,ps) && !DT->delete_elem[tri]){
        (*DT).delete_elem[tri] = true;
        (*DT).badtris[*nbad] = tri;
        (*nbad)++;
        eid = hfid2eid((*DT).sibhfs[tri][0]);
        if (eid!=0 && !(*DT).delete_elem[eid-1]){
            find_bad_tri_recursive(DT,eid-1,nbad,vid);
        }
        eid = hfid2eid((*DT).sibhfs[tri][1]);
        if (eid!=0 && !(*DT).delete_elem[eid-1]){
            find_bad_tri_recursive(DT,eid-1,nbad,vid);
        }
        eid = hfid2eid((*DT).sibhfs[tri][2]);
        if (eid!=0 && !(*DT).delete_elem[eid-1]){
            find_bad_tri_recursive(DT,eid-1,nbad,vid);
        }
    } else {
        DT->delete_elem[tri] = false;
        return;
    }
}

/**
 * @brief Find triangles which have point in their circumcircle (optimal recursive algorithm) and exit if the point is outside domain
 * 
 * @param DT Trangulation passed by reference
 * @param tri triangle to be focused on
 * @param nbad number of bad triangle passed by reference
 * @param vid point in question
 * @param exitf Flag if the point lies outside the domain
 */
void find_bad_tri_recursive(Mesh* DT, int tri, int *nbad, int vid, bool* exitf){
    if (!*exitf) {
        vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
        double d;
        int j;
        vector<double> u(2);
        vector<double> v(2);
        vector<double> ps = (*DT).coords[vid];
        xs[0] = (*DT).coords[(*DT).elems[tri][0]];
        xs[1] = (*DT).coords[(*DT).elems[tri][1]];
        xs[2] = (*DT).coords[(*DT).elems[tri][2]];
        int eid,lid;
        if (inside_circumtri(xs,ps) && !DT->delete_elem[tri]){
            (*DT).delete_elem[tri] = true;
            (*DT).badtris[*nbad] = tri;
            (*nbad)++;

            for(j = 0; j<3; j++){
            eid = hfid2eid((*DT).sibhfs[tri][j]);
            if (eid!=0 && !(*DT).delete_elem[eid-1]){
                find_bad_tri_recursive(DT,eid-1,nbad,vid,exitf);
            } else if (!(*DT).delete_elem[eid-1]){
                u =  DT->coords[DT->elems[tri][(j+1)%3]] - DT->coords[DT->elems[tri][j]];
                v =  DT->coords[vid] - DT->coords[DT->elems[tri][j]];
                *exitf = u[0]*v[1] - u[1]*v[0] < 0;
                if (*exitf) {
                    *nbad = elids2hfid(tri+1,j+1);
                    (*DT).delete_elem[tri] = false;
                    return;
                } 
            }
            }
        } else {
            DT->delete_elem[tri] = false;
            return;
        } 
    } else {
        return;
    }
}

/**
 * @brief Find triangles which have point in their circumcircle (unoptimal slow algorithm)
 * 
 * @param DT Trangulation passed by reference
 * @param tri triangle to be focused on
 * @param nbad number of bad triangle passed by reference
 * @param vid point in question
 */
void find_bad_tri(Mesh* DT, int tri, int *nbad, int vid){
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    vector<double> ps = (*DT).coords[vid];
    for (int i = 0; i<DT->nelems; i++){
        if (!DT->delete_elem[i]){
            xs[0] = (*DT).coords[(*DT).elems[i][0]];
            xs[1] = (*DT).coords[(*DT).elems[i][1]];
            xs[2] = (*DT).coords[(*DT).elems[i][2]];
            if (inside_circumtri(xs,ps)){
                (*DT).delete_elem[i] = true;
                (*DT).badtris[*nbad] = i;
                (*nbad)++;
            }
        }
    }
}

/**
 * @brief Delete triangles and reorganize data in Mesh DT
 * 
 * @param DT Mesh DT passed by reference
 */
void delete_tris(Mesh* DT){
    int i,j;
    int nelems = 0;
    int sz = (*DT).sibhfs[0].size();
    vector<int> idx((*DT).nelems);
    vector<int> idx_rev((*DT).nelems);
    fill(idx_rev.begin(), idx_rev.end(),-1);

    // delete triangles to be deleted
    for (i = 0; i<(*DT).nelems; i++){
        if (!(*DT).delete_elem[i]){
            (*DT).elems[nelems] = (*DT).elems[i];
            (*DT).sibhfs[nelems] = (*DT).sibhfs[i];
            (*DT).delete_elem[nelems] = false;
            idx[nelems] = i;
            nelems++;
        } else {
            DT->delete_elem[i] = false;
        }
    }
    (*DT).nelems = nelems;

    for (i = 0; i<(*DT).nelems; i++){
        idx_rev[idx[i]] = i;
    }

    int nside;
    int hfid, eid, lid;
    for (i = 0; i<nelems; i++){
        nside = 0;
        for (j = 0; j < sz; j++){
            hfid = (*DT).sibhfs[i][j];
            if (!hfid == 0){
                eid = hfid2eid(hfid);
                lid = hfid2lid(hfid);
                if (!DT->delete_elem[eid-1]){
                    if (idx_rev[eid-1] == -1){
                        DT->sibhfs[i][j] = 0;
                    } else {
                        (*DT).sibhfs[i][j] = elids2hfid(idx_rev[eid-1]+1,lid);
                    }
                } else {
                    DT->sibhfs[i][j] = 0;
                }
                nside++;
            } else {
                DT->sibhfs[i][j] = 0;
            }
        }
        if (nside == sz){
            (*DT).on_boundary[i] = false;
        } else {
            (*DT).on_boundary[i] = true;
        }
        (*DT).delete_elem[i] = false;
    }
    DT->elems.resize(nelems);
    DT->sibhfs.resize(nelems);
}

/**
 * @brief Delete triangles and reorganize data in Mesh DT and keep track of specific triangle tri
 * 
 * @param DT Mesh DT passed by reference
 * @param tri triangle pased by reference
 */
void delete_tris(Mesh* DT, int* tri){
    int i,j;
    int nelems = 0;
    int sz = (*DT).sibhfs[0].size();
    vector<int> idx((*DT).nelems);
    vector<int> idx_rev((*DT).nelems);

    // delete triangles to be deleted
    for (i = 0; i<(*DT).nelems; i++){
        if (!(*DT).delete_elem[i]){
            (*DT).elems[nelems] = (*DT).elems[i];
            (*DT).sibhfs[nelems] = (*DT).sibhfs[i];
            (*DT).delete_elem[nelems] = false;
            idx[nelems] = i;
            nelems++;
        } else {
            DT->delete_elem[i] = false;
        }
    }
    (*DT).nelems = nelems;

    for (i = 0; i<(*DT).nelems; i++){
        idx_rev[idx[i]] = i;
    }

    int nside;
    int hfid, eid, lid;
    for (i = 0; i<nelems; i++){
        nside = 0;
        for (j = 0; j < sz; j++){
            hfid = (*DT).sibhfs[i][j];
            if (!hfid == 0){
                eid = hfid2eid(hfid);
                lid = hfid2lid(hfid);
                if (!DT->delete_elem[eid-1]){
                    if (idx_rev[eid-1] == -1){
                        DT->sibhfs[i][j] = 0;
                    } else {
                        (*DT).sibhfs[i][j] = elids2hfid(idx_rev[eid-1]+1,lid);
                    }
                } else {
                    DT->sibhfs[i][j] = 0;
                }
                nside++;
            } else {
                DT->sibhfs[i][j] = 0;
            }
        }
        if (nside == sz){
            (*DT).on_boundary[i] = false;
        } else {
            (*DT).on_boundary[i] = true;
        }
        (*DT).delete_elem[i] = false;
    }
    *tri = idx_rev[*tri];
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
static vector<double> off_circumcenter(const vector<vector<double>> xs){
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
    return c1;
}

/// find whether point is inside the circumcircle of a triangle
static bool inside_circumtri(const vector<vector<double>> xs, const vector<double> ps){
    vector<double> C = circumcenter(xs);
    double R = pow(xs[0][0]-C[0],2) + pow(xs[0][1] - C[1],2);
    bool D = pow(ps[0]-C[0],2) + pow(ps[1] - C[1],2) < R;
    return (D);
}
/// find whether or not point lies inside diametral circle of line segment
static bool inside_diametral(Mesh* DT, int hfid, int vid){
    int eid = hfid2eid(hfid) - 1;
    int lid = hfid2lid(hfid) - 1;
    double v1[2] = {DT->coords[DT->elems[eid][lid]][0],DT->coords[DT->elems[eid][lid]][1]};
    double v2[2] = {DT->coords[DT->elems[eid][(lid+1)%3]][0],DT->coords[DT->elems[eid][(lid+1)%3]][1]};
    double r = sqrt(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2))/2;
    double p[2] = {(v1[0]+v2[0])/2, (v1[1]+v2[1])/2}; 
    double dist = sqrt(pow(DT->coords[vid][0]-p[0],2) + pow(DT->coords[vid][1]-p[1],2));
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
/// Find unique facets and reorder them for bowyer watson algorithm
vector<int> facet_reorder(Mesh* DT, int *nsegs){
    int nsegs2 = 0;
    int i,j;
    vector<vector<int>> segs2 = Zeros<int>(*nsegs,2);
    vector<int> vedge2(*nsegs);
    for (i=0;i<(*nsegs);i++){
        (*DT).bwork[i] = true;
    }
    for (i=0; i<(*nsegs); i++){
        if ((*DT).bwork[i]){
            for (j=0; j<(*nsegs); j++){
                if ((((*DT).facets[i][0] == (*DT).facets[j][0] && (*DT).facets[i][1] == (*DT).facets[j][1]) || \
                ((*DT).facets[i][0] == (*DT).facets[j][1] && (*DT).facets[i][1] == (*DT).facets[j][0])) && i!=j){
                    (*DT).bwork[i] = false;
                    (*DT).bwork[j] = false;
                }
            }
        }

        if ((*DT).bwork[i]){
            segs2[nsegs2][0] = (*DT).facets[i][0];
            segs2[nsegs2][1] = (*DT).facets[i][1];
            vedge2[nsegs2] = (*DT).vedge[i];
            nsegs2++;
        }
    }

    if (nsegs2 == 0){
        for (i=0; i<(*nsegs)/3; i++){
            cout << DT->badtris[i] << " ";
        }
        cout << endl;
        for (i=0; i<(*nsegs); i++){
            cout << (*DT).facets[i][0] << " " << (*DT).facets[i][1] << endl;
        }
    }

    *nsegs = nsegs2;
    for (i=0; i<(*nsegs); i++){
        (*DT).facets[i][0] = segs2[i][0];
        (*DT).facets[i][1] = segs2[i][1];
        (*DT).vedge[i] = vedge2[i];
        (*DT).bwork[i] = false;
    }
    vector<int> order(*nsegs);

    order[0] = 0;
    int vid;
    bool exitf;
    for (i=1; i<(*nsegs); i++){
        vid = (*DT).facets[order[i-1]][1];
        exitf = false;
        j=0;
        while ((j<(*nsegs)) && !exitf){
            if (vid == (*DT).facets[j][0]){
                order[i] = j;
                exitf = true;
            } else {
                j++;
            }
        }
    }
    return order;
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

/// reorder pointset by angle
int myrandom (int i) { return std::rand()%i;}
void reorder(vector<vector<double>> &xs){
    srand ( unsigned ( time(0) ) );
    int nv = xs.size();
    int ndims = xs[0].size();
    assert(ndims == 2);
    vector<double> quantity(nv);
    random_shuffle(xs.begin(), xs.end(),myrandom);
    return;

    vector<double> center = find_center(xs);

    // measuring some quantity for each point
    vector<double> e = {1,0};
    vector<double> u(ndims);
    for (int n = 0; n<nv; n++){
        u[0] = xs[n][0] - center[0];
        u[1] = xs[n][1] - center[1];
        quantity[n] = acos(inner(u,e)/norm(u));
    }

    // basic sorting algo for order points based on some quantity
    int min_idx; 
    int j,n;
    for (n = 0; n<nv-1; n++){
        min_idx = n;
        for (j=n+1; j<nv; j++){
            if (quantity[j] < quantity[min_idx]){
                min_idx = j;
            }

            if (min_idx != n){
                swap(quantity[min_idx],quantity[n]);
                swap(xs[min_idx],xs[n]);
            }
        }
    }
}
/// find area of a triangle
static double area_tri(const vector<vector<double>> &xs){
    vector<double> u(2);
    vector<double> v(2);
    u =  xs[1]-xs[0];
    v =  xs[2]-xs[0];
    return abs(u[0]*v[1] - u[1]*v[0])/2;
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