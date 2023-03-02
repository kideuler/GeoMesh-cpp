#include "GeoMesh.hpp"
using namespace std;

vector<vector<int>> find_boundary(Mesh* msh, bool findloop){
    int nv = msh->coords.size();
    vector<vector<int>> bdy = Zeros<int>(nv,2);
    int nb=0;
    for (int i = 0; i<msh->nelems; i++){
        for (int j = 0; j<3; j++){
            if (msh->sibhfs[i][j] == 0){
                bdy[nb][0] = msh->elems[i][j];
                bdy[nb][1] = msh->elems[i][(j+1)%3];
                nb++;
            }
        }
    }
    bdy.resize(nb);

    int j, v;
    if (findloop){
        for (int i = 0; i<nb; i++){
            j = i+1;
            v = bdy[i][1];
            while (j < nb){
                if (bdy[j][0] == v){
                    swap(bdy[i+1], bdy[j]);
                    break;
                }
                j++;
            }

        }
    }
    return bdy;
}

// Sanity functions
/**
 * @brief Check whether Mesh has valid half-faces
 * 
 * @param DT Mesh passed by reference
 * @return true 
 * @return false 
 */
bool check_sibhfs(Mesh* DT){
    const vector<vector<int>> edges = {{0,1},{1,2},{2,0}};
    int nelems = DT->nelems;
    int hfid,eid,lid;
    bool check = true;
    for (int i = 0; i<nelems; i++){
        for (int j = 0; j<3; j++){
            hfid = DT->sibhfs[i][j];
            if (hfid != 0){
            eid = hfid2eid(hfid);
            lid = hfid2lid(hfid);
            if (eid > nelems || lid > 3 || eid < 0 || DT->delete_elem[eid-1]){
                cout << "sibhfs is wrong at elem: " << i << " face: " << j << " oppeid: " << eid-1 << " opplid: " << lid-1 << endl;
                check = false;
            } 
            if (DT->elems[i][edges[j][0]] != DT->elems[eid-1][edges[lid-1][1]] || DT->elems[i][edges[j][1]] != DT->elems[eid-1][edges[lid-1][0]]){
                cout << "sides dont match is wrong at elem: " << i << " face: " << j << " oppeid: " << eid-1 << " opplid: " << lid-1 << " sibhfs: " << hfid << endl;
                check = false;
            }   
            }
        }
    }
    return check;
} 
/**
 * @brief Find whether Mesh has valid Jacobian determinants for all triangles
 * 
 * @param DT Mesh passed by reference
 * @return true 
 * @return false 
 */
bool check_jacobians(Mesh* DT){
    const vector<vector<double>> dphi = {{-1,-1},{1,0},{0,1}};
    vector<vector<double>> ps = {{0,0},{0,0},{0,0}};
    vector<vector<double>> J;
    double detJ;
    bool check = true;
    for (int i = 0; i<DT->nelems; i++){
        ps[0] = DT->coords[DT->elems[i][0]];
        ps[1] = DT->coords[DT->elems[i][1]];
        ps[2] = DT->coords[DT->elems[i][2]];
        J = Transpose(ps)*dphi;
        detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (detJ < 0){
            cout << "negative jacobian at eid: " << i << endl;
            check = false;
        }
    }
    return false;
}

void WrtieVtk_tri(const Mesh &msh, string filename){
    FILE *fid;
    fid = fopen(filename.c_str(),"w");
    fprintf(fid,"# vtk DataFile Version 3.0\n");
    fprintf(fid,"This file was written using writevtk_unstr.m\n");
    fprintf(fid,"ASCII\n");

    int nv = msh.coords.size();
    int ndims = msh.coords[0].size();
    int nelems = msh.elems.size();
    int nelem_coeff, cell_type;
    bool is_quadratic = msh.degree==2;
    if (is_quadratic){
        nelem_coeff = 7;
        cell_type = 22;
    } else {
        nelem_coeff = 4;
        cell_type = 5;
    }

    // header for points
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid, "POINTS %i double",nv);

    // write out vertices
    if (ndims == 2){
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],0.0);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    }

    // write out connectivity header
    fprintf(fid,"\n\nCELLS %i %i", nelems, nelem_coeff*nelems);
    if (is_quadratic){
        for (int i = 0; i<nelems; i++){
            fprintf(fid,"\n%d %d %d %d %d %d %d",6,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2],msh.elems[i][3],msh.elems[i][4],msh.elems[i][5]);
        }
    } else {
        for (int i = 0; i<nelems; i++){
            fprintf(fid,"\n%d %d %d %d",3,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2]);
        }
    }

    // write out cell types
    fprintf(fid, "\n\nCELL_TYPES %i", nelems);
    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\n%i",cell_type);
    }
    fclose(fid);
}

void WrtieVtk_mixed(const Mesh &msh, string filename){
    FILE *fid;
    bool q = msh.elems[0].size() > 3;
    fid = fopen(filename.c_str(),"w");
    fprintf(fid,"# vtk DataFile Version 3.0\n");
    fprintf(fid,"This file was written using writevtk_unstr.m\n");
    fprintf(fid,"ASCII\n");

    int nv = msh.coords.size();
    int ndims = msh.coords[0].size();
    int nelems = msh.nelems;
    

    bool* quads = new bool[nelems];
    int ntotal=0;
    if (q){
        for (int i = 0; i<nelems; i++){
            if (msh.elems[i][3] >= 0 && msh.elems[i][3] < nv){
                quads[i] = true;
                ntotal += 5;
            } else {
                ntotal += 4;
                quads[i] = false;
            }
        }
    } else {
        for (int i = 0; i<nelems; i++){
            quads[i] = false;
        }
        ntotal = 4*nelems;
    }

    // header for points
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid, "POINTS %i double",nv);

    // write out vertices
    if (ndims == 2){
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],0.0);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    }

    // write out connectivity header
    fprintf(fid,"\n\nCELLS %i %i", nelems, ntotal);
    for (int i = 0; i<nelems; i++){
        if (!msh.delete_elem[i]){
            if (quads[i]){
                fprintf(fid,"\n%d %d %d %d %d",4,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2],msh.elems[i][3]);
            } else{
                fprintf(fid,"\n%d %d %d %d",3,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2]);
            }
        }
    }


    // write out cell types
    fprintf(fid, "\n\nCELL_TYPES %i", nelems);
    for (int i = 0; i<nelems; i++){
        if (!msh.delete_elem[i]){
            if (quads[i]){
                fprintf(fid,"\n%i",9);
            } else{
                fprintf(fid,"\n%i",5);
            }
        }
    }

    fclose(fid);
    delete quads;
}

void WrtieVtk_tri(const Mesh &msh, const vector<double> &data, string filename){
    FILE *fid;
    fid = fopen(filename.c_str(),"w");
    fprintf(fid,"# vtk DataFile Version 3.0\n");
    fprintf(fid,"This file was written using writevtk_unstr.m\n");
    fprintf(fid,"ASCII\n");

    int nv = msh.coords.size();
    int ndims = msh.coords[0].size();
    int nelems = msh.elems.size();
    int nelem_coeff, cell_type;
    bool is_quadratic = msh.degree==2;
    if (is_quadratic){
        nelem_coeff = 7;
        cell_type = 22;
    } else {
        nelem_coeff = 4;
        cell_type = 5;
    }

    // header for points
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid, "POINTS %i double",nv);

    // write out vertices
    if (ndims == 2){
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],0.0);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    }

    // write out connectivity header
    fprintf(fid,"\n\nCELLS %i %i", nelems, nelem_coeff*nelems);
    if (is_quadratic){
        for (int i = 0; i<nelems; i++){
            fprintf(fid,"\n%d %d %d %d %d %d %d",6,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2],msh.elems[i][3],msh.elems[i][4],msh.elems[i][5]);
        }
    } else {
        for (int i = 0; i<nelems; i++){
            fprintf(fid,"\n%d %d %d %d",3,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2]);
        }
    }

    // write out cell types
    fprintf(fid, "\n\nCELL_TYPES %i", nelems);
    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\n%i",cell_type);
    }

    fprintf(fid, "\n\nPOINT_DATA %i",nv);
    fprintf(fid, "\nSCALARS meshdata double\n");
    fprintf(fid, "LOOKUP_TABLE default\n");
    for (int i = 0; i<nv; i++){
        fprintf(fid,"%g\n",data[i]);
    }

    fclose(fid);
}

Mesh ReadObj_tri(string filename){
    Mesh msh;

    string line;
    char c,c2;
    int i, j, k;
    double x, y, z;
    string si, sj, sk;

    ifstream in(filename);
    bool normal;
    bool texture;
    while ( getline(in, line) ){
        normal = (line.find("vn")<line.length());
        texture = (line.find("vt")<line.length());
        if ( line.find_first_of( "vVfF" ) == string::npos ) continue; 

        istringstream ss( line );                           
        ss >> c;     
        switch ( c )
        {
            case 'v':                                         
            case 'V':   
                if (!normal && !texture){                                      
                ss >> x >> y >> z;  
                msh.coords.push_back({x, y, z});     
                }
                break;
            case 'f':                                         
            case 'F':                                         
                ss >> si >> sj >> sk;                        
                i = stoi(si);  j = stoi(sj);  k = stoi( sk );
                msh.elems.push_back({i-1, j-1, k-1});
                break;
      }
   }
   in.close();
   msh.nelems = msh.elems.size();
   msh.delete_elem.resize(msh.nelems);
   msh.normals.resize(0);
   return msh;
}