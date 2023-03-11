#include "GeoMesh.hpp"
using namespace std;

void Mesh::Graph_init() {
    int nfacets = sibhfs[0].size();
    ncolors = nfacets;
    if (GraphNodes.size() > 0 && GraphEdges.size() > 0){
        for (int ii = 0; ii<GraphEdges.size(); ii++){
            GraphEdges[ii].color = 0;
        }
        return;
    }

    // resizing arrays
    GraphNodes.assign(nelems, vector<int>(nfacets, -1));
    GraphEdges.resize(nfacets*nelems);

    int E = 0;
    int eid,lid,hfid;
    for (int ii = 0; ii<nelems; ii++){
        for (int jj = 0; jj<nfacets; jj++){
            if (GraphNodes[ii][jj] < 0){
                GraphEdges[E].color = 0;
                GraphEdges[E].hfids[0] = elids2hfid(ii+1,jj+1);
                GraphNodes[ii][jj] = E;

                hfid = sibhfs[ii][jj];
                eid = hfid2eid(hfid)-1;
                lid = hfid2lid(hfid)-1;
                if (eid >= 0){
                    GraphEdges[E].hfids[1] = hfid;
                    GraphNodes[eid][lid] = E;
                } else {
                    GraphEdges[E].hfids[1] = -1;
                }
                E++;
            }
        }
    }

    GraphEdges.resize(E);
    return;
}

queue<int> Mesh::Graph_Color_greedy(bool userand){
    int nedges = GraphEdges.size();
    int color;
    queue<int> Q;
    for (int ii = 0; ii<nedges; ii++){
        color = find_nonconflict_color(ii);
        if (color == -1){
            Q.push(ii);
        }
        GraphEdges[ii].color = color;
    }

    bwork.assign(nedges,false);
    return Q;
}


// necessary functions
void Mesh::resolve_conflicting_edges(queue<int> Q){

} // overall function

bool Mesh::check_swap(int edge_id){
    
    int hfid1, hfid2, eid1, eid2, lid1, lid2, c;
    
    bool* used = new bool[ncolors]; 
    for (int i = 0; i<ncolors; i++){used[i] = false;}

    hfid1 = GraphEdges[edge_id].hfids[0];
    hfid2 = GraphEdges[edge_id].hfids[1];
    eid1 = hfid2eid(hfid1)-1;
    lid1 = hfid2lid(hfid1)-1;
    eid2 = hfid2eid(hfid2)-1;
    lid2 = hfid2lid(hfid2)-1;

    int nfacets = sibhfs[eid1].size();
    for (int i = 1; i<nfacets; i++){
        c = GraphEdges[GraphNodes[eid1][(lid1+i)%nfacets]].color;
        if (c > 0){
            if (!used[c-1]) {used[c-1] = true;} else{ return false;}
        }
    }

    for (int i = 0; i<ncolors; i++){used[i] = false;}
    if (eid2 >= 0){
        int nfacets = sibhfs[eid2].size();
        for (int i = 1; i<nfacets; i++){
            c = GraphEdges[GraphNodes[eid2][(lid2+i)%nfacets]].color;
            if (c > 0){
                if (!used[c-1]) {used[c-1] = true;} else{ delete used; return false;}
            }
        }
    }
    delete used; return true;
}

//bool Mesh::resolve(int edge_id); // tries to resolve edge id false if no immidiate resolution (hard)

//bool Mesh::edge_swap(int edge_id); // swaps edge while not revisiting any other previously visited edge. false if out of options (medium)

//void Mesh::create_conflict(int edge_id, queue<int> Q); // creates new conflict outlined in paper (easy i think)


int Mesh::find_nonconflict_color(int edge_id){
    default_random_engine re;

    int hfid1, hfid2, eid1, eid2, lid1, lid2;
    
    bool* used = new bool[ncolors]; 
    for (int i = 0; i<ncolors; i++){used[i] = false;}
    int* C = new int[ncolors];
    int count = 0;
    int c;

    hfid1 = GraphEdges[edge_id].hfids[0];
    hfid2 = GraphEdges[edge_id].hfids[1];
    eid1 = hfid2eid(hfid1)-1;
    lid1 = hfid2lid(hfid1)-1;
    eid2 = hfid2eid(hfid2)-1;
    lid2 = hfid2lid(hfid2)-1;

    int nfacets = sibhfs[eid1].size();
    for (int i = 1; i<nfacets; i++){
        c = GraphEdges[GraphNodes[eid1][(lid1+i)%nfacets]].color;
        if (c > 0){
            used[c-1] = true;
        }
    }

    if (eid2 >= 0){
        int nfacets = sibhfs[eid2].size();
        for (int i = 1; i<nfacets; i++){
            c = GraphEdges[GraphNodes[eid2][(lid2+i)%nfacets]].color;
            if (c > 0){
                used[c-1] = true;
            }
        }
    }

    for (int i = 0; i<ncolors; i++){
        if (!used[i]){ C[count] = i+1; count++;}
    }

    int val;
    if (count == 0){
        val = -1;
    } else {
        uniform_int_distribution<int> unif(0, count-1);
        val = C[unif(re)];
    }

    delete used;
    delete C;
    return val;
}

bool check_graph(const Mesh &msh){
    bool pass = true;
    int c;
    bool *used = new bool[msh.ncolors+1];
    for (int ii = 0; ii<msh.nelems; ii++){
        for (int jj=0; jj<msh.ncolors+1; jj++){used[jj]=false;}
        for (int jj=0; jj<3; jj++){
            c = msh.GraphEdges[msh.GraphNodes[ii][jj]].color;
            if (c >= 0){
                if (used[c]){
                    pass = false;
                    cout << "element: " << ii+1 << " is miscolored with color " << c << endl;
                } else {
                    used[c] = true;
                }
            }
        }
    }

    delete used;
    return pass;
}