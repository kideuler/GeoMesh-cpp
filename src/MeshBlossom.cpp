#include "GeoMesh.hpp"
#include <queue>
using namespace std;

void Tris2quads_blossom(Mesh *msh, Blossom* B){
    int nelems = msh->nelems;
    bool *processed = new bool[nelems];
    vector<vector<int>> quads(nelems, vector<int>(4,-1));
    int oppn[3] = {2,0,1};
    int eid1, eid2, lid1, lid2, hfid;
    int kk = 0;
    for (int i = 0; i<nelems; i++){
        eid1 = i;
        if (B->mate[i] > 0 && !processed[eid1]){
            hfid = B->mate[i];
            eid2 = hfid2eid(hfid)-1;
            lid2 = hfid2lid(hfid)-1;
            lid1 = hfid2lid(msh->sibhfs[eid2][lid2])-1;
            quads[kk] = {msh->elems[eid1][oppn[lid1]], msh->elems[eid1][lid1], msh->elems[eid2][oppn[lid2]], msh->elems[eid1][(lid1+1)%3]};
            processed[eid1] = true;
            processed[eid2] = true;
            kk++;
        } else if(!processed[eid1]) {
            quads[kk] = msh->elems[eid1];
            kk++;
        }
    }
    quads.resize(kk);
    msh->nelems = kk;

    for (int i = 0; i<kk; i++){
        for (int j = 0; j<4; j++){
            cout << quads[i][j] << " ";
        }
        cout << endl;
    }
    msh->elems = quads;
    delete processed;
    return;
}

Blossom Mesh2Graph(Mesh* msh){
    int hfid,eid,lid;
    Blossom B(msh->nelems);
    for (int i = 0; i<msh->nelems; i++){
        for (int j = 0; j<3; j++){
            hfid = msh->sibhfs[i][j];
            eid = hfid2eid(hfid)-1;
            lid = hfid2lid(hfid)-1;
            if (eid > -1){
                B.add_edge(i,eid,lid);
            }
        }
    }
    return B;
}


// Code here is based off of codeforce post https://codeforces.com/blog/entry/92339 on perfect matching blossom algorithm
void Blossom::add_edge(int u, int v) {
        G[u][v] = u;
        G[v][u] = v;
    }
void Blossom::add_edge(int u, int v, int lid) {
        G[u][v] = u;
        G[v][u] = v;
        E[u][v] = lid;
    }
void Blossom::match(int u, int v) {
        G[u][v] = G[v][u] = -1;
        mate[u] = elids2hfid(v+1, E[u][v]+1);
        mate[v] = elids2hfid(u+1, E[v][u]+1);
    }
vector<int> Blossom::trace(int x) {
    vector<int> vx;
    vx.reserve(50);
    while(true) {
        while(bl[x] != x) x = bl[x];
        if(!vx.empty() && vx.back() == x) break;
        vx.push_back(x);
        x = p[x];
    }
    return vx;
}
void Blossom::contract(int c, int x, int y, vector<int> &vx, vector<int> &vy) {
    b[c].clear();
    int r = vx.back();
    while(!vx.empty() && !vy.empty() && vx.back() == vy.back()) {
        r = vx.back();
        vx.pop_back();
        vy.pop_back();
    }
    b[c].push_back(r);
    b[c].insert(b[c].end(), vx.rbegin(), vx.rend());
    b[c].insert(b[c].end(), vy.begin(), vy.end());
    for(int i = 0; i <= c; i++) {
        G[c][i] = G[i][c] = -1;
    }
    for(int z : b[c]) {
        bl[z] = c;
        for(int i = 0; i < c; i++) {
            if(G[z][i] != -1) {
                G[c][i] = z;
                G[i][c] = G[i][z];
            }
        }
    }
}
vector<int> Blossom::lift(vector<int> &vx) {
    vector<int> A;
    while(vx.size() >= 2) {
        int z = vx.back(); vx.pop_back();
        if(z < n) {
            A.push_back(z);
            continue;
        }
        int w = vx.back();
        int i = (A.size() % 2 == 0 ? find(b[z].begin(), b[z].end(), G[z][w]) - b[z].begin() : 0);
        int j = (A.size() % 2 == 1 ? find(b[z].begin(), b[z].end(), G[z][A.back()]) - b[z].begin() : 0);
        int k = b[z].size();
        int dif = (A.size() % 2 == 0 ? i % 2 == 1 : j % 2 == 0) ? 1 : k - 1;
        while(i != j) {
            vx.push_back(b[z][i]);
            i = (i + dif) % k;
        }
        vx.push_back(b[z][i]);
    }
    return A;
}
int Blossom::solve() {
    int ans = 0;
    while(true) {
        fill(d.begin(), d.end(), 0);
        queue<int> Q;
        for(int i = 0; i < m; i++) bl[i] = i;
        for(int i = 0; i < n; i++) {
            if(mate[i] == -1) {
                Q.push(i);
                p[i] = i;
                d[i] = 1;
            }
        }
        int c = n;
        bool aug = false;
        while(!Q.empty() && !aug) {
            int x = Q.front(); Q.pop();
            if(bl[x] != x) continue;
            for(int y = 0; y < c; y++) {
                if(bl[y] == y && G[x][y] != -1) {
                    if(d[y] == 0) {
                        p[y] = x;
                        d[y] = 2;
                        p[hfid2eid(mate[y])-1] = y;
                        d[hfid2eid(mate[y])-1] = 1;
                        Q.push(hfid2eid(mate[y])-1);
                    }else if(d[y] == 1) {
                        vector<int> vx = trace(x);
                        vector<int> vy = trace(y);
                        if(vx.back() == vy.back()) {
                            contract(c, x, y, vx, vy);
                            Q.push(c);
                            p[c] = p[b[c][0]];
                            d[c] = 1;
                            c++;
                        }else {
                            aug = true;
                            vx.insert(vx.begin(), y);
                            vy.insert(vy.begin(), x);
                            vector<int> A = lift(vx);
                            vector<int> B = lift(vy);
                            A.insert(A.end(), B.rbegin(), B.rend());
                            for(int i = 0; i < (int) A.size(); i += 2) {
                                match(A[i], A[i + 1]);
                                if(i + 2 < (int) A.size()) add_edge(A[i + 1], A[i + 2]);
                            }
                        }
                        break;
                    }
                }
            }
        }
        if(!aug) return ans;
        ans++;
    }
}