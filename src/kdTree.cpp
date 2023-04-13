#include "GeoMesh.hpp"

kdNode* create_kdTree(const vector<vector<double>> &coords){
    kdNode* root = new kdNode;

    int nnodes = coords.size();
    int ndims = coords[0].size();
    for (int i = 0; i<nnodes; i++){
        root->insert_node(i, 0, coords, ndims);
    }

    return root;
}

void kdNode::insert_node(int nid, int node_depth, const vector<vector<double>> &coords, int ndims){
    if (vid == -1){
        vid = nid;
        depth = node_depth;
        left = new kdNode;
        right = new kdNode;
        printnode();
        return;
    } else {
        int axis = depth % ndims;
        if (coords[nid][axis] <= coords[vid][axis]){
            left->insert_node(nid, node_depth+1, coords, ndims); 
        } else {
            right->insert_node(nid, node_depth+1, coords, ndims);
        }
        return;
    }
}

void kdNode::printnode(){
    cout << "depth: " << depth << " node: " << vid << endl;
    return;
}