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

int kdNode::find_nearest_node(const vector<double> &xs, const vector<vector<double>> &coords){
    int ndims = xs.size();
    assert(ndims == coords[0].size());
    double best_distance = 1e308;
    int ans = -1;
    find_nearest_node_kernel(xs,coords,ndims,&ans,&best_distance);
    return ans;
}

void kdNode::find_nearest_node_kernel(const vector<double> &xs, const vector<vector<double>> &coords, int ndims, int* curr_node, double* best_distance){
    if (vid < 0){ // reached end of the tree
        return;
    }
    if (depth == 0){
        assert(xs.size() == ndims);
        *curr_node = vid;
        *best_distance = norm<double>(coords[*curr_node]-xs);
    } else {
        double d = norm<double>(coords[vid]-xs);
        if (d < *best_distance){
            *best_distance = d;
            *curr_node = vid;
        }
    }

    int axis = depth % ndims;
    double axis_vec = xs[axis] - coords[vid][axis];
    if (xs[axis] < coords[vid][axis]){
        left->find_nearest_node_kernel(xs,coords,ndims,curr_node,best_distance);
        if (xs[axis] + *best_distance >= coords[vid][axis]){
            right->find_nearest_node_kernel(xs,coords,ndims,curr_node,best_distance);
        }
    } else {
        right->find_nearest_node_kernel(xs,coords,ndims,curr_node,best_distance);
        if (xs[axis] - *best_distance <= coords[vid][axis]){
            left->find_nearest_node_kernel(xs,coords,ndims,curr_node,best_distance);
        }
    }
}

void kdNode::printnode(){
    cout << "depth: " << depth << " node: " << vid << endl;
    return;
}

// sorting algorithms for parallel tests
int quicksort_partition(int *x, int lo, int hi){
    int pivot = x[hi];

    int i = lo - 1;
    int temp;
    for (int j = lo; j<hi; j++){
        if (x[j] <= pivot){
            i = i+1;
            temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }

    i=i+1;
    temp = x[i];
    x[i] = x[hi];
    x[hi] = temp;
    return i;
}


void quicksort_kernel(int *x, int lo, int hi){
    if (lo >= hi || lo < 0){
        return;
    }

    int p = quicksort_partition(x,lo,hi);
    quicksort_kernel(x,lo,p-1);
    quicksort_kernel(x,p+1,hi);
    return;
}

void quicksort(int *x, int sz){
    quicksort_kernel(x,0,sz-1);
}

void mergesort(int *x, int sz1, int *y, int sz2, int* xy){
    int i=0,j=0,k=0;
    while (k < (sz1+sz2)){
        if (i >= sz1){
            xy[k] = y[j];
            j++;
        } else if(j >= sz2){
            xy[k] = x[i];
            i++;
        } else {
            if (x[i]<=y[j]){
                xy[k] = x[i];
                i++;
            } else {
                xy[k] = y[j];
                j++;
            }
        }
        k++;
    }
    return;
}