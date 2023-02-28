#include <vector>
#include <iostream>
#include <random>
#include <cassert>
using namespace std;

#ifndef LINALG_H
#define LINALG_H

/// Matrix operators

// matrix matrix addition
template <typename T>
vector<vector<T>> operator+(const vector<vector<T>> &A, const vector<vector<T>> &B){
    int m1 = A.size();
    int m2 = B.size();
    int n1 = A[0].size();
    int n2 = B[0].size();

    assert(n1==n1 && m1==m2);
    vector<vector<T>> C;
    C.assign(m1, vector<T>(n1, 0));
    for (int i = 0; i<m1; i++){
        for (int j = 0; j<n1; j++){
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

// matrix matrix subtraction
template <typename T>
vector<vector<T>> operator-(const vector<vector<T>> &A, const vector<vector<T>> &B){
    int m1 = A.size();
    int m2 = B.size();
    int n1 = A[0].size();
    int n2 = B[0].size();

    assert(n1==n1 && m1==m2);
    vector<vector<T>> C;
    C.assign(m1, vector<T>(n1, 0));
    for (int i = 0; i<m1; i++){
        for (int j = 0; j<n1; j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

// matrix addition
template <typename T>
void operator +=(vector<vector<T>> &lhs, const vector<vector<T>> &rhs){
    int m1 = lhs.size();
    int m2 = rhs.size();
    int n1 = lhs[0].size();
    int n2 = rhs[0].size();

    assert(n1==n1 && m1==m2);
    
    for (int i = 0; i<m1; i++){
        for (int j = 0; j<n1; j++){
            lhs[i][j] += rhs[i][j];
        }
    }
    return;
} 

// matrix subtraction
template <typename T>
void operator -=(vector<vector<T>> &lhs, const vector<vector<T>> &rhs){
    int m1 = lhs.size();
    int m2 = rhs.size();
    int n1 = lhs[0].size();
    int n2 = rhs[0].size();

    assert(n1==n1 && m1==m2);
    
    for (int i = 0; i<m1; i++){
        for (int j = 0; j<n1; j++){
            lhs[i][j] -= rhs[i][j];
        }
    }
    return;
} 

// matrix matrix multiplication
template <typename T>
vector<vector<T>> operator*(const vector<vector<T>> &A, const vector<vector<T>> &B){
    int i,j,k;
    int m1 = A.size();
    int m2 = B.size();
    int n1 = A[0].size();
    int n2 = B[0].size();

    assert(n1==m2);
    vector<vector<T>> prod;
    prod.assign(m1, vector<T>(n1, 0));
    for (i = 0; i<m1;i++){
        prod[i].resize(n2);
        for (j = 0; j<n2; j++){
            prod[i][j] = 0.0;
            for (k = 0; k<n1; k++){
                prod[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return prod;
}

// matrix vector multiplication
template <typename T>
vector<T> operator*(const vector<vector<T>> &A, const vector<T> &b){
    int m = A.size();
    int n = A[0].size();
    int x = b.size();

    vector<T> prod(m);
    assert(x==n);
    for(int i=0; i<m; i++){
        prod[i] = 0;
        for(int j=0; j<n;j++){
            prod[i] += A[i][j]*b[j];
        }
    }
    return prod;
}

// matrix scalar multiplication
template <typename T>
vector<vector<T>> operator*(const vector<vector<T>> &A, T a){
    int m = A.size();
    int n = A[0].size();

    vector<vector<T>> prod;
    prod.assign(m, vector<T>(n, 0));
    for(int i=0; i<m; i++){
        for(int j=0; j<n;j++){
            prod[i][j] = A[i][j]*a;
        }
    }
    return prod;
}

/// Vector operators

// vector scalar multiplication
template <typename T>
vector<T>  operator*(const vector<T>  &u, T a){
    int n = u.size();
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = u[i]*a;
    }
    return prod;
}
template <typename T>
vector<T>  operator*(T a, const vector<T>  &u){
    int n = u.size();
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = a*u[i];
    }
    return prod;
}
template <typename T>
void operator*=(vector<T>  &u, T a){
    int n = u.size();
    for (int i = 0; i<u.size(); i++){
        u[i] = u[i]*a;
    }
    return;
}

// vector scalar division
template <typename T>
vector<T>  operator/(const vector<T>  &u, T a){
    int n = u.size();
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = u[i]/a;
    }
    return prod;
}
template <typename T>
void operator/=(vector<T>  &u, T a){
    int n = u.size();
    for (int i = 0; i<u.size(); i++){
        u[i] = u[i]/a;
    }
    return;
}

// vector scalar addition
template <typename T>
vector<T>  operator+(const vector<T>  &u, T a){
    int n = u.size();
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = u[i]+a;
    }
    return prod;
}
template <typename T>
void operator+=(vector<T>  &u, T a){
    int n = u.size();
    for (int i = 0; i<u.size(); i++){
        u[i] = u[i]+a;
    }
    return;
}

// vector scalar subtraction
template <typename T>
vector<T>  operator-(const vector<T>  &u, T a){
    int n = u.size();
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = u[i]+a;
    }
    return prod;
}
template <typename T>
void operator-=(vector<T>  &u, T a){
    int n = u.size();
    for (int i = 0; i<u.size(); i++){
        u[i] = u[i]+a;
    }
    return;
}


// vector vector multiplication
template <typename T>
vector<T>  operator*(const vector<T>  &u, const vector<T>  &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = u[i]*v[i];
    }
    return prod;
}
template <typename T>
void operator*=(vector<T>  &u, const vector<T>  &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    for (int i = 0; i<u.size(); i++){
        u[i] = u[i]*v[i];
    }
    return;
}

// vector vector division
template <typename T>
vector<T>  operator/(const vector<T>  &u, const vector<T>  &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = u[i]/v[i];
    }
    return prod;
}
template <typename T>
void operator/=(vector<T>  &u, const vector<T>  &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    for (int i = 0; i<u.size(); i++){
        u[i] = u[i]/v[i];
    }
    return;
}

// vector vector addition
template <typename T>
vector<T>  operator+(const vector<T>  &u, const vector<T>  &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = u[i]+v[i];
    }
    return prod;
}
template <typename T>
void operator+=(vector<T>  &u, const vector<T>  &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    for (int i = 0; i<u.size(); i++){
        u[i] = u[i]+v[i];
    }
    return;
}

// vector vector subtraction
template <typename T>
vector<T>  operator-(const vector<T>  &u, const vector<T>  &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    vector<T> prod(n);
    for (int i = 0; i<u.size(); i++){
        prod[i] = u[i]-v[i];
    }
    return prod;
}
template <typename T>
void operator-=(vector<T>  &u, const vector<T>  &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    for (int i = 0; i<u.size(); i++){
        u[i] = u[i]-v[i];
    }
    return;
}


// matrix and vector initializers
template <typename T>
vector<T> randvec(int n, T lower=0, T upper=1){
    default_random_engine re(time(0));
    vector<T> rand(n);
    uniform_real_distribution<double> unif(lower, upper);
    for(int i = 0; i < n; i++){
        rand[i] = unif(re);
    }
    return rand;
}

template <typename T> 
vector<vector<T>> Zeros(int m, int n){
    vector<vector<T>> A;
    A.assign(m, vector<T>(n, 0));
    return A;
}

template <typename T>
vector<vector<T>> randMatrix(int m, int n, T lower=0, T upper=1){
    default_random_engine re(time(0));
    vector<vector<T>> rand = Zeros<T>(m,n);
    uniform_real_distribution<double> unif(lower, upper);
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            rand[i][j] = unif(re);
        }   
    }
    return rand;
}

// Matrix and vector utilities
template <typename T> 
void printvec(const vector<T> &A){
    for (int i = 0; i<A.size(); i++){
        cout << A[i] << endl;
    }
    cout << endl;
}

template <typename T> 
void printMat(const vector<vector<T>> &A){
    cout << endl;
    for (int i = 0; i<A.size(); i++){
        for (int j = 0; j<A[i].size(); j++){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

template <typename T>
T inner(const vector<T> &u, const vector<T> &v){
    int n = u.size();
    int n2 = v.size();
    assert(n==n2);
    T val = u[0]*v[0];
    for(int i = 1; i<n; i++){
        val += u[i]*v[i];
    }
    return val;
}

template <typename T>
T norm(const vector<T> &u){
    return sqrt(inner(u,u));
}

template <typename T>
vector<vector<T>> Transpose(vector<vector<T>> &A){
    int m = A.size();
    int n = A[0].size();
    double val;

   
    vector<vector<T>> B = Zeros<T>(n,m);
    for (int i = 0; i<m; i++){
        for (int j = 0; j<n; j++){
            B[j][i] = A[i][j];
        }
    }
    return B;
}

#endif