#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat join_diag(mat A, mat B){
    return join_vert( join_horiz(A, zeros(A.n_rows, B.n_cols)), join_horiz(zeros(B.n_rows, A.n_cols), B) );
}

vec join_diag(mat **lista, uint n){
    cout << -1 << endl;
    mat aux = *lista[0];
    cout << 0 << endl;
    cout << aux << endl;
    for(uint i = 1; i<n; i++){
        aux = join_diag(aux, *lista[i]);
        cout << i << endl;
        cout << aux << endl;
    }
    return *lista[0];
}

vec join_vert(field<vec> F){
    vec aux = F(0);
    for(uint i = 1; i< F.n_rows; i++)
        aux = join_vert(aux, F(i));
    return aux;
}

vec join_vert(vec *lista, uint n){
    vec aux = lista[0];
    for(uint i = 1; i<n; i++)
        aux = join_vert(aux, lista[i]);
    return aux;
}

vec join_vert(vec **lista, uint n){
    vec aux = *lista[0];
    for(uint i = 1; i<n; i++)
        aux = join_vert(aux, *lista[i]);
    return aux;
}

int main(void){
	cout << "Hello! This is a C++ program." << endl;
	vec A = randn(4);
    vec B = randn(2);
    vec B2 = randn(3);
    
    field<vec> F(3);
    F(0) = A;
    F(1) = B;
    F(2) = B2;

    
    //F.print("F:");
    cout << F << endl;
    cout << F(0) << endl;
    cout << F(1) << endl;
    mat C = randn(4,4);
    mat D = randn(2,2);
    cout << C << endl;
    cout << D << endl;
    cout << join_diag( join_diag(C,D), C ) << endl;

    vec *lista;
    lista = new vec[3];
    lista[0] = {1,2};
    lista[1] = {3,4,6};
    lista[2] = {7,8};
    cout << lista[0] << endl;
    cout << lista[1] << endl;
    cout << lista[2] << endl;

    vec lista2[] = {lista[0], lista[1], lista[2]};

    cout << join_vert(F) << endl;
    cout << join_vert(lista,3) << endl;
    cout << join_vert(lista2,3) << endl;

    vec **lista3 = new vec* [3];
    lista3[0] = &A;
    lista3[1] = &B;
    lista3[2] = &B2;

    mat **lista4 = new mat* [3];
    lista4[0] = &C;
    lista4[1] = &D;
    lista4[2] = &C;

    cout << *lista3[0] << endl;
    cout << join_vert(*lista3[0], *lista3[1]) << endl;
    cout << join_vert(lista3, 3) << endl;

    //mat aux = *lista4[0];
    //cout << aux << endl;
    //for(uint i = 1; i<3; i++){
    //    aux = join_diag(aux, *lista4[i]);
    //    cout << aux << endl;
    //}
    //cout << aux << endl;

    mat M = join_diag(lista4, 2);


    return 0;
}