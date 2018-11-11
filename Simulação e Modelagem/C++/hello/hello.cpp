#include <iostream>
#include <armadillo>
#define PI 3.1415926535897932384626433832795028841971693993751

using namespace std;
using namespace arma;

class Filter {
public:
    Filter(int size, mat AB_);
    ~Filter();
    vec Doit(vec u_);

    int order;
    int size;
    mat u__;
    mat y__;
    vec a_;
    vec b_;
};

Filter::Filter(int size, mat AB_){
    this->order = AB_.n_cols-1;
    this->size = size;
    this->a_ = (AB_.row(0)/AB_(0,0)).t();
    this->b_ = (AB_.row(1)/AB_(0,0)).t();
    u__.zeros(size,order+1);
    y__.zeros(size,order+1);
}

Filter::~Filter(){
    a_.clear();
    b_.clear();
    u__.clear();
    y__.clear();
}

vec Filter::Doit(vec u_){
    u__.col(0) = u_;
    y__.col(0) = b_(0)*u__.col(0);
    for(int i = 1; i <= order; i++)
        y__.col(0)  += b_(i)*u__.col(i) - a_(i)*y__.col(i);

    
    for(int i = order; i > 0; i--){
        u__.col(i) = u__.col(i-1);
        y__.col(i) = y__.col(i-1);
    }

    //u__ = shitf(u__,+1);
    //y__ = shitf(u__,+1);

    return y__.col(0);
}

mat Tustin(double T, double w0, mat ABs_){
    int order = ABs_.n_cols-1;
    mat AB_ = ABs_;
    double c;
    if (w0 == 0.0) c = 2/T;
    else c = w0/tan(0.5*w0*T);

    switch(order){
        case 1:
            AB_.col(0) = ABs_.col(0) + c*ABs_.col(1);
            AB_.col(1) = ABs_.col(0) - c*ABs_.col(1);
            break;
        case 2:
            AB_.col(0) =   ABs_.col(0) + c*ABs_.col(1) + c*c*ABs_.col(2);
            AB_.col(1) = 2*ABs_.col(0)               - 2*c*c*ABs_.col(2);
            AB_.col(2) =   ABs_.col(0) - c*ABs_.col(1) + c*c*ABs_.col(2);
            break;
        case 3:
            AB_.col(0) =   ABs_.col(0) + c*ABs_.col(1) + c*c*ABs_.col(2) +   c*c*c*ABs_.col(3);
            AB_.col(1) = 3*ABs_.col(0) + c*ABs_.col(1) - c*c*ABs_.col(2) - 3*c*c*c*ABs_.col(3);
            AB_.col(2) = 3*ABs_.col(0) - c*ABs_.col(1) - c*c*ABs_.col(2) + 3*c*c*c*ABs_.col(3);
            AB_.col(3) =   ABs_.col(0) - c*ABs_.col(1) + c*c*ABs_.col(2) -   c*c*c*ABs_.col(3);
            break;
    }


    return AB_;
}



double interpolacao(double x, double x0, double xf, vec v)
{
    int N = v.n_rows;
    if(x < x0) return v(0);
    else if(x > xf) return v(N-1);
    else{
        double dx = (xf-x0)/(N-1);
        double index = (x - x0)/dx;
        int index_esq = floor(index);
        int index_dir = ceil(index);
        double y0 = v(index_esq);
        double yf = v(index_dir);
        double y = y0 + (yf - y0)*(index - index_esq);
        return y;   
    }
}

double interpolacao_2D(double x, double y,  double x0, double y0, double xf, double yf, mat M)
{
    if( x >= x0 && x <= xf && yf >= y0 && y <= yf ){
        int Nx = M.n_rows;
        int Ny = M.n_cols;
        double dx = (xf-x0)/(Nx-1);
        double dy = (yf-y0)/(Ny-1);
        double index_x = (x - x0)/dx;
        int index_x1 = floor(index_x);
        int index_x2 = ceil(index_x);
        double index_y = (y - y0)/dy;
        int index_y1 = floor(index_y);
        int index_y2 = ceil(index_y);
        double z11 = M(index_x1, index_y1);
        double z12 = M(index_x1, index_y2);
        double z21 = M(index_x2, index_y1);
        double z22 = M(index_x2, index_y2);
        double z1 = z11 + (z21 - z11)*(index_x - index_x1);
        double z2 = z12 + (z22 - z12)*(index_x - index_x1);
        double z = z1 + (z2 - z1)*(index_y - index_y1);
        return z;
    }
    else return 0; 
}

mat interpolacao_2D_mat(double x, double y,  double x0, double y0, double xf, double yf, field<mat> M)
{
    if( x >= x0 && x <= xf && yf >= y0 && y <= yf ){
        int Nx = M.n_rows;
        int Ny = M.n_cols;
        double dx = (xf-x0)/(Nx-1);
        double dy = (yf-y0)/(Ny-1);
        double index_x = (x - x0)/dx;
        int index_x1 = floor(index_x);
        int index_x2 = ceil(index_x);
        double index_y = (y - y0)/dy;
        int index_y1 = floor(index_y);
        int index_y2 = ceil(index_y);
        mat z11 = M(index_x1, index_y1);
        mat z12 = M(index_x1, index_y2);
        mat z21 = M(index_x2, index_y1);
        mat z22 = M(index_x2, index_y2);
        mat z1 = z11 + (z21 - z11)*(index_x - index_x1);
        mat z2 = z12 + (z22 - z12)*(index_x - index_x1);
        mat z = z1 + (z2 - z1)*(index_y - index_y1);
        return z;
    }
    else return 0;
}

mat CrossMat(vec x)
{
    mat M;
    M << 0     << -x(2) <<  x(1) << endr
      << x(2)  << 0     << -x(0) << endr
      << -x(1) << x(0)  <<  0 << endr;
    return M;
}

int main(void){
	//cout << "Hello! This is a C++ program." << endl;

    /*
    mat A; A.zeros(2,2);
    A << 1 << 2 << endr
      << 3 << 4 << endr;
    vec b = {1, 2};
    /*

    //vec x = solve(A,b);
    /*A.print("A = ");
    b.print("b = ");
    x.print("x = ");
    cout << endl;*/

    /*
    int n = 20;
    vec v(n, fill::zeros);
    for(int i = 0; i<n; i++) v(i) = sin(2*PI*i/(n-1));
    cout << v << endl;

    double x0 = 0;
    double xf = 2*PI;
    double x = 3*2*PI/(n-1)+0.01;

    cout << interpolacao(x, x0, xf, v) << endl;
    cout << sin(x) << endl;
    */

    /*
    int n = 400;
    int nx = n;
    int ny = n;
    double x0 = -5.0;
    double y0 = -6.0;
    double xf = 4;
    double yf = 3;
    double dx = (xf-x0)/(nx-1);
    double dy = (yf-y0)/(ny-1);
    double x, y;

    field<mat> M(nx, ny);
    for(int i = 0; i<nx; i++){
        x = x0 + i*dx;
        for(int j = 0; j<ny; j++){
            y = y0 + j*dy;
            M(i,j) = zeros(2,2);
            M(i,j)(0,0) = sin( x*y ) ;
            M(i,j)(0,1) = cos( x*y ) ;
            M(i,j)(1,0) = sin( x*x ) ;
            M(i,j)(1,1) = cos( y*y ) ;
        }
    }
    //cout << M << endl;
    //M.save("mat_field");



    double xis = -1;
    double ypi = -2;
    cout << interpolacao_2D_mat(xis,ypi, x0, y0, xf, yf, M) << endl;
    //cout << interpolacao(xis, x0, xf, M.col(0)) << endl;
    //cout << sin(xis*ypi) << endl;
    */

    /*
    //Plotar área de trabalho
    uint nx = 96.0;
    uint ny = 56.0;
    double lx = 0.24;
    double ly = 0.28;
    double xi = -lx;
    double xf = lx;
    double yi = 0.0;
    double yf = ly;
    field<mat> fZ_(ny,nx); fZ_.load("fZ_field");
    field<mat> fMh_(ny,nx); fMh_.load("fMh_field");
    //field<vec> fgh_(ny,nx); fgh_.load("fgh_field");
    field<vec> fa1_(ny,nx); fa1_.load("fa1_field");
    field<vec> fa2_(ny,nx); fa2_.load("fa2_field");
    field<vec> fa12_(ny,nx); fa12_.load("fa12_field");

    fMh_(20,40).print("Mh_ = ");
    fZ_(20,40).t().print("Z^T = ");
    //cout << fMh_.n_rows << endl;
    //cout << fMh_.n_cols << endl;
    cout << endl << "Mh_ = " << endl <<interpolacao_2D_mat(-0.05, 0.20, xi, yi, xf, yf, fMh_) << endl;
    */

    /*vec vetor = {1,2,3};
    vec vetor2 = {4,5,6};
    cout << CrossMat(vetor) << endl;
    cout << CrossMat(vetor)*vetor2 << endl;
    cout << cross(vetor,vetor2)  << endl;
    */

    //y[k] = -a1*y[k-1] + b0*u[k] + b1*u[k-1]

    int n = 30;
    int size = 2;
    mat u__ = ones(2,n);
    mat y__ = zeros(2,n);
    vec as_ = {15,26.3351,18.4943,5.41166};
    vec bs_ = {15,0,0,0};
    mat ABs_ = join_horiz(as_,bs_).t();
    mat AB_ = Tustin(0.5,0.0,ABs_);

    cout << AB_ << endl;

    Filter *F1; F1 = new Filter(size,AB_);
    for(int i = 0; i < n; i++){
        y__.col(i) = F1->Doit(u__.col(i));
    }
    cout << y__.t() << endl;
    cout << F1->u__ << endl;


    return 0;
}