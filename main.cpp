#include <iostream>
#include <cmath>

using namespace std;

double *mult (int n, int m, double ** X, double * Y){
    double* Z = new double [n];
    for (int i = 0; i < n; ++i)
    {
        Z[i] = 0;
    }
    for (int i = 0; i < n; ++i)
    {
        for (int k = 0; k < m; ++k)
        {
            Z[i] += (X[i][k]) * (Y[k]);
        }
    }
    return Z;
};
void swapRows (int n, double **A, int a, int b){
    for (int i = 0; i < n; ++i)
    {
        swap (A[a][i], A[b][i]);
    }
};

double** LUP (int n, double ** A, double ** P){
    double** W = new double* [n];
    for (int i = 0; i < n; ++ i)
        W[i] = new double [n];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            W[i][j] = A[i][j];
            P[i][j] = 0;
        }
        P[i][i] = 1;
    }
    for (int i = 0; i < n; ++i)
    {
        //Ищем опорный элемент
        double sup = 0;
        int jsup = -1;
        for (int j = i; j < n; ++j)
        {
            if (abs(W[j][i]) > sup)
            {
                sup = abs(W[j][i]);
                jsup = j;
            }
        }
        if (sup != 0)
        {
            //Меняем местами i-ю строку и строку с опорным элементом
            swapRows(n, P, jsup, i);
            swapRows(n, W, jsup, i);

            for (int j = i + 1; j < n; ++j) {
                W[j][i] /= W[i][i];
                for (int k = i + 1; k < n; k++) {
                    W[j][k] -= W[j][i] * W[i][k];
                }
            }
        }
    }
    return W;
};
double* SLAU (int n, double ** LU, double * b, double ** P){
    double* X = new double [n];
    for (int i = 0; i < n; ++i)
    {
        X[i] = 0;
    }

    b = mult(n, n, P, b);

    double* Y = new double [n];
    for (int i = 0; i < n; ++i)
    {
        int k = i;
        Y[i] = b[i];
        while (k > 0)
        {
            Y[i] -= Y[k-1]*LU[i][k-1];
            --k;
        }
    }
    for (int i = n-1; i >= 0; --i)
    {
        int k = i;
        X[i] = Y[i];
        while (k < n-1)
        {
            X[i] -= X[k+1] * LU[i][k+1];
            ++k;
        }
        X[i] /= LU[i][i];
    }
    return X;
};

double f (double x) {
    return 2*cos(3.5*x)*exp(5*x/3) + 3*sin(1.5*x)*exp(-4*x)+3;
};
double* moments_n (double a, double b) {
    double* mu = new double [3];
    mu[0] = 5*(pow((b - 1.5), 0.8) - pow((a - 1.5), 0.8))/4;
    mu[1] = 5*(pow((b - 1.5), 0.8)*(8*b + 15) - pow((a - 1.5), 0.8)*(8*a + 15))/72;
    mu[2] = 5*(pow((b - 1.5), 0.8)*(24*b*b + 40*b + 75) - pow((a - 1.5), 0.8)*(24*a*a + 40*a + 75))/336;
    return mu;
};
double* moments_2n (double a, double b) {
    double* mu = new double [6];
    mu[0] = 5*(pow((b - 1.5), 0.8) - pow((a - 1.5), 0.8))/4;
    mu[1] = 5*(pow((b - 1.5), 0.8)*(8*b + 15) - pow((a - 1.5), 0.8)*(8*a + 15))/72;
    mu[2] = 5*(pow((b - 1.5), 0.8)*(24*b*b + 40*b + 75) - pow((a - 1.5), 0.8)*(24*a*a + 40*a + 75))/336;
    mu[3] = 5*(pow((b - 1.5), 0.8)*(224*b*b*b + 360*b*b + 600*b + 1125) - pow((a - 1.5), 0.8)*(224*a*a*a + 360*a*a + 600*a + 1125))/4256;
    mu[4] = 5*(pow((b - 1.5), 0.8)*(2128*b*b*b*b + 3360*b*b*b + 5400*b*b + 9000*b + 16875) - pow((a - 1.5), 0.8)*(2128*a*a*a*a + 3360*a*a*a + 5400*a*a + 9000*a + 16875))/51072;
    mu[5] = 5*(pow((b - 1.5), 0.8)*(34048*b*b*b*b*b + 53200*b*b*b*b + 84000*b*b*b + 135000*b*b + 225000*b + 421875) - pow((a - 1.5), 0.8)*(34048*a*a*a*a*a + 53200*a*a*a*a + 84000*a*a*a + 135000*a*a + 225000*a + 421875))/987392;
    return mu;
};

double IQF_NC (double a, double b){
    double ab = (a+b)/2;
    double* mu = new double [3]; // моменты весовой функции
    mu = moments_n(a, b);

    double** X = new double* [3];
    for (int i = 0; i < 3; ++i)
        X[i] = new double [3];
    X[0][0] = 1;
    X[0][1] = 1;
    X[0][2] = 1;
    X[1][0] = a;
    X[1][1] = ab;
    X[1][2] = b;
    X[2][0] = a*a;
    X[2][1] = ab*ab;
    X[2][2] = b*b;
    double* A = new double [3];
    double** P = new double* [3];
    for (int i = 0; i < 3; ++i)
        P[i] = new double [3];
    double** LU = new double* [3];
    for (int i = 0; i < 3; ++i)
        LU[i] = new double [3];
    LU = LUP(3, X, P);
    A = SLAU(3, LU, mu, P);
    double Sn = A[0]*f(a) + A[1]*f(ab) + A[2]*f(b);
    return Sn;
};
void QF_NC (double a, double b){
    double EXACT_VALUE = 32.21951452884234295708696008290380201405;
    double Sh1, Sh2 = 0, Sh3 = 0, Rn;
    int k, m = 0;
    double h = (b-a);
    int steps = 0;
    do
    {
        ++steps;
        Sh1 = Sh2;
        Sh2 = Sh3;
        Sh3 = 0;
        k = (b-a) / h;
        double* z = new double[k+1];

        for (int i = 0; i < k+1; ++i) // делим отрезок на k частей точками z[i]
        {
            z[i] = a + h*i;
        }
        for (int i = 0; i < k; ++i)
        {
            Sh3 += IQF_NC(z[i], z[i+1]);
        }
        m = -log((Sh3 - Sh2)/(Sh2 - Sh1))/log(2); // вычисляем m по Эйткена
        Rn = abs(EXACT_VALUE - (Sh2 + (Sh2 - Sh1)/(pow(2, m) - 1))); // вычисляем погрешность по Ричардсону
        h /= 2;
    } while(Rn > 0.00001);
    cout << "Steps: " << steps << endl;
    cout << "Sn = " << Sh3 << endl;
};
void QF_NC_hopt (double a, double b){
    double EXACT_VALUE = 32.21951452884234295708696008290380201405;
    double eps = 0.00001;
    double Sh1, Sh2 = 0, Sh3 = 0, Rn = 1;
    int k, m = 0;
    double h = (b-a), hopt;
    int steps = 0;

    do {
        ++steps;
        Sh1 = Sh2;
        Sh2 = Sh3;
        Sh3 = 0;
        k = (b-a) / h; // k - количество шагов

        double *z = new double[k + 1];
        for (int i = 0; i < k + 1; ++i) // делим отрезок на k частей точками z[i]
        {
            z[i] = a + h * i;
        }
        for (int i = 0; i < k; ++i) {
            Sh3 += IQF_NC(z[i], z[i + 1]);
        }

        if (Sh1 != 0)
        {
            m = trunc(-log((Sh3 - Sh2) / (Sh2 - Sh1)) / log(2)); // вычисляем m по Эйткенy
            Rn = abs(EXACT_VALUE - (Sh2 + (Sh2 - Sh1) / (pow(2, m) - 1))); // вычисляем погрешность по Ричардсону
        }
        if (k == 4)
        {
            hopt = 0.95 * h * pow((eps * (1 - pow(2, (-m))) / abs(Sh3 - Sh2)), (1. / m));
            h = (b-a)/ceil((b-a)/hopt);
        }
        else
            h /= 2;
    } while(Rn > eps);
    cout << "Steps: " << steps << endl;
    cout << "Sn = " << Sh3 << endl;
};

double* FCardano (double* A){
    double* x = new double[3];
    double p = (3*A[1] - A[2]*A[2])/3;
    double q = (2*A[2]*A[2]*A[2] - 9*A[2]*A[1] + 27*A[0])/27;
    double Q = p*p*p/27 + q*q/4;
//    if (Q < 0)
//    {
    double fi = acos(-q*pow((-3/p),1.5)/2);
    x[0] = 2*pow((-p/3),0.5)*cos(fi/3) - A[2]/3;
    x[1] = 2*pow((-p/3),0.5)*cos(fi/3 + 2*acos(-1.0)/3) - A[2]/3;
    x[2] = 2*pow((-p/3),0.5)*cos(fi/3 - 2*acos(-1.0)/3) - A[2]/3;
//    }
    return x;
};

double IQF_G (double a, double b){
    double ab = (a+b)/2;
    double* mu = new double [6]; // моменты весовой функции
    mu = moments_2n(a, b);

    double** MU = new double* [3];
    for (int i = 0; i < 3; ++i)
        MU[i] = new double [3];
    MU[0][0] = mu[0];
    MU[0][1] = mu[1];
    MU[0][2] = mu[2];
    MU[1][0] = mu[1];
    MU[1][1] = mu[2];
    MU[1][2] = mu[3];
    MU[2][0] = mu[2];
    MU[2][1] = mu[3];
    MU[2][2] = mu[4];
    double* a_ = new double [3];
    double** P = new double* [3];
    for (int i = 0; i < 3; ++i)
        P[i] = new double [3];
    double** LU = new double* [3];
    for (int i = 0; i < 3; ++i)
        LU[i] = new double [3];
    LU = LUP(3, MU, P);
    double* _mu = new double [3];
    _mu[0] = -mu[3];
    _mu[1] = -mu[4];
    _mu[2] = -mu[5];
    a_ = SLAU(3, LU, _mu, P);

    double* x = new double [3];
    x = FCardano (a_);

    double* A = new double [3];
    double** X = new double* [3];
    for (int i = 0; i < 3; ++i)
        X[i] = new double [3];
    X[0][0] = 1;
    X[0][1] = 1;
    X[0][2] = 1;
    X[1][0] = x[0];
    X[1][1] = x[1];
    X[1][2] = x[2];
    X[2][0] = x[0]*x[0];
    X[2][1] = x[1]*x[1];
    X[2][2] = x[2]*x[2];
    LU = LUP(3, X, P);
    _mu[0] = mu[0];
    _mu[1] = mu[1];
    _mu[2] = mu[2];
    A = SLAU(3, LU, _mu, P);

    double Sn = A[0]*f(x[0]) + A[1]*f(x[1]) + A[2]*f(x[2]);
    return Sn;
};
void QF_G (double a, double b){
    double EXACT_VALUE = 32.21951452884234295708696008290380201405;
    double eps = 0.00001;
    double Sh1, Sh2 = 0, Sh3 = 0, Rn = 1;
    int k, m = 0;
    double h = (b-a);
    int steps = 0;

    do
    {
        ++steps;
        Sh1 = Sh2;
        Sh2 = Sh3;
        Sh3 = 0;
        k = (b-a) / h;
        double* z = new double[k+1];
        for (int i = 0; i < k+1; ++i) // делим отрезок на k частей точками z[i]
        {
            z[i] = a + h*i;
        }
        for (int i = 0; i < k; ++i)
        {
            Sh3 += IQF_G(z[i], z[i+1]);
        }
        if (Sh1 != 0)
        {
            m = trunc(-log((Sh3 - Sh2) / (Sh2 - Sh1)) / log(2)); // вычисляем m по Эйткена
            Rn = abs(EXACT_VALUE - (Sh2 + (Sh2 - Sh1) / (pow(2, m) - 1))); // вычисляем погрешность по Ричардсону
        }
        h /= 2;
    } while(Rn > eps);
    cout << "Steps: " << steps << endl;
    cout << "Sn = " << Sh3 << endl;
};
void QF_G_hopt (double a, double b){
    double EXACT_VALUE = 32.21951452884234295708696008290380201405;
    double eps = 0.00001;

    double Sh1 = 0, Sh2 = 0, Sh3 = 0, Rn = 1;
    int k, m = 0;
    double h = (b-a), hopt;
    int steps = 0;

    do
    {
        ++steps;
        Sh1 = Sh2;
        Sh2 = Sh3;
        Sh3 = 0;
        k = (b-a) / h;
        double* z = new double[k+1];
        for (int i = 0; i < k+1; ++i) // делим отрезок на k частей точками z[i]
        {
            z[i] = a + h*i;
        }
        for (int i = 0; i < k; ++i)
        {
            Sh3 += IQF_G(z[i], z[i+1]);
        }
        if (Sh1 != 0)
        {
            m = trunc(-log((Sh3 - Sh2) / (Sh2 - Sh1)) / log(2)); // вычисляем m по Эйткена
            Rn = abs(EXACT_VALUE - (Sh2 + (Sh2 - Sh1) / (pow(2, m) - 1))); // вычисляем погрешность по Ричардсону
        }
        if (k == 4)
        {
            hopt = 0.95 * h * pow((eps * (1 - pow(2, (-m))) / abs(Sh3 - Sh2)), (1. / m));
            h = (b-a)/ceil((b-a)/hopt);
        }
        else
            h /= 2;
    } while(Rn > eps);
    cout << "Steps: " << steps << endl;
    cout << "Sn = " << Sh3 << endl;
};


int main() {
    double EXACT_VALUE = 32.21951452884234295708696008290380201405;
    double a = 1.5;
    double b = 2.3;

    cout.precision(9);
    cout << "Exact value: " << EXACT_VALUE << endl;
    cout << "EPS = " << 0.00001 << endl << endl;

    double Sn;

    cout << "_____NEWTON-COTES VARIANT_____" << endl;
    cout << "1.IQF" << endl;
    Sn = IQF_NC(a, b);
    cout << "Sn = " << Sn << endl;
    cout << "Methodical error: " << 0.01638766269011003 * 2275.1 / 6 << endl; // Mn = 2275.1; Rn = 2275.1/6 * int from {1.5} to {2.3} {|(x-1.5)^(4/5)*(x-1.9)*(x-2.3)|}
    cout << "Exact error: " << abs(EXACT_VALUE - Sn) << endl;
    cout << endl;

    cout << "2.CQF" << endl;
    QF_NC(a, b);
    cout << endl;

    cout << "3.Hopt" << endl;
    QF_NC_hopt (a,b);
    cout << endl;


    cout << "_________GAUSS VARIANT________" << endl;
    cout << "1.IQF" << endl;
    Sn = IQF_G(a, b);
    cout << "Sn = " << Sn << endl;
    cout << "Methodical error: " << 0.01638766269011003 * 2275.1 / 6 << endl; // Mn = 2275.1; Rn = 2275.1/6 * int from {1.5} to {2.3} {|(x-1.5)^(4/5)*(x-1.9)*(x-2.3)|}
    cout << "Exact error: " << abs(EXACT_VALUE - Sn) << endl;
    cout << endl;

    cout << "2.CQF" << endl;
    QF_G(a, b);
    cout << endl;

    cout << "3.Hopt" << endl;
    QF_G_hopt (a,b);
    cout << endl;

    return 0;
}