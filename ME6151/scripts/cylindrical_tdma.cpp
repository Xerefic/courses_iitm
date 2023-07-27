#include <bits/stdc++.h>
using namespace std;

const double PI = 3.14159265358979323846;
const double g = 9.80665;

const int ITER = 1e3;
const double ERROR = 1e-6;

double max(double a, double b)
{
    return (a>b)?a:b;
}
double min(double a, double b)
{
    return (a<b)?a:b;
}

double* initialize1d(int m, double value=0)
{
    double* x = new double[m];
    for (int i=0; i<m; i++)
    {
        x[i] = value;
    }
    return x;
}
double** initialize2d(int n, int m, double value=0)
{
    double** x = new double*[n];
    for (int i=0; i<n; i++)
    {
        x[i] = new double[m];
        for (int j=0; j<m; j++)
        {
            x[i][j] = value;
        }
    }
    return x;
}
double*** initialize3d(int l, int n, int m, double value=0)
{
    double*** x = new double**[l];
    for (int i=0; i<l; i++)
    {
        x[i] = new double*[n];
        for (int j=0; j<n; j++)
        {
            x[i][j] = new double[m];
            for (int k=0; k<m; k++)
            {
                x[i][j][k] = value;
            }
        }
    }
    return x;
}

double* copy1d(double* x, int m)
{
    double* x_ = initialize1d(m);
    for (int i=0; i<m; i++)
    {
        x_[i] = x[i];
    }
    return x_;
}
double** copy2d(double** x, int n, int m)
{
    double** x_ = initialize2d(n, m);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            x_[i][j] = x[i][j];
        }
    }
    return x_;
}
double*** copy3d(double*** x, int l, int n, int m)
{
    double*** x_ = initialize3d(l, n, m);
    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                x_[i][j][k] = x[i][j][k];
            }
        }
    }
    return x_;
}
double*** steady3d(double*** x, int l, int n, int m)
{
    double*** x_ = initialize3d(l, n, m);
    double temp;
    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            temp = 0;
            for (int k=0; k<m; k++)
            {
                temp += x[i][j][k]/m;
            }
            for (int k=0; k<m; k++)
            {
                x_[i][j][k] = temp;
            }
        }
    }
    return x_;
}

void print1d(double* x, int m)
{
    for (int i=0; i<m; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
}
void print2d(double** x, int n, int m)
{
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            cout << x[i][j] << " ";
        }
        cout << endl;
    }
}
void print3d(double*** x, int l, int n, int m)
{
    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                cout << x[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}
void print_paraview2d(double** x, int n, int m, double r_in, double r_out)
{
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;
    double r, theta;

    ofstream file("paraview_2d.txt");
    file << "X, Y, Z, T" << endl;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {   
            r = r_in + (i+0.5)*del_r;
            theta = (j+0.5)*del_theta;
            file << r*cos(theta) << ", " << r*sin(theta) << ", " << 0 << ", " << x[i][j] << endl;
        }
    }
    file.close();
}
void print_paraview3d(double*** x, int l, int n, int m, double r_in, double r_out, double t)
{
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;

    double r, theta;
    double temp;

    ofstream file("paraview_3d.txt");
    file << "X, Y, Z, T" << endl;
    double in_r = r_in/2;
    for (int i=0; i<l; i++)
    {
        for (int j=0; j<2; j++)
        {
            for (int k=0; k<m; k++)
            {
                r = 0 + (j+0.5)*in_r;
                theta = (k+0.5)*del_theta;
                file << r*cos(theta) << ", " << r*sin(theta) << ", " << (i+0.5)*del_z << ", " << 0 << endl;
            }
        }
    }
    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {   
            for (int k=0; k<m; k++)
            {
                r = r_in + (j+0.5)*del_r;
                theta = (k+0.5)*del_theta;
                file << r*cos(theta) << ", " << r*sin(theta) << ", " << (i+0.5)*del_z << ", " << x[i][j][k] << endl;
            }
        }
    }

    file.close();
}

double*** read_paraview3d(int l, int n, int m, double r_in, double r_out, double t)
{
    double*** x = initialize3d(l, n, m);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;

    double r, theta;

    ifstream file("paraview_3d.txt");
    string line;
    getline(file, line);
    for (int i=0; i<l; i++)
    {
        for (int j=0; j<2; j++)
        {
            for (int k=0; k<m; k++)
            {
                getline(file, line);
            }
        }
    }
    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                getline(file, line);
                stringstream ss(line);
                string token;
                getline(ss, token, ',');
                getline(ss, token, ',');
                getline(ss, token, ',');
                getline(ss, token, ',');
                x[i][j][k] = stod(token);
            }
        }
    }

    file.close();
    return x;
}


double error1d(double* x, double* x_, int m)
{
    double error = 0;
    for (int i=0; i<m; i++)
    {
        error += (x[i]-x_[i])*(x[i]-x_[i]);
    }
    return sqrt(error);
}
double error2d(double** x, double** x_, int n, int m)
{
    double error = 0;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            error += (x[i][j]-x_[i][j])*(x[i][j]-x_[i][j]);
        }
    }
    return sqrt(error);
}
double error3d(double*** x, double*** x_, int l, int n, int m)
{
    double error = 0;
    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                error += (x[i][j][k]-x_[i][j][k])*(x[i][j][k]-x_[i][j][k]);
            }
        }
    }
    return sqrt(error);
}

// TP - Theta Positive; TN - Theta Negative
double* tdma(double* P, double* TP, double* TN, double* B, int m)
{
    double* x = initialize1d(m);

    for (int i=0; i<m; i++)
    {
        if (i==0)
        {
            B[i] = B[i]/P[i];
            TP[i] = TP[i]/P[i];
            P[i] = 1;
        }
        else if (i==m-1)
        {
            B[i] = (B[i] + TN[i]*B[i-1])/(P[i] - TN[i]*TP[i-1]);
            TN[i] = 0;
            P[i] = 1;
        }
        else
        {
            B[i] = (B[i] + TN[i]*B[i-1])/(P[i] - TN[i]*TP[i-1]);
            TP[i] = TP[i]/(P[i] - TN[i]*TP[i-1]);
            TN[i] = 0;
            P[i] = 1;
        }
    }

    for (int i=m-1; i>=0; i--)
    {
        if (i==m-1)
        {
            x[i] = B[i];
        }
        else
        {
            x[i] = B[i] + TP[i]*x[i+1];
        }
    }

    return x;
}
double* circular_tdma(double* x, double* P, double* TP, double* TN, double* B, int m)
{
    double* tp = initialize1d(m);
    double* tn = initialize1d(m);
    double* p = initialize1d(m);
    double* b = initialize1d(m);

    for (int i=0; i<m; i+=int(m/4))
    {
        int end = (i+m-1)%m;
        int start = i;
        int iter = 1;

        tp[0] = TP[start];
        tn[0] = TN[start];
        p[0] = P[start];
        b[0] = B[start] + TN[start]*x[end];
        for (int j=(start+1)%m; iter<m; j=(j+1)%m)
        {
            tp[iter] = TP[j];
            tn[iter] = TN[j];
            p[iter] = P[j];
            b[iter] = B[j];

            iter += 1;
        }

        x = tdma(p, tp, tn, b, m);
    }
    double eta = 0, epsilon = 1;
    double* y = initialize1d(m);
    double* u = initialize1d(m);
    double* v = initialize1d(m);
    double* q = initialize1d(m);



    // P[0] = P[0] - epsilon;
    // P[m-1] = P[m-1] - TP[m-1]*TN[0]/epsilon;
    // u[0] = epsilon;
    // u[m-1] = TP[m-1];
    // v[0] = 1;
    // v[m-1] = TN[0]/epsilon;

    // for (int i=1; i<m; i++)
    // {
    //     P[i] = P[i] - TP[i-1]*TN[i]/P[i-1];
    //     B[i] = B[i] - B[i-1]*TN[i]/P[i-1];
    //     TP[i] = TP[i] - u[i-1]*TN[i]/P[i-1];
    // }

    // y[m-1] = B[m-1]/P[m-1];
    // q[m-1] = u[m-1]/P[m-1];
 
    // for (int i=m-2; i>=0; i--)
    // {
    //     y[i] = (B[i]-TP[i]*y[i+1])/P[i];
    //     q[i] = (u[i]-TP[i]*q[i+1])/P[i];
    // }
    // eta = (v[0]*y[0] + v[m-1]*y[m-1])/(1 + v[0]*q[0] + v[m-1]*q[m-1]);
    // for(int i=0; i<m; i++)
    // {
	// 	x[i] = y[i] - q[i]*eta;
	// }


    // P[0] = P[0] - TP[m-1]*TN[0]/P[m-1];
    // TP[0] = TP[0] - TN[m-1]*TN[0]/P[m-1];
    // B[0] = B[0] - B[m-1]*TN[0]/P[m-1];

    // for (int i=1; i<m; i++)
    // {
    //     P[i] = P[i] - TP[i-1]*TN[i]/P[i-1];
    //     TP[i] = TP[i] - TP[i-1]*TN[i]/P[i-1];
    //     B[i] = B[i] - B[i-1]*TN[i]/P[i-1];
    // }

    // x[m-1] = B[m-1]/P[m-1];
    // for (int i=m-2; i>=0; i--)
    // {
    //     x[i] = max((B[i]-TP[i]*x[i+1])/P[i], 0);
    // }

    return x;
}
double* converge_circular_tdma(double* x, double* P, double* TP, double* TN, double* B, double alpha, int m)
{
    double* x_ = initialize1d(m);

    double error = 1;
    int iter = 0;
    do
    {   
        for (int i=0; i<m; i++)
        {
            P[i] = P[i]/alpha;
            B[i] = B[i] + (1-alpha)*P[i]*x[i];
        }
        x = circular_tdma(x, P, TP, TN, B, m);
        error = error1d(x, x_, m);
        x_ = copy1d(x, m);
        iter += 1;
    } while ((error>ERROR) && (iter<ITER));

    return x;
}

// TP - Theta Positive; TN - Theta Negative
// RP - R Positive; RN - R Negative
double** polar_tdma(double** x, double** P, double** TP, double** TN, double** RP, double** RN, double** B, double alpha, int n, int m)
{
    double* tp = initialize1d(m);
    double* tn = initialize1d(m);
    double* p = initialize1d(m);
    double* b = initialize1d(m);

    for (int i=0; i<n; i++)
    {
        if (i==0)
            {
                for (int j=0; j<m; j++)
                {
                    tp[j] = TP[i][j];
                    tn[j] = TN[i][j];
                    p[j] = P[i][j];
                    b[j] = B[i][j] + RP[i][j]*x[i+1][j];
                }
                x[i] = converge_circular_tdma(x[i], p, tp, tn, b, alpha, m);
            }
            else if (i==n-1)
            {
                for (int j=0; j<m; j++)
                {
                    tp[j] = TP[i][j];
                    tn[j] = TN[i][j];
                    p[j] = P[i][j];
                    b[j] = B[i][j] + RN[i][j]*x[i-1][j];
                }
                x[i] = converge_circular_tdma(x[i], p, tp, tn, b, alpha, m);
            }
            else
            {
                for (int j=0; j<m; j++)
                {
                    tp[j] = TP[i][j];
                    tn[j] = TN[i][j];
                    p[j] = P[i][j];
                    b[j] = B[i][j] + RN[i][j]*x[i-1][j] + RP[i][j]*x[i+1][j];
                }
                x[i] = converge_circular_tdma(x[i], p, tp, tn, b, alpha, m);
            }
    }

    return x;
}
double** converge_polar_tdma(double** x, double** P, double** TP, double** TN, double** RP, double** RN, double** B, double alpha, int n, int m)
{
    double** x_ = initialize2d(n, m);

    double error = 1;
    int iter = 0;
    do
    {
        x = polar_tdma(x, P, TP, TN, RP, RN, B, alpha, n, m);
        error = error2d(x, x_, n, m);
        x_ = copy2d(x, n, m);
        iter += 1;
    } while ((error>ERROR) && (iter<ITER));

    return x;
}

// TP - Theta Positive; TN - Theta Negative
// RP - R Positive; RN - R Negative
// ZP - Z Positive; ZN - Z Negative
double*** cylindrical_tdma(double*** x, double*** P, double*** TP, double*** TN, double*** RP, double*** RN, double*** ZP, double*** ZN, double*** B, double alpha, int l, int n, int m)
{
    double** tp = initialize2d(n, m);
    double** tn = initialize2d(n, m);
    double** rp = initialize2d(n, m);
    double** rn = initialize2d(n, m);
    double** p = initialize2d(n, m);
    double** b = initialize2d(n, m);

    for (int i=0; i<l; i++)
    {   
        if (i==0)
        {
            for(int j=0; j<n; j++)
            {
                for (int k=0; k<m; k++)
                {
                    tp[j][k] = TP[i][j][k];
                    tn[j][k] = TN[i][j][k];
                    rp[j][k] = RP[i][j][k];
                    rn[j][k] = RN[i][j][k];
                    p[j][k] = P[i][j][k];
                    b[j][k] = B[i][j][k] + ZP[i][j][k]*x[i+1][j][k];
                }
            }
            x[i] = converge_polar_tdma(x[i], p, tp, tn, rp, rn, b, alpha, n, m);
        }
        else if (i==l-1)
        {
            for(int j=0; j<n; j++)
            {
                for (int k=0; k<m; k++)
                {
                    tp[j][k] = TP[i][j][k];
                    tn[j][k] = TN[i][j][k];
                    rp[j][k] = RP[i][j][k];
                    rn[j][k] = RN[i][j][k];
                    p[j][k] = P[i][j][k];
                    b[j][k] = B[i][j][k] + ZN[i][j][k]*x[i-1][j][k];
                }
            }
            x[i] = converge_polar_tdma(x[i], p, tp, tn, rp, rn, b, alpha, n, m);
        }
        else
        {
            for(int j=0; j<n; j++)
            {
                for (int k=0; k<m; k++)
                {
                    tp[j][k] = TP[i][j][k];
                    tn[j][k] = TN[i][j][k];
                    rp[j][k] = RP[i][j][k];
                    rn[j][k] = RN[i][j][k];
                    p[j][k] = P[i][j][k];
                    b[j][k] = B[i][j][k] + ZN[i][j][k]*x[i-1][j][k] + ZP[i][j][k]*x[i+1][j][k];
                }
            }
            x[i] = converge_polar_tdma(x[i], p, tp, tn, rp, rn, b, alpha, n, m);
        }
    }

    return x;
}
double*** converge_cylindrical_tdma(double*** x, double*** P, double*** TP, double*** TN, double*** RP, double*** RN, double*** ZP, double*** ZN, double*** B, double alpha, int l, int n, int m, double r_in, double r_out, double t)
{
    // double*** x = initialize3d(l, n, m);
    double*** x_ = initialize3d(l, n, m);

    double error = 1;
    int iter = 0;
    do
    {
        x = cylindrical_tdma(x, P, TP, TN, RP, RN, ZP, ZN, B, alpha, l, n, m);
        x = steady3d(x, l, n, m);
        error = error3d(x, x_, l, n, m);
        x_ = copy3d(x, l, n, m);
        iter += 1;
        print_paraview3d(x, l, n, m, r_in, r_out, t);
        cout << "Iteration [" << iter << "/" << ITER <<  "] Error: " << error << endl;
    } while ((error>ERROR) && (iter<ITER));

    return x;
}


double calculate_SP(int i, int j, int k, int l, int n, int m)
{
    return 0;
}
double calculate_SC(int i, int j, int k, int l, int n, int m)
{   
    return 0;
}

double*** coefficients_TP(int l, int n, int m, double gamma, double h, double r_in, double r_out, double t)
{
    double*** TP = initialize3d(l, n, m);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                TP[i][j][k] = gamma*del_z*del_r/((r_in + (j+1./2)*del_r) * del_theta);
            }
        }
    }

    return TP;
}
double*** coefficients_TN(int l, int n, int m, double gamma, double h, double r_in, double r_out, double t)
{
    double*** TN = initialize3d(l, n, m);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                TN[i][j][k] = gamma*del_z*del_r/((r_in + (j+1./2)*del_r) * del_theta);
            }
        }
    }

    return TN;
}
double*** coefficients_RP(int l, int n, int m, double gamma, double h, double r_in, double r_out, double t)
{
    double*** RP = initialize3d(l, n, m);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;
    double k;

    for (int i=0; i<l; i++)
    { 
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                if (j==n-1)
                {
                    k = 2*gamma/del_r;
                    RP[i][j][k] = h*k/(h+k) * del_z * (r_in + (j+1)*del_r)*del_theta;
                    // RP[i][j][k] = 2*gamma*del_z*(r_in + (j+1)*del_r)*del_theta/del_r;
                }
                else
                {
                    RP[i][j][k] = gamma*del_z*(r_in + (j+1)*del_r)*del_theta/del_r;
                }
            }
        }
    }

    return RP;
}
double*** coefficients_RN(int l, int n, int m, double gamma, double h, double r_in, double r_out, double t)
{
    double*** RN = initialize3d(l, n, m);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                if (j==0)
                {
                    RN[i][j][k] = 2*gamma*del_z*(r_in + j*del_r)*del_theta/del_r;
                }
                else
                {
                    RN[i][j][k] = gamma*del_z*(r_in + j*del_r)*del_theta/del_r;
                }
                    
            }
        }
    }

    return RN;
}
double*** coefficients_ZP(int l, int n, int m, double gamma, double h, double r_in, double r_out, double t)
{
    double*** ZP = initialize3d(l, n, m);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;
    double area, k;

    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                area = (r_in + (j+1)*del_r)*(r_in + (j+1)*del_r)*del_theta/2 - (r_in + j*del_r)*(r_in + j*del_r)*del_theta/2;
                if (i==l-1)
                {
                    k = 2*gamma/del_z;
                    ZP[i][j][k] = h*k/(h+k) * area;
                    // ZP[i][j][k] = 2*gamma*area/del_z;
                }
                else
                {
                    ZP[i][j][k] = gamma*area/del_z;
                }
            }
        }
    }

    return ZP;
}
double*** coefficients_ZN(int l, int n, int m, double gamma, double h, double r_in, double r_out, double t)
{
    double*** ZN = initialize3d(l, n, m);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;
    double area, k;

    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                area = (r_in + (j+1)*del_r)*(r_in + (j+1)*del_r)*del_theta/2 - (r_in + j*del_r)*(r_in + j*del_r)*del_theta/2;
                if (i==0)
                {
                    k = 2*gamma/del_z;
                    ZN[i][j][k] = h*k/(h+k) * area;
                    // ZN[i][j][k] = 2*gamma*area/del_z;
                }
                else
                {
                    ZN[i][j][k] = gamma*area/del_z;
                }
            }
        }
    }

    return ZN;
}
double*** coefficients_P(int l, int n, int m, double gamma, double h, double r_in, double r_out, double t)
{
    double*** P = initialize3d(l, n, m);
    double*** TP = coefficients_TP(l, n, m, gamma, h, r_in, r_out, t);
    double*** TN = coefficients_TN(l, n, m, gamma, h, r_in, r_out, t);
    double*** RP = coefficients_RP(l, n, m, gamma, h, r_in, r_out, t);
    double*** RN = coefficients_RN(l, n, m, gamma, h, r_in, r_out, t);
    double*** ZP = coefficients_ZP(l, n, m, gamma, h, r_in, r_out, t);
    double*** ZN = coefficients_ZN(l, n, m, gamma, h, r_in, r_out, t);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;
    double volume;

    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                volume = ((r_in + (j+1)*del_r)*(r_in + (j+1)*del_r)*del_theta/2 - (r_in + j*del_r)*(r_in + j*del_r)*del_theta/2)*del_z;
                P[i][j][k] = TP[i][j][k] + TN[i][j][k] + RP[i][j][k] + RN[i][j][k] + ZP[i][j][k] + ZN[i][j][k] - calculate_SP(i, j, k, l, n, m)*volume;
            }
        }
    }

    return P;
}
double*** coefficients_B(int l, int n, int m, double phi_r_in, double phi_r_out, double phi_inf_up, double phi_inf_down, double gamma, double h, double r_in, double r_out, double t)
{
    double*** B = initialize3d(l, n, m);
    double*** TP = coefficients_TP(l, n, m, gamma, h, r_in, r_out, t);
    double*** TN = coefficients_TN(l, n, m, gamma, h, r_in, r_out, t);
    double*** RP = coefficients_RP(l, n, m, gamma, h, r_in, r_out, t);
    double*** RN = coefficients_RN(l, n, m, gamma, h, r_in, r_out, t);
    double*** ZP = coefficients_ZP(l, n, m, gamma, h, r_in, r_out, t);
    double*** ZN = coefficients_ZN(l, n, m, gamma, h, r_in, r_out, t);
    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;
    double volume;

    for (int i=0; i<l; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<m; k++)
            {
                volume = ((r_in + (j+1)*del_r)*(r_in + (j+1)*del_r)*del_theta/2 - (r_in + j*del_r)*(r_in + j*del_r)*del_theta/2)*del_z;
                B[i][j][k] = calculate_SP(i, j, k, l, n, m)*volume;
                if (i==0)
                {
                    B[i][j][k] += ZN[i][j][k] * phi_inf_down;
                }
                if (i==l-1)
                {
                    B[i][j][k] += ZP[i][j][k] * phi_inf_up;
                }
                if (j==0)
                {
                    B[i][j][k] += RN[i][j][k] * phi_r_in;
                }
                if (j==n-1)
                {
                    B[i][j][k] += RP[i][j][k] * phi_r_out;
                }
            }
        }
    }

    return B;
}


int main()
{
    int l = 3, n = 32, m = 64;
    double gamma = 25, h = 10;
    double phi_r_in = 1000, phi_r_out = 40, phi_inf_up = 40, phi_inf_down = 40;
    double alpha = 1;
    double r_in = 1, r_out = 4, t = 1.6;

    double*** TP = coefficients_TP(l, n, m, gamma, h, r_in, r_out, t);
    double*** TN = coefficients_TN(l, n, m, gamma, h, r_in, r_out, t);
    double*** RP = coefficients_RP(l, n, m, gamma, h, r_in, r_out, t);
    double*** RN = coefficients_RN(l, n, m, gamma, h, r_in, r_out, t);
    double*** ZP = coefficients_ZP(l, n, m, gamma, h, r_in, r_out, t);
    double*** ZN = coefficients_ZN(l, n, m, gamma, h, r_in, r_out, t);
    double*** P = coefficients_P(l, n, m, gamma, h, r_in, r_out, t);
    double*** B = coefficients_B(l, n, m, phi_r_in, phi_r_out, phi_inf_up, phi_inf_down, gamma, h, r_in, r_out, t);

    double*** x = initialize3d(l, n, m);

    ifstream file;
    file.open("paraview_3d.txt");
    if (file)
    {
        x = read_paraview3d(l, n, m, r_in, r_out, t);
        file.close();
    }

    x = converge_cylindrical_tdma(x, P, TP, TN, RP, RN, ZP, ZN, B, alpha, l, n, m, r_in, r_out, t);
    print_paraview3d(x, l, n, m, r_in, r_out, t);
    cout << endl;
    // print3d(x, l, n, m);

    double del_z = t/l;
    double del_r = (r_out - r_in)/n;
    double del_theta = 2*PI/m;
    double area;

    double q=0, qmax=0;
    for (int j=0; j<n; j++)
    {
        for (int k=0; k<m; k++)
        {   
            area = (r_in + (j+1)*del_r)*(r_in + (j+1)*del_r)*del_theta/2 - (r_in + j*del_r)*(r_in + j*del_r)*del_theta/2;
            q += - gamma * area * (x[l-1][j][k]-phi_inf_up) - gamma * area * (x[0][j][k]-phi_inf_down);
            qmax += - gamma * area * (phi_r_in-phi_inf_up) - gamma * area * (phi_r_in-phi_inf_down);
        }
    }
    cout << "Thickness = " << t << " | "  << "Efficiency = " << q/qmax << endl;

    // for (int j=0; j<n; j++)
    // {
    //     cout << x[0][j][0] << ", ";
    // }

    return 0;
}