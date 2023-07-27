#include <bits/stdc++.h>
using namespace std;

const double PI = 3.14159265358979323846;
const double g = 9.80665;

double max(double a, double b)
{
    return (a>b)?a:b;
}

double* initialize1d(int n)
{
    double* x = new double[n];
    for (int i=0; i<n; i++)
    {
        x[i] = 0;
    }
    return x;
}
double** initialize2d(int n, int m)
{
    double** x = new double*[n];
    for (int i=0; i<n; i++)
    {
        x[i] = new double[m];
        for (int j=0; j<m; j++)
        {
            x[i][j] = 0;
        }
    }
    return x;
}

double* store1d(double* x, int n)
{
    double* x_ = initialize1d(n);
    for (int i=0; i<n; i++)
    {
        x_[i] = x[i];
    }
    return x_;
}

double** store2d(double** x, int n, int m)
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

void print1d(double* x, int n)
{
    for (int i=0; i<n; i++)
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
void print2darray(double** x, int n, int m)
{
    ofstream file("out.txt");
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            file << x[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}
void print_paraview2d(double** x, int n, int m, double rin, double rout)
{
    double dr = (rout-rin)/n;
    double dtheta = 2*PI/m;
    double r, theta;

    ofstream file("paraview_2d.txt");
    file << "X, Y, Z, T" << endl;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {   
            r = rin + (i+0.5)*dr;
            theta = (j+0.5)*dtheta;
            file << r*cos(theta) << ", " << r*sin(theta) << ", " << 0 << ", " << x[i][j] << endl;
        }
    }
    file.close();
}

double error1d(double* x, double* x_, int n)
{
    double err = 0;
    for (int i=0; i<n; i++)
    {
        err += pow(x[i] - x_[i], 2);
    }
    return sqrt(err);
}

double error2d(double** x, double** x_, int n, int m)
{
    double err = 0;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            err += pow(x[i][j] - x_[i][j], 2);
        }
    }
    return sqrt(err);
}

// TP - Theta Positive; TN - Theta Negative
double* tdma(double* P, double* TP, double* TN, double* B, int n)
{
    double* x = initialize1d(n);

    for (int i=0; i<n; i++)
    {
        if (i==0)
        {
            B[i] = B[i]/P[i];
            TP[i] = TP[i]/P[i];
            P[i] = 1;
        }
        else if (i==n-1)
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

    for (int i=n-1; i>=0; i--)
    {
        if (i==n-1)
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
double* circular_tdma(double* x, double* P, double* TP, double* TN, double* B, int n)
{
    double* tp = initialize1d(n);
    double* tn = initialize1d(n);
    double* p = initialize1d(n);
    double* b = initialize1d(n);

    for (int i=0; i<n; i++)
    {
        int end = (i+n-1)%n;
        int start = i;
        int iter = 1;

        tp[0] = TP[start];
        tn[0] = TN[start];
        p[0] = P[start];
        b[0] = B[start] + TN[start]*x[end];
        for (int j=(start+1)%n; iter<n; j=(j+1)%n)
        {
            tp[iter] = TP[j];
            tn[iter] = TN[j];
            p[iter] = P[j];
            b[iter] = B[j];

            iter += 1;
        }

        x = tdma(p, tp, tn, b, n);
    }

    return x;
}
double* converge_circular_tdma(double* x, double* P, double* TP, double* TN, double* B, double alpha, int n)
{
    double* x_ = initialize1d(n);

    double error = 1;
    do
    {   
        for (int i=0; i<n; i++)
        {
            P[i] = P[i]/alpha;
            B[i] = B[i] + (1-alpha)*P[i]*x[i];
        }
        x = circular_tdma(x, P, TP, TN, B, n);
        error = error1d(x, x_, n);
        x_ = store1d(x, n);
    } while (error>1e-3);

    return x;
}

// TP - Theta Positive; TN - Theta Negative
// RP - R Positive; RN - R Negative
double** polar_tdma(double** x, double**P, double** TP, double** TN, double** RP, double** RN, double** B, double alpha, int n, int m)
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
                    b[j] = B[i][j] + RN[i][j]*x[i-1][j];
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
                    b[j] = B[i][j] + RP[i][j]*x[i+1][j];
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
                    b[j] = B[i][j] + RP[i][j]*x[i+1][j] + RN[i][j]*x[i-1][j];
                }
                x[i] = converge_circular_tdma(x[i], p, tp, tn, b, alpha, m);
            }
    }

    return x;
}
double** converge_polar_tdma(double**P, double** TP, double** TN, double** RP, double** RN, double** B, double alpha, int n, int m)
{
    double** x = initialize2d(n, m);
    double** x_ = initialize2d(n, m);

    double error = 1;
    do
    {
        x = polar_tdma(x, P, TP, TN, RP, RN, B, alpha, n, m);
        error = error2d(x, x_, n, m);
        // cout << "Iteration Error: " << error << endl;
        x_ = store2d(x, n, m);
    } while (error>1e-3);

    return x;
}

double calculate_SP(int i, int j, int n, int m)
{
    return 0;
}
double calculate_SC(int i, int j, int n, int m)
{   
    return 0;
}


double** coefficients_TP(int n, int m, double gamma, double rin, double rout, double t)
{
    double** TP = initialize2d(n, m);
    double del_r = (rout - rin)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            TP[i][j] = gamma*t*del_r/((rin + (i+1./2)*del_r) * del_theta);
        }
    }

    return TP;
}
double** coefficients_TN(int n, int m, double gamma, double rin, double rout, double t)
{
    double** TN = initialize2d(n, m);
    double del_r = (rout - rin)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            TN[i][j] = gamma*t*del_r/((rin + (i+1./2)*del_r) * del_theta);
        }
    }

    return TN;
}
double** coefficients_RP(int n, int m, double gamma, double rin, double rout, double t)
{
    double** RP = initialize2d(n, m);
    double del_r = (rout - rin)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {   
            if (i==n-1)
            {
                RP[i][j] = 2*gamma*t*(rin + (i+1)*del_r)*del_theta/del_r;
            }
            else
            {
                RP[i][j] = gamma*t*(rin + (i+1)*del_r)*del_theta/del_r;
            }
        }
    }

    return RP;
}
double** coefficients_RN(int n, int m, double gamma, double rin, double rout, double t)
{
    double** RN = initialize2d(n, m);
    double del_r = (rout - rin)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            if (i==0)
            {
                RN[i][j] = 2*gamma*t*(rin + i*del_r)*del_theta/del_r;
            }
            else
            {
                RN[i][j] = gamma*t*(rin + i*del_r)*del_theta/del_r;
            }
        }
    }

    return RN;
}
double** coefficients_P(double** x, int n, int m, double gamma, double mu, double Pr, double rin, double rout, double t, double phi_infup, double phi_infdown)
{
    double** P = initialize2d(n, m);
    double** TP = coefficients_TP(n, m, gamma, rin, rout, t);
    double** TN = coefficients_TN(n, m, gamma, rin, rout, t);
    double** RP = coefficients_RP(n, m, gamma, rin, rout, t);
    double** RN = coefficients_RN(n, m, gamma, rin, rout, t);
    double del_r = (rout - rin)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {   
            double area = (rin + (i+1)*del_r)*(rin + (i+1)*del_r)*del_theta/2 - (rin + i*del_r)*(rin + i*del_r)*del_theta/2;
            double L = (PI*rout*rout-PI*rin*rin)/(2*PI*rout-2*PI*rin);
            double U = U_up(x[i][j], phi_infup, gamma, mu, Pr, L, t) + U_down(x[i][j], phi_infdown, gamma, mu, Pr, L, t);
            P[i][j] = TP[i][j] + TN[i][j] + RP[i][j] + RN[i][j] + U*area - calculate_SP(i, j, n, m)*area*t;
        }
    }
    

    return P;
}
double** coefficients_B(double** x, int n, int m, double gamma, double mu, double Pr, double rin, double rout, double t, double phi_infup, double phi_infdown, double phi_rin, double phi_rout)
{
    double** B = initialize2d(n, m);
    // double** TP = coefficients_TP(n, m, gamma, rin, rout, t);
    // double** TN = coefficients_TN(n, m, gamma, rin, rout, t);
    double** RP = coefficients_RP(n, m, gamma, rin, rout, t);
    double** RN = coefficients_RN(n, m, gamma, rin, rout, t);
    double del_r = (rout - rin)/n;
    double del_theta = 2*PI/m;

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {   
            double area = (rin + (i+1)*del_r)*(rin + (i+1)*del_r)*del_theta/2 - (rin + i*del_r)*(rin + i*del_r)*del_theta/2;
            double L = (PI*rout*rout-PI*rin*rin)/(2*PI*rout-2*PI*rin);
            double Uphi = U_up(x[i][j], phi_infup, gamma, mu, Pr, L, t)*phi_infup + U_down(x[i][j], phi_infdown, gamma, mu, Pr, L, t)*phi_infdown;
            if (i==0)
            {
                B[i][j] = RN[i][j]*phi_rin + Uphi*area + calculate_SC(i, j, n, m)*area*t;
            }
            else if (i==n-1)
            {
                B[i][j] = RP[i][j]*phi_rout + Uphi*area + calculate_SC(i, j, n, m)*area*t;
            }
            else
            {
                B[i][j] = Uphi*area + calculate_SC(i, j, n, m)*area*t;
            }
        }
    }

    return B;
}

int main()
{
    int n = 25, m = 8;
    double gamma = 25.32, mu = 15.89e-6, Pr = 0.707;
    double phi_rin = 500, phi_rout = 40, phi_infup = 40, phi_infdown = 40;
    double alpha = 0.8;
    double rin = 0.5, rout = 4, t = 0.001;

    double** x = initialize2d(n, m);
    double** x_ = initialize2d(n, m);

    double** TP = coefficients_TP(n, m, gamma, rin, rout, t);
    double** TN = coefficients_TN(n, m, gamma, rin, rout, t);
    double** RP = coefficients_RP(n, m, gamma, rin, rout, t);
    double** RN = coefficients_RN(n, m, gamma, rin, rout, t);

    double error = 1;
    do
    {
        double** P = coefficients_P(x, n, m, gamma, mu, Pr, rin, rout, t, phi_infup, phi_infdown);
        double** B = coefficients_B(x, n, m, gamma, mu, Pr, rin, rout, t, phi_infup, phi_infdown, phi_rin, phi_rout);

        x = converge_polar_tdma(P, TP, TN, RP, RN, B, alpha, n, m);
        error = error2d(x, x_, n, m);
        x_ = store2d(x, n, m);
        cout << "Iteration Error " << error << endl;
    } while (error>1e-6);
    
    
    print2darray(x, n, m);
    print_paraview2d(x, n, m, rin, rout);
    
    // cout << calculate_SP(100, phi_infup, phi_infdown, gamma, mu, Pr, t) << endl;

    return 0;
}