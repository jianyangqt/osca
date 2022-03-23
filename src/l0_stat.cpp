//
//  l0_stat.cpp
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l0_stat.h"


double betacf(double a, double b, double x)
{
    const int MAXIT=300;
    const double Eps=1.0e-08;
    const double FPMIN=(std::numeric_limits<double>::min)()/numeric_limits<double>::epsilon();
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if(abs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for(m=1;m<=MAXIT;m++)
    {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if(abs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if(abs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if(abs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if(abs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if(abs(del-1.0) <= Eps) break;
    }
    return h;
}


double gammln(double xx)
{
    int j;
    double x,y,tmp,ser;
    static const double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
        -0.5395239384953e-5};
    
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for(j=0;j<6;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

double betai(double a, double b, double x)
{
    double bt;    
    if(x < 0.0 || x > 1.0)
    {
        sprintf(logbuf, "Bad x in routine betai!.\n");
        logprintb();
        TERMINATE();
    }
    if(x == 0.0 || x == 1.0) bt=0.0;
    else bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
    
    if(x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
    else return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double F_prob(double df_1, double df_2, double F_value) {
    return betai(df_2 * 0.5, df_1 * 0.5, df_2 / (df_2 + df_1 * F_value));
}

double t_prob(double df, double t_value, bool two_tail)
{
    double p=betai( df*0.5, 0.5, df/(df+(t_value*t_value)) );
    
    if(two_tail) return p;
    
    return p*0.5;
}

double pchisq(double x, double df) {
    if (x < 0) return -9;
    
    double p, q;
    int st = 0; // error variable
    int w = 1; // function variable
    double bnd = 1; // boundary function
    
    // NCP is set to 0
    cdfchi(&w, &p, &q, &x, &df, &st, &bnd);
    
    // Check status
    if (st != 0) return -9;
    
    // Return p-value
    return q;
}
// q is not good to be less than 1e-161
double qchisq(double q, double df) {
    if (q < 0) return -9;
    else if (q >= 1) return 0;
    
    double x;
    double p = 1 - q;
    int st = 0; // error variable
    int w = 2; // function variable
    double bnd = 1; // boundary function
    
    // NCP is set to 0
    cdfchi(&w, &p, &q, &x, &df, &st, &bnd);
    
    // Check status
    if (st != 0) return -9;
    if(q<1e-161){
        double tmp= pchisq(x, 1);
        while(q/tmp > 10)
        {
            x-=10;
            tmp= pchisq(x, 1);
        }
        if(q/tmp>1)
        {
            while(q/tmp>1.01){
                x-=0.1;
                tmp= pchisq(x, 1);
            }
        } else if(q/tmp<1) {
            while(q/tmp<1.01){
                x+=0.1;
                tmp= pchisq(x, 1);
            }
        }
    }
    // Return p-value
    return x;
}


double median(vector<double> vec)
{
    // Computes approxiamte median of elements in a vector
    //  Input:
    //   vec: Vector with values to compute the median.
    //        The vector is copied since it should not be modified
    vector<double>::size_type mid(vec.size()/2); // compute half the size of the matrix (rounded to zero)
    nth_element (vec.begin(), vec.begin()+mid, vec.end()); // Find median element
    return vec[mid]; // Return median element
}
// Default: upper-tail
double pnorm(double x) {
    double z = 0.0;
    if(x>0) z = -1.0 * x;
    else z = x;
    double sqrt2pi = 2.50662827463;
    double t0, z1, p0;
    t0 = 1 / (1 + 0.2316419 * fabs(z));
    z1 = exp(-0.5 * z * z) / sqrt2pi;
    p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
    return x >= 0 ? p0 : 1.0 - p0;
}
double qnorm_sub(double x, double y) {
    return (y + 0.5 * x * pow(y, 2)+(2 * pow(x, 2) + 1) * pow(y, 3) / 6 + (6 * pow(x, 3) + 7 * x) * pow(y, 4) / 12);
}

double dnorm(double x) {
    return (0.39894228 * exp(-0.5 * pow(x, 2)));
}
// Default: lower-tail
double qnorm(double p, bool upper) {
    double x = 0;
    if (upper) p = 1 - p;
    for (int i = 0; i < 4; i++) x = x + qnorm_sub(x, (pnorm(x) - p) / dnorm(x));
    return (x);
}
void gcf(double &gammcf, const double a, const double x, double &gln) {
    const int ITMAX = 100;
    const double EPS = numeric_limits<double>::epsilon();
    const double FPMIN = (std::numeric_limits<double>::min)() / EPS;
    int i;
    double an, b, c, d, del, h;
    
    gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1; i <= ITMAX; i++) {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (Abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (Abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = d*c;
        h *= del;
        if (Abs(del - 1.0) <= EPS) break;
    }
    if (i > ITMAX) throw ("a too large, ITMAX too small in gcf");
    gammcf = exp(-x + a * log(x) - gln) * h;
}
void gser(double &gamser, const double a, const double x, double &gln) {
    const int ITMAX = 500;
    const double Eps = 1.0e-08;
    int n;
    double sum, del, ap;
    
    gln = gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) throw ("x less than 0 in routine gser");
        gamser = 0.0;
        return;
    } else {
        ap = a;
        del = sum = 1.0 / a;
        for (n = 0; n < ITMAX; n++) {
            ++ap;
            del *= x / ap;
            sum += del;
            if (Abs(del) < Abs(sum) * Eps) {
                gamser = sum * exp(-x + a * log(x) - gln);
                return;
            }
        }
        throw ("a too large, ITMAX too small in routine gser");
        return;
    }
}

double gammp(const double a, const double x) {
    double gamser, gammcf, gln;
    
    if (x < 0.0 || a <= 0.0) throw ("Invalid arguments in routine gammp");
    
    if (x < a + 1.0) {
        gser(gamser, a, x, gln);
        return gamser;
    } else {
        gcf(gammcf, a, x, gln);
        return 1.0 - gammcf;
    }
}

double chi_prob(double df, double chi_sqr_val) {
    return 1 - gammp(df * 0.5, chi_sqr_val * 0.5);
}
double ran1(int &idum) {
    const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32;
    const int NDIV = (1 + (IM - 1) / NTAB);
    const double EPS = 3.0e-16, AM = 1.0 / IM, RNMX = (1.0 - EPS);
    static int iy = 0;
    static vector<int> iv(NTAB);
    int j, k;
    double temp;
    
    if (idum <= 0 || !iy) {
        if (-idum < 1) idum = 1;
        else idum = -idum;
        for (j = NTAB + 7; j >= 0; j--) {
            k = idum / IQ;
            idum = IA * (idum - k * IQ) - IR*k;
            if (idum < 0) idum += IM;
            if (j < NTAB) iv[j] = idum;
        }
        iy = iv[0];
    }
    k = idum / IQ;
    idum = IA * (idum - k * IQ) - IR*k;
    if (idum < 0) idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = idum;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}


double gasdev(int &idum) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;
    
    if (idum < 0) iset = 0;
    if (iset == 0) {
        do {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            rsq = v1 * v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset = 0;
        return gset;
    }
}
// plink 1
bool realnum(double d)
{
    double zero = 0;
    if (d != d || d == 1/zero || d == -1/zero)
        return false;
    else
        return true;
}

double pT(double T, double df)
{
    
    if ( ! realnum(T) )
        return -9;
    
    T = abs(T);
    
    double p, q;
    int st = 0;      // error variable
    int w = 1;       // function variable
    double bnd = 1;  // boundary function
    
    // NCP is set to 0
    cdft(&w,&p,&q,&T,&df,&st,&bnd);
    
    // Check status
    if (st != 0 ) return -9;
    
    // Return two-sided p-value
    return 2*q;
    
}


double K(double zeta, VectorXd &lambda) {
    return -0.5*(1.0 - 2.0 * zeta * lambda.array()).log().sum();
}

double Kp(double zeta, VectorXd &lambda) {
    return (lambda.array() / (1.0 - 2.0 * zeta * lambda.array())).sum();
}

double Kpp(double zeta, VectorXd &lambda) {
    return 2.0 * (lambda.array().square() / (1.0 - 2.0 * zeta * lambda.array()).array().square()).sum();
}

double Kp_min_x(double zeta, VectorXd &lambda, double x) {
    return Kp(zeta, lambda) - x;
}

double Brents_Kp_min_x(VectorXd &lambda, double x, double lowerLimit, double upperLimit, double errorTol) {
    double a = lowerLimit;
    double b = upperLimit;
    double c = 0;
    double d = 1.7976931348623157E+308;
    
    double fa = Kp_min_x(a, lambda, x);
    double fb = Kp_min_x(b, lambda, x);
    
    double fc = 0;
    double s = 0;
    double fs = 0;
    
    // if f(a) f(b) >= 0 then error-exit
    if (fa * fb >= 0) {
        if (fa < fb)
            return a;
        else
            return b;
    }
    
    // if |f(a)| < |f(b)| then swap (a,b) end if
    if (fabs(fa) < fabs(fb)) {
        double tmp = a;
        a = b;
        b = tmp;
        tmp = fa;
        fa = fb;
        fb = tmp;
    }
    
    c = a;
    fc = fa;
    bool mflag = true;
    int i = 0;
    
    while (!(fb == 0) && (fabs(a - b) > errorTol)) {
        if ((fa != fc) && (fb != fc))
            // Inverse quadratic interpolation
            s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
        else
            // Secant Rule
            s = b - fb * (b - a) / (fb - fa);
        
        double tmp2 = (3 * a + b) / 4;
        if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) || (mflag && (fabs(s - b) >= (fabs(b - c) / 2))) || (!mflag && (fabs(s - b) >= (fabs(c - d) / 2)))) {
            s = (a + b) / 2;
            mflag = true;
        } else {
            if ((mflag && (fabs(b - c) < errorTol)) || (!mflag && (fabs(c - d) < errorTol))) {
                s = (a + b) / 2;
                mflag = true;
            } else
                mflag = false;
        }
        fs = Kp_min_x(s, lambda, x);
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        
        // if |f(a)| < |f(b)| then swap (a,b) end if
        if (fabs(fa) < fabs(fb)) {
            double tmp = a;
            a = b;
            b = tmp;
            tmp = fa;
            fa = fb;
            fb = tmp;
        }
        i++;
        if (i > 1000) return upperLimit+10;
    }
    return b;
}




double psatt(double x, VectorXd lambda) {
    double sum=lambda.sum();
    if(FloatEqual(sum,0.0)) return 2.0;
    
    double sq_sum=lambda.dot(lambda);
    double sum_sq=sum*sum;
    double a=sq_sum/sum;
    double b=sum_sq/sq_sum;
    
    return pchisq(x/a, b);
}

double psadd(double x, VectorXd lambda) {
    double d = lambda.maxCoeff();
    if (d <= 0.0) return 2.0;
    lambda = lambda.array() / d;
    x = x / d;
    
    double lmin = 0.0;
    double m = lambda.minCoeff();
    if (m < 0.0) lmin = 0.499995 / m;
    else if (x > lambda.sum()) lmin = -0.01;
    else lmin = -0.5 * (double) lambda.size() / x;
    double lmax = 0.499995 / lambda.maxCoeff();
    
    double hatzeta = Brents_Kp_min_x(lambda, x, lmin, lmax, 1e-08);
    if(hatzeta > lmax + 9) return 2.0;
    double sign = (hatzeta < 0.0) ? -1.0 : 1.0;
    double w = sign * sqrt(2 * (hatzeta * x - K(hatzeta, lambda)));
    double v = hatzeta * sqrt(Kpp(hatzeta, lambda));
    
    // debug
    //cout<<"hatzeta = "<<hatzeta<<endl;
    //cout<<"w = "<<w<<endl;
    //cout<<"v = "<<v<<endl;
    
    
    if (fabs(hatzeta) < 1e-04) return 2.0;
    else return pnorm(w + log(v / w) / w);
}
double pchisqsum(double x, VectorXd lambda) {
    double pval = psadd(x, lambda);
    if (pval > 1.0) pval = psatt(x, lambda);
    return pval;
}
