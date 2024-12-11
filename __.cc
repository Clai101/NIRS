int factorial(int n) {
    if (n == 0) return 1;
    else return n*factorial(n-1);
}

double poisson(double l, int n) {
  return pow(l, n)*exp(-l)/factorial(n);
}


double g_g_g(double *x, double *par) {

    double sqrt2pi = 2.50662827;

    double mean0   = par[2];
    double mean01   = par[2] + par[7];
    double mean1   = mean0 + par[4];
    double mean11   = mean0 + par[4] + par[7];
    
    double mean2 = mean0 + 2*par[4];
    double mean21 = mean0 + 2*par[4]+par[7];
    
    double mean3 = mean0 + 3*par[4];
    double mean31 = mean0 + 3*par[4] + par[7];
    
    double mean4 = mean0 + 4*par[4];
    double mean41 = mean0 + 4*par[4] + par[7];
    
    double mean5 = mean0 + 5*par[4];
    double mean51 = mean0 + 5*par[4] + par[7];

    double mean6 = mean0 + 6*par[4];
    double mean61 = mean0 + 6*par[4] + par[7];
    
    double mean7 = mean0 + 7*par[4];
    double mean71 = mean0 + 7*par[4] + par[7];
    
    double mean8 = mean0 + 8*par[4];
    double mean81 = mean0 + 8*par[4] + par[7];
    
    double sigmaG0 = par[3];
    double sigmaG01 = par[3]*par[8];
    
    double sigmaG1 = sqrt(par[3]*par[3]+par[5]*par[5]);
    double sigmaG11 = sqrt(par[3]*par[3]*par[8]*par[8]+par[5]*par[5]);
    
    double sigmaG2 = sqrt(par[3]*par[3]+2*par[5]*par[5]);
    double sigmaG21 = sqrt(par[3]*par[3]*par[8]*par[8]+2*par[5]*par[5]);
    
    double sigmaG3 = sqrt(par[3]*par[3]+3*par[5]*par[5]);
    double sigmaG31 = sqrt(par[3]*par[3]*par[8]*par[8]+3*par[5]*par[5]);
    
    double sigmaG4 = sqrt(par[3]*par[3]+4*par[5]*par[5]);
    double sigmaG41 = sqrt(par[3]*par[3]*par[8]*par[8]+4*par[5]*par[5]);
    
    double sigmaG5 = sqrt(par[3]*par[3]+5*par[5]*par[5]);
    double sigmaG51 = sqrt(par[3]*par[3]*par[8]*par[8]+5*par[5]*par[5]);
    
    double sigmaG6 = sqrt(par[3]*par[3]+6*par[5]*par[5]);
    double sigmaG61 = sqrt(par[3]*par[3]*par[8]*par[8]+6*par[5]*par[5]);
    
    double sigmaG7 = sqrt(par[3]*par[3]+7*par[5]*par[5]);
    double sigmaG71 = sqrt(par[3]*par[3]*par[8]*par[8]+7*par[5]*par[5]);
    
    double sigmaG8 = sqrt(par[3]*par[3]+8*par[5]*par[5]);
    double sigmaG81 = sqrt(par[3]*par[3]*par[8]*par[8]+8*par[5]*par[5]);

    double g00 = exp( -(*x - mean0)*(*x - mean0) / (2*sigmaG0*sigmaG0));
    g00 /= (sqrt2pi * sigmaG0);
    
    double g01 = exp( -(*x - mean01)*(*x - mean01) / (2*sigmaG01*sigmaG01));
    g01 /= (sqrt2pi * sigmaG01);
    
    double g10 = exp( -(*x - mean1)*(*x - mean1) / (2*sigmaG1*sigmaG1));
    g10 /= (sqrt2pi * sigmaG1);
    
    double g11 = exp( -(*x - mean11)*(*x - mean11) / (2*sigmaG11*sigmaG11));
    g11 /= (sqrt2pi * sigmaG11);
    
    double g20 = exp( -(*x - mean2)*(*x - mean2) / (2*sigmaG2*sigmaG2));
    g20 /= (sqrt2pi * sigmaG2);
    
    double g21 = exp( -(*x - mean21)*(*x - mean21) / (2*sigmaG21*sigmaG21));
    g21 /= (sqrt2pi * sigmaG21);
    
    double g30 = exp( -(*x - mean3)*(*x - mean3) / (2*sigmaG3*sigmaG3));
    g30 /= (sqrt2pi * sigmaG3);
    
    double g31 = exp( -(*x - mean31)*(*x - mean31) / (2*sigmaG31*sigmaG31));
    g31 /= (sqrt2pi * sigmaG31);
    
    double g40 = exp( -(*x - mean4)*(*x - mean4) / (2*sigmaG4*sigmaG4));
    g40 /= (sqrt2pi * sigmaG4);
    
    double g41 = exp( -(*x - mean41)*(*x - mean41) / (2*sigmaG41*sigmaG41));
    g41 /= (sqrt2pi * sigmaG41);
    
    double g50 = exp( -(*x - mean5)*(*x - mean5) / (2*sigmaG5*sigmaG5));
    g50 /= (sqrt2pi * sigmaG5);
    
    double g51 = exp( -(*x - mean51)*(*x - mean51) / (2*sigmaG51*sigmaG51));
    g51 /= (sqrt2pi * sigmaG51);
    
    double g60 = exp( -(*x - mean6)*(*x - mean6) / (2*sigmaG6*sigmaG6));
    g60 /= (sqrt2pi * sigmaG6);
    
    double g61 = exp( -(*x - mean61)*(*x - mean61) / (2*sigmaG61*sigmaG61));
    g61 /= (sqrt2pi * sigmaG61);
    
    double g70 = exp( -(*x - mean7)*(*x - mean7) / (2*sigmaG7*sigmaG7));
    g70 /= (sqrt2pi * sigmaG7);
    
    double g71 = exp( -(*x - mean71)*(*x - mean71) / (2*sigmaG71*sigmaG71));
    g71 /= (sqrt2pi * sigmaG71);
    
    double g80 = exp( -(*x - mean8)*(*x - mean8) / (2*sigmaG8*sigmaG8));
    g80 /= (sqrt2pi * sigmaG8);
    
    double g81 = exp( -(*x - mean81)*(*x - mean81) / (2*sigmaG81*sigmaG81));
    g81 /= (sqrt2pi * sigmaG81);
    
    return poisson(par[1],0)*((1-par[9])*g00 + par[9]*g01) + (1-par[6])*poisson(par[1],1)*((1-par[9])*g10 + par[9]*g11) \
        + (par[6]*poisson(par[1],1) + (1-2*par[6])*poisson(par[1],2))*((1-par[9])*g20 + par[9]*g21) \
        + (2*par[6]*poisson(par[1],2) + (1-3*par[6])*poisson(par[1],3))*((1-par[9])*g30 + par[9]*g31) \
        + (3*par[6]*poisson(par[1],3) + (1-4*par[6])*poisson(par[1],4))*((1-par[9])*g40 + par[9]*g41) \
        + (4*par[6]*poisson(par[1],4) + (1-5*par[6])*poisson(par[1],5))*((1-par[9])*g50 + par[9]*g51) \
        + (5*par[6]*poisson(par[1],5) + (1-6*par[6])*poisson(par[1],6))*((1-par[9])*g60 + par[9]*g61) \
        + (6*par[6]*poisson(par[1],6) + (1-7*par[6])*poisson(par[1],7))*((1-par[9])*g70 + par[9]*g71) \
        + (7*par[6]*poisson(par[1],7) + (1-8*par[6])*poisson(par[1],8))*((1-par[9])*g80 + par[9]*g81);
}