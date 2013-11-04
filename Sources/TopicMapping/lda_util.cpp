
// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)

// This file is part of LDA-C.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.




/*
 * given log(a) and log(b), return log(a + b)
 *
 */


double log_sum(double log_a, double log_b) {
    double v;
    
    if (log_a < log_b)
        {
        v = log_b+log(1 + exp(log_a-log_b));
        }
    else
        {
        v = log_a+log(1 + exp(log_b-log_a));
        }
    return(v);
}




double trigamma(double x)
{
    double p;
    int i;
    
    x=x+6;
    p=1/(x*x);
    p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)
         *p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;
    for (i=0; i<6 ;i++)
        {
        x=x-1;
        p=1/(x*x)+p;
        }
    return(p);
}


/*
 * taylor approximation of first derivative of the log gamma function
 *
 */

double digamma(double x)
{
    double p;
    x=x+6;
    p=1/(x*x);
    p=(((0.004166666666667*p-0.003968253986254)*p+
        0.008333333333333)*p-0.083333333333333)*p;
    p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
    return p;
}




double invert_digamma(double y) {
    
    // returns x suxh that digamma(x)=y
    
    double xold=0.;
    if(y>-2.22){
        xold=exp(y)+0.5;
    } else {
        xold= -1./(y -digamma(1));
    }
    
    double xnew;
    for(int iter=0; iter<5; iter++) {
        xnew= xold - (digamma(xold)-y)/(trigamma(xold));
        if(fabs(xnew-xold)<1e-6) {
            return xnew;
        }
    }
    return xnew;
}

