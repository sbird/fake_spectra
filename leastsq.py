# -*- coding: utf-8 -*-
"""
Module for computing the 5 least squares regression methods detailed in
        Linear regression in astronomy.
        Isobe, Takashi; Feigelson, Eric D.; Akritas, Michael G.; Babu, Gutti Jogesh
        Astrophysical Journal, Part 1 (ISSN 0004-637X), vol. 364, Nov. 20, 1990, p. 104-113
        http://adsabs.harvard.edu/abs/1990ApJ...364..104I

These methods are appropriate when the intrinsic scatter in the data is much larger
than the error on each data point.
"""

import numpy as np

def leastsq(x,y, method=5):
    """
       Compute the least squares fit to y = beta x + alpha,
       using one of the 5 methods outlined in
       http://adsabs.harvard.edu/abs/1990ApJ...364..104I
       Method 1 minimises distance from Y given X (ie, the standard least squares fit)
       Method 2 minimises distance from X given Y
       Method 3 (recommended) is the OLS bisector, which gives a line bisecting the above two.
       Method 4 (Orthogonal regression) minimises perpendicular distance from the line to points
       Method 5 is the geometric mean of the slopes from methods 1 and 2.
       Method 6 is the Theil-Sen estimator: the median of the pairwise slopes.
       (See Akritas 95,  http://www.tandfonline.com/doi/abs/10.1080/01621459.1995.10476499)
       Returns:
              (alpha, beta, bvar), the intercept slope and variance of the slope
    """
    #Define some sums
    xbar = np.mean(x)
    ybar = np.mean(y)
    xdif = x-xbar
    ydif = y-ybar
    sxx = np.sum(xdif**2)
    syy = np.sum(ydif**2)
    sxy = np.sum(ydif*xdif)

    #Check for zeros
    if sxx == 0 or syy == 0 or sxy == 0:
        raise ValueError("Least Squares ill-defined")
    if method > 6 or method < 1:
        raise ValueError("Method not recognised")

    #These formulas are taken from Table 1 of Isobe et al, page 3
    #Minimise distance from Y given X
    beta1 = sxy/sxx
    #Variance of b1
    bvar1 = np.sum(xdif**2*(ydif-beta1*xdif)**2)/sxx**2
    #Minimise distance from X given Y
    beta2 = syy/sxy
    #Variance of b2
    bvar2 = np.sum(ydif**2*(ydif-beta2*xdif)**2)/sxy**2
    #Covariance of b1 and b2
    covb12 = np.sum(xdif*ydif*(ydif-beta2*xdif)*(ydif-beta1*xdif))/(beta1*sxx**2)

    if method == 1:
        beta = beta1
        bvar = bvar1
    if method == 2:
        beta = beta2
        bvar = bvar2
    if method == 3:
        #OLS bisector: line that bisects the above two.
        beta1p1 = 1+beta1**2
        beta2p1 = 1+beta2**2
        beta = (beta1*beta2 - 1 + np.sqrt(beta1p1*beta2p1))/(beta1+beta2)
        #Variance
        prefac = beta**2 / ( (beta1 + beta2)**2 * beta1p1 * beta2p1)
        var = beta2p1**2 * bvar1 + 2 * beta1p1 * beta2p1 * covb12 + beta1p1**2 * bvar2
        bvar = prefac*var

    if method == 4:
        #Orthogonal: minimise perpendicular distance from line to points
        beta = 0.5*((beta2-1./beta1)+np.sign(sxy)*np.sqrt(4+(beta2-1./beta1)**2))
        prefac = beta**2 / (4*beta1**2 + (beta1*beta2 - 1)**2)
        bvar = prefac * ( bvar1/beta1**2 + 2*covb12 + beta1**2*bvar2 )

    if method == 5:
        #Reduced major axis:
        beta = np.sign(sxy)*np.sqrt(beta1*beta2)
        bvar = 0.25 * (beta2/beta1 * bvar1 + 2*covb12 + beta1/beta2 * bvar2)

    if method == 6:
        #Theil-Sen estimator for uncensored data: the median of the slopes.
        yy = np.subtract.outer(y,y)
        xx = np.subtract.outer(x,x)
        ind = np.where(xx != 0)
        beta = np.median(yy[ind]/xx[ind])
        #Can't find a formula for the variance
        bvar = 0

    #The intercept
    alpha = ybar - beta*xbar

    return (alpha, beta, bvar)

def pearson(x,y,alpha, beta, method=3):
    """Find the residual squares divided by total squares.
       ie, the Pearson correlation coefficient between the fit and the data.

       If method == 1, return
          r_y^2 = Σ ( y_i - f_i(x) )^2 / Σ ( y_i - ybar )^2
       For method == 2, substitute x for y above.
       For method == 3, symmetrise by taking the geometric mean of method 1 and 2,
          r_xy = (r_y*r_x)**2
    """
    xbar = np.mean(x)
    ybar = np.mean(y)
    #Vector of expected y from fit
    fity = beta*x + alpha
    #Vector of expected x from fit
    fitx = (y - alpha) / beta
    #Scatter from y axis: method 1 minimises this.
    yscat = np.sum((y-fity)**2)
    yvar = np.sum((y-ybar)**2)

    if method == 1:
        return np.sqrt( yscat / yvar )

    #Scatter from x axis: method 2 minimises this.
    xscat = np.sum((x-fitx)**2)
    xvar = np.sum((x-xbar)**2)

    if method == 2:
        return np.sqrt( xscat / xvar )

    if method == 3:
        return np.sqrt( np.sqrt( yscat / yvar ) * np.sqrt( xscat / xvar ) )
