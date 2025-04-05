import numpy as np

def analytic(a1, a2, w, ux0, n1, n2, es, e, dd):
    """
    Calculate analytic solutions of a viscously anisotropic layer subjected to simple shear.
    
    Parameters:
        a1, a2: Lower and upper bound depths of the anisotropic layer.
        w: Total depth of the model.
        ux0: Horizontal velocity on the top surface.
        n1, n2: Normal vectors of the weak anisotropy.
        es, e: Weak and strong (isotropic) anisotropic viscosities, respectively.
        dd: Grid size of the depth profile.
    
    Returns:
        d: Depth grid array.
        sig11, sig12, sig22: Arrays of three stress components (sxx, sxy, syy).
        str11, str12, str22: Arrays of three strain rate components.
        u1: Array for horizontal velocity profile.
        p: Array for pressure profile.
    """
    d_w = (a2 - a1) / w
    ux0_w = ux0 / w
    
    denom = 1 - (1 - es / e) * (1 - 4 * n1**2 * n2**2) * (1 - d_w)
    nom = 1 - (1 - es / e) * (1 - 4 * n1**2 * n2**2)

    s1 = ux0_w * nom / denom
    s2 = ux0_w / denom

#    s2 = ux0_w / tmp  # p_u1/p_y
#    s1 = (ux0_w - s2 * d_w) / (1 - d_w)
    
    d = np.arange(-w, 0 + dd, dd)  # Creating the depth array
    
    # Initialize arrays
    sig11 = np.zeros_like(d)
    sig12 = np.zeros_like(d)
    sig22 = np.zeros_like(d)
    str11 = np.zeros_like(d)
    str12 = np.zeros_like(d)
    str22 = np.zeros_like(d)
    p = np.zeros_like(d)
    u1 = np.zeros_like(d)
    
    # Calculate stress, strain rate, and pressure profiles
    for i in range(len(d)):
        if a1 <= d[i] <= a2:  # Anisotropic layer
            sig11[i] = -2 * (e - es) * (n1 * n2 - 2 * n1**3 * n2) * s2
            sig12[i] = e * s2 - (e - es) * (1 - 4 * n1**2 * n2**2) * s2
            sig22[i] = -2 * (e - es) * (n1 * n2 - 2 * n1 * n2**3) * s2
            str12[i] = s2 / 2
            p[i] = 2 * (e - es) * (n1 * n2 - 2 * n1 * n2**3) * s2
        else:  # Outside anisotropic layer
            sig12[i] = e * s1
            str12[i] = s1 / 2
    
    # Calculate the horizontal velocity profile
    for i in range(1, len(d)):
        if d[i] < a1:
            u1[i] = u1[i - 1] + dd * s1
        elif a1 <= d[i] < a2:
            u1[i] = u1[i - 1] + dd * s2
        else:
            u1[i] = u1[i - 1] + dd * s1
    
    return d, sig11, sig12, sig22, str11, str12, str22, u1, p

def calc_principal(sxx, syy, sxy):
    """
    Calculate 2D principal stresses and their orientations.
    
    Parameters:
        sxx, syy, sxy: Three components of a 2D stress tensor.
    
    Returns:
        smax, smin: Maximum and minimum principal stresses.
        nx0, ny0, nx1, ny1: Unit vectors for smax and smin orientations.
        J2: Second invariant of the stress tensor.
    """
    sig = np.array([
    [sxx, sxy],
    [sxy, syy]])
    
    t = np.degrees(np.arctan2(2 * sxy, (sxx - syy))) / 2

    Q = np.array([[np.cos(np.radians(t)), np.sin(np.radians(t))],
                  [-np.sin(np.radians(t)), np.cos(np.radians(t))]])
    
    sig = np.array([[sxx, sxy],
                    [sxy, syy]])
    
    princ_s = Q @ sig @ Q.T
    
    sig1 = princ_s[0, 0]
    sig2 = princ_s[1, 1]
    
    tol = 1e-6
    if abs(princ_s[0, 1]) > tol:
        _ = princ_s[0, 1]
    
    smax = sig1
    smin = sig2
    nx0 = np.cos(np.radians(t))
    ny0 = -np.sin(np.radians(t))
    nx1 = np.sin(np.radians(t))
    ny1 = np.cos(np.radians(t))
    
    # using numpy linalg to compute eig values and vecs. 

    eigvals, eigvecs = np.linalg.eigh(sig) 
    sig1, sig2 = eigvals
    dir1 = eigvecs[:,0]
    dir2 = eigvecs[:,1] 

    smax, smin = sig1, sig2 
    dirmax, dirmin = dir1, dir2

    J2 = sxx * syy - sxy**2
    
    return smax, smin, dirmax, dirmin, J2








