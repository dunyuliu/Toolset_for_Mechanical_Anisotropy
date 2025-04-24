import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
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
    sig = np.array([[sxx, sxy],
                    [sxy, syy]])
    
    # using numpy linalg to compute eig values and vecs. 
    eigvals, eigvecs = np.linalg.eigh(sig) 
    sig1, sig2 = eigvals
    dir1 = eigvecs[:,0]
    dir2 = eigvecs[:,1] 

    smax, smin = sig1, sig2 
    dirmax, dirmin = dir1, dir2

    J2 = np.abs(sxy) # eq (24a,b), only works for strain rate

    return smax, smin, dirmax, dirmin, J2

def plot_vectors_and_save(theta_n_deg,
                          theta_sigma_deg,
                          theta_eps_deg,
                          strain_enhance=1.0,
                          viscosity_contrast=1.0,
                          outfile="vectors_n_sigma_eps.png",
                          dpi=300):
    """
    Plot unit vectors n, sigma, ε given their CCW angles from +y axis,
    draw the corresponding angle arcs, and save the diagram with strain enhancement values.
    """

    # ---- helper -----------------------------------------------------------
    def unit(angle_deg):
        """Return (x,y) for a unit vector CCW from +y."""
        rad = np.deg2rad(angle_deg)
        return -np.sin(rad), np.cos(rad)
    def unit_mirror(angle_deg):
        """Return (x,y) for a unit vector that mirrors the original vector."""
        rad = np.deg2rad(angle_deg)
        return np.sin(rad), -np.cos(rad)

    def arc(angle_deg, r, npts=120):
        """Return (x,y) arrays for an arc from 90 to 90+angle_deg (°)."""
        # start at 0 rad (+y), sweep by sign(angle_deg)
        phi = np.linspace(np.deg2rad(90), np.deg2rad(90+angle_deg), npts)
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        return x, y
    # ----------------------------------------------------------------------

    # vector spec: label → (angle, color, arc_radius)
    specs = {
        'n'           : (theta_n_deg,     'tab:blue',   0.22),
        r'\sigma'     : (theta_sigma_deg, 'tab:green',  0.32),
        r'\varepsilon': (theta_eps_deg,   'tab:red',    0.42)
    }

    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_aspect('equal')

    # axes
    ax.axhline(0, color='k', lw=1)
    ax.axvline(0, color='k', lw=1)
    ax.text(1.2,0,'x',ha='left',va='center',fontsize=13)
    ax.text(0,1.2,'y',ha='center',va='bottom',fontsize=13)

    # +y dashed reference
    ax.plot([0,0],[0,1],'--',color='gray',lw=1)

    # loop over vectors
    for label,(ang,color,rad) in specs.items():
        vec = unit(ang)
        vec_mirror = unit_mirror(ang)
        # arrow
        ax.quiver(0,0,*vec,angles='xy',scale_units='xy',
                  scale=1,color=color,lw=2)
        ax.quiver(0,0,*vec_mirror,angles='xy',scale_units='xy',
                  scale=1,color=color,lw=2)
        ax.text(vec[0]*1.08, vec[1]*1.08,
                rf'$\mathbf{{{label}}}$',color=color,fontsize=13)

        # arc
        xarc,yarc = arc(ang, rad)
        ax.plot(xarc,yarc,color=color,lw=1.4)
        mid = np.deg2rad(ang)/2
        ax.text(-1.15*rad*np.sin(mid), 1.15*rad*np.cos(mid),
                rf'$\theta_{{{label}}}$', color=color, fontsize=12,
                ha='center', va='center')

    # limits / cosmetics
    lim = 1.35
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"Viscosity contrast{viscosity_contrast}; θ={theta_n_deg:.2f} CCW from y+); strain enhance={strain_enhance:.2f}",fontsize=12)
    plt.tight_layout()

    # ----- angle guides ----------------------------------------------------
    # radial spokes every 5°
    for ang in np.arange(0, 360, 5):
        a = np.deg2rad(ang)
        ax.plot([0, lim*np.sin(a)],  # x   (CCW from +y)
                [0, lim*np.cos(a)],  # y
                ls='--', lw=0.6, color='lightgray', zorder=0)

    # optional concentric circles every 0.25
    radii = np.arange(0.25, lim+0.25, 0.25)
    for r in radii:
        circle = plt.Circle((0,0), r, ls='--', lw=0.6,
                            color='lightgray', fill=False, zorder=0)
        ax.add_patch(circle)

    # Cartesian grid
    ax.set_xticks(np.arange(-lim, lim+0.1, 0.25))
    ax.set_yticks(np.arange(-lim, lim+0.1, 0.25))
    ax.grid(color='gainsboro', linestyle=':', linewidth=0.5)
    # -----------------------------------------------------------------------

    # save
    outfile = Path(outfile)
    fig.savefig(outfile,dpi=dpi,bbox_inches='tight')
    print(f"Diagram saved to {outfile.resolve()}")
    
    return fig






