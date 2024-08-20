#! /usr/bin/env python3
'''
This library contains functions used by Stokes_Mechanical_Anisotropy script. 
Last modified on 20231103. 

Function list:
  - D = D_isotropic(eta):
  - D = D_M_ti(R, eta, eta_s):
  
'''
import numpy as np
import math

np.set_printoptions(precision=2)

def D_isotropic(eta):
    # Create a 6x6 stiffness D matrix in Voigt form for isotropic anisotropy. 
    # Stress_vector = {Sxx, Syy, Szz, Sxy, Sxz, Syz}
    # Strain_rate_vector = {Exx, Eyy, Ezz, Rxy, Rxz, Ryz}, where Rxy=2Exy in the engineering strain convention. 
    
    D = np.zeros((6,6))
    D[0,0] = 2.*eta  # E_xx
    D[1,1] = 2.*eta  # E_yy
    D[2,2] = 2.*eta  # E_zz
    D[3,3] = 1.*eta  # R_yz, engineering strain rate.
    D[4,4] = 1.*eta  # R_xz
    D[5,5] = 1.*eta  # R_xy
    
    #print('Function D_isotropic is called. An isotropic Voigt stifness matrix D is created ... ...')
    return D

def D_M_ti(R, eta, eta_s):
    # The function D_M_ti returns a rotated anisotropic part of the D matrix based on Muhlhuas transverse isotropy.
    # It reads in the 3D rotation matrix R, 
    #   strong viscosity eta,
    #   weak viscosity eta_s,
    #   and returns the anisotropic part of rotated Voigt D matrix for viscosity.
    
    # The isotropic part is unchanged during rotation.
    
    # Obtain the isotropic D matrix for viscosity.
    D_iso  = D_isotropic(eta)   
    # First, create an empty 3D D matrix for Muhlhuas transverse isotropy.
    D_m_ti = np.zeros((6,6)) 
    # Second, fill in non-zero components with the normal director of the weak anisotropy aligned with y.
    # Therefore, the weak anisotropy eta_s works for R_yz and R_xy. 
    D_m_ti[0,0] = 2.*eta # E_xx
    D_m_ti[1,1] = 2.*eta # E_yy
    D_m_ti[2,2] = 2.*eta # E_zz
    D_m_ti[3,3] = eta_s  # R_yz, engineering strain rate.
    D_m_ti[4,4] = eta    # R_xz
    D_m_ti[5,5] = eta_s  # R_xy 
    
    D_m_ti_anisotropic = D_m_ti-D_iso # calculate the anisotropic part of D.
    # Calculate the rotated anisotropic part of D.
    D_m_ti_anisotropic_rot = D2Drot(D_m_ti_anisotropic,R) 
    
    return D_m_ti_anisotropic_rot
    
def D_HW_ti(R, eta, eta_s, eta_1):
    # The function D_HW_ti returns a rotated anisotropic part of the D matrix based on Han and Wahr (1997) transverse isotropy.
    # It reads in the 3D rotation matrix R,
    #   strong viscosity eta, the same as Muhlhaus formulation,
    #   weak viscosity eta_s, same as Muhlhaus formulation,
    #   eta1 for D[0,1] and D[1,0]
    #   and returns the anisotropic part of rotated Viogt D matrix for viscosity. 
    
    # Obtain the isotropic D matrix.
    D_iso  = D_isotropic(eta) 
    
    # First, create an empty 3D D matrix for Muhlhuas transverse isotropy.
    D_hw_ti = np.zeros((6,6)) 
    
    # Second, fill in non-zero components with the normal director of the weak anisotropy aligned with y.
    # Therefore, the weak anisotropy eta_s works for R_yz and R_xy. 
    D_hw_ti[0,0] = 2.*eta
    D_hw_ti[1,1] = 2.*eta
    D_hw_ti[2,2] = 2.*eta 

    # Adding eta_1
    D_hw_ti[0,0]+= eta_1
    D_hw_ti[2,2]+= eta_1
    D_hw_ti[1,1] = D_hw_ti[0,0] + eta_1
    D_hw_ti[2,0] = eta_1
    D_hw_ti[0,2] = D_hw_ti[2,0]

    # Define shear viscosity terms.
    D_hw_ti[3,3] = eta_s # R_yz, engineering strain rate.
    D_hw_ti[4,4] = eta   # R_xz
    D_hw_ti[5,5] = eta_s # R_xy
    
    D_hw_ti_anisotropic = D_hw_ti - D_iso # calculate the anisotropic part of D.
    print('Anisotropic part of D before rotation is ', D_hw_ti_anisotropic)

    # Calculate the rotated anisotropic part of D.
    D_hw_ti_anisotropic_rot = D2Drot(D_hw_ti_anisotropic,R)
    
    return D_hw_ti_anisotropic_rot
    
def getOrthorhombicD(R, eta, eta_s, eta_1, adjustedDiagNorm, adjustedOffNorm, adjustedShear):
    # The function D_ortho returns a rotated anisotropic part of the D matrix of an orthorhombic anisotropy, 
    #   which perturbs the D_HW_ti by wt*eta (wt is a percentage weight).
    
    # It reads in the 3D rotation matrix R,
    #   strong viscosity eta, the same as Muhlhaus formulation,
    #   weak viscosity eta_s, same as Muhlhaus formulation,
    #   eta1 for D[0,1] and D[1,0]
    #   and returns the anisotropic part of rotated Viogt D matrix for viscosity. 
    
    # Identity matrix. 
    R0 = np.zeros((3,3))
    R0[0,0] = R0[1,1] = R0[2,2] = 1.
    # Obtain the anisotropic part of HW TI D defined in the local coordinates with R0, i.e., no rotation.
    D_hw_ti_anisotropic = D_HW_ti(R0, eta, eta_s, eta_1)
    
    # Obtain the isotropic D.
    D_iso   = D_isotropic(eta)
    
    # wt, say wt=5%, the diagonal terms has 5% of eta perturbation.
    #wt_off  = wt*0.2 # The off-diagonal terms have 0.2*wt*eta perturbation. This is sort of arbitrary.   
    
    # Create an empty D matrix. 
    D_o_adjust = np.zeros((6,6))
    
    # the weak plane is xz. 
    # to avoid introducing tetragonal component, 
    # adjustments should be made to weak normal directions x, z.
    # Dxx, Dzz, D23, D21, Dyx(D66), Dyz(D44)
    D_o_adjust[0,0] -= adjustedDiagNorm # weaken sxx 
    #D_o_adjust[1,1] += 0.0 # no change to hexagonal axis normal component syy
    D_o_adjust[2,2] += adjustedDiagNorm # strengthen szz 
    D_o_adjust[0,1] += adjustedOffNorm
    D_o_adjust[1,0] = D_o_adjust[0,1]
    D_o_adjust[1,2] -= adjustedOffNorm 
    D_o_adjust[2,1] = D_o_adjust[1,2]
    D_o_adjust[3,3] -= adjustedShear 
    #D_o_adjust[4,4] no change
    D_o_adjust[5,5] += adjustedShear
    
    D_o_anisotropic = D_hw_ti_anisotropic + D_o_adjust
    
    print('Anisotropic part of D before rotation is ')
    print(D_o_anisotropic)
 
    D_o_anisotropic_rot = D2Drot(D_o_anisotropic,R)
    
    print('Rotated anisotropic part of D for the orthorhombic case is')
    print(D_o_anisotropic_rot)
    return D_o_anisotropic_rot

def M_ti_dd(eta, eta_s, norm1, norm2, norm3, dim):
    # Function M_ti_dd follows the formulation in Muhlhaus et al. (2002) to produce the 4th order tensor stiffness dd. 
    if dim == 2:
        n = np.zeros((1,dim))
        n[0][0] = norm1
        n[0][1] = norm3       
    elif dim == 3:
        n = np.zeros((1,dim))
        n[0][0] = norm1
        n[0][1] = norm2
        n[0][2] = norm3
    
    print(n)
    
    C_lambda = np.zeros((dim, dim, dim, dim))

    componentList = []
    for i in range(dim):
        componentList += [[],]
        for j in range(dim):
            componentList[i] += [[],]
            for k in range(dim):
                componentList[i][j] += [[],]
                for l in range(dim):
                    componentList[i][j][k] += [[],]    

    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                for l in range(dim):
                    a = 0
                    b = 0
                    c = 0
                    d = 0
                    e = n[0][i]*n[0][j]*n[0][k]*n[0][l]
                    f = 0
                    if l == j :
                        a = n[0][i]*n[0][k]
                    if i == l :
                        b = n[0][j]*n[0][k]
                    if k == j :
                        c = n[0][i]*n[0][l]
                    if i == k :
                        d = n[0][j]*n[0][l]
                        
                    if i==k and j==l:
                        f = 1
                        
                    componentList[i][j][k][l] = 2*(eta_s - eta) * ((a + b + c + d)/2 - 2*e)
                    C_lambda[i][j][k][l] = 2*(eta_s-eta)*((a + b + c + d)/2 - 2*e)
    #I4 = fourth_order_identity_tensor()
    #DI = np.eyes(6)*eta*2.
    #I = viogt2tensor(DI)
    C =  C_lambda
    # 4th order tensor of material property matrix C 
    # C = as_tensor(componentList) 
    #vtmp = as_tensor(epsilon(v))
    #i, j, k, l = ufl.indices(4)
    #C1 = ufl.as_tensor(C[i,j,k,l]*vtmp[k,l],(i,j))
    #C2 = 2*eta*sym(nabla_grad(v))
    #C3 = C1 + C2
    return C
    
def D2dd(D):
    # For 3D only.
    # The function D2dd converts stiffness matrix D(6x6) in Voigt form to the 4th order tensor dd(3x3x3x3)
    #   following Browaeys (2004) Table 1.
    
    # The mapping between the Voigt matrix indices and the tensor indices are:
    # voigt_map = {0: (0, 0), 1: (1, 1), 2: (2, 2), 3: (1, 2), 4: (0, 2), 5: (0, 1)}
    # voigt_map = {0: (x, x), 1: (y, y), 2: (z, z), 3: (y, z), 4: (x, z), 5: (x, y)}

    # # Fill in the tensor
    # for i in range(6):
        # for j in range(6):
            # p, q = voigt_map[i]
            # r, s = voigt_map[j]
            # tensor[p][q][r][s] = voigt_matrix[i][j]
            # tensor[q][p][r][s] = voigt_matrix[i][j]
            # tensor[p][q][s][r] = voigt_matrix[i][j]
            # tensor[q][p][s][r] = voigt_matrix[i][j]
    # return tensor
   dd = np.zeros((3,3,3,3))

   dd[0,0,0,0] = D[0,0] # c(1), c1111, C11
   dd[1,1,1,1] = D[1,1] # c(2), c2222, C22
   dd[2,2,2,2] = D[2,2] # c(3), c3333, C33
   dd[1,1,2,2] = D[1,2] # c(4), c2233=c3322, C23=C32
   dd[2,2,1,1] = D[1,2]
   dd[0,0,2,2] = D[0,2] # c(5), c1133=c3311, C13=C31
   dd[2,2,0,0] = D[0,2] 
   dd[0,0,1,1] = D[0,1] # c(6), c1122=c2211, C12=C21
   dd[1,1,0,0] = D[0,1]
   
   dd[1,2,1,2] = D[3,3] # c(7), c2323=c2332=c3223=c3232, C44
   dd[1,2,2,1] = D[3,3]
   dd[2,1,1,2] = D[3,3]
   dd[2,1,2,1] = D[3,3]
   
   dd[0,2,0,2] = D[4,4] # c(8), c1313=c1331=c3113=c3131, C55
   dd[0,2,2,0] = D[4,4]
   dd[2,0,0,2] = D[4,4]
   dd[2,0,2,0] = D[4,4]
   
   dd[0,1,0,1] = D[5,5] # c(9), c1212=c1221=c2112=c2121, C66 
   dd[0,1,1,0] = D[5,5]
   dd[1,0,0,1] = D[5,5] 
   dd[1,0,1,0] = D[5,5] 
   
   dd[0,0,1,2] = D[0,3] # c(10), c1123=c1132=c2311=c3211, C14=C41
   dd[0,0,2,1] = D[0,3]
   dd[1,2,0,0] = D[0,3]
   dd[2,1,0,0] = D[0,3]
   dd[1,1,0,2] = D[1,4] # c(11), c2213=c2231=c1322=c3122, C25=C52
   dd[1,1,2,0] = D[1,4]
   dd[0,2,1,1] = D[1,4]
   dd[2,0,1,1] = D[1,4]
   dd[2,2,0,1] = D[2,5] # c(12), c3312=c3321=c1233=c2133, C36=C63
   dd[2,2,1,0] = D[2,5]
   dd[0,1,2,2] = D[2,5]
   dd[1,0,2,2] = D[2,5]
   dd[2,2,1,2] = D[2,3] # c(13), c3323=c3332=c2333=c3233, C34=C43
   dd[2,2,2,1] = D[2,3]
   dd[1,2,2,2] = D[2,3]
   dd[2,1,2,2] = D[2,3] 
   dd[0,0,0,2] = D[0,4] # c(14), c1113=c1131=c1311=c3111, C15=C51
   dd[0,0,2,0] = D[0,4]
   dd[0,2,0,0] = D[0,4]
   dd[2,0,0,0] = D[0,4]
   dd[1,1,0,1] = D[1,5] # c(15), c2212=c2221=c1222=c2122, C26=C62
   dd[1,1,1,0] = D[1,5]
   dd[0,1,1,1] = D[1,5]
   dd[1,0,1,1] = D[1,5]
   dd[1,1,1,2] = D[1,3] # c(16), c2223=c2232=c2322=c3222, C24=C42
   dd[1,1,2,1] = D[1,3]
   dd[1,2,1,1] = D[1,3]
   dd[2,1,1,1] = D[1,3]
   dd[2,2,0,2] = D[2,4] # c(17), c3313=c3331=c1333=c3133, C35=C53
   dd[2,2,2,0] = D[2,4]
   dd[0,2,2,2] = D[2,4]
   dd[2,0,2,2] = D[2,4]
   dd[0,0,0,1] = D[0,5] # c(18), c1112=c1121=c1211=c2111, C16=C61
   dd[0,0,1,0] = D[0,5]
   dd[0,1,0,0] = D[0,5]
   dd[1,0,0,0] = D[0,5]
   dd[0,2,0,1] = D[4,5] # c(19), c1312=c1321=c3112=c3121=c1213=c1231=c2113=c2131, C56=C65
   dd[0,2,1,0] = D[4,5]
   dd[2,0,0,1] = D[4,5]
   dd[2,0,1,0] = D[4,5]
   dd[0,1,0,2] = D[4,5]
   dd[0,1,2,0] = D[4,5]
   dd[1,0,0,2] = D[4,5]
   dd[1,0,2,0] = D[4,5]
   dd[1,2,0,1] = D[3,5] # c(20), c2312=c2321=c3212=c3221=c1223=c1232=c2123=c2132, C46=C64
   dd[1,2,1,0] = D[3,5]
   dd[2,1,0,1] = D[3,5]
   dd[2,1,1,0] = D[3,5]
   dd[0,1,1,2] = D[3,5]
   dd[0,1,2,1] = D[3,5]
   dd[1,0,1,2] = D[3,5]
   dd[1,0,2,1] = D[3,5]
   dd[1,2,0,2] = D[3,4] # c(21), c2313=c2331=c3213=c3231=c1323=c1332=c3123=c3132, C45=C54
   dd[1,2,2,0] = D[3,4]
   dd[2,1,0,2] = D[3,4]
   dd[2,1,2,0] = D[3,4]
   dd[0,2,1,2] = D[3,4]
   dd[0,2,2,1] = D[3,4]
   dd[2,0,1,2] = D[3,4]
   dd[2,0,2,1] = D[3,4]  
   
   #print(dd)
   #print('Done converting Voigt matrix to 4th order tensor ... ...')
   return dd

def dd2D(dd):
    # For 3D only.
    # The function dd2D converts the 4th order tensor dd(3x3x3x3) to stiffness matrix D(6x6) in Voigt form 
    #   following Browaeys 2004 Table 1. 
    
    # The mapping between the Voigt matrix indices and the tensor indices are:
    # voigt_map = {0: (0, 0), 1: (1, 1), 2: (2, 2), 3: (1, 2), 4: (0, 2), 5: (0, 1)}
    # voigt_map = {0: (x, x), 1: (y, y), 2: (z, z), 3: (y, z), 4: (x, z), 5: (x, y)}
  
    # # Create an empty Voigt matrix
    # voigt_matrix = np.zeros((6, 6))

    # # Fill in the Voigt matrix
    # for i in range(3):
        # for j in range(3):
            # for k in range(3):
                # for l in range(3):
                    # voigt_matrix[voigt_map[(i, j)]][voigt_map[(k,l)]] = tensor[i][j][k][l]
    # return voigt_matrix
    #print('Converting 4th order tensor C to Voigt matrix ... ...')
    #print('Original C is', C)
    
    D = np.zeros((6,6))
    
    D[0,0] = dd[0][0][0][0] # c(1), c1111, C11
    D[1,1] = dd[1][1][1][1] # c(2), c2222, C22
    D[2,2] = dd[2][2][2][2] # c(3), c3333, C33
    D[1,2] = dd[1][1][2][2] # c(4), c2233=c3322, C23=C32
    D[2,1] = D[1,2] 
    D[0,2] = dd[0][0][2][2] # c(5), c1133=c3311, C13=C32
    D[2,0] = D[0,2]
    D[0,1] = dd[0][0][1][1] # c(6), c1122=c2211, C12=C21
    D[1,0] = D[0,1] 
    D[3,3] = dd[1][2][1][2] # c(7), c2323=c2332=c3223=c3232, C44
    D[4,4] = dd[0][2][0][2] # c(8), c1313=c1331=c3113=c3131, C55
    D[5,5] = dd[0][1][0][1] # c(9), c1212=c1221=c2112=c2121, C66
    D[0,3] = dd[0][0][1][2] # c(10), c1123=c1132=c2311=c3211, C14=c41
    D[3,0] = D[0,3]
    D[1,4] = dd[1][1][0][2] # c(11), c2213=c2231=c1322=c3122, C25=C52
    D[4,1] = D[1,4] 
    D[2,5] = dd[2][2][0][1] # c(12), c3312=c3321=c1233=c2133, C36=C63 
    D[5,2] = D[2,5]
    D[2,3] = dd[2][2][1][2] # c(13), c3323=c3332=c2333=c3233, C34=C43 
    D[3,2] = D[2,3]
    D[0,4] = dd[0][0][0][2] # c(14), c1113=c1131=c1311=c3111, C15=C51
    D[4,0] = D[0,4]
    D[1,5] = dd[1][1][0][1] # c(15), c2212=c2221=c1222=c2122, C26=C62
    D[5,1] = D[1,5]
    D[1,3] = dd[1][1][1][2] # c(16), c2223=c2232=c2322=c3222, C24=C42
    D[3,1] = D[1,3]
    D[2,4] = dd[2][2][0][2] # c(17), c3313=c3331=c1333=c3133, C35=C53
    D[4,2] = D[2,4]    
    D[0,5] = dd[0][0][0][1] # c(18), c1112=c1121=c1211=c2111, C16=C61
    D[5,0] = D[0,5]
    D[4,5] = dd[0][2][0][1] # c(19), c1312=c1321=c3112=c3121
                                     #=c1213=c1231=c2113=c2131, C56=C65
    D[5,4] = D[4,5]
    D[3,5] = dd[1][2][0][1] # c(20), c2312=c2321=c3212=c3221=c
                                     #=c1223=c1232=c2123=c2132, C46=C64
    D[5,3] = D[3,5]
    D[3,4] = dd[1][2][0][2] # c(21), c2313=c2331=c3213=c3231=c1323
                                     # = c1332=c3123=c3132,     C45=C54
    D[4,3] = D[3,4]
    #print('Done converting and Voigt matrix is', D)
    
    return D

def Rx(theta):
    # Rx returns a 3D rotation matrix, which will rotate points in xyz around axis x counterclockwisely through an angle of theta.
    # Right-hand rule should be applied.
    Rx = np.zeros((3,3))
    Rx[0,0] = 1.
    Rx[1,1] = np.cos(theta/180.*math.pi)
    Rx[2,2] = np.cos(theta/180.*math.pi)
    Rx[1,2] = -np.sin(theta/180.*math.pi)
    Rx[2,1] = np.sin(theta/180.*math.pi)
    return Rx
    
def Ry(theta):
    # Ry returns a 3D rotation matrix, which will rotate points in xyz around axis y ''''counterclockwisely''' through an angle of theta.
    # Or to rotate the coordinates xyz around the axis x '''clockwisely''' through the angle of theta.
    # Right-hand rule should be applied.
    
    Ry = np.zeros((3,3))
    Ry[0,0] = np.cos(theta/180.*math.pi)
    Ry[1,1] = 1.
    Ry[2,2] = np.cos(theta/180.*math.pi)
    Ry[0,2] = np.sin(theta/180.*math.pi)
    Ry[2,0] = -np.sin(theta/180.*math.pi)
    return Ry

def Rz(theta):
    # Rz returns a 3D rotation matrix, which will rotate points in xyz around axis z '''counterclockwisely''' through an angle of theta.
    # Or to rotate the coordinates xyz around the axis z '''clockwisely''' through the angle of theta.
    # Right-hand rule should be applied.
    Rz = np.zeros((3,3))
    Rz[0,0] = np.cos(theta/180.*math.pi)
    Rz[1,1] = np.cos(theta/180.*math.pi)
    Rz[2,2] = 1.
    Rz[0,1] = -np.sin(theta/180.*math.pi)
    Rz[1,0] = np.sin(theta/180.*math.pi)
    return Rz
    
def A2D(A,R):
    # For 3D only. 
    # Given the 6x6 compliance matrix A in Voigt form, rotation matrix R (3x3), 
    #   then return the 6x6 stiffness matrix D in Voigt form. 
        
    aa     = voigt2tensor(A)      # convert Voigt A to the 4th order tensor aa (3x3x3x3).
    aa_rot = rotate_4th_tensor(aa,R)   # rotate aa given R. aa_rot = RRRRaa.
    A_rot  = tensor2voigt(aa_rot) # convert rotated 4th order tensor aa_rot to 6x6 Voigt form A_rot. 
    D      = np.linalg.pinv(A_rot)# calculate the inverse of A_rot to get the stiffness matrix D.  
    return D

def D2Drot(D,R):
    # For 3D only.
    # Given the stiffness matrix D (6x6) in Voigt form, rotation matrix R (3x3), 
    #   then return the rotated stiffness matrix Drot in Voigt form. 
        
    dd     = D2dd(D)     # convert Voigt D to 4th order tensor dd. 
    dd_rot = rotate_4th_tensor(dd,R)  # rotate dd based on R. dd_rot = RRRRdd.
    Drot   = dd2D(dd_rot) # convert dd_rot back to Voigt form Drot. 
    
    return Drot
    
def rotate_4th_tensor(dd,R):
    # For 3D only, 
    # The function 'rotate_4th_tensor' rotate the 4th order tensor by the formulation
    #   dd_rot = RRRRdd.
    # dd(3x3x3x3) is the 4th order tensor of stiffness matrix/compliance matrix.
    # R(3x3) is the rotation matrix. 
        
    # for i1 in range(3):
        # for i2 in range(3):
            # for i3 in range(3):
                # for i4 in range(4):
    return np.einsum('ia,jb,kc,ld,abcd->ijkl',R,R,R,R,dd)
    
def check_tensor_symmetry(T):
    major_symmetry = True
    minor_symmetry = True

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if T[i][j][k][l] != T[k][l][i][j]:
                        major_symmetry = False
                    if T[i][j][k][l] != T[j][i][k][l] or T[i][j][k][l] != T[i][j][l][k]:
                        minor_symmetry = False

    full_symmetry = major_symmetry and minor_symmetry

    return major_symmetry, minor_symmetry, full_symmetry
'''
def D_HW_ti_30():
    ##################################
    # The following case is the rotated 30 degs anisotropic part of D.
    # Following Han and Wahr to add eta_1 and modified on top of the above
    # Maulhaus case. 
    # with the isotropic np.eye(6)*2. removed.
    # eta, eta_s = 1. and 0.1
    # eta_1 = 0.3
    # [Conclusion] This version works.
    anisoC = np.zeros((6,6))
    anisoC[0,0] = -0.375
    anisoC[1,1] = -0.375
    anisoC[2,2] = 0.6
    anisoC[0,1] = 0.975
    anisoC[1,0] = 0.975
    anisoC[0,5] = 0.39
    anisoC[5,0] = 0.39
    anisoC[1,5] = -0.39
    anisoC[5,1] = -0.39
    anisoC[3,3] = -0.675
    anisoC[4,4] = -0.225
    anisoC[5,5] = -0.225
    anisoC[3,4] = 0.39
    anisoC[4,3] = 0.39
    return anisoC
    
def D_M_ti_30():
    # The following case is the rotated 30 degs anisotropic part of D.
    # Following Maulhaus ti with eta = 1., and eta_s = 0.1
    # The isotropic part, np.eye(6)*2, is removed.
    # [Conclusion] This version works.
    anisoC = np.zeros((6,6))
    anisoC[0,0] = -0.675
    anisoC[1,1] = -0.675
    anisoC[2,2] = 0.0
    anisoC[0,1] = 0.675
    anisoC[1,0] = 0.675
    anisoC[0,5] = 0.39
    anisoC[5,0] = 0.39
    anisoC[1,5] = -0.39
    anisoC[5,1] = -0.39
    anisoC[3,3] = -0.675
    anisoC[4,4] = -0.225
    anisoC[5,5] = -0.225
    anisoC[3,4] = 0.39
    anisoC[4,3] = 0.39
    return anisoC
    
def D_ortho_30():
    ###################################
    # The following case is the rotated 30 degs anisotropic part of D.
    # Following the above Han and Wahr scenario but also adding 
    # an orthotropic adjustment. 
    #
    anisoC = np.zeros((6,6))
    anisoC[0,0] = -0.366
    anisoC[1,1] = -0.341
    anisoC[2,2] = 0.65
    anisoC[0,1] = 0.928
    anisoC[1,0] = anisoC[0,1]
    anisoC[0,2] = 0.005
    anisoC[2,0] = anisoC[0,2]
    anisoC[0,5] = 0.352
    anisoC[5,0] = anisoC[0,5]
    anisoC[1,2] = -0.005
    anisoC[2,1] = anisoC[1,2]
    anisoC[1,5] = -0.373
    anisoC[5,1] = anisoC[1,5]
    anisoC[2,5] = 0.009
    anisoC[5,2] = anisoC[2,5]
    anisoC[3,3] = -0.712
    anisoC[4,4] = -0.237
    anisoC[5,5] = -0.222
    anisoC[3,4] = 0.411
    anisoC[4,3] = anisoC[3,4]
    return anisoC
'''

