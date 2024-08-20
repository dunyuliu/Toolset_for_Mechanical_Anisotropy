
# In[ ]:


# Summary: the code is FEniCS-based implementation of Stokes and contiuum equations 
#   with mechanical anisotropy described in Maulhaus et al. (2002) and general anisotropy implementation by feeding in stiffness
#   D in Voigt form. 
# Part of the code is adopted from the FEniCS Stokes-iterative demo, whose link is provided below:
#    https://fenicsproject.org/olddocs/dolfin/1.3.0/python/demo/documented/stokes-iterative/python/documentation.html

# Now version 1.1.12.
# 20231103: improvement for case 35: solve with CG3+DG2 for 3D FEniCS mesh and project results to DG0 removes the wiggles in pressure and stress.
#   Also, Stokes_Mechanical_Anisotropy_Library.py is used to host key functions with full documentation.
# 20230919: adding case 34 with misoriented shear zone subjected to shortening for various viscous rheology. 
#   Adding case 35 for simple shearing with normal directors rotated - a 3D case for the case 26. 
# 20230125: using gmsh mesh for case 33. Cannot use MPI. Add a controlling parameter useMPI (True/False)
# updated on 01/20/2023: add the isotropically weak shear zone for case 33. Need additional parameter iso=True/False as input.
# Created on 08/02/2021. Last modfied on 06/29/2021. 
# Author: Dunyu Liu (dliu@ig.utexas.edu) 


# In[ ]:


from dolfin import *
import numpy as np
from lib import *
from prepare_case import *
import os, time, sys
from Stokes_Mechanical_Anisotropy_Library import * 

np.set_printoptions(precision=2)

#useMPI = False

st  = time.time()
#if useMPI == True:
#  parameters["ghost_mode"] = "shared_facet"

#os.environ['OMP_NUM_THREADS'] = '1'

#get_ipython().run_line_magic('matplotlib', 'inline')
#sys.setrecursionlimit(1500)

# Test for PETSc or Epetra
if not has_linear_algebra_backend("PETSc") and not has_linear_algebra_backend("Epetra"):
    info("DOLFIN has not been configured with Trilinos or PETSc. Exiting.")
    exit()

if not has_krylov_solver_preconditioner("amg"):
    info("Sorry, this code is only available when DOLFIN is compiled with AMG "
	 "preconditioner, Hypre or ML.");
    exit()


# In[ ]:


case = 35
case_info = { 23:'2D fossile shear zone.',
              26:'2D rectangle with a horizontal anisotropic layer. Same as the analytic solution.',
              30:'3D box with a horziontal anisotropic layer.',
              31:'3D box with vertical SAF fault zone parallel to boundaries.',
              32:'3D Leech River Schist above the Cascadia subduction zone.',
              33:'3D fossil mantle shear zone.',
              34:'3D fossil mantle shear zone with general anisotropy.',
              35:'3D vertical shear zone with normal director rotated and subject to simple shearing.'}
# case 26 for 2D rectangle with a horizontal anisotropic layer. Same as the analytic solution.
# case 30 for 3D box with a horziontal anisotropic layer.
# case 31 for 3D box with vertical SAF fault zone parallel to boundaries.
# case 32 for 3D Leech River Schist above the Cascadia subduction zone. 
# case 33 for 3D fossil mantle shear zone. 
# case 34 for 3D fossil mantle shear zone with general anisotropy.
# case 35 for 3D vertical shear zone with normal director rotated and subject to simple shearing. 

# date_stamp = '20221014_C100_n30_GaussianWeakViscosity' # Date stamp is used in the output folder name.
date_stamp = '20221025_ModelSize211_C10_n20_FossilMantleShearZone_FreeNorth' # Date stamp is used in the output folder name.


# In[ ]:


def stokes(delta, alpha, eta_strong, eta_weak, msh_path = None, msh_name = None, n = 30, args = None):
    
    # Use function prepare_case to generate key model specific data.
    mesh, boundaries, mf, bcs, f, bound_name_list, dim, W, norm1, norm2, norm3 = prepare_case(case, msh_path, msh_name, n, delta, alpha, args)
    #------- Case specific modifications end here --------
    
    anis_domain = bound_name_list['anis_domain']
    iso_domain = bound_name_list['iso_domain']

    g = Expression(('0', '-5', '0'), degree=2) # Neumann boundary condition
    
    ds = Measure("ds")(domain=mesh, subdomain_data = boundaries)
    dx = Measure('dx')(domain=mesh, subdomain_data = mf)

    if eta_strong>eta_weak:
        hete = 1
        print(hete)
	#iso = False
    elif eta_strong == eta_weak: 
        hete = 0
        print(hete)
        #iso = True

    # Define eta and eta_s.
    # If iso==False, eta is the strong viscosity in the anisotropic cases or the viscosity of isotropic case.
    #                eta_s is the weak viscosity in the anisotropic cases.
    # If iso==True,  eta is the viscosity of both anisotropic and isotropic domains but eta_weak will be assigned to the 'anisotropic' layer.
    if iso == True:
        eta = K(mf, eta_weak, eta_strong, anis_domain, iso_domain, degree=0)
    else:
        eta = eta_strong
    #eta_s = K(mf, eta_weak, eta_strong, anis_domain, iso_domain, degree=0) 
    
    # Assign stress_coef to be 1 in the anis_domain and 0 in the iso_domain
    # When calculating stress, only the anis_domain will need the anisotropic adjustment.
    stress_coef = K(mf, 1, 0, anis_domain, iso_domain, degree=0)
    
    if case==26 and guassian==True:
        eta_s = Expression('1-(eta0-eta1)*exp(-(pow((x[1]-zc),2)/(pow(std,2))))', eta0=eta_strong, eta1=eta_weak, zc=0.7, std=0.1, degree=5)
    else:
        eta_s = K(mf, eta_weak, eta_strong, anis_domain, iso_domain, degree=0)
       
    # Define functions -----------------------------------
    def epsilon(v):
        return sym(nabla_grad(v))
    # Define stress
    def sigma(v,p,eta):
        return 2*eta*sym(nabla_grad(v))-(-p)*Identity(dim)
    def J2(eps):
        I1     = eps[0,0] + eps[1,1] + eps[2,2]
        I2     = eps[0,0]*eps[1,1] + eps[1,1]*eps[2,2] + \
                 eps[2,2]*eps[0,0] - eps[0,1]*eps[0,1] - \
                 eps[1,2]*eps[1,2] - eps[2,0]*eps[2,0]
        J2     = 1/3*I1*I1 - I2
        J2     = sqrt(J2)
        return J2
    # https://www.continuummechanics.org/principalstressesandstrains.html
    # calc_princ is derived from the above link.
    def calc_princ(sigma):
        thetap = atan(2*sigma[0,1]/(sigma[0,0] - sigma[1,1]))
        thetap = thetap/math.pi*180/2

        sigma1 = (sigma[0,0] + sigma[1,1])/2 + sqrt(pow((sigma[0,0]/2 - sigma[1,1]/2),2) + sigma[0,1]*sigma[0,1])
        sigma3 = (sigma[0,0] + sigma[1,1])/2 - sqrt(pow((sigma[0,0]/2 - sigma[1,1]/2),2) + sigma[0,1]*sigma[0,1])

        return thetap, sigma1, sigma3
    
    # Define anisotropic viscosity constitutive relation
    # Based on Maulhaus et al. (2002). 
    def sigma_anisotropic(v, eta, eta_s, dim):
        if dim == 2:
            n = np.zeros((1,dim))
            n[0][0] = norm1
            n[0][1] = norm3       
        elif dim == 3:
            n = np.zeros((1,dim))
            n[0][0] = norm1
            n[0][1] = norm2
            n[0][2] = norm3
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
                        if l == j :
                            a = n[0][i]*n[0][k]
                        if i == l :
                            b = n[0][j]*n[0][k]
                        if k == j :
                            c = n[0][i]*n[0][l]
                        if i == k :
                            d = n[0][j]*n[0][l]
                        componentList[i][j][k][l] = 2*(eta_s - eta) * ((a + b + c + d)/2 - 2*e)
                        #C_lambda[i][j][k][l] = (a + b + c + d)/2 - 2*e

        C = as_tensor(componentList)
        vtmp = as_tensor(epsilon(v))
        i, j, k, l = ufl.indices(4)
        C1 = ufl.as_tensor(C[i,j,k,l]*vtmp[k,l],(i,j))
        C2 = 2*eta*sym(nabla_grad(v))
        C3 = C1 + C2
        return C3
        
    def sigma_general_anisotropy(v,D):
        # sigma_general_anisotropy calculates the stress given velocity vector u and 
        #   material 'stiffness' matrix C after rotation to the current x-y-z
        #   coordinate system.
        
        dd = as_tensor(D2dd(D))
        vtmp = as_tensor(epsilon(v)) 
        i, j, k, l = ufl.indices(4)
        res  = ufl.as_tensor(dd[i,j,k,l]*vtmp[k,l],(i,j))
        return res
        
    # Define variational problem
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    
    if iso==True:
        theta2    = -10000
        nametag   = '_isotropic' # Appended in the result files.
        # The following line is from the original FEniCS stokes tutorial.
        # However, the symmetry of grad(u) is a little confusing.
        #a         = inner(grad(u), eta*grad(v))*dx + div(v)*p*dx + q*div(u)*dx
        
        a         = inner(sym(nabla_grad(u))*2.0*eta, grad(v))*dx + div(v)*p*dx + q*div(u)*dx
        L         = inner(f, v)*dx
        # b is used to build the preconditioner matrix
        b         = inner(grad(u), eta*grad(v))*dx + p*q*dx
    elif iso==False:
        nametag   = '_anisotropic'
        a1        = inner(sigma_anisotropic(u,eta,eta_s,dim), grad(v))*dx
        a2        = div(v)*p*dx + q*div(u)*dx
        a         = a1 + a2 
        L         = inner(f, v)*dx# + dot(g, v) * ds(bound_name_list['north']) + dot(g, v) * ds(bound_name_list['south'])
        # b is used to build the preconditioner matrix
        b         = inner(grad(u), 2*eta*grad(v))*dx + p*q*dx    
    elif iso == 'general_anisotropic':
        print('Assembling general_anisotropic case with customized D matrix.')
        a1        = inner(sym(nabla_grad(u))*2.0, grad(v))*dx(iso_domain) # this line works.
        a2        = inner(sigma_general_anisotropy(u,anisoC) + sym(nabla_grad(u))*2.0, grad(v))*dx(anis_domain)
        a3        = div(v)*p*dx + q*div(u)*dx
        a         = a1 + a2 + a3
        L         = inner(f, v)*dx# + dot(g, v) * ds(bound_name_list['north']) + dot(g, v) * ds(bound_name_list['south'])
        # b is used to build the preconditioner matrix
        b         = inner(grad(u), 2.0*grad(v))*dx + p*q*dx   
    # Assemble the system
    A, bb         = assemble_system(a, L, bcs)
    # Assemble the preconditioner system
    P, btmp       = assemble_system(b, L, bcs)
    print('Assembling the system and preconditioner ...')
    
    if sol==1:
        # Create Krylov solver and AMG preconditioner
        #solver        = KrylovSolver("minres", "amg") # best.
        #solver = KrylovSolver("gmres", "amg") # second best.
        #solver = KrylovSolver("tfqmr", "amg") # least
        
        solver  = KrylovSolver(solverid, precond)

        solver.parameters["relative_tolerance"]      = rtol
        solver.parameters["absolute_tolerance"]      = 1.0e-15
        solver.parameters["divergence_limit"]        = 1.0e4
        solver.parameters["maximum_iterations"]      = 10000
        solver.parameters["error_on_nonconvergence"] = True
        solver.parameters["nonzero_initial_guess"]   = False
        solver.parameters["report"]                  = True
        solver.parameters["monitor_convergence"]     = True

        # Associate operator (A) and preconditioner matrix (P)
        solver.set_operators(A, P)    
        # Solve
        U             = Function(W)
        solver.solve(U.vector(), bb)
        print('Solving the system using Krylov solver, be patient...')
    elif sol==2:
        # Solve
        U             = Function(W)
        solve(a==L, U, bcs, solver_parameters={'linear_solver': 'mumps', 'preconditioner': 'amg'})
    
    # Get sub-functions
    u, p          = U.split()

    # Calculate stress and strain-rate tensors 
    if tensorFunc=='DG0':
        TensorFunc    = TensorFunctionSpace(mesh, "DG", 0) # Set DG0 for stress & strain rate
        ScalarFunc = FunctionSpace(mesh, "DG", 0)
    elif tensorFunc=='CG1':
        TensorFunc    = TensorFunctionSpace(mesh,"CG",1)
        ScalarFunc = FunctionSpace(mesh, "CG", 1)
    elif tensorFunc=='DG2':
        TensorFunc    = TensorFunctionSpace(mesh,"DG",2)
        ScalarFunc = FunctionSpace(mesh, "DG", 2)
    elif tensorFunc=='CG2':
        TensorFunc    = TensorFunctionSpace(mesh,"CG",2)
        ScalarFunc = FunctionSpace(mesh, "CG", 2)

    sig           = Function(TensorFunc, name="Stress")
    sig_iso       = Function(TensorFunc, name="Stress")
    sig_aniso     = Function(TensorFunc, name="Stress")
    strain_rate   = Function(TensorFunc, name="Strain")
    pressure      = Function(ScalarFunc)
    
    print(iso)
    # Project stresses to Vsig space.
    if iso==True:
        print('Case ',iso, ', using isotropic sig func')
        sig.assign(project(sigma(u,p,eta), TensorFunc,solver_type='cg', preconditioner_type='hypre_amg'))
    elif iso==False: # original Mulhuas implementation: 
        print('Case ',iso, ', using anisotropic sig func')
        sig.assign(project(sigma_anisotropic(u,eta,eta_s,dim)-(-p)*Identity(dim), TensorFunc, solver_type='cg', preconditioner_type='hypre_amg'))
    elif iso=='general_anisotropic': # general anisotropic function implementation
        print('Calculating stress from velocity fields ...')
        print('Case is ', case, ' with general anisotropic D matrix implementation ...')
        sig.assign(project(stress_coef*sigma_general_anisotropy(u,anisoC) + sym(nabla_grad(u))*2.0-(-p)*Identity(dim), TensorFunc, solver_type='cg', preconditioner_type='hypre_amg'))
        
    strain_rate.assign(project(epsilon(u), TensorFunc, solver_type='cg', preconditioner_type='hypre_amg'))
    DevStrainRate = epsilon(u) - (1./3)*tr(epsilon(u))*Identity(dim)
    J2StrainRate  = sqrt(3./2*inner(DevStrainRate, DevStrainRate))
    J2StrainRate  = project(J2StrainRate, ScalarFunc, solver_type='cg', preconditioner_type='hypre_amg')
    pressure      = project(p, ScalarFunc, solver_type='cg', preconditioner_type='hypre_amg')
    
    print('Writing out the results ...')
    # Save solution in xdmf format
    outpath   = "../res/case" + str(case) + "/" + str(date_stamp) + "/"
    prefix    = ""
    suffix    = "_theta" + "{:.1f}".format(delta) + "_hetero_" + str(hete) + ".xdmf"
    serial_write_xdmf(mesh, u, "velocity", outpath, prefix, suffix)
    serial_write_xdmf(mesh, p, "p", outpath, prefix, suffix)
    serial_write_xdmf(mesh, pressure, "pressure", outpath, prefix, suffix)
    serial_write_xdmf(mesh, sig, "stress", outpath, prefix, suffix)
    serial_write_xdmf(mesh, strain_rate,  "strain_rate",  outpath, prefix, suffix)
    serial_write_xdmf(mesh, J2StrainRate, "J2StrainRate", outpath, prefix, suffix)      
    print('Finished ...')


# In[ ]:


# theta is the angle of the layer normal to the horizontal axis in the square model.
# theta should be >90 and <180.
# delta = theta - 90.

print(case_info[case])
if case == 26:
    n= 40
    #useMPI=False
    iso=False
    guassian=True
    sol=1
    tensorFunc='DG0'
    eta_strong, eta_weak = 1, 0.01
    msh_path = None
    msh_name = None
    alpha = -999
    #date_stamp = '20230406_case26_StressCG1_C'+str(eta_strong/eta_weak) +'_iso_' + str(iso)
    date_stamp = '20230418_case26_C'+str(eta_strong/eta_weak)+'_n'+str(n)+'_iso_'+str(iso) + '_gaussian_'+str(guassian) #date_stamp for Guassian case.
    #for delta in np.arange(20,25,5):	
    for delta in np.arange(0,90+22.5,22.5):
        stokes(delta, alpha, eta_strong, eta_weak, msh_path, msh_name, n)
        
if case == 30: 
    # Figure 4b & 6.
    n = 10
    eta_strong, eta_weak = 1, 0.1
    msh_path = None
    msh_name = None
    delta = 10
    alpha = 45
    stokes(delta, alpha, eta_strong, eta_weak, msh_path, msh_name, n)
        
if case == 31:
    # Figure 4a & 5.
    n = 10
    eta_strong, eta_weak = 1, 0.1
    msh_path = None
    msh_name = None
    delta = 12.5
    stokes(delta, alpha, eta_strong, eta_weak, msh_path, msh_name, n)
        
if case == 32:
    # Figure 7 & 8.
    msh_path = "../msh/Cascadia_schist_mesh_case32/"
    msh_name = "cascadia_thrust"
    eta_strong, eta_weak = 1, 0.1
    delta = 30
    alpha = -999
    n = -999
    stokes(delta, alpha, eta_strong, eta_weak, msh_path, msh_name, n)

if case == 33:
#    useMPI     = False
    iso        = False
    n          = 50 
    eta_strong = 1
    eta_weak   = 0.1
    msh = 'case33_gmsh_0.02'
    msh_path   = './msh/'+msh +'/'
    sol = 1
    alpha      = -999
    tensorFunc = 'DG0' # DG0 or CG1
    #modeList = ['PureShear']
    #modeList   = ['FreeNorthSouth']
    modeList   = ['FreeNorthSouthFSTop']
    #modeList   = ['FreeNorth']
    #modeList   = ['FreeSlip']
    #modeList   = ['FreeNorthSouth','PureShear']
    
    #m    = sys.argv[1]
    #iso  = sys.argv[1]
    #n1   = int(sys.argv[2])
    
    #modeList = [m]
    
    for mode in modeList:
        
        date_stamp = '20230413_'+msh+'_C'+str(eta_strong/eta_weak) +'_' + mode + '_iso_' + str(iso) + '_solver_' + str(sol) + '_tensorFunc_'+tensorFunc
        for delta in np.arange(25,95,5):
        #for delta in np.arange(25,95,5):
            print(mode, iso, delta)
            msh_name   = "case33_" + str(delta) 
            st = time.time()
            stokes(delta, alpha, eta_strong, eta_weak, msh_path, msh_name, n, mode)
            et = time.time()
            print(et - st)

if case == 34: # fossil shear zone model with Hill viscous anisotropy inside the shear zone 
    iso        = 'general_anisotropic'
    n          = 50 
    eta_strong = 1
    eta_weak   = 0.1
    msh = 'case33_gmsh_0.02'
    msh_path   = './msh/'+msh +'/'
    sol = 1
    alpha      = -999
    tensorFunc = 'DG0' # DG0 or CG1
    #modeList   = ['FreeNorthSouth']
    modeList   = ['FreeNorthSouthFSTop']
    #modeList   = ['FreeNorth']
    #modeList   = ['FreeSlip']
    #modeList   = ['FreeNorthSouth','PureShear']
    
    rheology = 'HW' # 'M'/'HW'/'ortho'
    
    F,G,H,L,M,N = 0.5, 0.5, 0.5, 1.5, 1.5, 1.5
    A0 = HillA(F,G,H,L,M,N) # isotropic A
    F,G,H,L,M,N = 0.19,0.37,0.51,2.85,1.59,1.67
    A1 = HillA(F,G,H,L,M,N)
    F,G,H,L,M,N = 0.02, 0.23,0.37,8.92,2.15,2.21
    A2 = HillA(F,G,H,L,M,N) # A case2
    np.set_printoptions(precision=3)
    
    for mode in modeList:
        #date_stamp = '20230413_'+msh+'_C'+str(eta_strong/eta_weak) +'_' + mode + '_iso_' + str(iso) + '_solver_' + str(sol) + '_tensorFunc_'+tensorFunc
        date_stamp = '20230413_'+msh+'_' + rheology + '_' + mode + '_iso_' + str(iso) + '_solver_' + str(sol) + '_tensorFunc_'+tensorFunc
        
        for delta in np.arange(45,95,5):
        #for delta in np.arange(25,95,5):
            print(mode, iso, delta)
            t1 = -delta
            t2 = -delta+90.
            R = np.zeros((3,3))
            R[0,0] = np.cos(t1/180.*math.pi) 
            R[0,1] = np.sin(t1/180.*math.pi)
            R[1,0] = np.cos(t2/180.*math.pi)
            R[1,1] = np.sin(t2/180.*math.pi)
            R[2,2] = 1.
            # R0 = np.zeros((3,3))
            # R0[0,0] = 1.
            # R0[1,1] = 1.
            # R0[2,2] = 1.
            print('Rotation matrix R', R)
            
            if rheology == 'M':
                anisoC = D_M_ti_rot(R)
            elif rheology == 'HW':
                anisoC = D_HW_ti_rot(R)
            elif rheology == 'ortho':
                anisoC = D_ortho_rot(R)
            
            print('The general anisotropic rheology used is ')
            print(rheology)
            print(anisoC)
            
            msh_name   = "case33_" + str(delta) 
            st = time.time()
            stokes(delta, alpha, eta_strong, eta_weak, msh_path, msh_name, n, mode)
            et = time.time()
            print(et - st)   

if case == 35: 
    print('Case 35: 3D vertical shear zone with normal director rotated and subject to simple shearing.')
    iso        = 'general_anisotropic'
    #n          = 10 
    eta_strong = 1
    eta_weak   = 0.1
    eta_1      = 0.3
    #wt         = 0.05 # 5%
    OrthoType = 1
    if OrthoType==1:
        adjustedDiagNorm, adjustedOffNorm, adjustedShear = 0.6, 0.6, 0. 
    elif OrthoType==2:
        adjustedDiagNorm, adjustedOffNorm, adjustedShear = 0.0, 0.0, 0.3 #1., 0.3, 0.
    msh = 'case35GMSH001Pbc'
    msh_path   = './msh/'+msh +'/'
    msh_name   = msh
    sol        = 1 
    alpha      = -999
    tensorFunc = 'DG0' # DG0/CG1/CG2/DG2
    mode       = None
    print('=====================================================================================')
    print('Input parameters:')    
    print('For case35, please input delta1, delta2, beta, rheology, precond, solver, rtol, ny_or_nz')
    print('sys.argv is ', sys.argv)
    delta1 = float(sys.argv[1])
    delta2 = float(sys.argv[2])
    beta   = float(sys.argv[3])
    rheology = sys.argv[4]
    precond  = sys.argv[5]
    solverid   = sys.argv[6]
    rtol     = float(sys.argv[7])
    n        = int(sys.argv[8])

    #rheology = 'M' # 'M'/'HW'/'ortho'
    
    date_stamp = '20240401MPIC10FEniCSMshCG3DG2Ny'+str(n*4)+rheology+'Sol'+str(sol)+tensorFunc+str(rtol)+'beta'+str(int(beta))+precond+solverid
    # Use GMSH mesh instead.
    #date_stamp = '20231116'+msh+'CG3DG2'+str(iso)+rheology+ '_sol_' + str(sol)+tensorFunc+str(rtol)+'beta'+str(beta)+precond+solverid
    
    ddelta = 10
    for delta in np.arange(delta1,delta2+ddelta,ddelta):
    #for delta in np.arange(70,90+ddelta,ddelta):
        print('Case 35: mode, iso, delta are ', mode, iso, delta)
        print('Case 35: delta is defined as the angle from y+ to the normal director of weak anisotropy in Maulhaus TI counterclockwise.')
        print('Case 35: rotation matrix is defined by angles -delta and -delta+90 ... ...')
        
        R = np.dot(Rz(delta),Ry(beta))
        print('Case 35: Rotation matrix R is ', R)
        
        if rheology == 'M':
            anisoC = D_M_ti(R, eta_strong, eta_weak)
        elif rheology == 'HW':
            anisoC = D_HW_ti(R, eta_strong, eta_weak, eta_1)
        elif rheology == 'ortho':
            anisoC = getOrthorhombicD(R, eta_strong, eta_weak, eta_1, adjustedDiagNorm, adjustedOffNorm, adjustedShear)
        
        print('Case 35: The general anisotropic rheology used is ')
        print(rheology)
        print('Case 35: The anisotropic D matrix is ')
        print(anisoC)
        print('Case 35: The viscosity  D matrix is ')
        print(anisoC + D_isotropic(eta_strong) )
        
        st = time.time()
        stokes(delta, alpha, eta_strong, eta_weak, msh_path, msh_name, n, mode)
        et = time.time()
        print('Case 35: Total computational time consumed in seconds are ', et - st)   
            
if case == 23: # 2D fossile mantle shear zone models
#    useMPI     = False
    iso        = False
    n          = 50 
    eta_strong = 1
    eta_weak   = 0.1
    msh = 'case23_gmsh_0.005'
    msh_path   = './msh/'+msh +'/'
    sol = 2
    alpha      = -999
    tensorFunc = 'DG0' # DG0 or CG1
    #modeList = ['PureShear', 'FreeNorthSouth']
    #modeList   = ['FreeNorthSouth']
    modeList   = ['PureShear']
    
    #m    = sys.argv[1]
    #iso  = sys.argv[1]
    #n1   = int(sys.argv[2])
    
    #modeList = [m]
    
    for mode in modeList:
        
        date_stamp = '20230413_'+msh+'_C'+str(eta_strong/eta_weak) +'_' + mode + '_iso_' + str(iso) + '_solver_' + str(sol) + '_tensorFunc_'+tensorFunc
        for delta in np.arange(45,95,5):
        #for delta in np.arange(25,95,5):
            print(mode, iso, delta)
            msh_name   = "case23_" + str(delta) 
            st = time.time()
            stokes(delta, alpha, eta_strong, eta_weak, msh_path, msh_name, n, mode)
            et = time.time()
            print(et - st)
# In[ ]:




