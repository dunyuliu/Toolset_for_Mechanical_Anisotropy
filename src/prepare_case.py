from dolfin import *
import math

# This script hosts the function prepare_case that helps customize the mesh, function space, boundary conditions for different models.
# 
def prepare_case(case, path, name, n, delta, alpha):
    
    
    if case == 30:
        print('Simulating Case 30: 3D box model with a horizontal anisotropic layer.')
        
        theta = delta
        norm1 = cos(theta/180*math.pi)*cos(alpha/180*math.pi)
        norm2 = cos(theta/180*math.pi)*sin(alpha/180*math.pi)
        norm3 = sin(theta/180*math.pi)    
        print('Normal vectors for anisotrpy are ', norm1, norm2, norm3, 'in this scenario ...')
    
        mesh = UnitCubeMesh(n, n, n)
        dim = 3
        len1 = 1
        len2 = 1
        len3 = 1
        # Create periodic boundary condition
        # Sub domain for Periodic boundary condition
        class PeriodicBoundary(SubDomain):
            # Left and front boundarires are "target domain" G
            def inside(self, x, on_boundary):
                return bool((near(x[0],0) or near(x[1],0)) and 
                            (not((near(x[0],len1) and near(x[1],0)) or 
                                (near(x[0],0) and near(x[1],len2)))) and on_boundary)
            # Map right boundary (H) to left boundary (G)
            def map(self, x, y):
                if near(x[0],len1) and near(x[1],len2):
                    y[0] = x[0] - len1
                    y[1] = x[1] - len2
                    y[2] = x[2]
                elif near(x[0],len1):
                    y[0] = x[0] - len1
                    y[1] = x[1] 
                    y[2] = x[2] 
                elif near(x[1],len2):
                    y[0] = x[0] 
                    y[1] = x[1] - len2
                    y[2] = x[2]     
                else:
                    y[0] = -1000
                    y[1] = -1000
                    y[2] = -1000
                    
        pbc = PeriodicBoundary()#PeriodicBoundary()  

        # Pin the DOF on the bottom-left corner of the mesh, here 0, 0
        # Note that the choice of location is arbitrary.
        class PinPoint(SubDomain):
            def inside(self, x, on_boundary):
                return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS     

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)

        #W = V * Q # This expression is outdated
        New_element = MixedElement([V.ufl_element(), Q.ufl_element()])
        W = FunctionSpace(mesh, New_element, constrained_domain=pbc) # For the function space with periodic boundary condition 

        # Boundaries
        boundaries = MeshFunction("size_t", mesh, dim-1)

        # Boundaries
        class Top_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[2], len3, DOLFIN_EPS) and on_boundary
        class Bottom_boundary(SubDomain):
            def inside(self, x, on_boundary):
                 return near(x[2], 0.0, DOLFIN_EPS) and on_boundary
        class Left_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], 0, DOLFIN_EPS) and on_boundary
        class Right_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], len1, DOLFIN_EPS) and on_boundary
        class Front_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], 0, DOLFIN_EPS) and on_boundary
        class Back_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], len2, DOLFIN_EPS) and on_boundary
            
        boundaries.set_all(0)
        Top_boundary().mark(boundaries, 1)
        Bottom_boundary().mark(boundaries, 2)
        Left_boundary().mark(boundaries, 3)
        Right_boundary().mark(boundaries, 4)
        Front_boundary().mark(boundaries, 5)
        Back_boundary().mark(boundaries, 6)
        # Rename boundaries
        top = 1
        bottom = 2
        west = 3
        east = 4
        south = 5
        north = 6
            
        class Omega_0(SubDomain):
            def inside(self, x, on_boundary):
                return x[2] >= 0.5 and x[2] <= 0.9 # along depth
            
        mf = MeshFunction("size_t", mesh, dim)

        subdomain0 = Omega_0() # horizontal layer
        mf.set_all(0)
        subdomain0.mark(mf,1) # 1, the middle anisotrpic layer.

        bound_name_list = {'bottom':bottom,
                           'top':top,
                           'west':west,
                           'east':east,
                           'north':north,
                           'south':south, 
                           'anis_domain':1,
                           'iso_domain':0}
        
        bcs = [DirichletBC(W.sub(0), Constant((1.0, 0.0, 0.0)), boundaries, top),
                DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), boundaries, bottom)]
               
        print('Defining boundary conditions ...')

        f = Constant((0.0, 0.0, 0.0))
        
    elif case == 31: 
        print('Simulating Case 31: 3D box model with a vertical SAF fault zone')
        
        theta = 90 + delta
        norm1 = cos(theta/180*math.pi)
        norm2 = sin(theta/180*math.pi)
        norm3 = 0    
        print('Normal vectors for anisotrpy are ', norm1, norm2, norm3, 'in this scenario ...')
    
        mesh = UnitCubeMesh(n, n, n)
        dim = 3
        len1 = 1
        len2 = 1
        len3 = 1
        
        # Periodic boundary along x direction. 
        class PeriodicBoundary(SubDomain):
            # Left boundary is "target domain" G
            def inside(self, x, on_boundary):
                return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)
                #return bool(near(x[0],0)) and (not(near(x[0],len1)))
            # Map right boundary (H) to left boundary (G)
            def map(self, x, y):
                y[0] = x[0] - len1
                y[1] = x[1]            
                y[2] = x[2]
        pbc = PeriodicBoundary()#PeriodicBoundary()  

        # Pin the DOF on the bottom-left corner of the mesh, here 0, 0
        # Note that the choice of location is arbitrary.
        class PinPoint(SubDomain):
            def inside(self, x, on_boundary):
                return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS     

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1) 

        #W = V * Q # This expression is outdated
        New_element = MixedElement([V.ufl_element(), Q.ufl_element()])
        W = FunctionSpace(mesh, New_element, constrained_domain=pbc) # For the function space with periodic boundary condition 

        # Boundaries
        boundaries = MeshFunction("size_t", mesh, dim-1)

        # Boundaries
        class Top_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[2], len3, DOLFIN_EPS) and on_boundary
        class Bottom_boundary(SubDomain):
            def inside(self, x, on_boundary):
                 return near(x[2], 0.0, DOLFIN_EPS) and on_boundary
        class Left_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], 0, DOLFIN_EPS) and on_boundary
        class Right_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], len1, DOLFIN_EPS) and on_boundary
        class Front_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], 0, DOLFIN_EPS) and on_boundary
        class Back_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], len2, DOLFIN_EPS) and on_boundary
            
        boundaries.set_all(0)
        Top_boundary().mark(boundaries, 1)
        Bottom_boundary().mark(boundaries, 2)
        Left_boundary().mark(boundaries, 3)
        Right_boundary().mark(boundaries, 4)
        Front_boundary().mark(boundaries, 5)
        Back_boundary().mark(boundaries, 6)
        # Rename boundaries
        top = 1
        bottom = 2
        west = 3
        east = 4
        south = 5
        north = 6

        # The fault zone strike is along x and its width extends along y. 
        fault_zone_width = 0.2
        class Omega_0(SubDomain):
            def inside(self, x, on_boundary):
                return x[1] >= len2/2-fault_zone_width and x[1] <= len2/2+fault_zone_width # along depth
            
        mf = MeshFunction("size_t", mesh, dim)

        subdomain0 = Omega_0() # horizontal layer
        mf.set_all(0)
        subdomain0.mark(mf,1) # 1, the middle anisotrpic layer.

        bound_name_list = {'bottom':bottom,
                               'top':top,
                               'west':west,
                               'east':east,
                               'north':north,
                               'south':south, 
                               'anis_domain':1,
                               'iso_domain':0}
        
        bcs = [DirichletBC(W.sub(0), Constant((1.0, 0.0, 0.0)), boundaries, south),
               DirichletBC(W.sub(0), Constant((-1.0, 0.0, 0.0)), boundaries, north),
                DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, bottom)]
        f = Constant((0.0, 0.0, 0.0))
        
    if case == 32: # Cascadia Leech River Schist    
        print('Simulating Case 32: Simplified 3D box model with the Leech River Schist above Cascadia Subduction Zone')

        theta = 90 + delta
        norm1 = cos(theta/180*math.pi)
        norm2 = sin(theta/180*math.pi)
        norm3 = 0    
        print('Normal vectors for anisotrpy are ', norm1, norm2, norm3, 'in this scenario ...')
        
        # Load mesh
        dim = 3
        mesh = Mesh(path + name + '.xml')

        boundaries = MeshFunction("size_t", mesh, path + name + '_facet_region.xml')
        mf = MeshFunction("size_t", mesh, path + name + '_physical_region.xml')
        
        # Rename boundaries, labeled in Gmsh.
        bottom = 30
        west = 28
        east = 29
        north = 27
        south = 26
        schist = 31 
        other = 32

        bound_name_list = {'bottom':bottom,
                           'west':west,
                           'east':east,
                           'north':north,
                           'south':south,
                           'anis_domain':schist,
                          'iso_domain':other}
        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)

        #W = V * Q # This expression is outdated
        New_element = MixedElement([V.ufl_element(), Q.ufl_element()])
        W = FunctionSpace(mesh, New_element) # No need for periodic boundary condition in this case.
        
        bcs = [DirichletBC(W.sub(0), Constant((1.0,0.0,-1.0)), boundaries, west),
            DirichletBC(W.sub(0).sub(0), Constant((0.0)), boundaries, east),
            DirichletBC(W.sub(0).sub(1), Constant((0.0)), boundaries, north),
            DirichletBC(W.sub(0).sub(1), Constant((0.0)), boundaries, south),
            DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, bottom)]
        
        f = Constant((0.0, 0.0, 0.0)) # body force. Set to zero in this case.
    return mesh, boundaries, mf, bcs, f, bound_name_list, dim, W, norm1, norm2, norm3