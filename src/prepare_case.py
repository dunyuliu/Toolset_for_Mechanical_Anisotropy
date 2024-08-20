from dolfin import *
import math

# This script hosts the function prepare_case that helps customize the mesh, function space, boundary conditions for different models.
# 
def prepare_case(case, path, name, n, delta, alpha, args):
    
    if case == 26:   
        print('Simulating Case 26: Same as 2D analytic solution. Horizontal anisotropic layer subjected to simple shear.')
        
        # Defining norms. For 2D cases, norm1 and norm3 are used.
        theta = 90 + delta
        norm1, norm2, norm3 = cos(theta/180*math.pi), 0, sin(theta/180*math.pi)       

        dim = 2

        len1, len2 = 5, 1
        nx = 5*n
        ny = n
        de = 1/n
        mesh = RectangleMesh(Point(0, 0), Point(len1, len2), nx, ny, "crossed")

        # Sub domain for Periodic boundary condition
        class PeriodicBoundary(SubDomain):
            # Left boundary is "target domain" G
            def inside(self, x, on_boundary):
                return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)
                #return bool(near(x[0],0)) and (not(near(x[0],len1)))
            # Map right boundary (H) to left boundary (G)
            def map(self, x, y):
                y[0] = x[0] - len1
                y[1] = x[1]

        # Create periodic boundary condition
        pbc = PeriodicBoundary()#PeriodicBoundary()  

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)  
        Q = FunctionSpace(mesh, "CG", 1)

        #W = V * Q # This expression is outdated
        New_element = MixedElement([V.ufl_element(), Q.ufl_element()])
        W = FunctionSpace(mesh, New_element, constrained_domain=pbc)

        # Boundaries
        boundaries = MeshFunction("size_t", mesh, dim-1)

        class Top_boundary(SubDomain):
            def inside(self, x, on_boundary):
                #return x[1] > (len2 - DOLFIN_EPS) and on_boundary
                return near(x[1], len2, DOLFIN_EPS) and on_boundary

        class Bottom_boundary(SubDomain):
            def inside(self, x, on_boundary):
                #return x[1] < DOLFIN_EPS and on_boundary
                return near(x[1], 0.0, DOLFIN_EPS) and on_boundary
        class Left_boundary(SubDomain):
            def inside(self, x, on_boundary):
                #return x[0] < DOLFIN_EPS and on_boundary
                return near(x[0], 0, DOLFIN_EPS) and on_boundary

        class Right_boundary(SubDomain):
            def inside(self, x, on_boundary):
                #return x[0] > (len1 - DOLFIN_EPS) and on_boundary    
                return near(x[0], len1, DOLFIN_EPS) and on_boundary
        boundaries.set_all(0)
        Top_boundary().mark(boundaries, 1)
        Bottom_boundary().mark(boundaries, 2)
        Left_boundary().mark(boundaries, 3)
        Right_boundary().mark(boundaries, 4)
        # Rename boundaries
        top = 1
        bottom = 2
        west = 3
        east = 4

        class Omega_0(SubDomain):
            def inside(self, x, on_boundary):
                return x[1] >= 0.5 and x[1] <= 0.9 # Anisotropic layer depth between 0.5 and 0.9.
            
        mf = MeshFunction("size_t", mesh, 2)

        subdomain0 = Omega_0()
        mf.set_all(0)
        
        subdomain0.mark(mf,1) # 1, the middle anisotrpic layer.
        bound_name_list = {'bottom':bottom,
                           'top':top,
                           'west':west,
                           'east':east,
                           'anis_domain':1,
                           'iso_domain':0}

        bcs = [DirichletBC(W.sub(0), Constant((1.0, 0.0)), boundaries, top),
               DirichletBC(W.sub(0), Constant((0.0,0.0)), boundaries, bottom)]
        f = Constant((0.0, 0.0))    
        
    elif case == 30:
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
    if case == 33 or case == 34: # Fossil shear zone (case=33) or Hill Anisotropy (case=34)  
        if case == 33:
            print('Simulating Case 33: Fossil shear zone - gmsh mesh ...')
        elif case == 34:
            print('Simulating Case 34: Fossil shear zone with Hill viscous anisotropy - gmsh mesh ...')
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
        north  = 29
        south  = 30
        bottom = 31
        west   = 32
        east   = 33
        top    = 34
        surrounding = 36 
        shearzone = 35

        bound_name_list = {'bottom':bottom,
                           'top': top,
                           'west':west,
                           'east':east,
                           'north':north,
                           'south':south,
                           'anis_domain':shearzone,
                           'iso_domain':surrounding}
        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)

        #W = V * Q # This expression is outdated
        New_element = MixedElement([V.ufl_element(), Q.ufl_element()])
        W = FunctionSpace(mesh, New_element) # No need for periodic boundary condition in this case.
        
        if args     == 'FreeSlip':
            bcs = [DirichletBC(W.sub(0).sub(0), Constant((0.0)), boundaries, east),
               DirichletBC(W.sub(0).sub(0), Constant((1.0)), boundaries, west),
               DirichletBC(W.sub(1), Constant((0.0)), boundaries, south),
               DirichletBC(W.sub(1), Constant((0.0)), boundaries, north),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, bottom),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, top)]
        if args     == 'PureShear':
            bcs = [DirichletBC(W.sub(0).sub(0), Constant((0.0)), boundaries, east),
               DirichletBC(W.sub(0).sub(0), Constant((1.0)),  boundaries, west),
               DirichletBC(W.sub(0).sub(1), Constant((-0.2)), boundaries, south),
               DirichletBC(W.sub(0).sub(1), Constant((0.2)),  boundaries, north),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)),  boundaries, bottom)]
               #DirichletBC(W.sub(0).sub(2), Constant((0.0)),  boundaries, top)]
        elif args   == 'FreeNorth':
            bcs = [DirichletBC(W.sub(0).sub(0), Constant((0.0)), boundaries, east),
               DirichletBC(W.sub(0).sub(0), Constant((1.0)), boundaries, west),
               DirichletBC(W.sub(0).sub(1), Constant((0.0)), boundaries, south),
               #DirichletBC(W.sub(0).sub(1), Constant((0.0)), boundaries, north),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, bottom),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, top)]
        elif args   == 'FreeNorthSouth':
            bcs = [DirichletBC(W.sub(0).sub(0), Constant((0.0)), boundaries, east),
               DirichletBC(W.sub(0).sub(0), Constant((1.0)), boundaries, west),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, bottom),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, top)]
       	elif args   == 'FreeNorthSouthFSTop':
            # free surface
            bcs = [DirichletBC(W.sub(0).sub(0), Constant((0.0)), boundaries, east),
               DirichletBC(W.sub(0).sub(0), Constant((1.0)), boundaries, west),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, bottom)]
 
        f = Constant((0.0, 0.0, 0.0)) # body force. Set to zero in this case.
        
    elif case == 35: 
        print('Case 35: 3D box model with a vertical shear zone subjected to simple shearing ... ...')
        print('Case 35: the normal director is rotated.')
        print('Case 35: it is 3D version of case 26, the analytic solution.')
        print('Case 35: please choose which mesh to use - ')
        use_msh = 'built-in' # 'built-in', 'gmsh', 'built-in-refined'
        
        theta = 90 + delta
        norm1 = cos(theta/180*math.pi)
        norm2 = sin(theta/180*math.pi)
        norm3 = 0    
        print('Mesh used is ', use_msh)
        print('Normal vectors (x,y,z) for the anisotropy is ', norm1, norm2, norm3)
        
        class PinPoint(SubDomain):
            def inside(self, x, on_boundary):
                return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS and x[2] < DOLFIN_EPS
        
        tol = DOLFIN_EPS
        class PeriodicBoundary(SubDomain):
            # Left boundary is "target domain" G
            def inside(self, x, on_boundary):
                return bool(abs(x[0]) < tol)
                
                #return bool(near(x[0],0)) and (not(near(x[0],len1)))
            # Map right boundary (H) to left boundary (G)
            def map(self, x, y):
                y[0] = x[0] - len1
                y[1] = x[1]            
                y[2] = x[2]
        pbc = PeriodicBoundary(tol)#PeriodicBoundary()  
        
        # Boundaries
        class Top_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[2], len3, DOLFIN_EPS) and on_boundary
        class Bottom_boundary(SubDomain):
            def inside(self, x, on_boundary):
                 return near(x[2], 0.0, DOLFIN_EPS) and on_boundary
        class West_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], 0, DOLFIN_EPS) and on_boundary
        class East_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], len1, DOLFIN_EPS) and on_boundary
        class South_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], 0, DOLFIN_EPS) and on_boundary
        class North_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], len2, DOLFIN_EPS) and on_boundary
                
        fault_zone_width = 0.4
        yc = 0.5
        class Omega_0(SubDomain):
            def inside(self, x, on_boundary):
                return x[1] > yc-fault_zone_width/2 and x[1] < yc+fault_zone_width/2 #\
#			and abs(x[2]-len3/2.)<0.2 # along depth
        
        if use_msh == 'built-in':
            dim  = 3
            len1 = 1.
            len2 = 1.
            len3 = 1.
            #mesh = UnitCubeMesh(n, n, n)
            mesh = BoxMesh(Point(0,0,0),Point(len1,len2,len3), 10, int(4*n), int(4*n))
            # refine mesh
            iz = 2
            #for i in range(iz):
            #    cell_markers = MeshFunction("bool", mesh, 3)
            #    cell_markers.set_all(False)
            #    for cell in cells(mesh):
            #        if cell.midpoint()[2]>0.8 and cell.midpoint()[1]<0.4 and cell.midpoint()[1]>0.2 and cell.midpoint()[0]<0.2 and cell.midpoint()[0]>0.05:
            #            cell_markers[cell] = True
            #    mesh = refine(mesh, cell_markers)
            
            boundaries = MeshFunction("size_t", mesh, dim-1)
            
            print('Mesh used is FEniCS built-in')

            boundaries.set_all(0)
            Top_boundary().mark(boundaries, 1)
            Bottom_boundary().mark(boundaries, 2)
            West_boundary().mark(boundaries, 3)
            East_boundary().mark(boundaries, 4)
            South_boundary().mark(boundaries, 5)
            North_boundary().mark(boundaries, 6)
            PinPoint().mark(boundaries,7)
            # Rename boundaries
            top    = 1
            bottom = 2
            west   = 3
            east   = 4
            south  = 5
            north  = 6
            pinpoint = 7
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
        elif use_msh == 'built-in-refined':
            print('need implementation')
        elif use_msh == 'gmsh':
            ## Use FEniCS built-in mesh, which performs poorly with 40 grids along y. 
            ## 20231011. 
            # Use GMSH mesh instead. 
            # Load mesh
            print('Mesh used is GMSH mesh')
            dim  = 3
            len1 = 0.2 # used to define periodic bcs.

            mesh = Mesh(path + name + '.xml')

            boundaries = MeshFunction("size_t", mesh, path + name + '_facet_region.xml')
            mf = MeshFunction("size_t", mesh, path + name + '_physical_region.xml')
            
            # Rename boundaries, labeled in Gmsh.
            north  = 31
            south  = 32
            bottom = 36
            west   = 34
            east   = 33
            top    = 35
            surrounding = 29 
            shearzone = 30

            bound_name_list = {'bottom':bottom,
                               'top': top,
                               'west':west,
                               'east':east,
                               'north':north,
                               'south':south,
                               'anis_domain':shearzone,
                               'iso_domain':surrounding}
                           
        print('Case 35: Periodic boundary along x direction.')

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 3)
        Q = FunctionSpace(mesh, "DG", 2) 

        #W = V * Q # This expression is outdated
        New_element = MixedElement([V.ufl_element(), Q.ufl_element()])
        W = FunctionSpace(mesh, New_element, constrained_domain=pbc) # For the function space with periodic boundary condition 

        print('Case 35: Set the shear zone parallel to the simple shearing and located between 0.5 and 0.9') 
        print('Case 35: along fault-normal direction.') 
        

        
        bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), boundaries, south),
               DirichletBC(W.sub(0), Constant((1.0, 0.0, 0.0)), boundaries, north),
               #DirichletBC(W.sub(1), Constant((1.0)), boundaries, top),
               
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, top),
               DirichletBC(W.sub(0).sub(2), Constant((0.0)), boundaries, bottom)]
        f = Constant((0.0, 0.0, 0.0))
        
        print('Case 35: preparation done ... ...')
        print(' ')
        
    if case == 23: # 3D fossil mantle shear zone
        print('Simulating Case 23: 2D Fossil shear zone - gmsh mesh ...')

        theta = 90 + delta
        norm1 = cos(theta/180*math.pi)
        norm2 = 0.
        norm3 = sin(theta/180*math.pi)
            
        print('Normal vectors for anisotrpy are ', norm1, norm2, norm3, 'in this scenario ...')
        
        # Load mesh
        dim = 2
        mesh = Mesh(path + name + '.xml')

        boundaries = MeshFunction("size_t", mesh, path + name + '_facet_region.xml')
        mf = MeshFunction("size_t", mesh, path + name + '_physical_region.xml')
        
        # Rename boundaries, labeled in Gmsh.
        bottom = 14
        top   = 13
        left   = 11
        right    = 12
        surrounding = 15 
        shearzone = 16

        bound_name_list = {'bottom':bottom,
                           'top': top,
                           'left':left,
                           'right':right,
                           'anis_domain':shearzone,
                           'iso_domain':surrounding}
        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)

        #W = V * Q # This expression is outdated
        New_element = MixedElement([V.ufl_element(), Q.ufl_element()])
        W = FunctionSpace(mesh, New_element) # No need for periodic boundary condition in this case.
        
        if args     == 'PureShear':
            bcs = [DirichletBC(W.sub(0).sub(0), Constant((0.0)), boundaries, right),
               DirichletBC(W.sub(0).sub(0), Constant((1.0)),  boundaries, left),
               DirichletBC(W.sub(0).sub(1), Constant((-0.2)), boundaries, bottom),
               DirichletBC(W.sub(0).sub(1), Constant((0.2)),  boundaries, top)]
        elif args   == 'FreeNorthSouth':
            bcs = [DirichletBC(W.sub(0).sub(0), Constant((0.0)), boundaries, right),
               DirichletBC(W.sub(0).sub(0), Constant((1.0)), boundaries, left)]
 
        f = Constant((0.0, 0.0)) # body force. Set to zero in this case.
    return mesh, boundaries, mf, bcs, f, bound_name_list, dim, W, norm1, norm2, norm3
