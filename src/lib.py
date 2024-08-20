from dolfin import *
import numpy as np
import ufl as ufl
import sys
import math

# Define class K to set viscosity for two domains, anis_domain and iso_domain
class K(UserExpression):
    def __init__(self, subdomains, k_0, k_1, anis_domain, iso_domain, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.k_0 = k_0
        self.k_1 = k_1
        self.anis_domain = anis_domain
        self.iso_domain = iso_domain

    def eval_cell(self, values, x, cell):
        if self.subdomains[cell.index] == self.anis_domain: #faul zone/schist/or other anisotrpic domains.
            values[0] = self.k_0
        elif self.subdomains[cell.index] == self.iso_domain:
            values[0] = self.k_1
    def values_shape(self):
        return (1,)


def mpi_write_xdmf(mesh, var, varName, path, fileNamePrefix, fileNameSuffix):
    fullPath  = path + fileNamePrefix + varName + fileNameSuffix
    fid       = XDMFFile(mesh.mpi_comm(), fullPath)
    fid.parameters["functions_share_mesh"]  = True
    fid.parameters["rewrite_function_mesh"] = False
    fid.write(var)

def serial_write_xdmf(mesh, var, varName, path, fileNamePrefix, fileNameSuffix):
    fullPath  = path + fileNamePrefix + varName + fileNameSuffix
    fid       = XDMFFile(fullPath)
    fid.write(var)
        
