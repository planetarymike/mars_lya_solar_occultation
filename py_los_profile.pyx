#py_los_profile.pyx -- python wrapper for H LOS profile object
# distutils: language = c++

from libcpp cimport bool
from libcpp.string cimport string

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

cdef public api tonumpyarray(double* data, long long size) with gil:
    if not (data and size >= 0): raise ValueError
    cdef np.npy_intp dims = size
    #NOTE: it doesn't take ownership of `data`. You must free `data` yourself
    return np.PyArray_SimpleNewFromData(1, &dims, np.NPY_DOUBLE, <void*>data)

cdef extern from "nr3.h":
    cdef cppclass VecDoub:
        VecDoub() except +
        double* v
        int size()
        
cdef extern from "LOS_profile_object.h":
    cdef cppclass H_LOS_simulator:
        H_LOS_simulator() except +
        double nexo
        double Texo
        double lc
        double mpvel
        int nprof
        void read_scpos_from_file(string)
        void set_scpos(double*, double*,int)
        void set_scpos(VecDoub, VecDoub)
        void simulate(double, double)
        void write_to_file(string)
        VecDoub sc_alt_km
        VecDoub tanpt_alt_km
        VecDoub lb
        VecDoub sc_sza_deg
        VecDoub t0
        VecDoub* avec
        VecDoub* tvec
        
cdef class LOS_profile_sim:
    cdef H_LOS_simulator *thisptr #holds the reference to the cpp class
    def __cinit__(self):
        self.thisptr = new H_LOS_simulator()
    def __dealloc__(self):
        del self.thisptr
    def read_scpos_from_file(self, string infilename):
        self.thisptr.read_scpos_from_file(infilename)
    def write_to_file(self, string outfilename):
        self.thisptr.write_to_file(outfilename)
    def simulate(self, double nexo, double Texo):
        self.thisptr.simulate(nexo,Texo)
    def set_scpos(self, np.ndarray[np.double_t, ndim=1] tanpt_alt, np.ndarray[np.double_t, ndim=1] sc_sza, int nprof):
        self.thisptr.set_scpos(&tanpt_alt[0], &sc_sza[0], nprof)
    def get_avec(self):
        return tonumpyarray(self.thisptr.avec[0].v,self.thisptr.avec[0].size())
    def get_tvec(self,int i):
        return tonumpyarray(self.thisptr.tvec[i].v,self.thisptr.tvec[i].size())
    def get_alltvec(self):
        Tvec = np.ndarray(shape=(self.thisptr.nprof,self.thisptr.tvec[0].size()))
        for i in range(self.thisptr.nprof):
            Tvec[i]=tonumpyarray(self.thisptr.tvec[i].v,self.thisptr.tvec[i].size())
        return Tvec
    def get_sc_alt_km(self):
        return tonumpyarray(self.thisptr.sc_alt_km.v,self.thisptr.sc_alt_km.size())
    def get_tanpt_alt_km(self):
        return tonumpyarray(self.thisptr.tanpt_alt_km.v,self.thisptr.tanpt_alt_km.size())
    def get_lb(self):
        return tonumpyarray(self.thisptr.lb.v,self.thisptr.lb.size())
    def get_sc_sza_deg(self):
        return tonumpyarray(self.thisptr.sc_sza_deg.v,self.thisptr.sc_sza_deg.size())
    def get_t0(self):
        return tonumpyarray(self.thisptr.t0.v,self.thisptr.t0.size())
    def get_nexo(self):
        return self.thisptr.nexo
    def get_Texo(self):
        return self.thisptr.Texo
    def get_lc(self):
        return self.thisptr.lc
    def get_mpvel(self):
        return self.thisptr.mpvel
    def get_nprof(self):
        return self.thisptr.nprof

