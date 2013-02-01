#include <Python.h>
#include "numpy/arrayobject.h"
#include "global_vars.h"

/*Wraps the flux_extractor into a python module called spectra_priv. Don't call this directly, call the python wrapper.*/

/*****************************************************************************/
/*Interface for SPH interpolation*/
PyObject * Py_SPH_Interpolation(PyObject *self, PyObject *args)
{
    //Things which should be from input
    int nbins, NumLos, nspecies;
    long long Npart;
    double box100;
    npy_intp size[3];
    //Input variables in np format
    PyArrayObject *pos, *vel, *mass, *u, *ne, *h, *fractions;
    PyArrayObject *xx, *yy, *zz, *axis;
    PyArrayObject * rho_out, *temp_out, *vel_out;

    //For storing output
    interp species, metal_spec;
    interp* metal_ptr=&metal_spec;

    //Temp variables
    int nxx=0, i;
    los *los_table=NULL;
    sort_los *sort_los_table=NULL;
    struct particle_data P;
    //Get our input
    if(!PyArg_ParseTuple(args, "idO!O!O!O!O!O!O!O!O!O!O!O!",&nbins, &box100,  &PyArray_Type, &pos, &PyArray_Type, &vel, &PyArray_Type, &mass, &PyArray_Type, &u, &PyArray_Type, &ne, &PyArray_Type, &fractions, &PyArray_Type, &h, &PyArray_Type, &axis, &PyArray_Type, &xx, &PyArray_Type, &yy, &PyArray_Type, &zz) )
      return NULL;

    NumLos = PyArray_DIM(xx,0);
    Npart = PyArray_DIM(pos,0);
    //NOTE if nspecies == 1, fractions must have shape [1,N], rather than [N]
    nspecies = PyArray_DIM(fractions, 0);
    //Malloc stuff
    los_table=malloc(NumLos*sizeof(los));
    sort_los_table=malloc(NumLos*sizeof(sort_los));
    InitLOSMemory(&metal_spec,NumLos, nbins*nspecies);

    //Initialise P from the data in the input numpy arrays.
    P.Pos = (float *) PyArray_GETPTR2(pos,0,0);
    P.Vel = (float *) PyArray_GETPTR2(vel,0,0);
    P.Mass = (float *) PyArray_GETPTR1(mass,0);
    P.U = (float *) PyArray_GETPTR1(u,0);
    P.Ne = (float *) PyArray_GETPTR1(ne,0);
    P.h = (float *) PyArray_GETPTR1(h,0);
    P.fraction = (float *) PyArray_GETPTR2(fractions,0,0);

    //Initialise los_table from input
    for(i=0; i< NumLos; i++){
        los_table[i].axis = *(int *) PyArray_GETPTR1(axis,i);
        los_table[i].xx = *(float *) PyArray_GETPTR1(xx,i);
        los_table[i].yy = *(float *) PyArray_GETPTR1(yy,i);
        los_table[i].zz = *(float *) PyArray_GETPTR1(zz,i);
    }

    //Do the work
    populate_sort_los_table(los_table, NumLos, sort_los_table, &nxx);
    SPH_Interpolation(NULL,&species, nspecies, nbins, Npart, NumLos, box100, los_table,sort_los_table,nxx, &P);

    size[0] = NumLos;
    size[1] = nbins;
    //Number of metal species
    size[2] = nspecies;
    /*Is there a better way to do this?*/
    rho_out = (PyArrayObject *) PyArray_SimpleNew(3, size, NPY_FLOAT);
    vel_out = (PyArrayObject *) PyArray_SimpleNew(3, size, NPY_FLOAT);
    temp_out = (PyArrayObject *) PyArray_SimpleNew(3, size, NPY_FLOAT);
    for(i=0; i< NumLos; i++){
        for(int j=0; j< nbins; j++){
            for(int k=0; k < nspecies; k++){
                *(float *) PyArray_GETPTR3(rho_out,i,j,k) = metal_spec.rho[i*nbins+k*nspecies+k];
                *(float *) PyArray_GETPTR3(vel_out,i,j,k) = metal_spec.veloc[i*nbins+j*nspecies+k];
                *(float *) PyArray_GETPTR3(temp_out,i,j,k) = metal_spec.temp[i*nbins+j*nspecies+k];
            }
        }
    }
    //Build a tuple from the interp struct
	PyObject * for_return = Py_BuildValue("OOO",rho_out, vel_out, temp_out);
    //Free
    FreeLOSMemory(&species);
    free(los_table);
    free(sort_los_table);

    return for_return;
}


static PyMethodDef spectrae[] = {
  {"_SPH_Interpolate", Py_SPH_Interpolation, METH_VARARGS,
   "Find LOS density by SPH interpolation: "
   "    Arguments: nbins, box100, pos, vel, mass, u, ne, fractions, h, axis array, xx, yy, zz"
   "    "},
  {NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC
init_spectra_priv(void)
{
  Py_InitModule("_spectra_priv", spectrae);
  import_array();
}
