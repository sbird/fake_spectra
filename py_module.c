#include <Python.h>
#include "numpy/arrayobject.h"
#include "global_vars.h"

/*Wraps the flux_extractor into a python module called spectra_priv. Don't call this directly, call the python wrapper.*/

/*****************************************************************************/
/*Interface for SPH interpolation*/
PyObject * Py_SPH_Interpolation(PyObject *self, PyObject *args)
{
    //Things which should be from input
    int nbins, NumLos, nspecies, rho_H;
    long long Npart;
    double box100;
    double * rho_H_data=NULL;
    npy_intp size[3];
    //Input variables in np format
    PyArrayObject *pos, *vel, *mass, *u, *ne, *h, *fractions;
    PyArrayObject *cofm, *axis;
    PyArrayObject *rho_H_out=NULL, *rho_out, *temp_out, *vel_out;

    //For storing output
    interp species;

    //Temp variables
    int nxx=0, i;
    los *los_table=NULL;
    sort_los *sort_los_table=NULL;
    struct particle_data P;
    //Get our input
    if(!PyArg_ParseTuple(args, "iidO!O!O!O!O!O!O!O!O!",&rho_H, &nbins, &box100,  &PyArray_Type, &pos, &PyArray_Type, &vel, &PyArray_Type, &mass, &PyArray_Type, &u, &PyArray_Type, &ne, &PyArray_Type, &fractions, &PyArray_Type, &h, &PyArray_Type, &axis, &PyArray_Type, &cofm) )
      return NULL;

    NumLos = PyArray_DIM(cofm,0);
    Npart = PyArray_DIM(pos,0);
    //NOTE if nspecies == 1, fractions must have shape [N,1], rather than [N]
    nspecies = PyArray_DIM(fractions, 1);
    //Malloc stuff
    los_table=malloc(NumLos*sizeof(los));
    sort_los_table=malloc(NumLos*sizeof(sort_los));
    size[0] = NumLos;
    size[1] = nbins;
    //Number of metal species
    size[2] = nspecies;
    
    /* Allocate array space. This is (I hope) contiguous.
     * Note: for an array of shape (a,b), element (i,j) can be accessed as
     * [i*b+j] */
    rho_out = (PyArrayObject *) PyArray_SimpleNew(3, size, NPY_DOUBLE);
    vel_out = (PyArrayObject *) PyArray_SimpleNew(3, size, NPY_DOUBLE);
    temp_out = (PyArrayObject *) PyArray_SimpleNew(3, size, NPY_DOUBLE);

    if ( !rho_out || !vel_out || !temp_out)
        return NULL;

    //Initialise output arrays to 0.
    PyArray_FILLWBYTE(rho_out, 0);
    PyArray_FILLWBYTE(vel_out, 0);
    PyArray_FILLWBYTE(temp_out, 0);

    if(rho_H){
        rho_H_out = (PyArrayObject *) PyArray_SimpleNew(2, size, NPY_DOUBLE);
        PyArray_FILLWBYTE(rho_H_out, 0);
        rho_H_data = PyArray_DATA(rho_H_out);
    }

    //Here comes the cheat
    species.rho = (double *) PyArray_DATA(rho_out);
    species.temp = (double *) PyArray_DATA(temp_out);
    species.veloc = (double *) PyArray_DATA(vel_out);

    //Initialise P from the data in the input numpy arrays. 
    //Note: better be sure they are float32 in the calling function.
    P.Pos =(float *) PyArray_DATA(PyArray_GETCONTIGUOUS(pos));
    P.Vel =(float *) PyArray_DATA(PyArray_GETCONTIGUOUS(vel));
    P.Mass =(float *) PyArray_DATA(PyArray_GETCONTIGUOUS(mass));
    P.U =(float *) PyArray_DATA(PyArray_GETCONTIGUOUS(u));
    P.Ne =(float *) PyArray_DATA(PyArray_GETCONTIGUOUS(ne));
    P.h =(float *) PyArray_DATA(PyArray_GETCONTIGUOUS(h));
    P.fraction =(float *) PyArray_DATA(PyArray_GETCONTIGUOUS(fractions));

    //Initialise los_table from input
    for(i=0; i< NumLos; i++){
        los_table[i].axis = *(int *) PyArray_GETPTR1(axis,i);
        los_table[i].xx = *(double *) PyArray_GETPTR2(cofm,i,0);
        los_table[i].yy = *(double *) PyArray_GETPTR2(cofm,i,1);
        los_table[i].zz = *(double *) PyArray_GETPTR2(cofm,i,2);
    }

    //Do the work
    populate_sort_los_table(los_table, NumLos, sort_los_table, &nxx);
    SPH_Interpolation(rho_H_data,&species, nspecies, nbins, Npart, NumLos, box100, los_table,sort_los_table,nxx, &P);

    //Build a tuple from the interp struct
    PyObject * for_return;
    if(rho_H)
	    for_return = Py_BuildValue("OOOO",rho_H_out, rho_out, vel_out, temp_out);
    else
	    for_return = Py_BuildValue("OOO", rho_out, vel_out, temp_out);

    //Free
    free(los_table);
    free(sort_los_table);

    return for_return;
}


static PyMethodDef spectrae[] = {
  {"_SPH_Interpolate", Py_SPH_Interpolation, METH_VARARGS,
   "Find LOS density by SPH interpolation: "
   "    Arguments: nbins, box100, pos, vel, mass, u, ne, fractions, h, axis, cofm"
   "    "},
  {NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC
init_spectra_priv(void)
{
  Py_InitModule("_spectra_priv", spectrae);
  import_array();
}
