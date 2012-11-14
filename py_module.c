#include <Python.h>
#include "numpy/arrayobject.h"
#include "global_vars.h"

/*Wraps the flux_extractor into a python module called spectra_priv. Don't call this directly, call the python wrapper.*/

PyObject * Py_SPH_Interpolation(PyObject *self, PyObject *args)
{
    //Things which should be from input
    int nbins, NumLos;
    long long Npart;
    double box100;
    npy_intp size;
    //Input variables in np format
    PyObject * data;
    PyArrayObject *pos, *vel, *mass, *u, *nh0, *ne, *h;
    PyArrayObject *xx, *yy, *zz, *axis;
    PyArrayObject * rho_out, *temp_out, *vel_out;

    //For storing output
    interp species;

    //Temp variables
    int nxx=0, i;
    los *los_table=NULL;
    sort_los *sort_los_table=NULL;
    struct particle_data P;
    //Get our input
    if(!PyArg_ParseTuple(args, "idO!O!O!O!O!O!O!O!O!O!O!",&nbins, &box100, &PyArray_Type, &pos, &PyArray_Type, &vel, &PyArray_Type, &mass, &PyArray_Type, &u, &PyArray_Type, &nh0, &PyArray_Type, &ne, &PyArray_Type, &h, &PyArray_Type, &axis, &PyArray_Type, &xx, &PyArray_Type, &yy, &PyArray_Type, &zz) )
      return NULL;

    NumLos = xx->dimensions[0];
    Npart = pos->dimensions[0];
    //Malloc stuff
    los_table=malloc(NumLos*sizeof(los));
    sort_los_table=malloc(NumLos*sizeof(sort_los));
    InitLOSMemory(&species,NumLos, nbins);
    //Initialise P from the data in the input numpy arrays.
    P.Pos = (float *) pos->data;
    P.Vel = (float *) vel->data;
    P.Mass = (float *) mass->data;
    P.U = (float *) u->data;
    P.NH0 = (float *) nh0->data;
    P.Ne = (float *) ne->data;
    P.h = (float *) h->data;

    //Initialise los_table from input
    for(i=0; i< NumLos; i++){
        los_table[i].axis = ((int *) axis->data)[i];
        los_table[i].xx = ((float *) xx->data)[i];
        los_table[i].yy = ((float *) yy->data)[i];
        los_table[i].zz = ((float *) zz->data)[i];
    }

    //Do the work
    populate_sort_los_table(los_table, NumLos, sort_los_table, &nxx);
    SPH_Interpolation(NULL,&species,NULL,nbins, Npart, NumLos, box100, los_table,sort_los_table,nxx, &P);

    size = NumLos*nbins;
    rho_out = (PyArrayObject *) PyArray_SimpleNew(1, &size, NPY_FLOAT);
    vel_out = (PyArrayObject *) PyArray_SimpleNew(1, &size, NPY_FLOAT);
    temp_out = (PyArrayObject *) PyArray_SimpleNew(1, &size, NPY_FLOAT);
    /*Is there a better way to do this?*/
    for(i=0; i< NumLos*nbins; i++){
        ((float *) rho_out->data)[i] = species.rho[i];
        ((float *) vel_out->data)[i] = species.veloc[i];
        ((float *) temp_out->data)[i] = species.temp[i];
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
   "    Arguments: nbins, pos, vel, mass, u, nh0, ne, h, axis array, xx, yy, zz"
   "    "},
  {NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC
init_spectra_priv(void)
{
  Py_InitModule("_spectra_priv", spectrae);
  import_array();
}
