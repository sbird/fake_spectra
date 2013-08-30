#include <Python.h>
#include "numpy/arrayobject.h"
#include "part_int.h"
#include <set>

/*Wraps the flux_extractor into a python module called spectra_priv. Don't call this directly, call the python wrapper.*/

/*Check whether the passed array has type typename. Returns 1 if it doesn't, 0 if it does.*/
int check_type(PyArrayObject * arr, int npy_typename)
{
  return !PyArray_EquivTypes(PyArray_DESCR(arr), PyArray_DescrFromType(npy_typename));
}

int check_float(PyArrayObject * arr)
{
  return check_type(arr, NPY_FLOAT);
}


/* When handed a list of particles,
 * return a list of bools with True for those nearby to a sightline*/
extern "C" PyObject * Py_near_lines(PyObject *self, PyObject *args)
{
    int NumLos;
    long long Npart;
    double box100;
    PyArrayObject *cofm, *axis, *pos, *hh, *is_a_line;
    PyObject *out;

    if(!PyArg_ParseTuple(args, "dO!O!O!O!",&box100,  &PyArray_Type, &pos, &PyArray_Type, &hh, &PyArray_Type, &axis, &PyArray_Type, &cofm) )
      return NULL;

    NumLos = PyArray_DIM(cofm,0);
    Npart = PyArray_DIM(pos,0);

    if(NumLos != PyArray_DIM(axis,0) || 3 != PyArray_DIM(cofm,1)){
      PyErr_SetString(PyExc_ValueError, "cofm must have dimensions (np.size(axis),3) \n");
      return NULL;
    }
    if(check_type(cofm, NPY_DOUBLE) || check_type(axis,NPY_INT)){
      PyErr_SetString(PyExc_ValueError, "cofm must have 64-bit float type and axis must be a 32-bit integer\n");
      return NULL;
    }
    if(check_float(pos) || check_float(hh)){
       PyErr_SetString(PyExc_TypeError, "pos and h must have 32-bit float type\n");
       return NULL;
    }

    //Setup los_tables
    //PyArray_GETCONTIGUOUS increments the reference count of the object,
    //so to avoid leaking we need to save the PyArrayObject pointer.
    cofm = PyArray_GETCONTIGUOUS(cofm);
    axis = PyArray_GETCONTIGUOUS(axis);
    double * Cofm =(double *) PyArray_DATA(cofm);
    int32_t * Axis =(int32_t *) PyArray_DATA(axis);
    IndexTable sort_los_table(Cofm, Axis, NumLos, box100);

    //Set of particles near a line
    std::set<int> near_lines;
    //find lists
    //DANGER: potentially huge allocation
    pos = PyArray_GETCONTIGUOUS(pos);
    hh = PyArray_GETCONTIGUOUS(hh);
    const float * Pos =(float *) PyArray_DATA(pos);
    const float * h = (float *) PyArray_DATA(hh);
    #pragma omp parallel for
    for(long long i=0; i < Npart; i++){
	    std::map<int, double> nearby=sort_los_table.get_near_lines(&(Pos[3*i]),h[i]);
        if(nearby.size()>0){
           #pragma omp critical
           {
              near_lines.insert(i);
           }
        }
    }
    //Copy data into python
    npy_intp size = near_lines.size();
    is_a_line = (PyArrayObject *) PyArray_SimpleNew(1, &size, NPY_INT);
    int i=0;
    for (std::set<int>::const_iterator it = near_lines.begin(); it != near_lines.end() && i < size; ++it, ++i){
            *(npy_int *)PyArray_GETPTR1(is_a_line,i) = (*it);
    }
    out = Py_BuildValue("O", is_a_line);
    Py_DECREF(is_a_line);
    //Because PyArray_GETCONTIGUOUS incremented the reference count,
    //and may have made an allocation, in which case this does not point to what it used to.
    Py_DECREF(pos);
    Py_DECREF(hh);
    Py_DECREF(cofm);
    Py_DECREF(axis);
    return out;
}

/*****************************************************************************/
/*Interface for SPH interpolation*/
extern "C" PyObject * Py_Particle_Interpolation(PyObject *self, PyObject *args)
{
    //Things which should be from input
    int nbins, NumLos, compute_tau;
    long long Npart;
    double box100, velfac, lambda, gamma, fosc, amumass, atime;
    npy_intp size[2];
    //Input variables in np format
    PyArrayObject *pos, *vel, *dens, *temp, *h;
    PyArrayObject *cofm, *axis;

    //Get our input
    if(!PyArg_ParseTuple(args, "iidddddddO!O!O!O!O!O!O!", &compute_tau, &nbins, &box100,  &velfac, &atime, &lambda, &gamma, &fosc, &amumass, &PyArray_Type, &pos, &PyArray_Type, &vel, &PyArray_Type, &dens, &PyArray_Type, &temp, &PyArray_Type, &h, &PyArray_Type, &axis, &PyArray_Type, &cofm) )
    {
      PyErr_SetString(PyExc_AttributeError, "Incorrect arguments: use compute_tau, nbins, boxsize, velfac, atime, lambda, gamma, fosc, species mass (amu), pos, vel, dens, temp, h, axis, cofm\n");
      return NULL;
    }

    //Check that our input has the right types
    if(check_float(pos) || check_float(vel) || check_float(dens) || check_float(temp) || check_float(h)){
       PyErr_SetString(PyExc_TypeError, "One of the data arrays does not have 32-bit float type\n");
       return NULL;
    }
    if(check_type(cofm,NPY_DOUBLE)){
      PyErr_SetString(PyExc_TypeError, "Sightline positions must have 64-bit float type\n");
      return NULL;
    }
    if(check_type(axis, NPY_INT32)){
      PyErr_SetString(PyExc_TypeError, "Axis must be a 32-bit integer\n");
      return NULL;
    }

    NumLos = PyArray_DIM(cofm,0);
    Npart = PyArray_DIM(pos,0);
    //Malloc stuff
    size[0] = NumLos;
    size[1] = nbins;

    if(NumLos != PyArray_DIM(axis,0) || 3 != PyArray_DIM(cofm,1))
    {
      PyErr_SetString(PyExc_ValueError, "cofm must have dimensions (np.size(axis),3) \n");
      return NULL;
    }


    /* Allocate array space. This is (I hope) contiguous.
     * Note: for an array of shape (a,b), element (i,j) can be accessed as
     * [i*b+j] */
    PyArrayObject * colden_out = (PyArrayObject *) PyArray_SimpleNew(2, size, NPY_DOUBLE);
    PyArrayObject * tau_out;
    double * tau;
    if (compute_tau){
        tau_out = (PyArrayObject *) PyArray_SimpleNew(2, size, NPY_DOUBLE);
        tau = (double *) PyArray_DATA(tau_out);
    }
    else{
        tau_out = (PyArrayObject *) PyArray_SimpleNew(0, size, NPY_DOUBLE);
        tau = NULL;
    }


    if ( !colden_out || !tau_out ){
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for output arrays\n");
        return NULL;
    }

    //Initialise output arrays to 0.
    PyArray_FILLWBYTE(colden_out, 0);
    PyArray_FILLWBYTE(tau_out, 0);

    //Here comes the cheat
    double * colden = (double *) PyArray_DATA(colden_out);

    //Initialise P from the data in the input numpy arrays.
    //Note: better be sure they are float32 in the calling function.
    //PyArray_GETCONTIGUOUS increments the reference count of the object,
    //so to avoid leaking we need to save the PyArrayObject pointer.
    pos = PyArray_GETCONTIGUOUS(pos);
    vel = PyArray_GETCONTIGUOUS(vel);
    dens = PyArray_GETCONTIGUOUS(dens);
    temp = PyArray_GETCONTIGUOUS(temp);
    h = PyArray_GETCONTIGUOUS(h);
    float * Pos =(float *) PyArray_DATA(pos);
    float * Hh= (float *) PyArray_DATA(h);
    float * Vel =(float *) PyArray_DATA(vel);
    float * Dens =(float *) PyArray_DATA(dens);
    float * Temp =(float *) PyArray_DATA(temp);

    cofm = PyArray_GETCONTIGUOUS(cofm);
    axis = PyArray_GETCONTIGUOUS(axis);
    double * Cofm =(double *) PyArray_DATA(cofm);
    int32_t * Axis =(int32_t *) PyArray_DATA(axis);
    if( !Pos || !Vel || !Dens || !Temp || !Hh || !Cofm || !Axis ){
        PyErr_SetString(PyExc_MemoryError, "Getting contiguous copies of input arrays failed\n");
        return NULL;
    }
    ParticleInterp pint(tau, colden, nbins, lambda, gamma, fosc, amumass, box100, velfac, atime, Cofm, Axis ,NumLos);
    //Do the work
    pint.do_work(Pos, Vel, Dens, Temp, Hh, Npart);

    //Build a tuple from the interp struct
    PyObject * for_return;
	for_return = Py_BuildValue("OO", tau_out, colden_out);

    //Decrement the refcount
    if(tau_out){
        Py_DECREF(tau_out);
    }
    Py_DECREF(colden_out);
    //Because PyArray_GETCONTIGUOUS incremented the reference count,
    //and may have made an allocation, in which case this does not point to what it used to.
    Py_DECREF(pos);
    Py_DECREF(vel);
    Py_DECREF(dens);
    Py_DECREF(temp);
    Py_DECREF(h);
    Py_DECREF(cofm);
    Py_DECREF(axis);

    return for_return;
}

static PyMethodDef spectrae[] = {
  {"_Particle_Interpolate", Py_Particle_Interpolation, METH_VARARGS,
   "Find absorption and column density by interpolating particles. "
   "    Arguments: compute_tau nbins, boxsize, velfac, atime, lambda, gamma, fosc, species mass (amu), pos, vel, dens, temp, h, axis, cofm"
   "    "},
  {"_near_lines", Py_near_lines,METH_VARARGS,
   "Give a list of particles and sightlines, "
   "return a list of booleans for those particles near "
   "a sightline."
   "   Arguments: box, pos, h, axis, cofm"},
  {NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC
init_spectra_priv(void)
{
  Py_InitModule("_spectra_priv", spectrae);
  import_array();
}
