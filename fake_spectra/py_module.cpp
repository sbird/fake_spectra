#include <Python.h>
#include "numpy/arrayobject.h"
#include "part_int.h"

/*Wraps the flux_extractor into a python module called spectra_priv. Don't call this directly, call the python wrapper.*/

/*Check whether the passed array has type typename. Returns 1 if it doesn't, 0 if it does.*/
int check_type(PyArrayObject * arr, int npy_typename)
{
  //PyArray_DescrFromType hands out a new reference, which we own.
  PyArray_Descr * descr = PyArray_DescrFromType(npy_typename);
  const int badtype = !PyArray_EquivTypes(PyArray_DESCR(arr), descr);
  Py_DECREF(descr);
  return badtype;
}

int check_float(PyArrayObject * arr)
{
  return check_type(arr, NPY_FLOAT);
}

/* Holds the contiguous version of an input array. PyArray_GETCONTIGUOUS takes
 * a reference, and may hand back a fresh array rather than the one passed in,
 * so the reference has to be dropped on every path out of the interpolation,
 * including the error ones. */
class ContiguousArray
{
  public:
    explicit ContiguousArray(PyArrayObject * arr): contig(PyArray_GETCONTIGUOUS(arr)) {}
    ~ContiguousArray()
    {
        Py_XDECREF(contig);
    }
    //The reference is owned by exactly one of these.
    ContiguousArray(const ContiguousArray&) = delete;
    ContiguousArray& operator=(const ContiguousArray&) = delete;
    //False if getting a contiguous copy failed, in which case data() is not usable.
    bool valid() const
    {
        return contig != NULL;
    }
    void * data() const
    {
        return PyArray_DATA(contig);
    }
  private:
    PyArrayObject * contig;
};


/*****************************************************************************/
/*Interface for SPH interpolation*/
extern "C" PyObject * Py_Particle_Interpolation(PyObject *self, PyObject *args)
{
    //Things which should be from input
    int nbins, NumLos, compute_tau, kernel;
    long long Npart;
    double box100, velfac, lambda, gamma, fosc, amumass, atime, tautail;
    npy_intp size[2];
    //Input variables in np format
    PyArrayObject *pos, *vel, *dens, *temp, *h;
    PyArrayObject *cofm, *axis;

    //Get our input
    if(!PyArg_ParseTuple(args, "iiiddddddddO!O!O!O!O!O!O!", &compute_tau, &nbins, &kernel, &box100,  &velfac, &atime, &lambda, &gamma, &fosc, &amumass, &tautail, &PyArray_Type, &pos, &PyArray_Type, &vel, &PyArray_Type, &dens, &PyArray_Type, &temp, &PyArray_Type, &h, &PyArray_Type, &axis, &PyArray_Type, &cofm) )
    {
      PyErr_SetString(PyExc_AttributeError, "Incorrect arguments: use compute_tau, nbins, boxsize, velfac, atime, lambda, gamma, fosc, species mass (amu), min. tau,  pos, vel, dens, temp, h, axis, cofm\n");
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

    /*Check the shapes before reading any dimension past the first: the
     * interpolation indexes these arrays directly, so one that is too short
     * is an out of bounds read rather than an exception. Note PyArray_DIM
     * itself reads past the shape for an array of too small a rank.*/
    if(PyArray_NDIM(pos) != 2 || PyArray_DIM(pos,1) != 3)
    {
      PyErr_SetString(PyExc_ValueError, "pos must have dimensions (npart,3)\n");
      return NULL;
    }
    if(PyArray_NDIM(cofm) != 2 || PyArray_DIM(cofm,1) != 3)
    {
      PyErr_SetString(PyExc_ValueError, "cofm must have dimensions (np.size(axis),3) \n");
      return NULL;
    }
    if(PyArray_NDIM(dens) != 1 || PyArray_NDIM(h) != 1 || PyArray_NDIM(axis) != 1)
    {
      PyErr_SetString(PyExc_ValueError, "dens, h and axis must be one-dimensional\n");
      return NULL;
    }

    NumLos = PyArray_DIM(cofm,0);
    Npart = PyArray_DIM(pos,0);
    //Malloc stuff
    size[0] = NumLos;
    size[1] = nbins;

     if(Npart != PyArray_DIM(dens,0) || Npart  != PyArray_DIM(h,0))
    {
      PyErr_SetString(PyExc_ValueError, " Dens, pos and h must have the same length\n");
      return NULL;
    }

    if(NumLos != PyArray_DIM(axis,0))
    {
      PyErr_SetString(PyExc_ValueError, "cofm must have dimensions (np.size(axis),3) \n");
      return NULL;
    }

    /*Vel and temp are indexed with the same particle indices as pos, but only
     * when we are computing tau: the caller is free to pass a dummy for them
     * when all we want is the column density.*/
    if(compute_tau)
    {
        if(PyArray_NDIM(vel) != 2 || PyArray_DIM(vel,0) != Npart || PyArray_DIM(vel,1) != 3)
        {
          PyErr_SetString(PyExc_ValueError, "vel must have dimensions (npart,3) to compute tau\n");
          return NULL;
        }
        if(PyArray_NDIM(temp) != 1 || PyArray_DIM(temp,0) != Npart)
        {
          PyErr_SetString(PyExc_ValueError, "temp must be the same length as pos to compute tau\n");
          return NULL;
        }
    }

    //Initialise P from the data in the input numpy arrays.
    //Note: better be sure they are float32 in the calling function.
    //These hold the references taken by PyArray_GETCONTIGUOUS until they go
    //out of scope, so returning early does not leak them.
    ContiguousArray pos_c(pos), dens_c(dens), h_c(h), cofm_c(cofm), axis_c(axis);
    if( !pos_c.valid() || !dens_c.valid() || !h_c.valid() || !cofm_c.valid() || !axis_c.valid() ){
        PyErr_SetString(PyExc_MemoryError, "Getting contiguous copies of input arrays failed\n");
        return NULL;
    }
    float * Pos =(float *) pos_c.data();
    float * Hh= (float *) h_c.data();
    float * Dens =(float *) dens_c.data();

    double * Cofm =(double *) cofm_c.data();
    int32_t * Axis =(int32_t *) axis_c.data();
    //The index table buckets each line by its axis, so an axis outside this
    //range would be an out of bounds access rather than a bad answer.
    for(int i = 0; i < NumLos; i++){
        if(Axis[i] < 1 || Axis[i] > 3){
            PyErr_SetString(PyExc_ValueError, "axis must be 1, 2 or 3\n");
            return NULL;
        }
    }
    ParticleInterp pint(nbins, lambda, gamma, fosc, amumass, box100, velfac, atime, Cofm, Axis ,NumLos, kernel, tautail);

    PyObject * for_return;
    /* Allocate array space. This is (I hope) contiguous.
     * Note: for an array of shape (a,b), element (i,j) can be accessed as
     * [i*b+j] */
    if (compute_tau){
        ContiguousArray vel_c(vel), temp_c(temp);
        if( !vel_c.valid() || !temp_c.valid() ){
          PyErr_SetString(PyExc_MemoryError, "Getting contiguous copies of Vel and Temp failed\n");
          return NULL;
        }
        float * Vel =(float *) vel_c.data();
        float * Temp =(float *) temp_c.data();

        PyArrayObject * tau_out = (PyArrayObject *) PyArray_SimpleNew(2, size, NPY_DOUBLE);
        if ( !tau_out ){
          PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for tau\n");
          return NULL;
        }
        double * tau = (double *) PyArray_DATA(tau_out);
        PyArray_FILLWBYTE(tau_out, 0);
        //Do the work
        pint.compute_tau(tau, Pos, Vel, Dens, Temp, Hh, Npart);

        //Build a tuple from the interp struct
        for_return = Py_BuildValue("O", tau_out);
        Py_DECREF(tau_out);
    }
    else{
        PyArrayObject * colden_out = (PyArrayObject *) PyArray_SimpleNew(2, size, NPY_DOUBLE);
        if ( !colden_out ){
          PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for colden\n");
          return NULL;
        }
        double * colden = (double *) PyArray_DATA(colden_out);
        //Initialise output arrays to 0.
        PyArray_FILLWBYTE(colden_out, 0);

        //Do the work
        pint.compute_colden(colden, Pos, Dens, Hh, Npart);

        //Build a tuple from the interp struct
        for_return = Py_BuildValue("O", colden_out);
        Py_DECREF(colden_out);
    }

    return for_return;
}

static PyMethodDef spectrae[] = {
  {"_Particle_Interpolate", Py_Particle_Interpolation, METH_VARARGS,
   "Find absorption or column density by interpolating particles. "
   "    Arguments: compute_tau nbins, boxsize, velfac, atime, lambda, gamma, fosc, species mass (amu), pos, vel, dens, temp, h, axis, cofm"
   "    "},
  {NULL, NULL, 0, NULL},
};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_spectra_priv", /* m_name */
  "C functions for accelerating spectral work",      /* m_doc */
  -1,                  /* m_size */
  spectrae,            /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC
PyInit__spectra_priv(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);
    import_array();
    if (m == NULL)
        return NULL;
    return m;
}
