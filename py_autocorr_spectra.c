/* Python module to calculate the 3D autocorrelation function along some spectra */

#include <Python.h>
#include "numpy/arrayobject.h"
#include <omp.h>

/*Check whether the passed array has type typename. Returns 1 if it doesn't, 0 if it does.*/
int check_type(PyArrayObject * arr, int npy_typename)
{
  return !PyArray_EquivTypes(PyArray_DESCR(arr), PyArray_DescrFromType(npy_typename));
}

//Compute square of distance between two sightlines, assuming same axis
double spec_distance2(const double * a, const double * b)
{
    const int dif[2] = {(*(a)-*(b)), (*(a+1)-*(b+1))};
    return dif[0]*dif[0]+dif[1]*dif[1];
}

/*Find the autocorrelation function from a list of spectra
   Spectra are assumed to be along the same axis.
   slist - list of quantity along spectra to autocorrelate. npix * nspectra
   spos -  positions of the spectra: 2x nspectra: (x, y).
   nbins - number of bins in output autocorrelation function
   pixsz - Size of a pixel in units of the cofm.
*/
PyObject * _autocorr_spectra(PyObject *self, PyObject *args)
{
    PyArrayObject *slist, *spos;
    int nbins;
    double pixsz;
    if(!PyArg_ParseTuple(args, "O!O!di",&PyArray_Type, &slist, &PyArray_Type, &spos, &pixsz, &nbins) )
    {
        PyErr_SetString(PyExc_AttributeError, "Incorrect arguments: use slist, spos, pixsz, nbins\n");
        return NULL;
    }
    if(check_type(slist, NPY_DOUBLE) || check_type(spos, NPY_DOUBLE))
    {
          PyErr_SetString(PyExc_AttributeError, "Input arrays are not float64.\n");
          return NULL;
    }
    npy_intp npix = PyArray_DIM(slist,0);
    npy_intp nspectra = PyArray_DIM(slist,1);
    npy_intp npnbins = nbins;
    //Bin autocorrelation, must cover sqrt(dims)*size
    //so each bin has size sqrt(dims)*size /nbins
    const int nproc = omp_get_num_procs();
    double autocorr_C[nproc][nbins];
    int modecount_C[nproc][nbins];
    memset(modecount_C,0,nproc*nbins*sizeof(int));
    memset(autocorr_C,0,nproc*nbins*sizeof(int));
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
	const double binsz = nbins / (npix * pixsz*sqrt(3));
        #pragma omp for
        for(int b=0; b<nspectra; b++){
            for(int a=0; a<nspectra; a++){
	      double sdist2 = spec_distance2((double *) PyArray_GETPTR2(spos,0, a), (double *) PyArray_GETPTR2(spos,0, b));
	      double * speca = (double *) PyArray_GETPTR2(slist,0, a);
	      double * specb = (double *) PyArray_GETPTR2(slist,0, b);
	      for(int bb=0; bb<npix; bb++){
		for(int aa=0; aa<npix; aa++){
                        double rr = sqrt(sdist2+ (bb-aa)*(bb-aa)*pixsz);
                        //Which bin to add this one to?
                        int cbin = floor(rr * binsz);
                        autocorr_C[tid][cbin]+=speca[aa]*specb[bb];
			modecount_C[tid][cbin]+=1;
		}
	      }
            }
        }
    }
    PyArrayObject *autocorr = (PyArrayObject *) PyArray_SimpleNew(1,&npnbins,NPY_DOUBLE);
    PyArray_FILLWBYTE(autocorr, 0);
    PyArrayObject *modecount = (PyArrayObject *) PyArray_SimpleNew(1,&npnbins,NPY_INT);
    PyArray_FILLWBYTE(modecount, 0);

    for(int tid=0; tid < nproc; tid++){
        for(int nn=0; nn< nbins; nn++){
            *(double *)PyArray_GETPTR1(autocorr,nn)+=autocorr_C[tid][nn];
	    *(int *)PyArray_GETPTR1(modecount,nn)+=modecount_C[tid][nn];
        }
    }
    return Py_BuildValue("OO", modecount, autocorr);
}


static PyMethodDef __autocorr[] = {
  {"autocorr_spectra", _autocorr_spectra, METH_VARARGS,
   "Find the autocorrelation function from a list of spectra"
   "Spectra are assumed to be along the same axis."
   "slist - list of quantity along spectra to autocorrelate. npix * nspectra"
   "spos -  positions of the spectra: 2x nspectra: (x, y). "
   "pixsz - Size of a pixel in units of the cofm."
   "nbins - number of bins in output autocorrelation function"
   "    "},
  {NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC
init_autocorr_spectra_priv(void)
{
  Py_InitModule("_autocorr_spectra_priv", __autocorr);
  import_array();
}
