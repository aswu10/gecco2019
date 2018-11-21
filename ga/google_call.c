/*
 * C code that calls Python function gmap_dist in file
 * google_func.py.
 * 
 * It creates C arrays of source and destination locations.  Each
 * element is a latitude or longitude.  It then turns them into
 * Python tuples to pass to the Python function, gmap_dist.
 * 
 * The result returned from gmap_dist is a list of two lists:
 *       a list of distances and
 *       a list of times
 * Each list contains one value for each element in the cross product
 * of sources with destinations
 * 
 * It unpacks these Python lists into C arrays: dists, and times
 *
 */

#include<stdio.h>
#include<unistd.h>
#include<Python.h>

PyObject *pName, *pModule, *pFunc, *pDict, *pValue, *pArgs, *dTuple, *sTuple;
PyObject *py_dists, *py_times;


double google_dist(double source[2], double dest[2])
{    
    // create the Python objects needed for function call and parameters
    // PyObject *pName, *pModule, *pFunc, *pDict, *pValue, *pArgs, *dTuple, *sTuple;

    // Because C doesn't provide a function to determine the number of elements
    // in an array, we use sizeof to dtermine this.
    // elt_size: number of bytes in a double
    int elt_size = sizeof(double);

    // source_num: number of elements in the source array
    //             each location consists of two elements: lat and lon
    // int source_num = sizeof(source)/elt_size;
    int source_num = 2;
    
    // dest_num: number of elements in the dest array
    //           each location consists of two elements: lat and lon
    // int dest_num = sizeof(dest)/elt_size;
    int dest_num = 2;
    
    // result_num: number of elements in the cross product
    //             divide by 2 because each location is two elements: lat and lon
    int result_num = (source_num/2) * (dest_num/2);
    
    // C arrays to hold the returned results after conversion from Python objects
    double dists[result_num];
    double times[result_num];
    
    // required by the C Python library
    // Py_Initialize();

    // pArgs is passed as argument to gmaps_dist
    // it is s tuple of two tuples:
    //      sTuple: the source locations (comma-separated lats and lons)
    //      dTuple: the destination locations (comma-separated lats and lons)
    pArgs = PyTuple_New(2);
    sTuple = PyTuple_New(source_num);
    dTuple = PyTuple_New(dest_num);
    
    // build sTuple from the array of sources
    // printf("Creating sTuple:\n");
    for(int i = 0; i < source_num; i++)
    {
        // printf("    source[%d]: %f\n", i, source[i]);
        pValue = PyFloat_FromDouble(source[i]);
        PyTuple_SetItem(sTuple, i, pValue);
    }
    
    // build dTuple from the array of destinations
    for(int i = 0; i < dest_num; i++)
    {
        pValue = PyFloat_FromDouble(dest[i]);
        PyTuple_SetItem(dTuple, i, pValue);
    }
    
    // add sTuple and dTuple to pArgs
    PyTuple_SetItem(pArgs, 0, sTuple);
    PyTuple_SetItem(pArgs, 1, dTuple);
    
    // create Python string with name of Python code file
    pName = PyString_FromString("google_func");

    // pModule is a reference to the imported source code file
    pModule = PyImport_Import(pName);
    
    // creates a dictionary of data about the imported module
    pDict = PyModule_GetDict(pModule);

    // get reference to the function gmap_dist from the dictionary
    pFunc = PyDict_GetItemString(pDict, "gmap_dist");
    
    // check if the function is callable
    if(PyCallable_Check(pFunc))
    {
        usleep(10000);
        // call the function and save the returned object in pValue
        // pValue is a Python list containing 2 Python lists
        pValue = PyObject_CallObject(pFunc, pArgs);
        
        if(pValue != NULL)
        {
            // get the two Python lists from pValue
            // PyObject *py_dists, *py_times;
            py_dists = PyList_GET_ITEM(pValue, 0);
            py_times = PyList_GET_ITEM(pValue, 1);

            printf("google.c: got dists and times from Python\n");
            // extract results from the Python lists and store in C arrays
            for(int i = 0; i < result_num; i++)
            {
                dists[i] = (double)PyInt_AsLong(PyList_GET_ITEM(py_dists, i));
                times[i] = (double)PyInt_AsLong(PyList_GET_ITEM(py_times, i));
            }
            
//            Py_DECREF(py_dists);
//            Py_DECREF(py_times);
        }
        else
            printf("pValue is NULL\n");
    }
    else
    {
        PyErr_Print();
    }

//    Py_DECREF(pModule);
//    Py_DECREF(pName);
//    Py_DECREF(pArgs);
//    Py_DECREF(sTuple);
//    Py_DECREF(dTuple);

    // Py_Finalize();
    
    return dists[0];
}


double google_dist_single(double source[2], double dest[2])
{    
    // create the Python objects needed for function call and parameters
    PyObject *pName, *pModule, *pFunc, *pDict, *pValue, *pArgs, *dTuple, *sTuple;

    // the lists of sources and destinations
    // double source[] = {43.912454,-91.181146};
    // double dest[] = {43.813145,-91.233841,43.907036,-91.236639};
    
    // Because C doesn't provide a function to determine the number of elements
    // in an array, we use sizeof to dtermine this.
    // elt_size: number of bytes in a double
    int elt_size = sizeof(double);
    // source_num: number of elements in the source array
    //             each location consists of two elements: lat and lon
    int source_num = 2;
    // dest_num: number of elements in the dest array
    //           each location consists of two elements: lat and lon
    int dest_num = 2;
    // result_num: number of elements in the cross product
    //             divide by 2 because each location is two elements: lat and lon
    int result_num = 1;
    
    // C arrays to hold the returned results after conversion from Python objects
    double dist;
    double time;
    
    // required by the C Python library
    Py_Initialize();

    // pArgs is passed as argument to gmaps_dist
    // it is s tuple of two tuples:
    //      sTuple: the source locations (comma-separated lats and lons)
    //      dTuple: the destination locations (comma-separated lats and lons)
    pArgs = PyTuple_New(2);
    sTuple = PyTuple_New(source_num);
    dTuple = PyTuple_New(dest_num);
    
    // build sTuple from the array of sources
    for(int i = 0; i < source_num; i++)
    {
        pValue = PyFloat_FromDouble(source[i]);
        PyTuple_SetItem(sTuple, i, pValue);
    }
    
    // build dTuple from the array of destinations
    for(int i = 0; i < dest_num; i++)
    {
        pValue = PyFloat_FromDouble(dest[i]);
        PyTuple_SetItem(dTuple, i, pValue);
    }
    
    // add sTuple and dTuple to pArgs
    PyTuple_SetItem(pArgs, 0, sTuple);
    PyTuple_SetItem(pArgs, 1, dTuple);
    
    // create Python string with name of Python code file
    pName = PyString_FromString("google_func");

    // pModule is a reference to the imported source code file
    pModule = PyImport_Import(pName);
    
    // creates a dictionary of data about the imported module
    pDict = PyModule_GetDict(pModule);

    // get reference to the function gmap_dist from the dictionary
    pFunc = PyDict_GetItemString(pDict, "gmap_dist");
    
    // check if the function is callable
    if(PyCallable_Check(pFunc))
    {
        // call the function and save the returned object in pValue
        // pValue is a Python list containing 2 Python lists
        pValue = PyObject_CallObject(pFunc, pArgs);
        
        if(pValue != NULL)
        {
            // get the two Python lists from pValue
            PyObject *py_dists, *py_times;
            py_dists = PyList_GET_ITEM(pValue, 0);
            py_times = PyList_GET_ITEM(pValue, 1);

            // extract results from the Python lists and store in C arrays
            dist = (double)PyInt_AsLong(PyList_GET_ITEM(py_dists, 0));
            printf("goole dist: %f\n", dist);
            time = (double)PyInt_AsLong(PyList_GET_ITEM(py_times, 0));
            printf("google time: %f\n", time);
            
            Py_DECREF(py_dists);
            Py_DECREF(py_times);
            printf("after decref of py_dists and py_times\n");
        }
        else
            printf("pValue is NULL\n");
    }
    else
    {
        PyErr_Print();
    }

    Py_DECREF(pModule);
    Py_DECREF(pName);
    Py_DECREF(pArgs);
    Py_DECREF(sTuple);
    Py_DECREF(dTuple);

    printf("about to finalize\n");
    Py_Finalize();

    printf("about to return dist\n");
    return dist;
}
