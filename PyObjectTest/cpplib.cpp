
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <thread>

#include "Python.h"

using std::cout;
using std::endl;

// compile as C code, but use C++ syntax for sanity reasons
extern "C"{

void mycall(PyObject *python_callable, int n){

    cout << "Hello from the lib: mycall() called from tasmanian_cpp_method()" << endl;

    for(int i=0; i<n; i++){ // just a simple loop to make multiple lambda calls
        double x = (double) (i+1);
        cout << "callback with x = " << x << endl;

        // passing parameters to Python means they have to be "packed"
        // make a simple python floating point value
        PyObject* val =  Py_BuildValue("f", x);
        // pack one floating point value
        PyObject* newt = PyTuple_Pack(1, val);

        // this is how to call a callable (lambda or function)
        PyObject *result = PyObject_CallObject(python_callable, newt);

        // convert the result to a floating point
        double res = PyFloat_AsDouble(result);
        cout << "Library got " << res << endl;

        // this is Python version of "delete" and "free()"
        Py_DECREF(result);
        Py_DECREF(newt);
        Py_DECREF(val);
    }
}
// pretend this is a wrapper to the Tasmanian method
static PyObject *tasmanian_cpp_method( PyObject *, PyObject *args ){

    // we are getting some args
    // parse the object to convert the values to C++ values
    int n = 0;
    PyObject *python_callable = nullptr;
    PyArg_ParseTuple(args, "O|i", &python_callable, &n); // read one object and one integer

    mycall(python_callable, n); // go into work mode

    return Py_BuildValue(""); // Python equivalent of void
}


// callable methods cannot deviate too much from this signature
// the second class is still there, I guess passing nullptr or empty typle
// because the function is registered as: METH_NOARGS
static PyObject *tasmanian_dummy_test( PyObject *, PyObject * ){

    cout << "the C++ code reports a call to tasmanian_dummy_test()" << endl;

    return Py_BuildValue(""); // Python equivalent of void
}

static void test_numpy(int dims){

    cout << "\n\n" << endl;
    cout << std::scientific; cout.precision(8);
    cout << " enter with dims = " << dims << endl;

    PyObject *tx = PyTuple_New(dims);

    std::vector<double> x(dims);
    for(int i=0; i<dims; i++) x[i] = (double)(i+1);

    for(int i=0; i<dims; i++){
        PyObject* val =  Py_BuildValue("f", x[i]);
        PyTuple_SetItem(tx, i, val);
    }

    PyObject *module = PyImport_ImportModule("pytestmod");
    std::string s = "model_proxy";
    PyObject *callable_proxy = PyObject_GetAttrString(module, s.c_str());

    PyObject* call_args = PyTuple_Pack(1, tx);
    PyObject* result = PyObject_CallObject(callable_proxy, call_args);

    if (result == nullptr) throw std::runtime_error("ouch");

    for(int i=0; i<dims; i++){
        double r = PyFloat_AsDouble(PyTuple_GET_ITEM(result, i));
        cout << " at i: " << i << "  " << r << "  expected: " << std::sin(x[i]) << endl;
    }


    //PyGILState_STATE gstate = PyGILState_Ensure();

    //std::vector<std::thread> works(2);
    //for(auto &w : works)
    //    w = std::thread([&]()->void{
    //        cout << "making new interp" << endl;
    //        //PyGILState_STATE gstate = PyGILState_Ensure();
    //        PyThreadState* pythread = Py_NewInterpreter();
    //        PyGILState_STATE gstate = PyGILState_Ensure();
    //
    //        PyObject* tresult = PyObject_CallObject(callable_proxy, call_args);
    //
    //        if (tresult == nullptr) throw std::runtime_error("ouch");
    //
    //        for(int i=0; i<dims; i++){
    //            double r = PyFloat_AsDouble(PyTuple_GET_ITEM(tresult, i));
    //            cout << " at i: " << i << "  " << r << "  expected: " << std::sin(x[i]) << endl;
    //        }
    //
    //        PyGILState_Release(gstate);
    //        Py_EndInterpreter(pythread);
    //        //PyGILState_Release(gstate);
    //    });
    //
    //for(auto &w : works) w.join();

    //PyGILState_Release(gstate);

    Py_DECREF(result);
    Py_DECREF(call_args);
    Py_DECREF(callable_proxy);
    Py_DECREF(module);
    Py_DECREF(tx);
}

static PyObject *tasmanian_array_test( PyObject *, PyObject *args ){

    // we are getting some args
    // parse the object to convert the values to C++ values
    int dims;
    PyObject *obj = nullptr;
    PyArg_ParseTuple(args, "i|O", &dims, &obj); // read one object and one integer

    size_t pntr = PyLong_AsLong(obj);

    std::vector<double> *v = reinterpret_cast<std::vector<double>*>(pntr);
    cout << v->size() << endl;
    cout << (*v)[0] << endl;

    try{
        test_numpy(dims); // go into work mode
    }catch(std::runtime_error &){
        return nullptr;
    }

    return Py_BuildValue(""); // Python equivalent of void
}

void* makevec(){ return (void*) new std::vector<double>(7, 11.0); }

// in order for Python to call a C++ method, the method has to be registered
// using function pointer and specify simple args
PyMethodDef methoddef{ "tasmanian_cpp_method", &tasmanian_cpp_method, METH_VARARGS, "no-docs"};

// test to register two methods
static
PyMethodDef methoddefs[4] = {{ "tasmanian_cpp_method", &tasmanian_cpp_method, METH_VARARGS, "no-docs"},
                             { "tasmanian_dummy_test", &tasmanian_dummy_test, METH_NOARGS, "no-docs"},
                             { "tasmanian_array_test", &tasmanian_array_test, METH_VARARGS, "no-docs"},
                             { nullptr, nullptr, METH_NOARGS, nullptr }};

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_libcpplib(void){
    // the extra nullptrs are for multi-phase init/deinit
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,  // so that Python knows what struc is this (all python objects are effectively nullptrs)
        "libcpplib",            // the name of the module
        "This is cpp code",     // documentation
        -1,                     // module has global state and does not support sub-interpreters
        //&methoddef,             // the array of methods, single method example
        methoddefs,             // this registers two methods
        nullptr, nullptr, nullptr, nullptr };

    return PyModule_Create(&moduledef);

    // deprecated style:
    // Py_InitModule3(name, methods, doc);
}
#else
void initlibcpplib(void){
    // 3 stands for the number of input parameters
    Py_InitModule3("libcpplib", methoddefs, "documentation");
}
#endif

} // closes extern "C"
