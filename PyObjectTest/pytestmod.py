# import the module, see inside for names and such
import numpy
import libcpplib
import ctypes

def callSample(f):
    # call a dummy "hello world" method
    libcpplib.tasmanian_dummy_test()

    #print("got f: ", f) # helps debug pointers

    # pass the callable object and an integer
    libcpplib.tasmanian_cpp_method(f, 4)


current_model = lambda x : 2 * x

def model_proxy(lX):
    print("calling proxy", lX)
    x = numpy.array(lX)
    y = current_model(x)
    print("result", tuple([y[i] for i in range(len(y))]))
    if (len(y) == 3):
        print("ERROR: size is 3")
        raise Exception("Test")
    return tuple([y[i] for i in range(len(y))])

def callNumpy(f):

    pLibTSG = ctypes.cdll.LoadLibrary("./libcpplib.so")
    pLibTSG.makevec.restype = ctypes.c_void_p

    p = pLibTSG.makevec()
    print(p)

    # make a call using numpy-arrays

    dims = 3
    x = numpy.zeros(dims)

    global current_model
    current_model = f;

    libcpplib.tasmanian_array_test(5, p)
    try:
        libcpplib.tasmanian_array_test(4, p)
    except:
        print("Caught the error")

