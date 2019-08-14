import numpy
import pytestmod

print("Hello Cruel World")

def f(x):
    print("python f(x) got an x = ", x)
    return 2 * x

# both work:
#pytestmod.callSample(lambda x : 2 * x)
pytestmod.callSample(f)


def g(x):
    return numpy.sin(x)

print(g)

pytestmod.callNumpy(g)
