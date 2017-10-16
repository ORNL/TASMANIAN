#!/usr/bin/python

# necessary import for every use of TASMANIAN
#
import TasmanianSG
import numpy as np
import math

from random import uniform

import matplotlib.pyplot as plt
import matplotlib.colors as cols
from mpl_toolkits.mplot3d import Axes3D

import scipy.linalg as linalg

grid = TasmanianSG.TasmanianSparseGrid()
grid.makeGlobalGrid(2, 1, 10, 'level', 'leja')
grid.write("tst", True)

print "Written"
grid.read("tst", True)


exit(0)

#grid = TasmanianSG.TasmanianSparseGrid()
#print("GPU Name:   {0:1s}".format(grid.getGPUname(0)))
#print("Number of GPUs: {0:1d}\n".format(grid.getNumGPUs()))
#
#for iI in range(grid.getNumGPUs()):
#    print("GPU Number: {0:1d}".format(iI))
#    print("GPU Name:   {0:1s}".format(grid.getGPUname(iI)))
#    print("GPU Memory: {0:1d}MB\n".format(grid.getGPUmemory(iI)))


lS = [[ -1.0000000000000000e+00, -1.0000000000000000e+00 ],
[ -9.4999999999999996e-01, -8.9308697941380522e-01 ],
[ -8.9999999999999991e-01, -8.0694437527264384e-01 ],
[ -8.4999999999999987e-01, -7.3431296933571744e-01 ],
[ -7.9999999999999982e-01, -6.7078258498244259e-01 ],
[ -7.4999999999999978e-01, -6.1361965332418866e-01 ],
[ -6.9999999999999973e-01, -5.6108228141437999e-01 ],
[ -6.4999999999999969e-01, -5.1202201970728334e-01 ],
[ -5.9999999999999964e-01, -4.6565226257165915e-01 ],
[ -5.4999999999999960e-01, -4.2141266500440605e-01 ],
[ -4.9999999999999961e-01, -3.7888865673137279e-01 ],
[ -4.4999999999999962e-01, -3.3776262937565893e-01 ],
[ -3.9999999999999963e-01, -2.9778349233483314e-01 ],
[ -3.4999999999999964e-01, -2.5874705517542806e-01 ],
[ -2.9999999999999966e-01, -2.2048293657180265e-01 ],
[ -2.4999999999999967e-01, -1.8284551312103192e-01 ],
[ -1.9999999999999968e-01, -1.4570743696295380e-01 ],
[ -1.4999999999999969e-01, -1.0895482545767013e-01 ],
[ -9.9999999999999686e-02, -7.2483556382354597e-02 ],
[ -4.9999999999999684e-02, -3.6196295687066143e-02 ],
[ 3.1918911957973251e-16, 2.3097292789814214e-16 ],
[ 5.0000000000000322e-02, 3.6196295687066608e-02 ],
[ 1.0000000000000032e-01, 7.2483556382355069e-02 ],
[ 1.5000000000000033e-01, 1.0895482545767059e-01 ],
[ 2.0000000000000034e-01, 1.4570743696295430e-01 ],
[ 2.5000000000000033e-01, 1.8284551312103245e-01 ],
[ 3.0000000000000032e-01, 2.2048293657180312e-01 ],
[ 3.5000000000000031e-01, 2.5874705517542856e-01 ],
[ 4.0000000000000030e-01, 2.9778349233483364e-01 ],
[ 4.5000000000000029e-01, 3.3776262937565943e-01 ],
[ 5.0000000000000033e-01, 3.7888865673137334e-01 ],
[ 5.5000000000000038e-01, 4.2141266500440672e-01 ],
[ 6.0000000000000042e-01, 4.6565226257165998e-01 ],
[ 6.5000000000000047e-01, 5.1202201970728389e-01 ],
[ 7.0000000000000051e-01, 5.6108228141438088e-01 ],
[ 7.5000000000000056e-01, 6.1361965332418955e-01 ],
[ 8.0000000000000060e-01, 6.7078258498244359e-01 ],
[ 8.5000000000000064e-01, 7.3431296933571877e-01 ],
[ 9.0000000000000069e-01, 8.0694437527264495e-01 ],
[ 9.5000000000000073e-01, 8.9308697941380666e-01 ],
[ 1.0000000000000000e+00, 1.0000000000000000e+00 ]]
aS = np.array( lS )

plt.figure(1)
plt.plot( aS[:,0], aS[:,1], '-' )
#plt.axis([0.9,5.1,0.9,5.1])
plt.axis([-1.05,1.05,-1.05,1.05])

plt.show()


exit(0)

grid = TasmanianSG.TasmanianSparseGrid()

try:
        grid.makeGlobalGrid(2, 1, 10, 'level', 'leja')
except TasmanianSG.TasmanianInputError as t:
        print( t.sVariable )
        print( t.sMessage )

exit(0)
        
grid.makeGlobalGrid(2, 1, 10, 'level', 'leja')
grid.loadNeededPoints(np.array([[0,0] for i in range(66)]))


iDim = 2
iSamples = 200

#lfX = np.array([uniform(-1.0,1.0) for i in range(iDim)])
#lfX = np.array([[0.0 for i in range(iDim)] for j in range(2)])
#lfX[1][1] = -1.0

llfX = np.array([[uniform(-1.0,1.0) for i in range(iDim)] for j in range(iSamples)])

#grid.makeGlobalGrid( iDim, 1, 3, "level", "leja" )
#grid.makeLocalPolynomialGrid(iDim, 1, 6, 2, "localp")
#grid.makeWaveletGrid( iDim, 1, 3, 1 )
#grid.makeSequenceGrid(iDim, 1, 10, 'level', 'leja')
grid.makeGlobalGrid(iDim, 1, 10, 'level', 'leja')

#print grid.getNumPoints()

aVanMatrix = grid.evalBatchHierarchicalFunctions( llfX )

aRHS = np.array( [[ math.exp( -x[0] -x[1] ) for x in llfX ] for j in range(1)] )
aC = linalg.lstsq( aVanMatrix, aRHS.T )[0]

#print aC

#try:
#        grid.setHierarchicalCoefficients( aC )
#except t:
#        print t.sMessage

grid.loadNeededPoints( np.array( [ [0.0] for i in range(grid.getNumNeeded()) ] ) )
grid.setHierarchicalCoefficients( aC )

#aC2 = grid.getSurpuses()

#print np.column_stack( [aC, aC2 ] )

llfX = np.array([[uniform(-1.0,1.0) for i in range(iDim)] for j in range(iSamples)])
aRHS = np.array( [ math.exp( -x[0] -x[1] ) for x in llfX ] )
lRes = grid.evaluateBatch( llfX )

fMax = 0.0
for i in range(iSamples):
        fMax = max( math.fabs(aRHS[i] - lRes[i][0]), fMax )
        
print("Projection error = {0:1.4e}".format( fMax ))


#grid.makeLocalPolynomialGrid(iDim, 1, 4, 2, "localp")
#grid.makeWaveletGrid( iDim, 1, 3, 1 )
#grid.makeSequenceGrid(iDim, 1, 10, 'level', 'leja')
grid.makeGlobalGrid(iDim, 1, 10, 'level', 'leja')

llfX = grid.getNeededPoints()
aV = np.array( [ [math.exp( -x[0] -x[1] )] for x in llfX ] )
grid.loadNeededPoints( aV )

#print grid.getSurpuses()

llfX = np.array([[uniform(-1.0,1.0) for i in range(iDim)] for j in range(iSamples)])
aRHS = np.array( [ math.exp( -x[0] -x[1] ) for x in llfX ] )
lRes = grid.evaluateBatch( llfX )

fMax = 0.0
for i in range(iSamples):
        fMax = max( math.fabs(aRHS[i] - lRes[i][0]), fMax )

print("Interpolation error = {0:1.4e}".format( fMax ))

#print grid.evalBatchHierarchicalFunctions( lfX )

exit(0)


grid = TasmanianSG.TasmanianSparseGrid()

iPlot = 2

if ( iPlot == 1 ):
        # plot PWC 1-D basis functions
        N = 90

        x = np.linspace(-0.999, 0.999, N ).reshape(N,1)

        grid.makeLocalPolynomialGrid( 1, 1, 0, 0 )
        grid.loadNeededPoints( np.array( [[ 1.0, ],] ) )

        f = grid.evaluateBatch( x )
        #print f

        #plt.figure( 1 )

        dum, ax = plt.subplots( 3, figsize=(10, 12) )

        ax[0].plot( x, f, 'b-', linewidth=3.0, markersize = 9.0 )
        #ax[0].plot( x, np.zeros( x.shape ), 'k-', linewidth=3.0, markersize = 9.0 )
        ax[0].plot( [-1.1,1.1], [0,0], 'k-', linewidth=3.0, markersize = 9.0 )
        ax[0].axis([-1.1,1.1,-0.1,1.1])


        grid.makeLocalPolynomialGrid( 1, 1, 1, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0, ], [1.0, ], [0.0, ]] ) )
        f = grid.evaluateBatch( x )
        ax[1].plot( x[0:N/3], f[0:N/3], 'r-', linewidth=3.0, markersize = 9.0 )

        grid.makeLocalPolynomialGrid( 1, 1, 1, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0, ], [0.0, ], [1.0, ]] ) )
        f = grid.evaluateBatch( x )
        ax[1].plot( x[2*N/3:N], f[2*N/3:N], 'g-', linewidth=3.0, markersize = 9.0 )

        #ax[1].plot( x, np.zeros( x.shape ), 'k-', linewidth=3.0, markersize = 9.0 )
        ax[1].plot( [-1.1,1.1], [0,0], 'k-', linewidth=3.0, markersize = 9.0 )
        ax[1].axis([-1.1,1.1,-0.1,1.1])


        grid.makeLocalPolynomialGrid( 1, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,], [0.0,], [0.0,],[ 1.0,], [0.0,], [0.0,],[ 0.0,], [0.0,], [0.0,]] ) )
        f = grid.evaluateBatch( x )
        ax[2].plot( x[0:N/9], f[0:N/9], 'b-', linewidth=3.0, markersize = 9.0 )

        grid.makeLocalPolynomialGrid( 1, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,], [0.0,], [0.0,],[ 0.0,], [1.0,], [0.0,],[ 0.0,], [0.0,], [0.0,]] ) )
        f = grid.evaluateBatch( x )
        ax[2].plot( x[2*N/9:N/3], f[2*N/9:N/3], 'r-', linewidth=3.0, markersize = 9.0 )

        grid.makeLocalPolynomialGrid( 1, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,], [0.0,], [0.0,],[ 0.0,], [0.0,], [1.0,],[ 0.0,], [0.0,], [0.0,]] ) )
        f = grid.evaluateBatch( x )
        ax[2].plot( x[N/3:4*N/9], f[N/3:4*N/9], 'g-', linewidth=3.0, markersize = 9.0 )

        grid.makeLocalPolynomialGrid( 1, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,], [0.0,], [0.0,],[ 0.0,], [0.0,], [0.0,],[ 1.0,], [0.0,], [0.0,]] ) )
        f = grid.evaluateBatch( x )
        ax[2].plot( x[5*N/9:6*N/9], f[5*N/9:6*N/9], 'c-', linewidth=3.0, markersize = 9.0 )

        grid.makeLocalPolynomialGrid( 1, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,], [0.0,], [0.0,],[ 0.0,], [0.0,], [0.0,],[ 0.0,], [1.0,], [0.0,]] ) )
        f = grid.evaluateBatch( x )
        ax[2].plot( x[6*N/9:7*N/9], f[6*N/9:7*N/9], 'm-', linewidth=3.0, markersize = 9.0 )

        grid.makeLocalPolynomialGrid( 1, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,], [0.0,], [0.0,],[ 0.0,], [0.0,], [0.0,],[ 0.0,], [0.0,], [1.0,]] ) )
        f = grid.evaluateBatch( x )
        ax[2].plot( x[8*N/9:N], f[8*N/9:N], 'y-', linewidth=3.0, markersize = 9.0 )

        #ax[2].plot( x, np.zeros( x.shape ), 'k-', linewidth=3.0, markersize = 9.0 )
        ax[2].plot( [-1.1,1.1], [0,0], 'k-', linewidth=3.0, markersize = 9.0 )
        ax[2].axis([-1.1,1.1,-0.1,1.1])

        for tick in ax[0].xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
        for tick in ax[1].xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
        for tick in ax[2].xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
        for tick in ax[0].yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
        for tick in ax[1].yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
        for tick in ax[2].yaxis.get_major_ticks():
                tick.label.set_fontsize(14)

        plt.show()

if ( iPlot == 2 ):
        iN = 90
        
        x = np.linspace(-0.999, 0.999, iN )
        y = np.linspace(-0.999, 0.999, iN )

        XX, YY = np.meshgrid( x, y )
        
        #dum, ax = plt.subplots( 3, 3, figsize=(12, 12), projection='3d' )
        fig = plt.figure( 1, figsize=(12, 12) )
        ax = fig.add_subplot(3, 3, 1, projection='3d')

        grid.makeLocalPolynomialGrid( 2, 1, 0, 0 )
        grid.loadNeededPoints( np.array( [[ 1.0, ],] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(0.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX, YY, ZZ, facecolors=fcolors )
        ax.set_xlim3d(-1.1,1.1)
        ax.set_ylim3d(-1.1,1.1)
        ax.set_zlim3d(-0.1,1.1)
        
        ax.set_xticklabels( () )
        ax.set_yticklabels( () )
        ax.set_zticklabels( () )
        
        
        ax = fig.add_subplot(3, 3, 2, projection='3d')
        
        grid.makeLocalPolynomialGrid( 2, 1, 1, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(0.8*np.ones( ZZ.shape ))
        ax.plot_surface( XX[0:iN/3,:], YY[0:iN/3,:], ZZ[0:iN/3,:], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 1, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(-0.8*np.ones( ZZ.shape ))
        ax.plot_surface( XX[2*iN/3:iN,:], YY[2*iN/3:iN,:], ZZ[2*iN/3:iN,:], facecolors=fcolors )
        
        ax.set_xlim3d(-1.1,1.1)
        ax.set_ylim3d(-1.1,1.1)
        ax.set_zlim3d(-0.1,1.1)
        ax.set_xticklabels( () )
        ax.set_yticklabels( () )
        ax.set_zticklabels( () )
        
        ax = fig.add_subplot(3, 3, 4, projection='3d')
        
        grid.makeLocalPolynomialGrid( 2, 1, 1, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(0.8*np.ones( ZZ.shape ))
        ax.plot_surface( XX[:,0:iN/3], YY[:,0:iN/3], ZZ[:,0:iN/3], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 1, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-0.8, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(-0.5*np.ones( ZZ.shape ))
        ax.plot_surface( XX[:,2*iN/3:iN], YY[:,2*iN/3:iN], ZZ[:,2*iN/3:iN], facecolors=fcolors )
        
        ax.set_xlim3d(-1.1,1.1)
        ax.set_ylim3d(-1.1,1.1)
        ax.set_zlim3d(-0.1,1.1)
        ax.set_xticklabels( () )
        ax.set_yticklabels( () )
        ax.set_zticklabels( () )
        
        
        ax = fig.add_subplot(3, 3, 3, projection='3d')
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='rainbow')
        m.set_array([])
        fcolors = m.to_rgba(-1.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[0:iN/9,:], YY[0:iN/9,:], ZZ[0:iN/9,:], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='gnuplot2')
        m.set_array([])
        fcolors = m.to_rgba(0.9*np.ones( ZZ.shape ))
        ax.plot_surface( XX[2*iN/9:3*iN/9,:], YY[2*iN/9:3*iN/9,:], ZZ[2*iN/9:3*iN/9,:], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(0.7*np.ones( ZZ.shape ))
        ax.plot_surface( XX[3*iN/9:4*iN/9,:], YY[3*iN/9:4*iN/9,:], ZZ[3*iN/9:4*iN/9,:], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(0.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[5*iN/9:6*iN/9,:], YY[5*iN/9:6*iN/9,:], ZZ[5*iN/9:6*iN/9,:], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='gnuplot2')
        m.set_array([])
        fcolors = m.to_rgba(-1.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[6*iN/9:7*iN/9,:], YY[6*iN/9:7*iN/9,:], ZZ[6*iN/9:7*iN/9,:], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='gist_rainbow')
        m.set_array([])
        fcolors = m.to_rgba(1.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[8*iN/9:9*iN/9,:], YY[8*iN/9:9*iN/9,:], ZZ[8*iN/9:9*iN/9,:], facecolors=fcolors )
        
        ax.set_xlim3d(-1.1,1.1)
        ax.set_ylim3d(-1.1,1.1)
        ax.set_zlim3d(-0.1,1.1)
        ax.set_xticklabels( () )
        ax.set_yticklabels( () )
        ax.set_zticklabels( () )
        
        
        
        ax = fig.add_subplot(3, 3, 7, projection='3d')
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='rainbow')
        m.set_array([])
        fcolors = m.to_rgba(-1.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[:,0:iN/9], YY[:,0:iN/9], ZZ[:,0:iN/9], facecolors=fcolors )
        #ax.plot_surface( XX, YY, ZZ, facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='gnuplot2')
        m.set_array([])
        fcolors = m.to_rgba(0.9*np.ones( ZZ.shape ))
        ax.plot_surface( XX[:,2*iN/9:3*iN/9], YY[:,2*iN/9:3*iN/9], ZZ[:,2*iN/9:3*iN/9], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(0.7*np.ones( ZZ.shape ))
        ax.plot_surface( XX[:,3*iN/9:4*iN/9], YY[:,3*iN/9:4*iN/9], ZZ[:,3*iN/9:4*iN/9], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
        m.set_array([])
        fcolors = m.to_rgba(-0.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[:,5*iN/9:6*iN/9], YY[:,5*iN/9:6*iN/9], ZZ[:,5*iN/9:6*iN/9], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='gnuplot2')
        m.set_array([])
        fcolors = m.to_rgba(-1.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[:,6*iN/9:7*iN/9], YY[:,6*iN/9:7*iN/9], ZZ[:,6*iN/9:7*iN/9], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[1.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='gist_rainbow')
        m.set_array([])
        fcolors = m.to_rgba(1.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[:,8*iN/9:9*iN/9], YY[:,8*iN/9:9*iN/9], ZZ[:,8*iN/9:9*iN/9], facecolors=fcolors )
        
        ax.set_xlim3d(-1.1,1.1)
        ax.set_ylim3d(-1.1,1.1)
        ax.set_zlim3d(-0.1,1.1)
        ax.set_xticklabels( () )
        ax.set_yticklabels( () )
        ax.set_zticklabels( () )
        
        
        ax = fig.add_subplot(3, 3, 5, projection='3d')
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='rainbow')
        m.set_array([])
        fcolors = m.to_rgba(1.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[0:iN/3,0:iN/3], YY[0:iN/3,0:iN/3], ZZ[0:iN/3,0:iN/3], facecolors=fcolors )
        #ax.plot_surface( XX, YY, ZZ, facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='rainbow')
        m.set_array([])
        fcolors = m.to_rgba(-1.0*np.ones( ZZ.shape ))
        ax.plot_surface( XX[2*iN/3:iN,0:iN/3], YY[2*iN/3:iN,0:iN/3], ZZ[2*iN/3:iN,0:iN/3], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='gnuplot2')
        m.set_array([])
        fcolors = m.to_rgba(-0.9*np.ones( ZZ.shape ))
        ax.plot_surface( XX[0:iN/3,2*iN/3:iN], YY[0:iN/3,2*iN/3:iN], ZZ[0:iN/3,2*iN/3:iN], facecolors=fcolors )
        
        grid.makeLocalPolynomialGrid( 2, 1, 2, 0 )
        grid.loadNeededPoints( np.array( [[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 1.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[ 0.0,],[0.0,]] ) )
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )

        norm = cols.Normalize(-1.0, 1.0)
        m = plt.cm.ScalarMappable(norm=norm, cmap='gist_rainbow')
        m.set_array([])
        fcolors = m.to_rgba(-0.1*np.ones( ZZ.shape ))
        ax.plot_surface( XX[2*iN/3:iN,2*iN/3:iN], YY[2*iN/3:iN,2*iN/3:iN], ZZ[2*iN/3:iN,2*iN/3:iN], facecolors=fcolors )
        
        ax.set_xlim3d(-1.1,1.1)
        ax.set_ylim3d(-1.1,1.1)
        ax.set_zlim3d(-0.9,1.1)
        ax.set_xticklabels( () )
        ax.set_yticklabels( () )
        ax.set_zticklabels( () )
        
        fig.subplots_adjust(hspace=0.01, wspace=0.01)
        
        plt.show()








exit(0)

#print TasmanianSG.bTsgPlotting

grid.makeLocalPolynomialGrid( 2, 1, 0, 0 )

#print grid.getNumPoints()

grid.loadNeededPoints( np.array( [[ 1.0, ],] ) )

iN = 100

x = np.linspace(-1.0, 1.0, iN )
y = np.linspace(-1.0, 1.0, iN )

XX, YY = np.meshgrid( x, y )

#print XX, YY

#print np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T.shape

try:
        ZZ = grid.evaluateBatch( np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T ).reshape( (iN,iN) )
except TasmanianSG.TasmanianInputError as e:
        #print np.vstack( [ XX.reshape( [iN*iN,] ), YY.reshape( [iN*iN,] ) ] ).T.shape
        #print e.sMessage
        pass
        
fig = plt.figure( 1 )
ax = fig.gca(projection='3d')
        
#Axes3D.plot_surface( XX, YY, ZZ )
ax.plot_surface( XX, YY, ZZ )
plt.show()


#grid.makeGlobalGrid( 2, 1, 3, "level", "clenshaw-curtis" )

#aTestP = np.array( [ [uniform(-1.0,1.0), uniform(-1.0,1.0)] for i in range(1000) ] )
#aTestVals = np.empty( (aTestP.shape[0],) )
#for iI in range(aTestP.shape[0]):
        #aTestVals[iI] = math.exp(- aTestP[iI,0] ** 2 - aTestP[iI,1] ** 2)

#lResult = [ [] for i in range(5) ]

#for iOrder in range(5):
        #iIter = 11
        #if ( iOrder == 0 ): iIter = 8
        #for iI in range( iIter ):
                ##iE = iOrder +1
                ##if ( iOrder == 4 ):iE = -1
                #iE = iOrder
                
                #grid.makeLocalPolynomialGrid( 2, 1, iI, iE, "localp" )

                #aPoints = grid.getNeededPoints()

                #aVals = np.empty([aPoints.shape[0], 1])

                #for iI in range(aPoints.shape[0]):
                        #aVals[iI] = math.exp(- aPoints[iI,0] ** 2 - aPoints[iI,1] ** 2)
                #grid.loadNeededPoints(aVals)
                
                #aRes = grid.evaluateBatch( aTestP )
                
                #fError = max( np.fabs( aTestVals - aRes[:,0] ) )
                #iNumP = grid.getNumPoints()
                #lResult[iOrder].append( [iNumP, fError, -math.log(fError) / math.log(iNumP+1)] )
                #print iOrder+1, iNumP, fError, -math.log(fError) / math.log((iNumP+1) ** 0.5)
                
                

#plt.figure( 1 )
#iOrder = 0
#p1, = plt.semilogy( [ l[0] for l in lResult[iOrder] ], [ l[1] for l in lResult[iOrder] ], label = "Linear", linewidth=3.0 )

#iOrder = 1
#p2, = plt.semilogy( [ l[0] for l in lResult[iOrder] ], [ l[1] for l in lResult[iOrder] ], label = "Quadratic", linewidth=3.0 )

#iOrder = 2
#p3, = plt.semilogy( [ l[0] for l in lResult[iOrder] ], [ l[1] for l in lResult[iOrder] ], label = "Cubic", linewidth=3.0 )

#iOrder = 3
#p4, = plt.semilogy( [ l[0] for l in lResult[iOrder] ], [ l[1] for l in lResult[iOrder] ], label = "Quartic", linewidth=3.0 )

#iOrder = 4
#p5, = plt.semilogy( [ l[0] for l in lResult[iOrder] ], [ l[1] for l in lResult[iOrder] ], label = "Maximal", linewidth=3.0 )

#plt.legend( handles = [ p1, p2, p3, p4, p5 ], loc = 1, fontsize=18 )

#plt.axis([1,8000,1.E-10,1.E+0])

#plt.xlabel("Number of points", fontsize=18)
#plt.ylabel("Error", fontsize=18)

#plt.show()

#print lResult

#exit(0)

#grid.removePointsBySurplus( 1.E-7, -1 )

#pP = grid.plotPoints2D( sStyle="ro", iNumFigure=1, iMarkerSize=3, bShow=False )
#pP.xlabel( "first dimension", fontsize=20 )
#pP.show()

#grid.removePointsBySurplus( 1.E-3, -1 )
#print "Kept: {0:1d}".format(grid.getNumPoints())
#grid.removePointsBySurplus( 15.E-3, -1 )

#pP = grid.plotPoints2D( sStyle="ro", iNumFigure=1, iMarkerSize=3, bShow=True )

#grid.plotResponse2D()

