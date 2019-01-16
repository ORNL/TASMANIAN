########################################################################
# Pre-configured for GNU Make.
# The information here will be omnipresent throughout the testing
# environment, which will reduce the need for configuring.
########################################################################

bEnableSyncTests = True

sLibPath = "./libtasmaniansparsegrid.so"

iGPUID = -1 # not used by the GNU Make

bHasBlas = False
bHasCuBlas = False
bHasCuda = False

bUsingMSVC = False

sGaussPattersonTableFile = "GaussPattersonRule.table"
