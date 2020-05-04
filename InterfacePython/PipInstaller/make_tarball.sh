
this_folder=`pwd`

if [ ! -d "$this_folder/SparseGrids/" ] || [ ! -d "$this_folder/InterfacePython/" ]; then
    echo "ERROR: must run this script from the Tasmanian source root folder"
    exit 1
fi

cp ./InterfacePython/PipInstaller/setup.py .
cp ./InterfacePython/PipInstaller/MANIFEST.in .
cp ./InterfacePython/PipInstaller/pyproject.toml .

python3 setup.py sdist

rm MANIFEST.in
rm MANIFEST
rm pyproject.toml
rm setup.py

echo "tarball build in $this_folder/dist/"
