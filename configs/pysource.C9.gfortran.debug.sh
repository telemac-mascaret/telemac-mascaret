# This file is a template for a Linux environement file
# runnning "source pysource.template.sh" will position all
# the necessary environement variable for telemac
# To adapt to your installation replace word <word> by their local value
###
### TELEMAC settings -----------------------------------------------------------
###
# Path to telemac root dir
export HOMETEL=$HOME/opentelemac/git/trunk
# Adding python scripts to PATH
export PATH=$HOMETEL/scripts/python27:.:$PATH
# Configuration file
export SYSTELCFG=$HOMETEL/configs/systel.edf.cfg
# Name of the configuration to use
export USETELCFG=C9.gfortran.debug
# Path to this file
export SOURCEFILE=$HOMETEL/configs/pysource.${USETELCFG}.sh
### Python
# To force python to flush its output
export PYTHONUNBUFFERED='true'
### API
export PYTHONPATH=$HOMETEL/scripts/python27:$PYTHONPATH
export LD_LIBRARY_PATH=$HOMETEL/builds/$USETELCFG/wrap_api/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$HOMETEL/builds/$USETELCFG/wrap_api/lib:$PYTHONPATH
###
### COMPILERS -----------------------------------------------------------
###
# Here are a few exemple for external libraries
SYSTEL=$HOME/opt

### MPI -----------------------------------------------------------
###
### EXTERNAL LIBRARIES -----------------------------------------------------------
###
### HDF5 -----------------------------------------------------------
export HDF5HOME=$SYSTEL/hdf5-1.8.14/arch/C9
export LD_LIBRARY_PATH=$HDF5HOME/lib:$LD_LIBRARY_PATH
### MED  -----------------------------------------------------------
export MEDHOME=$SYSTEL/med-3.2.0/arch/C9
export LD_LIBRARY_PATH=$MEDHOME/lib:$LD_LIBRARY_PATH
export PATH=$MEDHOME/bin:$PATH
### MUMPS -------------------------------------------------------------
export MUMPSHOME=$SYSTEL/mumps/gnu
export LD_LIBRARY_PATH=$MUMPSHOME/lib:$LD_LIBRARY_PATH
export SCALAPACKHOME=$SYSTEL/scalapack/gnu
export LD_LIBRARY_PATH=$SCALAPACKHOME/lib:$LD_LIBRARY_PATH
### METIS -------------------------------------------------------------
export METISHOME=$SYSTEL/metis-5.1.0/arch/C9
export LD_LIBRARY_PATH=$METISHOME/lib:$LD_LIBRARY_PATH
