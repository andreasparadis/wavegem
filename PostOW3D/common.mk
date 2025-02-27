# For inclusion in makefile for use with GNU make (gmake)
# 
# Purpose: Modify to local configuration by user.
#

# Program name
PROGNAME = OceanWave3D
LIBNAME  = libOceanWave3D_botp.so

# Installation directory
INSTALLDIR = $(HOME)/OceanWave3D
LIBINSTALLDIR = $(HOME)/lib

# Build directory where object files are stored 
BUILDDIR = $(PWD)/build

# The build environment is set either by the choice of a compiler 
# flag, or by creating a block for a specific $USER.  
# Choose the Fortran compiler on this system
# E.g. pathf90, f90, gfortran, gf90, ifort
FC = gfortran
#FC = gfortran44
#FC = gfortran-4.4
#FC = gf90

USER =

# First the blocks based on compiler name:  

ifeq ($(FC),gfortran44)
  # gbar linux machines
  LIBDIRS  = -L $(HOME)/lib/
  LINLIB   = -lharwell -lskit -llapack -lblas
  DBFLAGS  = -pg -g -O0 -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
  OPTFLAGS = -O3 -ffpe-trap=invalid,zero,overflow -ffree-line-length-none
endif

ifeq ($(FC),ifort)
  # hbb work machine with intel compiler
  LIBDIRS  = -L/usr/local/lib/ 
  LINLIB   = -lharwell_intel -lskit_intel -llapack_intel -lblas_intel
  DBFLAGS  = -g -CB -fpe0 -fpstkchk  -traceback 
  OPTFLAGS = -O3 -tpp7
#  OPTFLAGS = -g -fast
endif


ifeq ($(FC),gfortran)
  # MacOS, apek
  LIBDIRS  = -L /lib -L $(HOME)/lib/ -L $(HOME)/local/lib/  -L $(HOME)/local/include/ -L $(HOME)/Harwell/ -L $(HOME)/ThirdParty-HSL/ -L $(HOME)/ThirdParty-ASL/ -L $(HOME)/SPARSKIT2/ -L $(HOME)/VTKFortran/ -L$(HOME)/Desktop/OceanWave3D/hsl-galahad/
# -L/Users/apek/Documents/Fortran/Harwell/lib -L/usr/local/atlas/lib -L/Users/apek/Documents/Fortran/SPARSKIT2
 LINLIB   = -lharwell_gfortran -lskit  -llapack -Lib_VTK_IO -lblas
#  LINLIB   = -lharwell -lskit -latlas -l_VTK_IO -framework veclib
#   LINLIB   = -lharwell -lskit -latlas -framework veclib
  DBFLAGS  = -pg -g -O -fcheck=all -ffpe-trap=invalid,zero,overflow
  OPTFLAGS = -O2
endif

ifeq ($(FC),pathf90)
  # Niflheim cluster, DTU, apek
  LIBDIRS  = -L/opt/acml3.5.0/pathscale64/lib -L$(HOME)/lib
  LINLIB   = -lharwellFPATH90 -lskitPATHF90 -llinpackPATHF90 -lacml
  DBFLAGS  = -pg -g -O -static -woffoptions
  OPTFLAGS = -O2 -static -woffoptions -ffortran-bounds-check
endif

ifeq ($(FC),f90)
  # gbar, DTU, apek
  LIBDIRS  = -L$(HOME)/lib/ 
  LINLIB   = -lharwell -lskit -xlic_lib=sunperf
  DBFLAGS  = -pg -g -O0 # -static -woffoptions
  OPTFLAGS = -O -fast 
#  OPTFLAGS = -g -fast
endif

# Then the blocks for specific users (this clobbers the above info.)
ifeq ($(USER),hbb)
  # hbb machine, gfortran
  FC=gfortran	
  LIBDIRS  = -L $(HOME)/lib/ 
  LINLIB   = -lharwell -lskit -llapack -lblas
  DBFLAGS  = -pg -g -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none  -fno-automatic
  OPTFLAGS = -pg -O3 -ffree-line-length-none -fno-automatic -ffpe-trap=invalid,zero,overflow
#  OPTFLAGS = -g -fast
endif

ifeq ($(USER),botp-dev)
  # botp kubuntu, 10.04-64bit
  FC       = gfortran
  LIBDIRS  = -L$(HOME)/lib/ 
  LINLIB   = -lharwell_gfortran -ltmglib_gfortran -llapack_gfortran -lblas -lskit_gfortran
  DBFLAGS  = -pg -g -O0 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
  OPTFLAGS = -O3 -fPIC -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
  SHLIBFLAGS  = -shared -O2 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
endif


ifeq ($(USER),olli)
  # olli linux machine
  FC       = gfortran
  LIBDIRS  = -L$(HOME)/lib/ -Ldep/SPARSKIT2/ -Ldep/Harwell/
  LINLIB   = -lharwell -lskit -llapack -lblas 
  DBFLAGS  = -pg -g -O0 -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
  OPTFLAGS = -O3 -ffpe-trap=invalid,zero,overflow -ffree-line-length-none
endif


ifeq ($(USER),gudu)
  # gudu linux machine
  FC       = ifort
  LIBDIRS  = -L$(HOME)/lib/
  LINLIB   = -lharwell -lskit -llapack -lblas 
  DBFLAGS  = -O0 -traceback -check all -warn all 
  OPTFLAGS = -O3 -xHOST -ipo -ip
endif

ifeq ($(USER),apek)
  # MacOS, apek
  FC       = gfortran
  LIBDIRS  = -L $(HOME)/lib/
# -L/Users/apek/Documents/Fortran/Harwell/lib -L/usr/local/atlas/lib -L/Users/apek/Documents/Fortran/SPARSKIT2                                                                                   
#  LINLIB   = -lharwell -lskit -latlas -llapack -l_VTK_IO                                                                                                                                        
#  LINLIB   = -lharwell -lskit -latlas -l_VTK_IO -framework veclib                                                                                                                               
  LINLIB   = -lharwell -lskit -latlas -framework veclib
  DBFLAGS  = -pg -g -O -fcheck=all -ffpe-trap=invalid,zero,overflow
  OPTFLAGS = -O2
endif
