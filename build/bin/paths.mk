#Makefile include paths.mk

# RAMS root directory.

RAMS_ROOT=./../..

# Versions.

BRAMS_VERSION=5.0
JULES_VERSION=3.0
RAMS_VERSION=5.04
REVU_VERSION=2.3.1
UTILS_VERSION=2.0
HYPACT_VERSION=1.2.0
POST_VERSION=2.0

# Source directories.

# New paths like rams 6.x
BC=$(RAMS_ROOT)/src/brams/bc
CATT=$(RAMS_ROOT)/src/brams/catt
MODEL=$(RAMS_ROOT)/src/brams/model

#--(DMK-CCATT-INI)---------------------------------------
ADVC=$(RAMS_ROOT)/src/brams/advect
CCATT=$(RAMS_ROOT)/src/brams/ccatt
MODEL_CHEM=$(CCATT)/$(CHEM)
ISAN_CHEM=$(RAMS_ROOT)/src/brams/isan_chem
TUV=$(CCATT)/TUV
#--(DMK-CCATT-FIM)---------------------------------------

CUPARM=$(RAMS_ROOT)/src/brams/cuparm
FDDA=$(RAMS_ROOT)/src/brams/fdda
INIT=$(RAMS_ROOT)/src/brams/init
IO=$(RAMS_ROOT)/src/brams/io
ISAN=$(RAMS_ROOT)/src/brams/isan
MEMORY=$(RAMS_ROOT)/src/brams/memory
MICRO=$(RAMS_ROOT)/src/brams/micro
MKSFC=$(RAMS_ROOT)/src/brams/mksfc
MPI=$(RAMS_ROOT)/src/brams/mpi
NESTING=$(RAMS_ROOT)/src/brams/nesting
RADIATE=$(RAMS_ROOT)/src/brams/radiate
SIB=$(RAMS_ROOT)/src/brams/sib
SOIL_MOISTURE=$(RAMS_ROOT)/src/brams/soil_moisture
SURFACE=$(RAMS_ROOT)/src/brams/surface
JULES=$(RAMS_ROOT)/src/jules
TEB_SPM=$(RAMS_ROOT)/src/brams/teb_spm
TURB=$(RAMS_ROOT)/src/brams/turb
STILT=$(RAMS_ROOT)/src/brams/stilt
MATRIX=$(RAMS_ROOT)/src/brams/matrix

UTILS_LIB=$(RAMS_ROOT)/src/utils/lib
EFF=$(RAMS_ROOT)/src/utils/eff
NCARGD=$(RAMS_ROOT)/src/utils/ncargd
UTILS_INCS=$(RAMS_ROOT)/src/utils/include
UTILS_MODS=$(RAMS_ROOT)/src/utils/lib/modules
UTILS_DUMP=$(RAMS_ROOT)/src/utils/plot

ISAN=$(RAMS_ROOT)/src/brams/isan
ISAN_MODS=$(RAMS_ROOT)/src/brams/isan

POST_SRC = $(RAMS_ROOT)/src/post
POST_INCS = $(UTILS_INCS)

