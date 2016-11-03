#CXX ?= mpiFCCpx

########at52
CC = gcc
CXX= g++ 
# FC="mpifrtpx" 
CFLAGS= -O -Wall
CXXFLAGS= -fopenmp -O2 -Wall 
########k-computer
#CXX ?= mpiFCCpx
# CC = mpifccpx
# CXX= mpiFCCpx 
# FC="mpifrtpx" 
#CPPFLAGS = -I/opt/klocal/include
#CFLAGS=-Xg -O2 -KPIC -O -Wall
#FCFLAGS=-X9 -O2 -KPIC -O -fopenmp
#CXXFLAGS= -Xg -O2 -KPIC -O -Wall

LDFLAGS = -fopenmp
# LDFLAGS = -Wl,-rpath,/opt/klocal/lib -Wl,-rpath,/opt/klocal/lib \
# -L/opt/klocal/lib  -lpetsc -SCALAPACK -SSL2 -lpthread -ltrtmetcpp \
# -Wl,-rpath,/opt/FJSVpxtof/sparc64fx/lib64 -L/opt/FJSVpxtof/sparc64fx/lib64 \
# -Wl,-rpath,/opt/FJSVtclang/GM-1.2.0-15/lib64 -L/opt/FJSVtclang/GM-1.2.0-15/lib64 \
# -Wl,-rpath,//opt/FJSVtclang/GM-1.2.0-15/lib64 -L//opt/FJSVtclang/GM-1.2.0-15/lib64 \
# -Wl,-rpath,/opt/FJSVxosmmm/lib64 -L/opt/FJSVxosmmm/lib64 \
# -lmpi_cxx -lfjdemgl -lstd_mt -lpthread -lstdc++ -ltrtmetcpp -lmpi_cxx -lfjdemgl -lstd_mt \
# -lpthread -lstdc++ -ltrtmet -ltrtmet_c -ldl -lmpi -ltofucom -ltofutop -lnsl -lutil -ltrtfdb \
# -lfjrtcl -ltrtth -lmpg -lmpgpthread -lpapi -lrt -lelf -lgcc_s -ldl

# LDFLAGS = -Wl,-rpath,/opt/klocal/lib -Wl,-rpath,/opt/klocal/lib \
# -L/opt/klocal/lib  -lpetsc -SCALAPACK -SSL2 -lpthread -ltrtmetcpp \
# -Wl,-rpath,/opt/FJSVpxtof/sparc64fx/lib64 -L/opt/FJSVpxtof/sparc64fx/lib64 \
# -Wl,-rpath,/opt/FJSVxosmmm/lib64 -L/opt/FJSVxosmmm/lib64 \
# -lmpi_cxx -lfjdemgl -lstd_mt -lpthread -lstdc++ -ltrtmetcpp -lmpi_cxx -lfjdemgl -lstd_mt \
# -lpthread -lstdc++ -ltrtmet -ltrtmet_c -ldl -lmpi -ltofucom -ltofutop -lnsl -lutil -ltrtfdb \
# -lfjrtcl -ltrtth -lmpg -lmpgpthread -lpapi -lrt -lelf -lgcc_s -ldl

program_NAME := knl_omp_ell 
program_C_SRCS := $(wildcard ./*.c)
program_CXX_SRCS := $(wildcard ./*.cpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
#program_INCLUDE_DIRS := /usr/local/PETSc/3.4.3/include /usr/local/PETSc/3.4.3/include /opt/FJSVfxlang/GM-1.2.1-08/include/mpi/fujitsu
#program_LIBRARY_DIRS :=
#program_LIBRARIES :=



#CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir)) -std=gnu99 
#LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
#LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))
#LDFLAGS = -Wl,-rpath,/usr/local/PETSc/3.4.3/lib -Wl,-rpath,/usr/local/PETSc/3.4.3/lib \
#-L/usr/local/PETSc/3.4.3/lib -lpetsc -SSL2 -lpthread -ltrtmetcpp \
# -Wl,-rpath,/opt/FJSVpxtof/sparc64fx/lib64 -L/opt/FJSVpxtof/sparc64fx/lib64 \
# -Wl,-rpath,/opt/FJSVfxlang/GM-1.2.1-08/lib64 -L/opt/FJSVfxlang/GM-1.2.1-08/lib64 \
# -Wl,-rpath,//opt/FJSVfxlang/GM-1.2.1-08/lib64 -L//opt/FJSVfxlang/GM-1.2.1-08/lib64 \
# -Wl,-rpath,/opt/FJSVxosmmm/lib64 -L/opt/FJSVxosmmm/lib64 -lmpi_cxx -lfjdemgl -lstd_mt \
# -lpthread -lstdc++ -ltrtmetcpp -lmpi_cxx -lfjdemgl -lstd_mt -lpthread -lstdc++ -ltrtmet \
# -ltrtmet_c -lmpi_f77 -lmpi_f90 -ldl -lmpi -ltofucom -ltofutop -lnsl -lutil -ltrtfdb -lfj90i \
#-lfj90fmt -lfj90f -lfjcrt -lfjrtcl -ltrtth -lmpg -lmpgpthread -lpapi -lrt -lelf -lgcc_s -ldl

.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(CXX) $(program_OBJS) -o $(program_NAME) $(LDFLAGS)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)

distclean: clean

