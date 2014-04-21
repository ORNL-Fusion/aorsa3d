# Name:     makefile
# Function: Provides the make utility with instructions for
#           compiling the aorsa3d code into an executable module.
# Date:     10/25/2002
# Revised:

# Macro definitions

HOME =/dfs/home/$(LOGNAME)

CODENAME=xaorsa3d

CODEDIR = AORSA3D/04-15-03_AORSA3D_NEWLAB3_SCALAPACK

SOURCE_DIR = $(HOME)/$(CODEDIR)/src

INCLUDE_DIR = $(HOME)/$(CODEDIR)/src

OBJECT_DIR = $(HOME)/$(CODEDIR)/obj

OBJFILES = \
 $(OBJECT_DIR)/kron3_mod.o \
 $(OBJECT_DIR)/prof_mod.o \
 $(OBJECT_DIR)/aorsa3dMain.o \
 $(OBJECT_DIR)/aorsaSubs.o \
 $(OBJECT_DIR)/sigma.o \
 $(OBJECT_DIR)/zfunction.o \
 $(OBJECT_DIR)/ztable.o \
 $(OBJECT_DIR)/current.o \
 $(OBJECT_DIR)/wdot.o \
 $(OBJECT_DIR)/fourier.o \
 $(OBJECT_DIR)/assert.o \
 $(OBJECT_DIR)/setupblacs.o \
 $(OBJECT_DIR)/check.o

INCFILES = $(SOURCE_DIR)/kron3_mod.f
INCFILES = $(SOURCE_DIR)/prof_mod.f

F77 = mpxlf90_r -qfixed \
      -O3 \
      -qinitauto=00 \
      -qnolm -qarch=auto -qfixed -qnosmp\
      -qmaxmem=-1 -qmoddir=/tmp/jaegeref/mod -I/tmp/jaegeref/mod \
      -qsave -qrealsize=8 -q32
      
F90 = mpxlf90_r -qfixed \
      -O3 -qstrict \
      -qinitauto=00 \
      -qnolm -qarch=auto -qfree=f90 -qnosmp\
      -qmaxmem=-1 \
      -qsave -qrealsize=8 -q32
      
      

OPTIMIZATIONg =  -g -C -qextchk
OPTIMIZATIONO =  -g -O -qmaxmem=-1 -qextchk
OPTIMIZATION = $(OPTIMIZATIONO)  


F77FLAGS= -c $(OPTIMIZATION)  -I $(INCLUDE_DIR)
F77FLAGSg= -c $(OPTIMIZATIONg) -I $(INCLUDE_DIR)
F77FLAGS4= -c $(OPTIMIZATION4) -I $(INCLUDE_DIR)

COMPILE = $(F77) $(F77FLAGS)
COMPILE90 = $(F90) $(F77FLAGS)
COMPILEg = $(F77) $(F77FLAGSg)
COMPILE4 = $(F77) $(F77FLAGS4)



LOADFLAGS = -q32 -bmaxdata:2000000000 -bmaxstack:256000000 \
        -L/tmp/gpfs200a/efdazedo/LIB \
        -L/usr/local/lib/ -lscalapack \
        -lblacsF77init -lblacs \
        -lpessl -lessl


LOAD = $(F77) $(OPTIMIZATION) $(LOADFLAGS)

# Compile the program

xaorsa3d:          $(OBJFILES)
	          $(LOAD) -o ./$(CODENAME)  $(OBJFILES)

# Dependencies

$(OBJECT_DIR)/kron3_mod.o:      $(SOURCE_DIR)/kron3_mod.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/kron3_mod.o \
                                $(SOURCE_DIR)/kron3_mod.f
				
$(OBJECT_DIR)/prof_mod.o:       $(SOURCE_DIR)/prof_mod.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/prof_mod.o \
                                $(SOURCE_DIR)/prof_mod.f				

$(OBJECT_DIR)/aorsa3dMain.o:    $(SOURCE_DIR)/aorsa3dMain.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/aorsa3dMain.o \
                                $(SOURCE_DIR)/aorsa3dMain.f

$(OBJECT_DIR)/aorsaSubs.o:      $(SOURCE_DIR)/aorsaSubs.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/aorsaSubs.o \
                                $(SOURCE_DIR)/aorsaSubs.f

$(OBJECT_DIR)/sigma.o:          $(SOURCE_DIR)/sigma.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/sigma.o \
                                $(SOURCE_DIR)/sigma.f
                                
$(OBJECT_DIR)/zfunction.o:      $(SOURCE_DIR)/zfunction.f $(INCFILES)
	                        $(COMPILE) -o $(OBJECT_DIR)/zfunction.o \
                                $(SOURCE_DIR)/zfunction.f

$(OBJECT_DIR)/ztable.o:         $(SOURCE_DIR)/ztable.f $(INCFILES)
	                        $(COMPILE90) -o $(OBJECT_DIR)/ztable.o \
                                $(SOURCE_DIR)/ztable.f                                
                                
$(OBJECT_DIR)/current.o:        $(SOURCE_DIR)/current.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/current.o \
                                $(SOURCE_DIR)/current.f
                                
$(OBJECT_DIR)/wdot.o:           $(SOURCE_DIR)/wdot.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/wdot.o \
                                $(SOURCE_DIR)/wdot.f                                

$(OBJECT_DIR)/fourier.o:        $(SOURCE_DIR)/fourier.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/fourier.o \
                                $(SOURCE_DIR)/fourier.f
                                
$(OBJECT_DIR)/assert.o:         $(SOURCE_DIR)/assert.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/assert.o \
                                $(SOURCE_DIR)/assert.f
                                
$(OBJECT_DIR)/descinit.o:       $(SOURCE_DIR)/descinit.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/descinit.o \
                                $(SOURCE_DIR)/descinit.f
                                                                                                
$(OBJECT_DIR)/infog2l.o:        $(SOURCE_DIR)/infog2l.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/infog2l.o \
                                $(SOURCE_DIR)/infog2l.f
                                
$(OBJECT_DIR)/pxerbla.o:        $(SOURCE_DIR)/pxerbla.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/pxerbla.o \
                                $(SOURCE_DIR)/pxerbla.f
                                
$(OBJECT_DIR)/setupblacs.o:     $(SOURCE_DIR)/setupblacs.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/setupblacs.o \
                                $(SOURCE_DIR)/setupblacs.f                                                                                                

$(OBJECT_DIR)/check.o:          $(SOURCE_DIR)/check.f $(INCFILES)
	                        $(COMPILEg) -o $(OBJECT_DIR)/check.o \
                                $(SOURCE_DIR)/check.f



