# Makefile for LAKE model

 exec = lake.out
 #FC=ifort#mpif90
 FC=gfortran
 check_keys = # -check bounds -check pointers
 debug_keys = #-g # debugger

 ifeq ($(FC),ifort)
   opt_keys = -openmp #-O3
 endif
 ifeq ($(FC),gfortran)
   opt_keys = -fopenmp #-O3
 endif

 objfiles_path = ./objfiles/
 model_path = ./source/model/
 driver_path = ./source/driver/
 shared_path = ./source/shared/
 Flake_path = ./source/Flake/

#  Build the executable
 all :
	cd ./source && make all && cd ..
	$(FC) $(objfiles_path)*.o $(debug_keys) $(check_keys) $(opt_keys) -o $(exec)

 doc :
	cd ./docs/doxygen && doxygen mkdoc && cd latex && make && evince refman.pdf && cd .

# Clean all
clean :
	rm -f $(exec)
	rm -f $(objfiles_path)*.o
	rm -f ./source/*.mod

