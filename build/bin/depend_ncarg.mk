ncarg_dummy.o  : $(NCARGD)/ncarg_dummy.f90 
	 cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90)
	 rm -f $(<F:.f90=.f90)

