f90=ifort 
opt=
optf77=-O
lib=/home/fb2/wittenburg/F90_library/CDFLIB90/source/cdflib90_1.2/source
exe=fbayesWG

objects=kindnumber sharedTypes inout misc nuisance inverse estimate fastBayesWG
OBJS=$(addsuffix .o, $(objects)) $(exe).o

.SUFFIXES: .f90 .f

.f90.o:
	$(f90) $(opt) -I$(lib) -c $<
.f.o:
	$(f90) $(opt) $(optf77) -c $<
$(exe): $(OBJS)
	$(f90) $(opt) -o $(exe) $(OBJS) $(lib)/libcdflib90.a 
	rm -f *.o *.mod
clean:
	rm -f *.o *.mod $(exe)
