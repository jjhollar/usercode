S1      = ilpair-cms-pp.f
INTER1  = lpair
#
LIBRARY = `cernlib jetset74 geant mathlib pdflib804`
LIBINT	= libilpairg.a
OBJ	= evertex.o
#
###############################################################################
#
# Rule for objects:

FC = g77

CFLAGS = -g -fno-automatic -w -o

RM = /bin/rm  
RMFLAGS = -rf

# Make the executable

clean: 
	$(RM) $(RMFLAGS) ${INTER1}.exe last.kumac* paw.metafile *.log
Clean:
	$(RM) $(RMFLAGS) ${INTER1}.exe last.kumac* paw.metafile *.log *ps *.out ${INTER1}*.hbook OUT/*.* source/*.o
#
run:
	./lpair.exe <ilpair-cms-pp.card > OUT/ilapir-cms-pp.out &
#
# library
#
lib:
	cd ./source;pwd;g77  -fno-automatic -c *.f;/usr/bin/ar cf libilpairg.a *.o;pwd;mv libilpairg.a ..;cd ..
#
# Make the executable
#
$(INTER1):
	$(FC) $(CFLAGS) $@.exe $(S1) $(LIBINT) $(OBJ) $(LIBRARY)
#
