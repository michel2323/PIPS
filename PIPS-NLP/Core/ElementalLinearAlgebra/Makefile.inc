EMTLLINALGDIR = $(srcdir)/ElementalLinearAlgebra

EMTLOBJECTS =  $(EMTLLINALGDIR)/EmtlGenIndefSolver.o \
  $(EMTLLINALGDIR)/EmtlStochSymIndefSolver.o \
	$(EMTLLINALGDIR)/EmtlSymPSDSolver.o \
	$(EMTLLINALGDIR)/EmtlDenSymMatrix.o \
	$(EMTLLINALGDIR)/EmtlDenGenMatrix.o \
	$(EMTLLINALGDIR)/EmtlLinearAlgebraPackage.o \
	$(EMTLLINALGDIR)/EmtlVector.o \
  $(EMTLLINALGDIR)/EmtlContext.o

EMTLHEADERS = $(patsubst $(EMTLLINALGDIR)/%,$(localincludedir)/%,\
	$(wildcard $(EMTLLINALGDIR)/*.h))

libooqpemtl_STATIC = $(locallibdir)/libooqpemtl.a
libooqpemtl_OBJ    = $(EMTLOBJECTS)
