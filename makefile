# Makefile for s2fil


# ======== OPTIONS ========

USEPGPLOT = no
USEFFTW   = yes


# ======== COMPILER ========

FC      = gfortran
#FC      = f95
#FC      = g95

ifneq ($(USEPGPLOT),yes)
  OPTPGPLOT     = -DNO_PGPLOT
endif
OPT = $(OPTPGPLOT) -m64 -O3 -DS2FIL_VERSION=\"1.0b2\" -DS2FIL_BUILD=\"`svnversion -n .`\" 


# ======== LINKS ========

PROGDIR = ..

HPIXDIR = $(PROGDIR)/Healpix
HPIXLIB = $(HPIXDIR)/lib
HPIXLIBNM= healpix
HPIXINC = $(HPIXDIR)/include

S2DIR  = $(PROGDIR)/s2
S2LIB  = $(S2DIR)/lib
S2LIBNM= s2
S2INC  = $(S2DIR)/include
S2DOC  = $(S2DIR)/doc

CSWTDIR  = $(PROGDIR)/fastcswt
CSWTLIB  = $(CSWTDIR)/lib
CSWTLIBNM= fastcswt
CSWTINC  = $(CSWTDIR)/include
CSWTDOC  = $(CSWTDIR)/doc

COMBDIR  = $(PROGDIR)/comb
COMBLIB  = $(COMBDIR)/lib
COMBLIBNM= comb
COMBINC  = $(COMBDIR)/include
COMBDOC  = $(COMBDIR)/doc

S2FILDIR  = $(PROGDIR)/s2fil
S2FILLIB  = $(S2FILDIR)/lib
S2FILLIBNM= s2fil
S2FILINC  = $(S2FILDIR)/include
S2FILSRC  = $(S2FILDIR)/src/mod
S2FILPROG = $(S2FILDIR)/src/prog
S2FILBIN  = $(S2FILDIR)/bin
S2FILDOC  = $(S2FILDIR)/doc

CFITSIOLIB   = $(PROGDIR)/cfitsio/lib
CFITSIOLIBNM = cfitsio

PGPLOTLIB    = $(PROGDIR)/pgplot
PGPLOTLIBNM  = pgplot
X11LIB       = /usr/X11R6/lib
X11LIBNM     = X11

FFTWLIB      = $(PROGDIR)/fftw-2.1.5/lib
FFTWLIBNM    = fftw
RFFTWLIBNM   = rfftw


# ======== FFFLAGS ========

FFLAGS  = -I$(HPIXINC) -I$(S2INC) -I$(CSWTINC) -I$(COMBINC) -I$(S2FILINC)


# ======== LDFLAGS ========

ifeq ($(USEPGPLOT),yes)
  LDFLAGSPGPLOT = -L$(PGPLOTLIB) -L$(X11LIB) \
                  -l$(PGPLOTLIBNM) -l$(X11LIBNM)
endif

ifeq ($(USEFFTW),yes)
  LDFLAGSFFTW   = -L$(FFTWLIB) \
                  -l$(RFFTWLIBNM) -l$(FFTWLIBNM) 
endif

LDFLAGS =  -L$(S2FILLIB) -l$(S2FILLIBNM) \
	   -L$(CSWTLIB) -l$(CSWTLIBNM) \
	   -L$(COMBLIB) -l$(COMBLIBNM) \
	   -L$(S2LIB) -l$(S2LIBNM) \
           -L$(HPIXLIB) -l$(HPIXLIBNM) \
           -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) \
           $(LDFLAGSFFTW) $(LDFLAGSPGPLOT)


# ======== PPFLAGS ========

ifeq ($(FC),f95)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),g95)
  PPFLAGS = -cpp $(OPT)
else ifeq ($(FC),gfortran)
  PPFLAGS = -x f95-cpp-input $(OPT)
endif

# ======== OBJECT FILES TO MAKE ========

S2FILOBJ = $(S2FILINC)/s2fil_types_mod.o  \
           $(S2FILINC)/s2fil_error_mod.o  \
           $(S2FILINC)/s2fil_filter_mod.o \
           $(S2FILINC)/s2fil_field_mod.o


# ======== MAKE RULES ========

default: all

all:     lib prog

lib:	 $(S2FILLIB)/lib$(S2FILLIBNM).a

prog:    $(S2FILBIN)/s2fil_filter_construct   \
         $(S2FILBIN)/s2fil_field_construct    \
	 $(S2FILBIN)/s2fil_localisation_thres \
         $(S2FILBIN)/s2fil_draw_dots          \
         $(S2FILBIN)/s2fil_draw_dots_only     \
         $(S2FILBIN)/s2fil_about

$(S2FILINC)/%.o: $(S2FILSRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(S2FILINC)

$(S2FILINC)/%.o: $(S2FILPROG)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 


# Library

$(S2FILLIB)/lib$(S2FILLIBNM).a: $(S2FILOBJ)
	ar -r $(S2FILLIB)/lib$(S2FILLIBNM).a $(S2FILOBJ)


# Documentation

docs:
	f90doc $(S2FILSRC)/*.f90
	f90doc $(S2FILPROG)/*.f90
	ln_multi $(S2DOC)/s2_*.html
	ln_multi $(S2DOC)/index_s2.html
	ln_multi $(COMBDOC)/comb_*.html
	ln_multi $(COMBDOC)/index_comb.html
	ln_multi $(CSWTDOC)/cswt_*.html
	ln_multi $(CSWTDOC)/index_cswt.html
	mv *.html $(S2FILDOC)/.
	addstyle $(S2FILDOC)/s2fil_*.html

cleandocs:
	rm -f $(S2FILDOC)/s2fil_*.html
	rm -f $(S2FILDOC)/s2_*.html
	rm -f $(S2FILDOC)/index_s2.html
	rm -f $(S2FILDOC)/comb_*.html
	rm -f $(S2FILDOC)/index_comb.html
	rm -f $(S2FILDOC)/cswt_*.html
	rm -f $(S2FILDOC)/index_cswt.html


# Clean up

clean:	tidy
	rm -f $(S2FILINC)/*.mod
	rm -f $(S2FILINC)/*.o
	rm -f $(S2FILLIB)/lib$(S2FILLIBNM).a
	rm -f $(S2FILBIN)/*

tidy:
	rm -f *.mod
	rm -f $(S2FILSRC)/*~
	rm -f $(S2FILPROG)/*~


# Module dependencies

$(S2FILINC)/s2fil_types_mod.o: $(S2FILSRC)/s2fil_types_mod.f90
$(S2FILINC)/s2fil_error_mod.o: $(S2FILSRC)/s2fil_error_mod.f90  \
                                 $(S2FILINC)/s2fil_types_mod.o
$(S2FILINC)/s2fil_filter_mod.o:$(S2FILSRC)/s2fil_filter_mod.f90 \
                                 $(S2FILINC)/s2fil_types_mod.o  \
                                 $(S2FILINC)/s2fil_error_mod.o
$(S2FILINC)/s2fil_field_mod.o: $(S2FILSRC)/s2fil_field_mod.f90  \
                                 $(S2FILINC)/s2fil_filter_mod.o \
                                 $(S2FILINC)/s2fil_types_mod.o  \
                                 $(S2FILINC)/s2fil_error_mod.o


# Program dependencies and compilation

$(S2FILINC)/s2fil_test.o:	$(S2FILPROG)/s2fil_test.f90 lib
$(S2FILBIN)/s2fil_test:	$(S2FILINC)/s2fil_test.o
	$(FC) -o $(S2FILBIN)/s2fil_test $(S2FILINC)/s2fil_test.o $(LDFLAGS) $(PPFLAGS)

$(S2FILINC)/s2fil_filter_construct.o: $(S2FILPROG)/s2fil_filter_construct.f90 lib
$(S2FILBIN)/s2fil_filter_construct:   $(S2FILINC)/s2fil_filter_construct.o
	$(FC) -o $(S2FILBIN)/s2fil_filter_construct $(S2FILINC)/s2fil_filter_construct.o \
        $(LDFLAGS) $(PPFLAGS)

$(S2FILINC)/s2fil_field_construct.o: $(S2FILPROG)/s2fil_field_construct.f90 lib
$(S2FILBIN)/s2fil_field_construct:   $(S2FILINC)/s2fil_field_construct.o
	$(FC) -o $(S2FILBIN)/s2fil_field_construct $(S2FILINC)/s2fil_field_construct.o \
        $(LDFLAGS) $(PPFLAGS)

$(S2FILINC)/s2fil_localisation_thres.o: $(S2FILPROG)/s2fil_localisation_thres.f90 lib
$(S2FILBIN)/s2fil_localisation_thres:   $(S2FILINC)/s2fil_localisation_thres.o
	$(FC) -o $(S2FILBIN)/s2fil_localisation_thres \
        $(S2FILINC)/s2fil_localisation_thres.o $(LDFLAGS) $(PPFLAGS)

$(S2FILINC)/s2fil_draw_dots.o: $(S2FILPROG)/s2fil_draw_dots.f90 lib
$(S2FILBIN)/s2fil_draw_dots:   $(S2FILINC)/s2fil_draw_dots.o
	$(FC) -o $(S2FILBIN)/s2fil_draw_dots \
        $(S2FILINC)/s2fil_draw_dots.o $(LDFLAGS) $(PPFLAGS)

$(S2FILINC)/s2fil_draw_dots_only.o: $(S2FILPROG)/s2fil_draw_dots_only.f90 lib
$(S2FILBIN)/s2fil_draw_dots_only:   $(S2FILINC)/s2fil_draw_dots_only.o
	$(FC) -o $(S2FILBIN)/s2fil_draw_dots_only \
        $(S2FILINC)/s2fil_draw_dots_only.o $(LDFLAGS) $(PPFLAGS)

$(S2FILINC)/s2fil_about.o: $(S2FILPROG)/s2fil_about.f90 lib
$(S2FILBIN)/s2fil_about:   $(S2FILINC)/s2fil_about.o
	$(FC) -o $(S2FILBIN)/s2fil_about \
        $(S2FILINC)/s2fil_about.o $(LDFLAGS) $(PPFLAGS)

