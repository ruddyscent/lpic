bin_PROGRAMS = lpic

if PVM
PVMLIB = -lpvm3
PVMINC =
#PVMINC = $(PVM_ROOT)/include
else
PVMLIB =
PVMINC =
endif

INCLUDES = $(PVMINC) 
LDFLAGS = $(PVMLIB)

dist_data_DATA = Makefile.old plain.h parallel.h parallel_or_plain.h

lpic_SOURCES = \
	error.C \
	parameter.C \
	readfile.C \
	box.C \
	domain.C \
	pulse.C \
	diagnostic_stepper.C \
	diagnostic_trace.C \
	diagnostic_spacetime.C \
	diagnostic_energy.C \
	diagnostic_reflex.C \
	diagnostic_flux.C \
	diagnostic_poisson.C \
	diagnostic_phasespace.C \
	diagnostic_snapshot.C \
	diagnostic_velocity.C \
	diagnostic.C \
	propagate.C \
	propagate_fields.C \
	propagate_particles.C \
	stack.C \
	matrix.C \
	uhr.C \
	main.C \
	network.C

include_HEADERS = \
	box.h \
	cell.h \
	common.h \
	debug.h \
	diagnostic.h \
	diagnostic_stepper.h \
	diagnostic_trace.h \
	diagnostic_spacetime.h \
	diagnostic_energy.h \
	diagnostic_reflex.h \
	diagnostic_flux.h \
	diagnostic_poisson.h \
	diagnostic_phasespace.h \
	diagnostic_snapshot.h \
	diagnostic_velocity.h \
	domain.h \
	error.h \
	main.h \
	matrix.h \
	network.h \
	parameter.h \
	particle.h \
	propagate.h \
	pulse.h \
	readfile.h \
	stack.h \
	uhr.h \
	units.h \
	network.h
