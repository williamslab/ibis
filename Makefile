#CPPSRCS = TransposeBSQUID_parallel.cc
#CPPSRCS = TransposeBSQUID_no_leniency_2.cc


CPPSRCS = IBIS.cc

CSRCS= 
OBJS= $(patsubst %.cc,%.o,$(CPPSRCS)) $(patsubst %.c,%.o,$(CSRCS))

#EXEC = IBISv0.99_no_leniency_2
#EXEC = IBISv0.99_parallel_3

EXEC = ibis
#EXEC = IBISv1.01


GPP = g++
GCC = gcc
DEFINES= 
CFLAGS = -Wall $(DEFINES) -fopenmp
CPPFLAGS = -std=c++11 $(CFLAGS)
ifdef DEBUG           # to use run `make DEBUG=1`
  CFLAGS += -g
else
  CFLAGS += -O2
endif

# profiling:
ifdef PROFILE       # to use run `make PROFILE=1
  CFLAGS += -pg
endif

LIBS= -lgenetio -lz -lpthread

# dependency variables / commands
DEPDIR = .deps
df = $(DEPDIR)/$(*F)

all: $(EXEC)

$(EXEC): $(OBJS) $(HEADERS)
	$(GPP) -o $(EXEC) $(OBJS) $(CFLAGS) $(LIBS)

# This way of building dependencies (per-file) described at
# http://make.paulandlesley.org/autodep.html

.c.o:
	@mkdir -p $(DEPDIR)
	$(GCC) -MMD $(CFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d

.cc.o:
	@mkdir -p $(DEPDIR)
	$(GPP) -MMD $(CPPFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d


# include the .P dependency files, but don't warn if they don't exist (the -)
-include $(CPPSRCS:%.cc=$(DEPDIR)/%.P)
-include $(CSRCS:%.c=$(DEPDIR)/%.P)
# The following applies if we don't use a dependency directory:
#-include $(SRCS:.cc=.P)

tags: $(SRCS) *.h
	ctags --language-force=c++ --extra=+q --fields=+i --excmd=n *.cc *.h

clean:
	rm -f $(EXEC) $(OBJS)

clean-deps:
	rm -f $(DEPDIR)/*.P
