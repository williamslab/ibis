GENETIOSRCS = genetio/personio.cc genetio/marker.cc genetio/superperson.cc genetio/personloopdata.cc genetio/personbulk.cc genetio/nuclearfamily.cc genetio/util.cc
CPPSRCS = IBIS.cc $(GENETIOSRCS)
CSRCS= 
OBJS= $(patsubst %.cc,%.o,$(CPPSRCS)) $(patsubst %.c,%.o,$(CSRCS))


EXEC = ibis

GPP = g++
GCC = gcc
DEFINES= 
CFLAGS = -I. -Wall $(DEFINES) -fopenmp -mpopcnt
# -march=native
CPPFLAGS = -std=c++11 $(CFLAGS)
ifdef DEBUG           # to use run `make DEBUG=1`
  CFLAGS += -g
else
  CFLAGS += -O3
endif

# profiling:
ifdef PROFILE       # to use run `make PROFILE=1
  CFLAGS += -pg
endif

LIBS= -lz
ifdef RELEASE
  LIBS += -static-libstdc++ -static-libgcc
endif
# dependency variables / commands
DEPDIR = .deps
df = $(DEPDIR)/$(*F)

all: $(EXEC) bseg2seg seg2coef

$(EXEC): $(OBJS) $(HEADERS)
	$(GPP) -o $(EXEC) $(OBJS) $(CFLAGS) $(LIBS)

bseg2seg: bseg2seg.cc
	$(GPP) $(CPPFLAGS) -o $@ $^

seg2coef: seg2coef.cc
	$(GPP) $(CPPFLAGS) -o $@ $^

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
	rm -f $(EXEC) $(OBJS) bseg2seg seg2coef

clean-deps:
	rm -f $(DEPDIR)/*.P
