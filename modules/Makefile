ROOTDIR = $(abspath ..)
include $(ROOTDIR)/Make.config

DEFAULT_MATLAB_DIR_LINUX = "/usr/local/MATLAB/$(DEFAULT_MATLAB_VERSION)"
DEFAULT_MATLAB_DIR_DARWIN = "/Applications/MATLAB_${DEFAULT_MATLAB_VERSION}.app"

UNAME = $(shell uname)

ifeq ($(strip $(MATLAB_DIR)),)
  ifeq ($(UNAME), Linux)
    MATLAB_DIR = "${DEFAULT_MATLAB_DIR_LINUX}"
  else
    ifeq ($(UNAME), Darwin)
      MATLAB_DIR = "${DEFAULT_MATLAB_DIR_DARWIN}"
    endif
  endif
endif
ifeq ($(strip $(MATLAB_DIR)),)
  $(error "MATLAB directory must be specified in the file Make.config for OS's other than Linux or Mac OS X")
endif

source_name1 = compute_true_frequencies
output_name1 = "compute_true_frequencies_mex"

source_name2 = read_alignment_fasta
output_name2 = "read_alignment_fasta_mex"

CFLAGS = -Wall -O3 -g -fPIC -lm -I${MATLAB_DIR}/extern/include -std=c99 -DPACKBITS=$(PACKBITS)

ifeq ($(OPENMP),1)
  CFLAGS += -fopenmp -DOMP
  LDFLAGS += -fopenmp -lgomp
endif

.PHONY: all matlab-module1 matlab-module2

allmodules = matlab-module1 matlab-module2

all: $(allmodules)

$(allmodules) :
	@test -d ${MATLAB_DIR} || { echo "MATLAB directory not found. You may set the correct path in the Make.config file" && exit 1; }

matlab-module1: .mlab1

.mlab1: ${source_name1}.c
	${MATLAB_DIR}/bin/mex -v -largeArrayDims "CFLAGS=\$$CFLAGS ${CFLAGS}" "LDFLAGS=\$$LDFLAGS ${LDFLAGS}" -lut ${source_name1}.c -output ${output_name1} && touch .mlab1

matlab-module2: .mlab2

.mlab2: ${source_name2}.c
	${MATLAB_DIR}/bin/mex -v -largeArrayDims "CFLAGS=\$$CFLAGS ${CFLAGS}" "LDFLAGS=\$$LDFLAGS ${LDFLAGS}" -lz -lut ${source_name2}.c -output ${output_name2} && touch .mlab2

clean:
	rm -f *.o .mlab* *.mex*
