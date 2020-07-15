CFLAGS=-Wall -O3
CXXFLAGS=$(CFLAGS) -std=c++17
LIBS=-lz
PROG=fress

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

KMCAPI=kmc_api/kmer_api.cpp kmc_api/kmc_file.cpp kmc_api/mmer.cpp

all:$(PROG)

fress:fress.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $(KMCAPI) $< $(LIBS)

clean:
	rm $(PROG)
	rm *.o
