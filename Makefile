all:executable

CFLAGS=-Wall -O3
CXXFLAGS=$(CFLAGS) -std=c++17
KMC_API_DIR = kmc_api
PROG=fress

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

## KMCAPI=kmc_api/kmer_api.cpp kmc_api/kmc_file.cpp kmc_api/mmer.cpp

KMC_API_OBJS = \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o

FRESS_OBJS = \
fresslib.o \
fress.o

$(FRESS_OBJS) $(KMC_API_OBJS): %.o: %.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

executable: $(KMC_API_OBJS) $(FRESS_OBJS)
	$(CXX) -o $(PROG) $(FRESS_OBJS) $(KMC_API_OBJS)

clean:
	rm $(PROG) 
	rm *.o
	rm $(KMC_API_OBJS)
