##
SOURCES := 	mcpsc.C \
			ce_align.C \
			pom/pdbutil.C  \
			pom/pom.C \
			pom/ipdb.C \
			cmp_util.C \
			pom/miscutil.C \
			pom/jdate.C \
			pom/ipdb_df.C \
			pom/linkedid.C \
			pom/derive.C \
			tmalign.C \
			usm.C

MODULES := 	ce_align.C \
			pom/pdbutil.C  \
			pom/pom.C \
			pom/ipdb.C \
			cmp_util.C \
			pom/miscutil.C \
			pom/jdate.C \
			pom/ipdb_df.C \
			pom/linkedid.C \
			pom/derive.C \
			tmalign.C \
			usm.C

PROGRAM := mcpsc
OBJECTS := $(SOURCES:.C=.o)
MOBJECTS := $(MODULES:.C=.o)

CC := g++
CE_CFLAGS := $(CFLAGS) -DFUNCPROTO -DREAD_WRITE -DSGI -O3 -fopenmp -std=c++11
CE_LFLAGS := $(LFLAGS) -lm -ffast-math -fopenmp -lboost_iostreams -lboost_system -lpthread -lz #-static 
#-lboost_iostreams -lboost_system -lpthread

test: $(MOBJECTS) 
	$(CC) $(CE_CFLAGS) -Wall -c mcpsc.C -o mcpsc.o -D ONLY_TMALIGN #-D ONLY_TMALIGN -D ONLY_CE -D ONLY_USM
	$(CC) -o mcpsc mcpsc.o $(MOBJECTS) $(CE_LFLAGS) 

$(PROGRAM): $(OBJECTS) $(RCKSKEL_ARCHIVE) $(RCCE_ARCHIVE)
	@echo [$(OBJECTS)] are the objects
	$(CC) -o $@ $(OBJECTS) $(CE_LFLAGS)

clean:
	rm -f *.o pom/*.o core $(PROGRAM)

.C.o: 
	$(CC) $(CE_CFLAGS) -Wall -c $< -o $@

##
