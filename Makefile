LIB_DIR = ./lib
CFITSIO_DIR = cfitsio
RNG_DIR = rngstreams-1.0.1

CC = gcc -Ofast -march=native -mfpmath=sse -flto -ftree-vectorize -pedantic -std=c11 -fopenmp -ggdb -fno-omit-frame-pointer

#for Apple clang uncomment the following line
#CC = clang -Ofast -march=native -mfpmath=sse -flto -ftree-vectorize -pedantic -std=c11  -ggdb -fno-omit-frame-pointer


## Desctription of targets

all:  ./bin/squeeze

./bin/squeeze:  ./src/read_fits.c ./src/free_fits.c ./src/squeeze.c ./src/squeeze.h ./src/regularizations.c ./src/models/modelcode.c ./src/extract_oifits.c $(LIB_DIR)/$(CFITSIO_DIR)/libcfitsio.a $(LIB_DIR)/$(RNG_DIR)/src/.libs/librngstreams.la

	@echo 'Building SQUEEZE...' &&  ${CC} ./src/squeeze.c ./src/read_fits.c ./src/free_fits.c -I$(LIB_DIR)/$(CFITSIO_DIR) -I$(LIB_DIR)/$(RNG_DIR)/src -o ./bin/squeeze -lm -L$(LIB_DIR)/$(RNG_DIR)/src/.libs/ -lrngstreams -L$(LIB_DIR)/${CFITSIO_DIR} -lcfitsio && echo 'done !' 

$(LIB_DIR)/$(RNG_DIR)/Makefile:
	@cd $(LIB_DIR)/$(RNG_DIR) && ./configure --silent 

$(LIB_DIR)/$(RNG_DIR)/src/.libs/librngstreams.la: $(LIB_DIR)/$(RNG_DIR)/Makefile
	@cd $(LIB_DIR)/$(RNG_DIR) && echo 'Building the random number generator library...please wait...' && make -s

$(LIB_DIR)/$(CFITSIO_DIR)/Makefile:
	@echo 'Configuring a reentrant local version of cfitsio' && cd $(LIB_DIR)/$(CFITSIO_DIR) && ./configure --silent --enable-reentrant --enable-ssse3

$(LIB_DIR)/$(CFITSIO_DIR)/libcfitsio.a: $(LIB_DIR)/$(CFITSIO_DIR)/Makefile
	@cd $(LIB_DIR)/$(CFITSIO_DIR) && echo 'Building the extracted reentrant local version of cfitsio... please wait...'  && make -s 

cleanall:
	rm -f ./bin/squeeze 
	rm -f ./src/*.o
	cd $(LIB_DIR)/$(CFITSIO_DIR) && make clean
	rm -f $(LIB_DIR)/$(CFITSIO_DIR)/Makefile
	cd $(LIB_DIR)/$(RNG_DIR) && make clean
	rm -f $(LIB_DIR)/$(RNG_DIR)/Makefile

clean:
	rm -f ./bin/squeeze
	rm -f ./src/*.o
