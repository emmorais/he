CC=g++
IDIR =include
CFLAGS=-I$(IDIR) -std=c++11
SRC=src
OBJP=obj
LIBS= -g -O2 -pthread -L/usr/local/lib -lntl -lgmp -lm  
TEST=test

# Dependencies:
_DEPS = ringgsw.h gsw.h ltv.h ltv_base.h mp.h plainlwe.h ringlwe.h mpsphf.h rregsphf.h pregsphf.h util.h sample.h functions.h encoding.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
_OBJ = ringgsw.o gsw.o ltv.o ltv_base.o mp.o plainlwe.o ringlwe.o mpsphf.o rregsphf.o pregsphf.o util.o sample.o functions.o encoding.o
OBJ = $(patsubst %,$(OBJP)/%,$(_OBJ))
_OTESTS = gswtest.o ltvtest.o mptest.o lwetest.o sphftest.o
OTESTS = $(patsubst %,$(TEST)/%,$(_TESTS))

$(OBJP)/%.o: $(SRC)/%.cpp 
	$(CC) -g -c -o $@ $< $(CFLAGS)

# Test cases:
test: test_gsw test_ltv test_lwe test_sphf

test_gsw:
	./test/gswtest

test_ltv:
	./test/ltvtest

test_lwe:
	./test/lwetest

test_sphf:
	./test/sphftest

test_mp:
	./test/mptest

all: mptest ltvtest gswtest lwetest sphftest

# Build:
gswtest: $(OBJ) $(OTESTS)
	$(CC) -g test/gswtest.cpp -o $(TEST)/$@ $^ $(CFLAGS) $(LIBS)

ltvtest: $(OBJ) $(OTESTS)
	$(CC) -g test/ltvtest.cpp -o $(TEST)/$@ $^ $(CFLAGS) $(LIBS)

mptest: $(OBJ) $(OTESTS)
	$(CC) -g test/mptest.cpp -o $(TEST)/$@ $^ $(CFLAGS) $(LIBS)

lwetest: $(OBJ) $(OTESTS)
	$(CC) -g test/lwetest.cpp -o $(TEST)/$@ $^ $(CFLAGS) $(LIBS)

sphftest: $(OBJ) $(OTESTS)
	$(CC) -g test/sphftest.cpp -o $(TEST)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: test clean

# Cleaning:
clean:
	rm -f $(OBJP)/*.o *~ core $(INCDIR)/*~ 
	rm -f $(LOBJP)/*.o $(TEST)/*~

