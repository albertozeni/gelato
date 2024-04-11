CXX = g++
BIN_DIR = bin
CXX_FLAGS = -O3 -std=c++17 -fopenmp
LIBS = -lpthread -lstdc++fs -lz 
PROG = $(BIN_DIR)/gelato
OBJS = src/main.o src/gelato.o src/graphexplore.o src/kmercount.o
INCLUDES = -Isrc/

ifdef RANDOM_REF
	CXX_FLAGS += -DRANDOM_REF
endif

.SUFFIXES:.cpp .o

%.o: %.cpp
	$(CXX) -c $(CXX_FLAGS) $(INCLUDES) $< -o $@ 

all:$(PROG)

 $(BIN_DIR)/gelato: $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXX_FLAGS) $^ -o $@ $(LIBS) 

clean:
	rm -rf bin/ src/*.o

src/main.o : src/*.hpp
src/gelato.o : src/*.hpp
src/graphexplore.o : src/*.hpp
src/kmercount.o : src/*.hpp
