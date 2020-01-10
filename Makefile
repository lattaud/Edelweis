BINS = Plot_Energy 

CXX = g++
CCFLAGS = -O0 -ggdb -std=c++11 


LD = g++
LDFLAGS = 

LIBS = $(shell root-config --libdir) 
INCS = $(shell root-config --incdir)  

.PHONY: clean rebuild

OBJ = main.o Plot_HEnergy_Voltage.o  
 
all: $(BINS) 

$(BINS): $(OBJ) 
	$(CXX) $^ -o $(BINS) $(CCFLAGS) $(shell root-config --libs) -lMLP  -lTreePlayer  -lXMLIO  -lMLP -lm -g


%.o: %.cpp
	$(CXX) $(CCFLAGS) $^ -c -o $@  $(addprefix -I, $(INCS))

clean: 
	rm -vf *.o */*.o 	$(BINS)

rebuild: clean $(BINS)
