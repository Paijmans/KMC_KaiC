# Configuration of the executable
TARGET = KMCKaiC
#ODIR=obj
RUNDIR = RUN
 

# Compiler configuration
CXX      = g++ 
#CXXFLAGS = -Wall -Werror -Wextra -Wshadow -g
CXXFLAGS = -Wall -Wextra -Wshadow -O3 
CXXFLAGSPG = -Wall -Wextra -Wshadow -O3 -pg
COMP = -c 
#INC = -I ~/include -L ~/lib
#INC_LIB = -lgsl -lgslcblas -lm  


OBJ = main.o \
      Monomer.o \
      Hexamer.o \
      PropensityContainer.o \
      propagate.o \
      random.o \
          
#OBJ = $(patsubst %,$(ODIR)/%,$(OBJS)) 
    
default: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ) $(INC) $(INC_LIB)   
	
new:
	@echo "clean ..."
	@$(RM) $(TARGET) *.o *~
	$(MAKE)

# clean
clean:
	@echo "clean ..."
	@$(RM) $(TARGET) *.o *~
	
#With profiling
pg: $(OBJ) 
	$(CXX) $(CXXFLAGSPG) -o $(TARGET) $(OBJ) $(INC) $(INC_LIB) 
  
