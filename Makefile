MOD_DIR = mod
OBJ_DIR = obj
SRC_DIR = src
OUT_DIR = out

FC = gfortran
FFLAGS = -J$(MOD_DIR) -Wall -g -O2 -ffpe-trap=invalid,zero -ffree-line-length-512 -fcheck=all

EXEC = run_IsochroneSplitting_QuickExample.x

SRC = ErrorManager.f90 \
      HamiltonianSplitting.f90 \
      Isochrone.f90 \
      IsochronePropagator.f90 \
      KeplerEquation.f90 \
      MathConstant.f90 \
      Plummer.f90 \
      run_IsochroneSymplecticIntegrator_QuickExample.f90 \
      SymplecticIntegrator.f90 
      
OBJ_NOPREFIX=$(SRC:.f90=.o)
OBJ=$(addprefix $(OBJ_DIR)/,$(OBJ_NOPREFIX))

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) $^ -o $@

$(OBJ_DIR)/ErrorManager.o: $(SRC_DIR)/ErrorManager.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/HamiltonianSplitting.o: $(SRC_DIR)/HamiltonianSplitting.f90 $(OBJ_DIR)/MathConstant.o $(OBJ_DIR)/SymplecticIntegrator.o $(OBJ_DIR)/IsochronePropagator.o $(OBJ_DIR)/Isochrone.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/Isochrone.o: $(SRC_DIR)/Isochrone.f90 $(OBJ_DIR)/MathConstant.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(OBJ_DIR)/IsochronePropagator.o: $(SRC_DIR)/IsochronePropagator.f90 $(OBJ_DIR)/MathConstant.o $(OBJ_DIR)/ErrorManager.o $(OBJ_DIR)/KeplerEquation.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(OBJ_DIR)/KeplerEquation.o: $(SRC_DIR)/KeplerEquation.f90 $(OBJ_DIR)/MathConstant.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(OBJ_DIR)/MathConstant.o: $(SRC_DIR)/MathConstant.f90
	$(FC) $(FFLAGS) -c $< -o $@
	
$(OBJ_DIR)/Plummer.o: $(SRC_DIR)/Plummer.f90 $(OBJ_DIR)/MathConstant.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/run_IsochroneSymplecticIntegrator_QuickExample.o: $(SRC_DIR)/run_IsochroneSymplecticIntegrator_QuickExample.f90 $(OBJ_DIR)/MathConstant.o $(OBJ_DIR)/HamiltonianSplitting.o $(OBJ_DIR)/Isochrone.o $(OBJ_DIR)/Plummer.o
	$(FC) $(FFLAGS) -c $< -o $@	

$(OBJ_DIR)/SymplecticIntegrator.o: $(SRC_DIR)/SymplecticIntegrator.f90 $(OBJ_DIR)/MathConstant.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(MOD_DIR):
	mkdir -p $(MOD_DIR)
	
$(OUT_DIR):
	mkdir -p $(OUT_DIR)
	
clean:
	rm -rf $(OBJ_DIR)/* $(MOD_DIR)/* $(EXEC)
