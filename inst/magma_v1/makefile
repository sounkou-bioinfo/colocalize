#C++ compiler
CXX ?= g++

#Flags for linker
LD_FLAGS=

#Flags for compiler
CXX_FLAGS=-w -O2 -std=c++11


###########################################################

CXX_FLAGS += -I src/eigen -I src/boost/boost_1_74_0

OBS=src/annotation.o src/baseinput.o src/baseinput_indiv.o src/baseinput_plinkinput.o src/baseinput_snps.o src/buffer.o src/data_genedata.o src/data_snpdata.o src/engineutils.o src/geneanalysis.o src/genecorrelations.o src/gene.o src/geneinput.o src/genemodelengine_base.o src/genemodelengine_compound.o src/genemodelengine_pcreg.o src/genemodelengine_snpwise.o src/geneoutput.o src/genestats_base.o src/genestats_main.o src/genestats_subset.o src/geneutils.o src/input.o src/magma.o src/mathutils.o src/output.o src/parse.o src/permutation.o src/setanalysis.o src/setinput.o src/setmodel.o src/setoutput.o src/setstats.o src/setstats_utils.o src/statutils.o src/utils.o src/utils_dataframe.o src/utils_dataframe_ids.o src/utils_dataframe_variables.o src/utils_dataframe_variables_genecorrelation.o src/utils_dataframe_variables_storage.o

magma: $(OBS)
	$(CXX) $(LD_FLAGS) -o magma $(OBS)

magma.exe: $(OBS)
	$(CXX) $(LD_FLAGS) -o magma.exe $(OBS)
%.o:
	$(CXX) $(CXX_FLAGS) -c $*.cpp -o $*.o

