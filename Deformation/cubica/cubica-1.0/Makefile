SHELL := /bin/bash -e

#BINS = ArpackLMA BasisViewer BdTreeTest ClonedViewer CollisionCubatureViewer CollisionViewer Constrainer CubatureGenerator CubatureViewer Embedder GlslBasicViewer GlslObjViewer GlslSceneViewer IsoStuffer MeshExerciser MeshViewer ModalDerivatives ModalDerivativesCubatureGenerator MultiArpackLMA MultiBdTreeTester MultiPCAFullPartitionedCubatureGenerator MultiPCAPartitionedCubatureGenerator MultiPartitionedCubatureGenerator MultiPartitionedEmbedder MultiPartitioner MultiTrainingSampleViewer MultiZeroEigsLMA OdeCapsulePlacer OdeCubatureGenerator OdeCubatureViewer OdeCubatureViewerAOV OdeDropTestCubatureViewer OdeIsoStuffer OdePartitioner OdePoseLMACubatureGenerator OdePoseSpaceCubatureViewer OdePoseSpaceSampler OdeSkeletonViewer OdeSkinnedPreview OdeTetAnimator OdeTetViewer OnlineBeamViewer OnlineMeshViewer PCACubatureGenerator Partitioner PinocchioAOV PinocchioCubatureGenerator PinocchioLMAPCACubatureGenerator PinocchioMultiAnimator PinocchioMultiAnimatorAOV PinocchioMultiMocapCubatureGenerator PinocchioPartitioner PinocchioPrettyMultiAnimator PinocchioTetAnimator PinocchioTetViewer PinocchioViewer ZeroEigsLMA
BINS = IsoStuffer ArpackLMA ModalDerivatives SkinnedViewer CubatureViewer ModalDerivativesCubatureGenerator OdeIsoStuffer OdeTetViewer OdePartitioner OdePoseSpaceSampler OdePoseLMACubatureGenerator OdePrecomputeSandwiches OdeCubatureViewer PrintVector
ODEBINS = OdeIsoStuffer OdeTetViewer ArpackLMA ModalDerivatives OdePartitioner OdePoseSpaceSampler OdePoseLMACubatureGenerator OdePrecomputeSandwiches OdeCubatureViewer

SRCS = integrator cubature geometry linearalgebra material rendering util 

all : 
	-for d in $(BINS); do (echo -e "\n==== Entering $$d ====\n";cd ./projects/$$d; make -f Makefile.osx;cd ../..); done

osx : 
	-for d in $(BINS); do (echo -e "\n==== Entering $$d ====\n";cd ./projects/$$d; make -f Makefile.osx;cd ../..); done

linux : 
	-for d in $(BINS); do (echo -e "\n==== Entering $$d ====\n";cd ./projects/$$d; make -j 8 -f Makefile.linux;cd ../..); done

odebins : 
	-for d in $(ODEBINS); do (echo -e "\n==== Entering $$d ====\n";cd ./projects/$$d; make -j 8 -f Makefile.osx;cd ../..); done

regressionTest: 
	-for d in $(BINS); do (echo -e "\n==== Entering $$d ====\n";cd ./projects/$$d; make -f Makefile.osx regression;cd ../..); done

clean: 
	-for d in $(BINS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./projects/$$d; make -f Makefile.osx clean;cd ../..); done 
	-for d in $(SRCS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./src/$$d; rm *.o; cd ../..); done

objclean:
	-for d in $(SRCS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./src/$$d; rm -f *.o; cd ../..); done

ctags:
	ctags -R *
