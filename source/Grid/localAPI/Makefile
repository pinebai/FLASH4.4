Grid += gr_ptInterface.o gr_ptInit.o gr_ptMove.o\
        gr_ptFinalize.o gr_ptMarkRefineDerefine.o \
	gr_ptMoveOffBlk.o\
        gr_ptFillBlkParticleInfo.o gr_ensureValidNeighborInfo.o \
        gr_solversInit.o gr_solversFinalize.o gr_solversTestInit.o \
	gr_bcInterface.o gr_bcHseInit.o gr_bcApplyToAllBlks.o gr_bcApplyToOneFace.o \
	gr_hgInit.o  gr_mpoleInit.o gr_isoMpoleInit.o gr_pfftInit.o \
	gr_hgFinalize.o  gr_mpoleFinalize.o gr_isoMpoleFinalize.o gr_pfftFinalize.o \
        gr_bcInit.o gr_bcFinalize.o gr_hgInterface.o gr_pfftInterface.o \
	gr_isoInterface.o gr_ptMapInit.o gr_hgPfftInit.o gr_hgPfftFinalize.o \
	gr_interface.o gr_interfaceTypeDecl.o gr_pfftInterfaceTypeDecl.o \
	gr_getInteriorBlkPtr.o gr_releaseInteriorBlkPtr.o gr_setBlockType.o \
	gr_mgInterface.o gr_mgInit.o gr_mgFinalize.o gr_mgPfftInit.o gr_mgPfftFinalize.o \
        gr_bicgInterface.o gr_bicgInit.o gr_bicgFinalize.o 

Grid += gr_hg_amr_1blk_bcset_work.o gr_pfftPoissonDirect.o

Grid +=  gr_hypreComputeB.o gr_hypreSetupSolver.o gr_hypreDestroySolver.o  \
         gr_hypreDestroyGrid.o gr_hypreSetupGrid.o gr_hypreCreateMatrix.o \
	 gr_hypreCreateMatrixAnisoCond.o gr_hypreCreateMatrixFcBAnisoCond.o \
         gr_hypreCreateMatrix1Blk.o \
         gr_hypreCreateMatrixFcB.o \
         gr_hypreSolve.o gr_hypreUpdateSoln.o gr_hypreSetIniGuess.o gr_hypreExchangeFacB.o \
         gr_hypreGridStatus.o gr_hypreInit.o gr_hypreFinalize.o gr_hypreApplyBcToFace.o \
         gr_hypreGetFaceBFcB.o gr_hypreExchangeFacBFcB.o \
         gr_hypreGetFaceB.o gr_hypreAddGraph.o

Grid += gr_mpoleAllocateRadialArrays.o   \
        gr_mpoleCen1Dspherical.o         \
        gr_mpoleCen2Dcylindrical.o       \
        gr_mpoleCen2Dspherical.o         \
        gr_mpoleCen3Dcartesian.o         \
        gr_mpoleCen3Dcylindrical.o       \
        gr_mpoleCenterOfExpansion.o      \
        gr_mpoleCollectMoments.o         \
        gr_mpoleDeallocateRadialArrays.o \
        gr_mpoleDumpMoments.o            \
        gr_mpoleFinalize.o               \
        gr_mpoleHeapsort.o               \
        gr_mpoleInit.o                   \
        gr_mpoleMom1Dspherical.o         \
        gr_mpoleMom2Dcylindrical.o       \
        gr_mpoleMom2Dspherical.o         \
        gr_mpoleMom3Dcartesian.o         \
        gr_mpoleMom3Dcylindrical.o       \
        gr_mpoleMomBins1Dspherical.o     \
        gr_mpoleMomBins2Dcylindrical.o   \
        gr_mpoleMomBins2Dspherical.o     \
        gr_mpoleMomBins3Dcartesian.o     \
        gr_mpoleMomBins3Dcylindrical.o   \
        gr_mpoleMoments.o                \
        gr_mpolePot1Dspherical.o         \
        gr_mpolePot2Dcylindrical.o       \
        gr_mpolePot2Dspherical.o         \
        gr_mpolePot3Dcartesian.o         \
        gr_mpolePot3Dcylindrical.o       \
        gr_mpolePotentials.o             \
        gr_mpolePrintRadialInfo.o        \
        gr_mpoleRad1Dspherical.o         \
        gr_mpoleRad2Dcylindrical.o       \
        gr_mpoleRad2Dspherical.o         \
        gr_mpoleRad3Dcartesian.o         \
        gr_mpoleRad3Dcylindrical.o       \
        gr_mpoleRadialSampling.o         \
        gr_mpoleSetInnerZoneGrid.o       \
        gr_mpoleSetOuterZoneGrid.o       \
        gr_mpoleSetRadialBinData.o

Grid += gr_sbInterface.o gr_sbInit.o gr_sbSendBoundBox.o\
	gr_sbSendParticles.o gr_sbCreateParticles.o gr_sbGetProcBlock.o gr_sbSendPosn.o gr_sbStoreParticlesPerProc.o \
	gr_sbSendParticleCount.o gr_sbSendParticles.o gr_sbUpdateForces.o gr_sbSendForces.o gr_sbFinalize.o

Grid += gr_bhInterface.o gr_bhInit.o gr_bhFinalize.o

Grid += gr_xyzToBlockLevel.o gr_xyzToBlock.o

gr_interface.o     : gr_interfaceTypeDecl.o
gr_pfftInterface.o : gr_pfftInterfaceTypeDecl.o
