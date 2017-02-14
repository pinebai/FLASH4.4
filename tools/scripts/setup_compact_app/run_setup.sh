#!/usr/bin/env bash

#Abort if any error
set -e

FLASH_ROOT="../../.."
REVISION=$(svn info | awk '/Revision/{print $2}')
COMPACT_APPS="compact_apps_r${REVISION}"
COMPACT_DIR="${FLASH_ROOT}/${COMPACT_APPS}"

if [ -d ${COMPACT_DIR} ]; then
    echo "Compact app directory already exists.  Remove before running this script"
    exit -1
fi

mkdir ${COMPACT_DIR}
mkdir ${COMPACT_DIR}/UG
mkdir ${COMPACT_DIR}/PM

# ********** EOS
#We don't gain anything by running EOS tests with AMR.
./setup_compact_app.sh unitTest/Eos/Multigamma -auto +nofbs -3d -portable
mv ${FLASH_ROOT}/Multigamma ${COMPACT_DIR}/UG/

./setup_compact_app.sh unitTest/Eos/Helmholtz -auto +nofbs -3d -portable
mv ${FLASH_ROOT}/Helmholtz ${COMPACT_DIR}/UG/

# ********** Parallel FFT and hybrid parallel FFT
./setup_compact_app.sh unitTest/PFFT_PoissonFD -auto -3d +nofbs -parfile=test_UG_4p_3d_128cube.par -portable
mv ${FLASH_ROOT}/PFFT_PoissonFD ${COMPACT_DIR}/UG/

./setup_compact_app.sh unitTest/Gravity/Poisson3 -auto -3d +pm4dev -maxblocks=600 -unit=Grid/GridSolvers/Multigrid/PfftTopLevelSolve PfftSolver=HomBcTrigSolver -portable
mv ${FLASH_ROOT}/Poisson3 ${COMPACT_DIR}/PM/Poisson3_multigrid

# ********** Multipole
./setup_compact_app.sh unitTest/Gravity/Poisson3 -auto -3d +nofbs -parfile=flash.par.ug +newMpole -portable
mv ${FLASH_ROOT}/Poisson3 ${COMPACT_DIR}/UG/Poisson3_multipole

./setup_compact_app.sh unitTest/Gravity/Poisson3 -auto -3d +pm4dev -maxblocks=600 +newMpole -portable
mv ${FLASH_ROOT}/Poisson3 ${COMPACT_DIR}/PM/Poisson3_multipole

# ********** Guardcell exchange
./setup_compact_app.sh unitTest/Grid/GCell -auto -3d +nofbs -parfile=test_3d.par -portable
mv ${FLASH_ROOT}/GCell ${COMPACT_DIR}/UG/

./setup_compact_app.sh unitTest/Grid/GCell -auto -3d +pm4dev -parfile=test_3d.par -portable
mv ${FLASH_ROOT}/GCell ${COMPACT_DIR}/PM/

# ********** I/O
./setup_compact_app.sh unitTest/IO -auto -3d +nofbs withParticles=True +parallelio -portable
mv ${FLASH_ROOT}/IO ${COMPACT_DIR}/UG/

./setup_compact_app.sh unitTest/IO/IOTypes -auto -3d -nofbs +ug -unit=source/IO/IOMain/hdf5/parallel/PM_argonne withParticles=True parallelIO=True -portable
mv ${FLASH_ROOT}/IOTypes ${COMPACT_DIR}/UG/

./setup_compact_app.sh unitTest/IO -auto -3d +pm4dev withParticles=True +parallelio -portable
mv ${FLASH_ROOT}/IO ${COMPACT_DIR}/PM/

./setup_compact_app.sh unitTest/IO/IOTypes -auto -3d +pm4dev -unit=source/IO/IOMain/hdf5/parallel/PM_argonne withParticles=True parallelIO=True -portable
mv ${FLASH_ROOT}/IOTypes ${COMPACT_DIR}/PM/

# ********** Hydrodynamics
./setup_compact_app.sh Sedov -auto -2d +nofbs +parallelio -parfile=test_UG_nofbs_2d.par -portable
mv ${FLASH_ROOT}/Sedov ${COMPACT_DIR}/UG/Sedov_split

./setup_compact_app.sh Sedov -auto -2d +nofbs +parallelio +unsplitHydro -parfile=test_UG_nofbs_2d.par -portable
mv ${FLASH_ROOT}/Sedov ${COMPACT_DIR}/UG/Sedov_unsplit

./setup_compact_app.sh Sedov -auto -2d +pm4dev +parallelio -parfile=coldstart_pm.par -portable
mv ${FLASH_ROOT}/Sedov ${COMPACT_DIR}/PM/Sedov_split

./setup_compact_app.sh Sedov -auto -2d +pm4dev +parallelio +unsplitHydro -parfile=coldstart_pm.par -portable
mv ${FLASH_ROOT}/Sedov ${COMPACT_DIR}/PM/Sedov_unsplit

# ********** Diffusion
# (The MGDStep applications below depend on Hypre)
./setup_compact_app.sh ConductionDeltaSaDiff -2d -auto +nofbs -unit=source/physics/Diffuse/DiffuseMain/CG -parfile=flash_comp.par -portable
mv ${FLASH_ROOT}/ConductionDeltaSaDiff ${COMPACT_DIR}/UG/

./setup_compact_app.sh ConductionDeltaSaDiff -2d -auto +pm4dev -unit=source/physics/Diffuse/DiffuseMain/CG -parfile=flash_comp.par -portable
mv ${FLASH_ROOT}/ConductionDeltaSaDiff ${COMPACT_DIR}/PM/

./setup_compact_app.sh MGDStep -2d -auto +parallelIO +nofbs -parfile=coldstart_ug_2p_2d.par mgd_meshgroups=4 -portable
mv ${FLASH_ROOT}/MGDStep ${COMPACT_DIR}/UG/

./setup_compact_app.sh MGDStep -2d -auto +parallelIO +pm4dev -parfile=flash_2d_polar.par mgd_meshgroups=4 -portable
mv ${FLASH_ROOT}/MGDStep ${COMPACT_DIR}/PM/


echo -e "\nCreating compressed tar archive named ${COMPACT_APPS}.tar.gz"
cd ${FLASH_ROOT}
tar -cf "${COMPACT_APPS}.tar" "${COMPACT_APPS}/"
gzip -9 "${COMPACT_APPS}.tar"
rm -rf "${COMPACT_APPS}/"
