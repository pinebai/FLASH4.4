!!****if* source/Particles/ParticlesInitialization/pt_setMask
!!
!! NAME
!!
!!  pt_setMask
!!
!! SYNOPSIS
!!
!!  call pt_setMask()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***

#include "Flash.h"
#include "Particles.h"
#include "constants.h"
subroutine pt_setMask()
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Particles_data, ONLY : pt_gcMaskForAdvance, pt_typeInfo, pt_gcMaskForWrite, pt_meshVar,&
       pt_meshMe, pt_numAttributes, pt_attributes
  use Simulation_interface, ONLY : Simulation_mapStrToInt,Simulation_mapParticlesVar
  use pt_interface, ONLY : pt_mapStringParamToInt, pt_picInit

  implicit none
  logical :: massive
  character (len=MAX_STRING_LENGTH) :: partAttrPrefix,partAttr
  integer, save :: attributes
  integer :: i,j,k

  ! Guard cell masks, with loads of ifdefs.  This first one is for
  ! use in advancing particles, depending upon the method used
 

  pt_gcMaskForAdvance = .FALSE.
#ifdef DENS_VAR
  ! this is needed if velocities are obtained from the mesh because
  ! the prolongation may depend upon conversion to conservative etc
  pt_gcMaskForAdvance(DENS_VAR) = .TRUE.
#endif

#ifdef VELX_VAR
  pt_gcMaskForAdvance(VELX_VAR) = .TRUE.
#endif

#if NDIM > 1
#ifdef VELY_VAR
  pt_gcMaskForAdvance(VELY_VAR) = .TRUE.
#endif
#endif
  
#if NDIM > 2
#ifdef VELZ_VAR
  pt_gcMaskForAdvance(VELZ_VAR) = .TRUE.
#endif
#endif

#ifdef GPOT_VAR
  pt_gcMaskForAdvance(GPOT_VAR) = .FALSE.
  do i = 1,NPART_TYPES
     massive=(pt_typeInfo(PART_ADVMETHOD,i) == LEAPFROG).or.&
             (pt_typeInfo(PART_ADVMETHOD,i) == LEAPFROG_COSMO).or.&
             (pt_typeInfo(PART_ADVMETHOD,i) == EULER_MAS)
     pt_gcMaskForAdvance(GPOT_VAR) = pt_gcMaskForAdvance(GPOT_VAR).or.massive
  end do
#endif

  !! this section sets the mask for guardcell fill before particles write. Here all
  !! the variable corresponding to particle attributes that are going to be output
  !! need to be included in the fill

  pt_gcMaskForWrite = .FALSE.
  pt_numAttributes=0
  partAttrPrefix='particle_attribute_'

  do i = 1,PT_MAX_ATTRIBUTES

     call pt_mapStringParamToInt(attributes,partAttrPrefix,MAPBLOCK_PART,i)
     if(attributes>0) then
        pt_numAttributes=pt_numAttributes+1
        j=pt_numAttributes
        pt_attributes(j) = attributes
        call Simulation_mapParticlesVar(pt_attributes(j),&
             pt_meshVar(PT_VAR,j), pt_meshVar(PT_MAP,j))

#ifdef DEBUG_PARTICLES
        print*,'  attribute #',pt_numAttributes,' =',pt_attributes(j),'->meshVar',pt_meshVar(PT_VAR,j),pt_meshVar(PT_MAP,j)
#endif
        select case (pt_meshVar(PT_MAP,j))
        case(PARTICLEMAP_UNK)
           pt_gcMaskForWrite(pt_meshVar(PT_VAR,j))=.TRUE.
#if NSCRATCH_GRID_VARS > 0
        case(PARTICLEMAP_SCRATCH)
           ! DO NOTHING
#endif
#if NFACE_VARS > 0
        case(PARTICLEMAP_FACEX)
           pt_gcMaskForWrite(NUNK_VARS+pt_meshVar(PT_VAR,j))=.TRUE.
#if NDIM > 1
        case(PARTICLEMAP_FACEY)
           pt_gcMaskForWrite(NUNK_VARS+NFACE_VARS+pt_meshVar(PT_VAR,j))=.TRUE.
#if NDIM > 2
        case(PARTICLEMAP_FACEZ)
           pt_gcMaskForWrite(NUNK_VARS+NFACE_VARS+NFACE_VARS+pt_meshVar(PT_VAR,j))=.TRUE.
#endif
#endif
#endif
        case default
           if(pt_meshMe==MASTER_PE) then
              call concatStringWithInt(partAttrPrefix,i,partAttr) !for diagonstic msg only
              print*,"WARNING: ",trim(partAttr)," does not map to any mesh variable, will be ignored: PARTICLEMAP "
99            format(1x,'Particles_init: particle property',i2,' from variable',i2,', map',i2,', attr #',i2,'.')
              print 99, pt_attributes(j), pt_meshVar(PT_VAR,j), pt_meshVar(PT_MAP,j), j
           end if
           pt_numAttributes = pt_numAttributes - 1 ! ignore this one.
        end select
     else
        if(pt_meshMe==MASTER_PE)then
           if (attributes .NE. 0) then
              call concatStringWithInt(partAttrPrefix,i,partAttr) !for diagonstic msg only
              print*,"WARNING: ",trim(partAttr)," not recognized as a particle property, will be ignored"
           end if
        end if
     end if
  end do
#ifdef DEBUG_PARTICLES
  print*,'pt_gcMaskForAdvance:',pt_gcMaskForAdvance
  print*,'pt_gcMaskForWrite:  ',pt_gcMaskForWrite
#endif
  call Logfile_stampVarMask(pt_gcMaskForAdvance, .FALSE., '[pt_setMask]', 'gcMaskForAdvance')
  call Logfile_stampVarMask(pt_gcMaskForWrite, .FALSE., '[pt_setMask]', 'gcMaskForWrite')


end subroutine pt_setMask
