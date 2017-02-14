#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "AMRLevelFlashFactory.h"

AMRLevelFlashFactory::AMRLevelFlashFactory()
{
}

AMRLevelFlashFactory::~AMRLevelFlashFactory()
{
}

void AMRLevelFlashFactory::define(const flash_amr_info_t& flashAMRInfo,
				  const mesh_info_t& meshInfo)
{
  m_flashAMRInfo = flashAMRInfo;
  m_meshInfo = meshInfo;
}

AMRLevel* AMRLevelFlashFactory::new_amrlevel() const
{
  AMRLevelFlash* amrLevelFlashPtr = new AMRLevelFlash(-1);
  amrLevelFlashPtr->defineParams(m_flashAMRInfo, m_meshInfo,false);
  return (static_cast <AMRLevel*> (amrLevelFlashPtr));
}
