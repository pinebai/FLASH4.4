#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.

5.2.4     Class AMRLevelFactory
The class AMRLevelFactory is a pure virtual base class, with only one
member function:

   • virtual AMRLevel* new_amrlevel() const = 0; This is the only
     member function of AMRLevelFactory, and it must be defined by the
     user in a derived class. The derived function will return a
     pointer to a physics-specific class derived from AMRLevel.

A pointer to an object of this class is passed to the define function
of an AMR object, which uses it to construct the various AMRLevel
objects that it requires.
*/
#endif

#ifndef _AMR_LEVEL_FLASH_FACTORY_H_
#define _AMR_LEVEL_FLASH_FACTORY_H_

#include "AMRLevelFactory.H"
#include "AMRLevelFlash.h"

class AMRLevelFlashFactory : public AMRLevelFactory
{

public:
  AMRLevelFlashFactory();

  virtual ~AMRLevelFlashFactory();

  virtual void define(const flash_amr_info_t& flashAMRInfo,
		      const mesh_info_t& meshInfo);

  virtual AMRLevel* new_amrlevel() const;

private:
  flash_amr_info_t m_flashAMRInfo;
  mesh_info_t m_meshInfo;


  void operator=(const AMRLevelFlashFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  AMRLevelFlashFactory(const AMRLevelFlashFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
