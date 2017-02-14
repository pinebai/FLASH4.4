"""Code to generate files, e.g. makefiles, Flash.h, mapVartoString and the like
some of them just calculate certain strings and use the template."""

import string
import time
import os
import re
import glob
import shutil
import sys
import itertools
import types

# local imports
# relative imports needed
from ..dsl.configuration import EOS_STATIC_LIST
from . import setup_globals
from .setup_globals import gvars, SetupError
from .utils import strictlyCaseSensitiveFilenames
from .template import Template
from .lib_union import LibUnion
from .lazy_file import LazyFile



def generate_flash_defines(config_info):
    """ Generate Flash.h using template """

    # sub routine
    def makeDefines(list, strname, num):
      # if dictionary then extract keys
      if type(list) == type({}): list = list.keys()
      list = map(string.upper,list)
      ans = []
      i = num
      for up_str in list:
        i += 1
        ans.append("#define %s_%s %d" % (up_str,strname,i))
      return string.join(ans,"\n")+"\n"



    fname = os.path.join(gvars.project_setup_dir, setup_globals.FLASH_DEFINES_FILENAME)
    tname = os.path.join(gvars.binDir, setup_globals.FLASH_DEFINES_TEMPLATE)
    gvars.out.put('generating %s' % setup_globals.FLASH_DEFINES_FILENAME,setup_globals.INFO)

    # compute necessary values
    tpl = Template(tname)

    # stuff from gvars
    tpl['ndim'] = gvars.dimension
    tpl['maxblocks'] = gvars.maxblocks
    tpl['nguard'] = config_info.guardcells
    if gvars.setup_vars.get("fixedBlockSize"):
       tpl['fixedblocksize'] = 1
    else: tpl['fixedblocksize'] = 0

    tpl['nxb'] = gvars.nxb
    tpl['nyb'] = gvars.nyb
    tpl['nzb'] = gvars.nzb

    tpl['gridGeometry'] = gvars.gridGeometry
    tpl['curvilinear']  = gvars.curvilinear

    tpl['strictParams'] = gvars.strictParams

    if gvars.setup_vars.get("npg"):
       tpl['npg'] = 1
    else: tpl['npg'] = 0

    # counts of various things
    if config_info.variable:
       tpl['nvars'] = len(config_info.variable)
    else: tpl['nvars'] = 1
    tpl['nspecies'] = config_info.nspecies
    tpl['nmass_scalars'] = config_info.nmassscalars
    tpl['nfacevars'] = len(config_info.facevar)
    tpl['nflux'] = len(config_info.flux)
    tpl['nrealprops'] = config_info.n_real_props
    tpl['nscratchvars'] = len(config_info.scratchvar)
    tpl['nscratchcentervars'] = len(config_info.scratchcentervar)
    tpl['nscratchfacexvars'] = len(config_info.scratchfacexvar)
    tpl['nscratchfaceyvars'] = len(config_info.scratchfaceyvar)
    tpl['nscratchfacezvars'] = len(config_info.scratchfacezvar)
    tpl['max_plot_vars'] = config_info.max_plot_vars

    # sets of # defines
    tpl["ppDefines"] = config_info.ppdefines_lines
    tpl["variableDefines"] = makeDefines(config_info.variable_names, "VAR", 0)
    tpl["speciesDefines"]  = makeDefines(config_info.species_names, "SPEC", tpl['nvars'])   
    tpl["massscalarsDefines"]  = makeDefines(config_info.massscalars, "MSCALAR", 
                                 tpl['nvars']+tpl['nspecies'])
    tpl["nmassscalar_groups"]  = len(config_info.mass_scalar_groups)
    tpl["facevarDefines"]  = makeDefines(config_info.facevar_names, "FACE_VAR", 0)
    tpl["fluxDefines"]     = makeDefines(config_info.flux_names, "FLUX", 0)
    tpl["scratchvarDefines"]  = makeDefines(config_info.scratchvar_names, "SCRATCH_GRID_VAR", 0)
    tpl["scratchcentervarDefines"]  = makeDefines(config_info.scratchcentervar_names, 
                                                  "SCRATCH_CENTER_VAR", 0)
    tpl["scratchfacexvarDefines"]  = makeDefines(config_info.scratchfacexvar_names, 
                                                 "SCRATCH_FACEX_VAR", 0)
    tpl["scratchfaceyvarDefines"]  = makeDefines(config_info.scratchfaceyvar_names, 
                                                 "SCRATCH_FACEY_VAR", 0)
    tpl["scratchfacezvarDefines"]  = makeDefines(config_info.scratchfacezvar_names, 
                                                 "SCRATCH_FACEZ_VAR", 0)

    if config_info.particleprop:
       tpl["partpropDefines"] = makeDefines(config_info.realproperty, "PART_PROP", 0)
    else: tpl["partpropDefines"] = ""

    tpl['nparticletypes'] = config_info.nparticletypes
    tpl["particletypeDefines"]  = makeDefines(config_info.particletype, "PART_TYPE", 0)

#  The use of mfrac_specDefine is a leftover from times when we used to enforce that there
#  was always one mass fraction variable in UNK (giving it the name MFRAC_SPEC if no
#  mass fractions were defined in any Config file), like this:
##     if len(config_info.species) == 0:
##         tpl["mfrac_specDefine"] = "#define MFRAC_SPEC %s" % str(tpl["nvars"]+1) 
    tpl["mfrac_specDefine"] = ""
    
    # non-replicated unk vars
    chain = itertools.chain
    def starts(xs):
        n = 0
        for x in xs:
            yield n; n += x
        yield n
    
    nonrep = config_info.nonrep
    tp2fmt = { 'MASS_SCALAR':'%s_MSCALAR', 'VARIABLE':'%s_VAR' }
    tpl["nonrepDefines"] = ''.join(chain([ \
        '#define %s_NONREP %d\n' % (name, i+1) + \
        '#define %s_NONREP_LOC2UNK(loc) (%s)\n' % (name, rec['nlocs'] > 0 and ('(loc)-1+'+(tp2fmt[rec['tp']] % (rec['locf'] % 1))) or '-1') + \
        '#define %s_NONREP_MAXLOCS %d\n' % (name, rec['nlocs']) + \
        '#define %s_NONREP_RPCOUNT "%s"\n' % (name, rec['rpcount']) \
        for i,(name,rec) in enumerate(nonrep.iteritems()) \
    ], [ \
        '#define NONREP_COUNT %d\n' % len(nonrep), \
        '#define NONREP_NAMEF_FLAT_LWR "%s"\n' % ''.join([ rec['namef'] for rec in nonrep.itervalues() ]).lower(), \
        '#define NONREP_NAMEF_START (/%s/)\n' % ','.join([ str(x+1) for x in starts([ len(rec['namef']) for rec in nonrep.itervalues() ]) ]), \
        '#define NONREP_MAXLOCS (/%s/)\n' % ','.join(['0'] + [ str(rec['nlocs']) for rec in nonrep.itervalues() ]), \
        '#define NONREP_LOCUNK1 (/%s/)\n' % ','.join(['0'] + [ rec['nlocs'] > 0 and (tp2fmt[rec['tp']] % (rec['locf'] % 1)) or '-1' for rec in nonrep.itervalues() ]), \
        '#define NONREP_RPCOUNT_FLAT "%s"\n' % ''.join([ rec['rpcount'] for rec in nonrep.itervalues() ]), \
        '#define NONREP_RPCOUNT_START (/%s/)\n' % ','.join([ str(x+1) for x in starts([ len(rec['rpcount']) for rec in nonrep.itervalues() ]) ]) \
    ]))
    
    # write the file if it really changed
    tpl.generate(fname)
    gvars.out.put('writing %s' % fname,setup_globals.INFO)


def generate_pm3rp_defines(config_info):
    """ Generate amr_runtime_parameters using template """

    fname = os.path.join(gvars.flash_src_dir, gvars.project_setup_dir, setup_globals.PM3_RP_FILENAME)
    # Use template in the object directory if it exists, otherwise the bin directory
    tname = os.path.join(gvars.flash_src_dir, gvars.project_setup_dir, setup_globals.PM3_RP_TEMPLATE)
    if not os.path.isfile(tname):
        tname = os.path.join(gvars.binDir, setup_globals.PM3_RP_TEMPLATE)
    gvars.out.put('generating %s' % setup_globals.PM3_RP_FILENAME,setup_globals.INFO)

    # compute necessary values
    tpl = Template(tname)

    # stuff from gvars
    tpl['ndim'] = gvars.dimension
    tpl['maxblocks'] = gvars.maxblocks
    tpl['nguard'] = config_info.guardcells
    if gvars.setup_vars.get("fixedBlockSize"):
       tpl['fixedblocksize'] = 1
    else: tpl['fixedblocksize'] = 0

    tpl['nxb'] = gvars.nxb
    tpl['nyb'] = gvars.nyb
    tpl['nzb'] = gvars.nzb

    if gvars.setup_vars.get("npg"):
       tpl['npg'] = 'T'
    else: tpl['npg'] = 'F'

    if gvars.curvilinear:
       tpl['curvilinear'] = 'T'
    else: 
       tpl['curvilinear'] = 'F'

    if (gvars.gridGeometry == setup_globals.GRID_GEOM_UNDEF):
       tpl['cartesian_eff'] = 'T'
       tpl['cartesian'] = 'F'
       tpl['noncartesian'] = 'F'
    elif (gvars.gridGeometry == setup_globals.GRID_GEOM_CARTESIAN):
       tpl['cartesian_eff'] = 'T'
       tpl['cartesian'] = 'T'
       tpl['noncartesian'] = 'F'
    else:
       tpl['cartesian_eff'] = 'F'
       tpl['cartesian'] = 'F'
       tpl['noncartesian'] = 'T'


    # counts of various things
    if config_info.variable_names:
       tpl['nvars'] = len(config_info.variable_names)
    else: tpl['nvars'] = 1
    tpl['nspecies'] = config_info.nspecies
    tpl['nmass_scalars'] = config_info.nmassscalars
    tpl['nfacevars'] = len(config_info.facevar)
    tpl['nflux'] = len(config_info.flux)

    # some defines derived from the previous ones
    tpl["nunk_vars"] = tpl["nvars"] + tpl["nspecies"] + tpl["nmass_scalars"]
    tpl["nfluxes"] = tpl["nflux"] + tpl["nspecies"] + tpl["nmass_scalars"]

    if config_info.ppdefines.has_key('NBOUNDARIES') and config_info.ppdefines['NBOUNDARIES']:
        tpl["nboundaries"] = config_info.ppdefines['NBOUNDARIES']
    else:
        tpl["nboundaries"] = 2*tpl["ndim"]

    # write the file if it really changed
    tpl.generate(fname)
    gvars.out.put('writing %s' % fname,setup_globals.INFO)


def write_simulation_files(config_info):
    """
    'config_info' is an instance of class FlashUnitUnion
    """
    ## This function generates all the Simulation_*.F90 files
    ## which setup needs to generate.

    ## mapIntToStr: maps a variable enumerated constant to string.
    ## mapStrToVar: reverse of above map 
    ## mapParticlesVar: returns the variable a particle property maps to
    ##
    ## Handles: UNK vars, Fluxes, Particle Properties
    ##
    ## These functions are useful with IO so we know
    ## which variables we are outputing in terms of a string
    ## rather than a number

    def handleblock(vars,keys,mapblock,dictkey):
        # vars and keys are modified in place
        x = 0
        tpl = "(%s * MAPBLOCKSIZE)+ %%2d" % mapblock
        for v in getattr(config_info, dictkey.lower()):
            x += 1
            vars.append(v.lower())
            keys.append(tpl % x)

    def handlemap (part_keys, var_keys, var_types, var_type, mapblock, dictkey):
        #creates the set of keys for mappings from the pairs in dictkey
        for v in getattr(config_info, dictkey.lower()):
            tmpkey = "%s_PART_PROP" % v[0]
            tmpvar = "%s_%s" %(v[1], var_type)
            part_keys.append(tmpkey.upper())
            var_keys.append(tmpvar.upper())
            var_types.append(mapblock.upper())


    #EOS functions:
    #---------------------------------------------------------------------------
    #Returns an EOS dictionary initialised to NONEXISTENT.
    def makeEosVaRDict():
        dEos = {}
        for eosKey in EOS_STATIC_LIST:
            dEos[eosKey.upper()]=(gvars.nonexistent,gvars.nonexistent)
        return dEos


    #Fill/update the EOS dictionary with values.  Can be called multiple
    #times to update the pre-existing values in the dictionary.
    def associateEosVars(eosVaRDict, zippedData, strSuffix):
        for (scratchVar, eosMapIn, eosMapOut) in zippedData:
            #scratchVarStr is the pre-processor definition.
            scratchVarStr = scratchVar.upper() + strSuffix

            i = eosMapIn.upper()
            if i in eosVaRDict:
                (dIn,dOut) = eosVaRDict[i]
                if dIn != scratchVarStr and dIn != "NONEXISTENT":
                    gvars.out.put('%s: Conflicting specifications for EOSMAPIN: %s - %s and %s,%s.' %
                                  ("WARNING",eosMapIn,dIn,scratchVarStr," ignoring the latter"),
                                  setup_globals.WARN)
                else:
                    eosVaRDict[i] = (scratchVarStr,dOut)

            i = eosMapOut.upper()
            if i in eosVaRDict:
                (dIn,dOut) = eosVaRDict[i]
                if dOut != scratchVarStr and dOut != "NONEXISTENT":
                    raise SetupError('Conflicting specifications for EOSMAPOUT: %s - %s and %s.' %
                                     (eosMapOut,dOut,scratchVarStr))
                eosVaRDict[i] = (dIn,scratchVarStr)


    #Returns a list of EOS Fortran code.
    def makeEosPrintableString(eosVaRDict, eosStr):
        tplList = []; tmp1 = []; tmp2 = []
        for (var, (eosIn, eosOut)) in eosVaRDict.items():
            fixedStr = eosStr + '(EOSMAP_' + var
            tmp1.append(fixedStr + ',EOSIN) = ' + eosIn)
            tmp2.append(fixedStr + ',EOSOUT) = ' + eosOut)
        tplList.append("!1. EOS input variables:")
        tplList.extend(tmp1)
        tplList.append("!2. EOS output variables:")
        tplList.extend(tmp2)
        return tplList
    #---------------------------------------------------------------------------


    gvars.out.put('generating Simulation mapping files',setup_globals.INFO)

    # All unk variables and their keys
    if config_info.variable_names:
       vars = [v.lower() for v in config_info.variable_names]
    else: vars = ["dummy"]
    vars.extend([v.lower() for v in config_info.species_names])
    vars.extend([v.lower() for v in config_info.massscalars])
    # compute their corresponding keys starting from 1
    keys = []
    for x in xrange(len(vars)):
        keys.append("(MAPBLOCK_UNK  * MAPBLOCKSIZE) + %2d" % (x+1))

    # add fluxes and their keys
    handleblock(vars,keys,"MAPBLOCK_FLUX","flux")
    # add particle properties and keys
    handleblock(vars,keys,"MAPBLOCK_PART","realproperty")
    # add scratch vars and their keys
    handleblock(vars,keys,"MAPBLOCK_SCRATCH","scratchvar")
    # add face vars and their keys
    handleblock(vars,keys,"MAPBLOCK_FACES","facevar")
    # add scratch center vars and their keys
    handleblock(vars,keys,"MAPBLOCK_SCRATCH_CENTER","scratchcentervar")
    # add scratch face vars and their keys
    handleblock(vars,keys,"MAPBLOCK_SCRATCH_FACEX","scratchfacexvar")
    handleblock(vars,keys,"MAPBLOCK_SCRATCH_FACEY","scratchfaceyvar")
    handleblock(vars,keys,"MAPBLOCK_SCRATCH_FACEZ","scratchfacezvar")



    part_keys = []
    var_keys = []
    var_types = []
    #add unk vars to particle maps
    handlemap(part_keys, var_keys, var_types, "VAR", "PARTICLEMAP_UNK", "particlemaps_variable")
    #add mass scalars to particle maps
    handlemap(part_keys, var_keys, var_types, "MSCALAR", "PARTICLEMAP_UNK", "particlemaps_mscalar")
    #add species to particle maps
    handlemap(part_keys, var_keys, var_types, "SPEC", "PARTICLEMAP_UNK", "particlemaps_species")
    #add scratchvars to particle maps
    handlemap(part_keys, var_keys, var_types, "SCRATCH_GRID_VAR", "PARTICLEMAP_SCRATCH", "particlemaps_scratchvar")
    #add x-facevars to particle maps
    handlemap(part_keys, var_keys, var_types, "FACE_VAR", "PARTICLEMAP_FACEX", "particlemaps_facex")
    #add y-facevars to particle maps
    handlemap(part_keys, var_keys, var_types, "FACE_VAR", "PARTICLEMAP_FACEY", "particlemaps_facey")
    #add z-facevars to particle maps
    handlemap(part_keys, var_keys, var_types, "FACE_VAR", "PARTICLEMAP_FACEZ", "particlemaps_facez")
    #add scratchcentervars to particle maps
    handlemap(part_keys, var_keys, var_types, "SCRATCH_CENTER_VAR", "PARTICLEMAP_SCRATCH_CTR", "particlemaps_scratchcentervar")
    #add scratchfacexvars to particle maps
    handlemap(part_keys, var_keys, var_types, "SCRATCH_FACEX_VAR", "PARTICLEMAP_SCRATCH_FACEX", "particlemaps_scratchfacexvar")
    #add scratchfaceyvars to particle maps
    handlemap(part_keys, var_keys, var_types, "SCRATCH_FACEY_VAR", "PARTICLEMAP_SCRATCH_FACEY", "particlemaps_scratchfaceyvar")
    #add scratchfacezvars to particle maps
    handlemap(part_keys, var_keys, var_types, "SCRATCH_FACEZ_VAR", "PARTICLEMAP_SCRATCH_FACEZ", "particlemaps_scratchfacezvar")

    # partition: given a list (u,v) returns a dictionary from u to a list of v
    # xs: iterator of pairs to partition
    def partition(xs):
        d = {}
        for u,v in xs:
            if d.has_key(u):
                d[u].append(v)
            else:
                d[u] = [ v ];
        return d
    chain = itertools.chain
    
    # generate the files
    tpl = Template(os.path.join(gvars.binDir,setup_globals.SIM_INT_TO_STR_TEMPLATE))
    tpl["values"] = vars
    tpl["keys"] = keys
    tpl["cases"] = [ 'case(%s); str="%s"' % (k,v) for k,v in zip(keys,vars) ]
    tpl.generate(os.path.join(gvars.flash_src_dir, gvars.project_setup_dir, setup_globals.SIM_INT_TO_STR_FILENAME))
    del tpl
    tpl = Template(os.path.join(gvars.binDir, setup_globals.SIM_STR_TO_INT_TEMPLATE))
    tpl["values"] = vars
    tpl["keys"] = keys
    tpl["select_key_from_strlwr_and_map"] = map( \
        lambda tab_txt: '    '*tab_txt[0] + tab_txt[1], \
        chain( \
            [ (0,"select case(strlwr)") ], \
            chain(*map( \
                (lambda v_ks: (lambda v,ks: chain( \
                    [ (0,'case("%s")' % v) ], \
                        [ (1,'select case(map)') ], \
                        [ (1,'case((%s)/MAPBLOCKSIZE); key = %s' % (k,k)) for k in ks ], \
                        [ (1,'end select') ] \
                ))(v_ks[0],v_ks[1])), \
                partition(zip(vars, keys)).iteritems() \
            )), \
            [ (0,'end select') ] \
        ) \
    )
    tpl.generate(os.path.join(gvars.flash_src_dir, gvars.project_setup_dir, setup_globals.SIM_STR_TO_INT_FILENAME))
    del tpl
    tpl = Template(os.path.join(gvars.binDir,setup_globals.SIM_PARTICLES_VAR_TEMPLATE))
    tpl["partkeys"] = part_keys
    tpl["varkeys"] = var_keys
    tpl["vartypes"] = var_types
    tpl.generate(os.path.join(gvars.flash_src_dir, gvars.project_setup_dir, setup_globals.SIM_PARTICLES_VAR_FILENAME))

    # now for renormGroup code
    tpl = Template(os.path.join(gvars.binDir,setup_globals.RENORM_GROUP_TEMPLATE))
    it = config_info.massscalars_map.items()
    tpl["mscalars"] = [x.upper()+"_MSCALAR" for (x,y) in it]
    tpl["groups"] = [y for (x,y) in it]
    # ideally the group_names formatting should be done in the template
    # but I did not build a way to escape the ! in the list specification
    # So this is a hack around the template scheme of things
    tpl["group_names"] = " \n   !! ".join(["%s -> %d" % (x,y) for (x,y) in \
                                           config_info.massscalars_group_map.items()])
    #tpl.generate(os.path.join(gvars.flash_src_dir,gvars.project_build_dir,setup_globals.RENORM_GROUP_FILENAME))
    tpl.generate(os.path.join(gvars.project_setup_dir, setup_globals.RENORM_GROUP_FILENAME))

    # now for variable type code
    tpl = Template(os.path.join(gvars.binDir,setup_globals.VARNAME_TYPE_TEMPLATE))
    tpl["varnames"] = [ "%s_VAR" % x.upper() for x in config_info.variable_names ]
    tpl["vartypes"] = [ "VARTYPE_%s" % x for x in config_info.var_types ]
    #tpl.generate(os.path.join(gvars.flash_src_dir,gvars.project_build_dir,setup_globals.VARNAME_TYPE_FILENAME))
    tpl.generate(os.path.join(gvars.project_setup_dir, setup_globals.VARNAME_TYPE_FILENAME))


    #We are adding a file which contains information about all particle types.
    #If we have not included the particles unit then we will have a value of 0
    #for nparticletypes - in this case we do not create the file!
    if (config_info.nparticletypes > 0):
        tpl = Template(os.path.join(gvars.binDir,setup_globals.SIM_PARTICLE_TYPE_TEMPLATE))

        #We can slide in alternative particle methods by honouring the
        #particlesmethod command.
        for particleType in config_info.particle_type:
            if gvars.particleMethods.has_key(particleType):
                config_infoIndex = config_info.particle_type.index(particleType)
                overwriteData = gvars.particleMethods.get(particleType)

                print "NOTE: particlemethods option given on command line!"
                print "  Editing particleType:", particleType
                for (name,value) in overwriteData:
                    if (name == "INIT"):
                        print "  Init method was:", config_info.init_method[config_infoIndex], \
                              "now:", value
                        config_info.init_method[config_infoIndex] = value
                    elif (name == "MAP"):
                        print "  Map method was:", config_info.mapmethod[config_infoIndex], \
                              "now:", value
                        config_info.mapmethod[config_infoIndex] = value
                    elif (name == "ADV"):
                        print "  adv method was:", config_info.advmethod[config_infoIndex], \
                              "now:", value
                        config_info.advmethod[config_infoIndex] = value
                    else:
                        print "  Ignoring unknown option:", name
        
        tpl["particleInitVect"] = config_info.init_method
        tpl["particleMapVect"] = config_info.mapmethod
        tpl["particleAdvVect"] = config_info.advmethod
        #tpl.generate(os.path.join(gvars.flash_src_dir, gvars.project_build_dir, setup_globals.SIM_PARTICLE_TYPE_FILENAME))
        tpl.generate(os.path.join(gvars.flash_src_dir, gvars.project_setup_dir, setup_globals.SIM_PARTICLE_TYPE_FILENAME))
    


    #Add a file to select Eos_map roles from the various FLASH data structures.
    tpl = Template(os.path.join(gvars.binDir,setup_globals.EOS_MAP_TEMPLATE))

    #Create a dictionary for each of the FLASH data structures.
    dEosMap = {}
    dEosMap["_VAR"] = makeEosVaRDict()
    dEosMap["_FACE_VAR"] = makeEosVaRDict()
    dEosMap["_SCRATCH_GRID_VAR"] = makeEosVaRDict()
    dEosMap["_SCRATCH_CENTER_VAR"] = makeEosVaRDict()
    dEosMap["_SCRATCH_FACEX_VAR"] = makeEosVaRDict()
    dEosMap["_SCRATCH_FACEY_VAR"] = makeEosVaRDict()
    dEosMap["_SCRATCH_FACEZ_VAR"] = makeEosVaRDict()

    #EOS variables are associated with the following KEYWORDS.
    unkType_VAR = zip(config_info.variable_names, 
                      config_info.eosmapin_unkvars,
                      config_info.eosmapout_unkvars)
    unkType_MSCALAR = zip(config_info.massscalars, 
                          config_info.eosmapin_ms,
                          config_info.eosmapout_ms)
    unkType_FACE = zip(config_info.facevar_names, 
                       config_info.eosmapin_facevars,
                       config_info.eosmapout_facevars)
    unkType_SCRATCH = zip(config_info.scratchvar_names, 
                          config_info.eosmapin_scratchvars,
                          config_info.eosmapout_scratchvars)
    unkType_SCRATCHCENTER = zip(config_info.scratchcentervar_names, 
                          config_info.eosmapin_scratchcentervars,
                          config_info.eosmapout_scratchcentervars)
    unkType_SCRATCHFACEX = zip(config_info.scratchfacexvar_names, 
                          config_info.eosmapin_scratchfacexvars,
                          config_info.eosmapout_scratchfacexvars)
    unkType_SCRATCHFACEY = zip(config_info.scratchfaceyvar_names, 
                          config_info.eosmapin_scratchfaceyvars,
                          config_info.eosmapout_scratchfaceyvars)
    unkType_SCRATCHFACEZ = zip(config_info.scratchfacezvar_names, 
                          config_info.eosmapin_scratchfacezvars,
                          config_info.eosmapout_scratchfacezvars)

    #Use KEYWORD data to fill the values in the appropriate dictionary.
    associateEosVars(dEosMap["_VAR"], 
                     unkType_VAR, "_VAR")
    associateEosVars(dEosMap["_VAR"], 
                     unkType_MSCALAR, "_MSCALAR")
    associateEosVars(dEosMap["_FACE_VAR"], 
                     unkType_FACE, "_FACE_VAR")
    associateEosVars(dEosMap["_SCRATCH_GRID_VAR"], 
                     unkType_SCRATCH, "_SCRATCH_GRID_VAR")
    associateEosVars(dEosMap["_SCRATCH_CENTER_VAR"], 
                     unkType_SCRATCHCENTER, "_SCRATCH_CENTER_VAR")
    associateEosVars(dEosMap["_SCRATCH_FACEX_VAR"], 
                     unkType_SCRATCHFACEX, "_SCRATCH_FACEX_VAR")
    associateEosVars(dEosMap["_SCRATCH_FACEY_VAR"], 
                     unkType_SCRATCHFACEY, "_SCRATCH_FACEY_VAR")
    associateEosVars(dEosMap["_SCRATCH_FACEZ_VAR"], 
                     unkType_SCRATCHFACEZ, "_SCRATCH_FACEZ_VAR")

    #Generate Fortran code to populate data structures named
    #eos_unk, eos_face, eos_scratch in eos_variableMap.F90:
    tpl["eos_unk"] = makeEosPrintableString(dEosMap["_VAR"],
                                            "eosmap_unk")
    tpl["eos_face"] = makeEosPrintableString(dEosMap["_FACE_VAR"],
                                             "eosmap_face")
    tpl["eos_scratch"] = makeEosPrintableString(dEosMap["_SCRATCH_GRID_VAR"],
                                                "eosmap_scratch")
    tpl["eos_scratch_ctr"] = makeEosPrintableString(dEosMap["_SCRATCH_CENTER_VAR"],
                                                "eosmap_scratch_ctr")
    tpl["eos_scratch_facexvar"] = makeEosPrintableString(dEosMap["_SCRATCH_FACEX_VAR"],
                                                "eosmap_scratch_facexvar")
    tpl["eos_scratch_faceyvar"] = makeEosPrintableString(dEosMap["_SCRATCH_FACEY_VAR"],
                                                "eosmap_scratch_faceyvar")
    tpl["eos_scratch_facezvar"] = makeEosPrintableString(dEosMap["_SCRATCH_FACEZ_VAR"],
                                                "eosmap_scratch_facezvar")

    #Get length of an EOS list.  NOTE: All 3 lists have the same length.
    #We subtract 2 because our lists contains 2 printed comments.
    #We divide by 2 because the single list contains IN and OUT.
    tpl["len_eos_lists"] = (len(tpl["eos_unk"])-2)/2
    #tpl.generate(os.path.join(gvars.flash_src_dir,gvars.project_build_dir,setup_globals.EOS_MAP_FILENAME))
    tpl.generate(os.path.join(gvars.project_setup_dir, setup_globals.EOS_MAP_FILENAME))
    del tpl



def generate_makefile(config_info, machdir):
    flagRedirect = setRedirectFlags(os.path.join(gvars.project_setup_dir, 'Makefile.h'), 
                                    gvars.buildFlag, config_info.lib_config_info)
    taumakefile= ''
    pdbvar = ''
    pdbrule = ''
    taumakeline = ''
    tauinserts = ''
    
    cDefines = " ".join(gvars.defines)
    fDefines = ' '.join(["$(MDEFS)" + item for item in gvars.defines])

    setupVaRDict = gvars.setup_vars.vars
    threadKeysFound = [x for x in setupVaRDict.keys()
                       if x in ["threadBlockList","threadWithinBlock","threadRayTrace"]]
    useopenmp = 0
    for key in threadKeysFound:
        #Set USEOPENMP = 1 in Makefile if at least one thread setup variable is True.
        if type(setupVaRDict[key]) == types.BooleanType and setupVaRDict[key]:
            useopenmp = 1
            break

    if gvars.build_tau:
        shutil.copyfile(os.path.join(gvars.tauDir,'select.tau'),
                        os.path.join(gvars.project_setup_dir,'select.tau'))
                        #os.path.join(gvars.flash_src_dir,gvars.project_build_dir,'select.tau'))
        taumakefile = gvars.build_tau
        taumakeline="TAU_MAKEFILE="+gvars.build_tau
        tauinserts = """
TAU_OPTIONS='-optVerbose -optTauSelectFile=./select.tau -optPdtGnuFortranParser'

FCOMP=tau_f90.sh -tau_makefile=$(TAU_MAKEFILE) -tau_options=$(TAU_OPTIONS)
LINK=tau_f90.sh -tau_makefile=$(TAU_MAKEFILE) -tau_options=$(TAU_OPTIONS)
CCOMP=tau_cc.sh -tau_makefile=$(TAU_MAKEFILE) -tau_options=$(TAU_OPTIONS)
CPPCOMP=tau_cxx.sh -tau_makefile=$(TAU_MAKEFILE) -tau_options=$(TAU_OPTIONS)"""
                                                    
    makefiles = glob.glob(os.path.join(gvars.project_setup_dir, 'Makefile.*'))
    makefiles = [os.path.relpath(mf, gvars.project_setup_dir) for mf in makefiles]

    #FIXME remove this once all those makefile dependencies are fixed.
    makefiles.sort()
    ", ".join(makefiles)

    try:
       makefiles.remove('Makefile.h')
    except ValueError:
       pass  # if not found dont bother 
    includeList = "".join(["include %s\n"%file for file in makefiles])
    includeMacros = "".join(["\n       $(%s) \\" % os.path.splitext(file)[1][1:] for file in makefiles])

    makedisplay = gvars.makedisplay
    indexReorder = 0
    if gvars.indexReorder:
        indexReorder = 1
    if strictlyCaseSensitiveFilenames():
        dependFlags = '--generateINTERMEDIATElines'
    else:
        dependFlags = ''
    tpl = Template(os.path.join(gvars.binDir, setup_globals.MAKEFILE_TEMPLATE))
    tpl.update(locals())

    target = os.path.join(gvars.project_setup_dir, 'Makefile')
    tpl.generate(target)

###########################

def generate_buildstamp_generator():
    #OUTFILE = os.path.join(gvars.flash_src_dir,gvars.project_build_dir,setup_globals.BUILD_STAMP_GEN_FILENAME)
    OUTFILE = os.path.join(gvars.project_setup_dir, setup_globals.BUILD_STAMP_GEN_FILENAME)
    tname = os.path.join(gvars.binDir, setup_globals.BUILD_STAMP_TEMPLATE)

    tpl = Template(tname)
    tpl["date"] = time.asctime(time.localtime(time.time()))
    tpl["uname"] = os.uname()
    tpl.generate(OUTFILE)
    os.chmod(OUTFILE, 0744)

################################################# Internally called code

def setRedirectFlags(makefile, buildFlag, lib_config_info):
    """FIXME: proper documentation
    
    If the Makefile.h uses the mechanism but hasn't defined it for the given
    buildFlag (i.e., it has FFLAGS_OPT but not FFLAGS_TEST) we default to _OPT
    for whatever compilers are missing the right flag. If even that is missing
    we look for just plain FFLAGS (without any _OPT)

    Also, we have added support for "internal" libraries. If a LIBRARY
    requirement is not found in Makefile.h, then the directory
    flash_src_dir/lib/name/object is searched, where name is the name of
    the LIBRARY as specified in the Config file. If that directory exists and
    contains a file named libname.a, then Makefile has that info added to
    its LIB macro. Otherwise, setup attempts to execute a file called build.csh
    in a directory flash_src_dir/lib/name/source. This file contains commands
    for building the library and placing the library libname.a in the
    lib/name/object dir.

    Some further things to consider are adding support for multiple libraries
    within a lib directory and specifying any library dependencies.
    """
    #import pdb; pdb.set_trace()

    #Find all macros defined in makefile
    if not os.path.isfile(makefile): 
       text = ''
       makefilename = ""
    else: 
       text = open(makefile).read()
       makefilename = os.path.abspath(makefile)
    ro = re.compile("^ *([A-Za-z0-9_]+) *=",re.M)
    Makefile_macros = ro.findall(text)

    lib_config_info.setMacros(Makefile_macros) # inform our libConfig class about the macros

    newDef = {}    
    for compiler in setup_globals.COMPILERS:
        newDef[compiler] = ['$(%s)'%compiler]
        for macro in [ '%s_%s'%(compiler,buildFlag), 
                       '%s_%s'%(compiler,setup_globals.DEFLTFLAG)]:
            if macro in Makefile_macros:
                newDef[compiler] = ['$(%s)'%macro]
                break

    # process the libraries 
    for lib, args in lib_config_info.libOrder:
        gvars.out.push()
        libFlags = lib_config_info.getLibFlags(lib, buildFlag=buildFlag,
                                             args=args, makefilename=makefilename)
        if not libFlags or not libFlags.has_key("LIB"):
            gvars.out.put('ERROR: A Config in your simulation requires the %s library\n'\
                          'but I cannot find any info about it. ' % lib,setup_globals.ERROR)
            gvars.out.put('If you automatically link in that library, create a variable '\
                          'LIB_%s in your Makefile.h and make it empty' % lib.upper(),setup_globals.ERROR)
            raise SetupError("Error getting info for library %s" % lib)
        else:
            # the fact that we prepend (not append) and the specific order of the libraries
            # should ensure that the ordering here is correct
            for key,val in libFlags.items(): 
                if newDef.has_key(key): newDef[key].insert(1,val)
        gvars.out.pop()

    outList = []
    for compiler in newDef.keys():
        if (len(newDef[compiler])==1) and \
           (newDef[compiler][0]=='$(%s)'%compiler): continue
        macros = newDef[compiler]
        #want "main" macros like LIB_OPT to be at the end of the line
        macros = macros[1:]+[macros[0]]
        outList.append('%s := %s\n' % (compiler, string.join(macros)))

    # Make sure that this is a str, not unicode!
    redirect_flags = str("".join(outList))
    return redirect_flags
