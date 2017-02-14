"""Bunch of small functions doing for useful stuff for flmake"""
import re
import sys
import os
import string
import glob
import socket
import re
import time
import shutil
import subprocess
import tempfile

# relative imports needed
from .. import _metadata, FLASH_SRC_DIR, FLASH_CLEAN_SRC_DIR
from ..utils import strip_comments, warning
from . import setup_globals
from .setup_globals import gvars, SetupError
from .lazy_file import LazyFile


TEMPDIR = tempfile.gettempdir()


#
# version contol functions
#


def git_info(d):
    """Return the git version information for a directory."""
    branch_cmd = ['git', 'branch', '--no-color']
    branch_out = subprocess.check_output(branch_cmd, cwd=d)
    branch = re.search('\* (.*?)\s*\n', branch_out).group(1)
    hash_cmd = ['git', 'rev-list', 'HEAD', '-1']
    hash_out = subprocess.check_output(hash_cmd, cwd=d)
    hash = hash_out.strip()
    return 'git', branch, hash


def hg_info(d):
    """Return the mecurial version information for a directory."""
    id_cmd = ['hg', 'id', '-i', '-b']
    id_out = subprocess.check_output(id_cmd, cwd=d)
    id, branch = id_out.split()
    return 'hg', branch, id
    

def svn_info(d):
    """Return the subversion information for a directory."""
    info_cmd = ['svn', 'info']
    info_out = subprocess.check_output(info_cmd, cwd=d)
    branch = re.search('Repository Root: (.*?)\s*\n', info_out).group(1)
    branch = os.path.split(branch)[-1]
    rev = re.search('Revision: (.*?)\s*\n', info_out).group(1)
    return 'svn', branch, rev


def release_info(d):
    """Returns the version information for a stable FLASH release."""
    with open(os.path.join(d, 'RELEASE')) as f:
        ver = f.read().strip()
    if os.path.exists(FLASH_CLEAN_SRC_DIR):
        rtn = 'release', '', ver
    else:
        rtn = tuple()
        print warning("Flash release used, but clean copy not found "
                      "at FLASH_CLEAN_SRC_DIR:"  + FLASH_CLEAN_SRC_DIR)
    return rtn


def vc_info(d):
    """Gets the version control information of a directory.  Returns either an 
    empty tuple when no repo is found or a 3-tuple of (repo-type, branch, version)
    when it is."""
    ls = [f for f in sorted(os.listdir(d)) if f in ['.git', '.hg', '.svn', 'RELEASE']]
    vc_type = None if 0 == len(ls) else ls[0]
    vc_switch = {None: lambda x: tuple(),
                 '.git': git_info,
                 '.hg': hg_info,
                 '.svn': svn_info,
                 'RELEASE': release_info,
                 }
    if hasattr(subprocess, 'check_output'):
        try:
            info = vc_switch[vc_type](d)
        except OSError:
            print warning("{0} data found, but {0} info failed".format(vc_type))
            info = tuple()
    else:
        print warning("{0} data found, but {0} info failed".format(vc_type))
        info = tuple()
    return info


def release_diff(src, clean_src):
    """Diffs a release directory with the clean version of the release dir."""
    clean_relpath =  os.path.relpath(clean_src, src)
    visible_files = [f for f in os.listdir(src) if not f.startswith('.')]
    diffs = ""
    for f in visible_files:
        cmd = ['diff', '-rupN', os.path.join(clean_relpath, f), f]
        diffs += subprocess.check_output(cmd, cwd=src) + "\n"
    return diffs


def vc_diff(d, vc_type=None):
    """Makes diff of the version controlled directory."""
    if vc_type is None:
        ls = [f for f in sorted(os.listdir(d)) if f in ['.git', '.hg', '.svn', 'RELEASE']]
        vc_type = None if 0 == len(ls) else ls[0].replace('.', '').lower()
    vc_func = lambda p, v: subprocess.check_output([v, 'diff'], cwd=p)
    vc_switch = {None: lambda p, v: "",
                 'git': vc_func,
                 'hg': vc_func,
                 'svn': vc_func,
                 'release': release_diff,
                 }

    # change results for release
    if vc_type == 'release':
        v = FLASH_CLEAN_SRC_DIR
        if not os.path.isdir(FLASH_CLEAN_SRC_DIR):
            vc_type = None
    else:
        v = vc_type

    if hasattr(subprocess, 'check_output'):
        try:
            diff = vc_switch[vc_type](d, v)
        except OSError:
            print warning("{0} data found, but {0} diff failed".format(vc_type))
            diff = ""
    else:
        print warning("{0} data found, but {0} diff failed".format(vc_type))
        diff = ""
    return diff


def vc_patch(d, diff, vc_type=None):
    """Applies a patch to a version controlled directory."""
    if vc_type is None:
        ls = [f for f in sorted(os.listdir(d)) if f in ['.git', '.hg', '.svn', 'RELEASE']]
        vc_type = None if 0 == len(ls) else ls[0].replace('.', '').lower()
    patchfile = os.path.join(TEMPDIR, "patchfile")
    with open(patchfile, 'w') as f:
        f.write(diff)
    patch_level = {'git': '-p1', 
                   'hg': '-p1',
                   'svn': '-p0',
                   'release': '',
                   }
    if not hasattr(subprocess, 'check_output'):
        return
    try:
        subprocess.check_output(" ".join(['patch', patch_level[vc_type], '<', patchfile]), 
                                cwd=d, shell=True)
    except OSError:
        print warning("patching {0} failed".format(d))
    finally:
        os.remove(patchfile)

        

def git_checkout(src, info, base=TEMPDIR, cmd=None):
    """Checks out out a copy of the git repo at info point in history in the 
    base directory."""
    bn = "{0}-{1}{2}".format(os.path.basename(src),
                             ("" if cmd is None else cmd + '-'),
                             "_".join(info))
    targ = os.path.join(base, bn)
    if os.path.exists(targ):
        return targ
    subprocess.check_output(['git', 'clone', src, targ])
    subprocess.check_output(['git', 'checkout', info[1]], cwd=targ)
    subprocess.check_output(['git', 'reset', '--hard', info[2]], cwd=targ)
    return targ


def hg_checkout(src, info, base=TEMPDIR, cmd=None):
    """Checks out out a copy of the hg repo at info point in history in the 
    base directory."""
    bn = "{0}-{1}{2}".format(os.path.basename(src),
                             ("" if cmd is None else cmd + '-'),
                             "_".join(info))
    targ = os.path.join(base, bn)
    if os.path.exists(targ):
        return targ
    subprocess.check_output(['hg', 'clone', src, targ])
    subprocess.check_output(['hg', 'update', '-C', info[1]], cwd=targ)
    subprocess.check_output(['hg', 'update', '-r', info[2]], cwd=targ)
    return targ


def svn_checkout(src, info, base=TEMPDIR, cmd=None):
    """Checks out out a copy of the svn repo at info point in history in the 
    base directory."""
    bn = "{0}-{1}{2}".format(os.path.basename(src),
                             ("" if cmd is None else cmd + '-'),
                             "_".join(info))
    targ = os.path.join(base, bn)
    if os.path.exists(targ):
        return targ
    src_info = svn_info(src)
    if info[1] != src_info[1]:
        msg = "branch {0} of svn repo {1} does not match {2}"
        msg = msg.format(src_info[1], src, info[1])
        raise RuntimeError(msg)
    shutil.copytree(src, targ)
    subprocess.check_output(['svn', 'revert', '-R', '.'], cwd=targ)
    revert_command = ['svn', 'merge', '-r{0}:{1}'.format(src_info[1], info[1]), info[2], '.']
    subprocess.check_output(revert_command, cwd=targ)
    return targ


def release_checkout(src, info, base=TEMPDIR, cmd=None):
    """Checks out out a copy of the  repo at info point in history in the 
    base directory."""
    bn = ("" if cmd is None else cmd + '-') + info[0] + '_' + info[2] 
    targ = os.path.join(base, bn)
    if os.path.exists(targ):
        return targ
    shutil.copytree(src, targ)    
    return targ


def vc_checkout(src, info, base=TEMPDIR, cmd=None):
    """Checks out out a copy of the source repo at info point in history in the 
    base directory."""
    vc_switch = {'git': git_checkout,
                 'hg': hg_checkout,
                 'svn': svn_checkout,
                 'release': release_checkout,
                 }
    src = FLASH_CLEAN_SRC_DIR if info[0] == 'release' else src
    if hasattr(subprocess, 'check_output'):
        try:
            targ = vc_switch[info[0]](src, info, base, cmd)
        except OSError:
            print warning("failed to checkout {0} via {1} at {2} {3}.".format(src, *info))
            targ = ""
    else:
        print warning("failed to checkout {0} via {1} at {2} {3}.".format(src, *info))
        targ = ""
    return targ



#
# Other utilites
#

def search_paths(paths, file_or_dir):
    """Searches the paths for a file (or directory).

    Parameters
    ----------
    paths : list of str
        Paths to search.
    file_or_dir : str
        File or directory name to find.

    Returns
    -------
    fullpath : str
        Full path of file or directory, '/path/to/file'."""
    fullpath = None
    for path in paths:
        tmpfullpath = os.path.join(path, file_or_dir)

        # Ensure valid path
        if not os.path.isdir(path):
            continue

        # Search path
        if os.path.exists(tmpfullpath):
            fullpath = tmpfullpath
            break
        else:
            subpaths = [os.path.abspath(os.path.join(path, f)) for f in os.listdir(path)]
            subsearch = search_paths(subpaths, file_or_dir)
            if subsearch is not None:
                fullpath = subsearch
                break

    return fullpath


def desc_cmd_metadata(project_dir='.'):
    """Gathers desc metadata common to all flmake commands."""
    desc_meta = {}
    desc_meta['env'] = dict(os.environ)
    desc_meta['metadata'] = _metadata
    desc_meta['timestamp'] = time.time()
    desc_meta['command'] = sys.argv
    desc_meta['source_version'] = vc_info(FLASH_SRC_DIR)
    desc_meta['source_diff'] = vc_diff(FLASH_SRC_DIR)
    desc_meta['project_version'] = vc_info(project_dir)
    desc_meta['project_diff'] = vc_diff(project_dir)
    rcfilename = os.path.join(project_dir, 'flashrc.py')
    if os.path.isfile(rcfilename):
        with open(rcfilename) as f:
            rcfile = f.read()
    else:
        rcfile = ""
    desc_meta['rcfile'] = rcfile
    return desc_meta


def hash_list_to_dict(hashes, d=None):
    """Converts a list of hashes to a recursive dictionary."""
    if d is None:
        d = {}

    curr_d = d
    for h in hashes:
        if h not in curr_d:
            curr_d[h] = {}
        curr_d = curr_d[h]

    return d


def hash_dict_to_str(hash_dict, labels=None, indent_str=""):
    """Converts a hash dictionary to a str.  Nice for printing."""
    if labels is None:
        labels = {}

    s = ""
    hd_len = len(hash_dict)

    if 0 < hd_len:
        for i, h in enumerate(sorted(hash_dict.keys())):
            next_indent_str = indent_str + ("| "  if 1 < hd_len and i != hd_len - 1 else "  ")
            s += indent_str 
            s += "+-" + ("{0} ({1})".format(h, labels[h]) if h in labels else h)
            s += "\n" + hash_dict_to_str(hash_dict[h], labels, next_indent_str)

    return s
    


def get_rel_path(filename, basedir):
    """Return the relative path to FILENAME from basedir"""
    sep = os.sep
    srcdir = os.path.abspath(basedir) + sep
    tgtdir = os.path.abspath(os.path.dirname(filename)) + sep

    # most common prefix (will end with "/" or contain extra chars
    cp = os.path.dirname(os.path.commonprefix([srcdir, tgtdir]))

    # Handles the special case when common prefix is "/"
    if cp[-1] != sep: 
        cp = cp + sep
    src = srcdir[len(cp):]
    tgt = tgtdir[len(cp):]

    # src and tgt contains name relative to cp
    c = src.count(sep) # how many levels up to reach common dir
    if c == 0:
      prefix = "." + sep
    else: 
        prefix = c*(".." + sep)
    return os.path.join(prefix + tgt, os.path.basename(filename))


def dir_glob(pathname):
    """Takes a pattern (absolute or relative to current directory)
    and returns a list of directories matching pattern. The match is 
    made case-insensitive"""
    # maps a -> [aA] but "/" -> "/", "1" -> "1" 
    mapfn = lambda x: (x.upper() != x.lower() and "[{0}{1}]".format(x.lower(), x.upper())) or x
    globstr = "".join(map(mapfn, pathname)) # concatenate

    # find files which match pathname (except for case)
    files = glob.glob(globstr)

    # return only those of which are directories
    return [name for name in files if os.path.isdir(name)]  


def get_os_type(prototypes_dir):
    ostype = sys.platform.lower()
    if '-' in ostype:
        ostype = ostype[:ostype.find('-')]

    for proto in os.listdir(prototypes_dir):
        if ostype.count(proto.lower()):
            return proto

    return ostype


def getHostName(sitesDir):
    """Returns the hostname to use."""

    #Change made by sam 09/01/2011. A bug occurs when trying to get the 
    #host by address for machines without dns entries in a name server 
    #or on a machine without a connection to a dns server
    
    tempHostName = socket.gethostname()

    try:	
    	temp = socket.gethostbyaddr(tempHostName)
    except:
	temp = (socket.getfqdn(),[])
        
	
    fallback = temp[0]
    namesToTry = [temp[0]]
    namesToTry.append(socket.gethostname())
    namesToTry.extend(temp[1]) # list of addl names for the current host

    # Read the alias file into memory
    aliasLines = []
    try:
       aliasFile = open(os.path.join(sitesDir,'Aliases'))
       gvars.out.put('checking sites Aliases file', setup_globals.IMPINFO)
       for line in aliasFile.readlines():
           line = strip_comments(line, '#','"')
           line = string.strip(line)
           if line:
              parts = string.split(line)
              if len(parts) <> 2:
                 gvars.out.put("Ignoring bad Aliases file line '%s'" % string.strip(line),globals.WARN)
              else: aliasLines.append(parts)
       aliasFile.close()
    except IOError:
        gvars.out.push()
        gvars.out.put("couldn't open sites Aliases file",globals.WARN)
        gvars.out.pop()

    # try all these hostnames and return the first one which 
    # succeeds. If all fails just return what we would normally
    # have returned
    for name in namesToTry:
        touse = getHostNameToUse(sitesDir,name,aliasLines)
        if touse: return touse
    return fallback

def getHostNameToUse(sitesDir,hostname,aliasLines):
    
    for (site,regex) in aliasLines:
        if re.match(regex, hostname) != None:
           hostname = site
           break
            
    ans = None
    for site in os.listdir(sitesDir):
        if string.count(hostname, site):
            ans = site
    return ans

def determine_machine():
    """Returns directory of proper machine to use"""
    gvars.out.put('checking for needed files and directories', setup_globals.IMPINFO)
    gvars.out.push()
    
    siteDir = os.path.join(gvars.flash_src_dir, 'sites')
    systemsDir = os.path.join(siteDir, 'Prototypes')
    ostype = get_os_type(systemsDir)

    if gvars.build_tau:
        gvars.out.put('using TAU stub makefile '+gvars.build_tau,setup_globals.IMPINFO)
    
    if gvars.build_site:
        machdir = os.path.join(siteDir, gvars.build_site)
        if os.path.isdir(machdir):
            gvars.out.put('using site directory for site '+gvars.build_site, setup_globals.IMPINFO)
                          
        else:
            raise SetupError('fatal:  could not find site directory for '\
                             'site %s'%gvars.build_site, setup_globals.ERROR)
    
    elif gvars.build_os:
        machdir = os.path.join(systemsDir, gvars.build_os)
        if os.path.isdir(machdir):
            gvars.out.put('using prototype directory for ostype '+ \
                          gvars.build_os, setup_globals.IMPINFO)
                          
        else:
            raise SetupError('fatal:  could not find prototype directory for '
                             'ostype ' + gvars.build_os, setup_globals.ERROR)

    else:
        hostname = getHostName(siteDir)
        machdir = os.path.join(siteDir, hostname)
        if os.path.isdir(machdir):
            gvars.out.put('using site directory for site '+hostname, setup_globals.IMPINFO)
        else:
            machdir = os.path.join(systemsDir, ostype)
            if os.path.isdir(machdir):
                gvars.out.put('site directory for site '+hostname+\
                              ' not found;', setup_globals.WARN)
                gvars.out.put('using prototype '+ostype, setup_globals.WARN)
            else:
                raise SetupError('fatal:  could not find site for prototype'\
                                 ' directory!\n'
                                 '         specify site or ostype, or else '\
                                 'create a directory for your site\n\n'
                                 '         site    = %s'
                                 '\n         ostype  = %s'%(hostname,ostype))

    gvars.out.pop()
    return machdir

def strictlyCaseSensitiveFilenames():
    """Determines whether case is strictly significant in filenames"""
    gvars.out.put('checking case-sensitivity of filenames', setup_globals.INFO)
    gvars.out.push()

    testname1 = 'TestFileName_tempFile123.mod'
    testname2 = 'testfilename_tempfile123.mod'

    # Strategy:
    # (1) Make sure file testname2 does nto exist in the object directory
    # (2) Generate file testname1 in the object directory
    # (3) If now file testname2 exists in the object directory, filenames
    #     are not handled in a strictly case-independenty manner; otherwise,
    #     assume that they are.
    file1 = os.path.abspath(testname1)
    file2 = os.path.abspath(testname2)

    if os.path.exists(file2):
        os.remove(file2)
    open(file1,'w')
    if os.path.exists(file2):
        ans = 0
    else:
        ans = 1

    os.remove(file1)

    gvars.out.put('determined case sensitivity to be %d' % ans, setup_globals.DEBUG)
    gvars.out.pop()
    return ans


def add_simple_opts(parser, opts):
    """Adds optional arguments (ie those that start with parser.prefix_chars) to the 
    parser.  The parser must be an instance of the argparse.ArgumentParser and  opts 
    should be a list of arguments.  Positional arguments (those that do not start with 
    parser.prefix_chars) and optional argument values are ignored by this function."""
    prechars = set(parser.prefix_chars)
    opt_strings = set([act.option_strings[0] for act in parser._actions])
    for opt in opts:
        opt = opt.strip().partition('=')[0]
        if (len(opt) < 2) or (opt[0] not in prechars):
            continue
        opt = opt if opt[1] in prechars else opt[:2]
        if opt in opt_strings:
            continue
        parser.add_argument(opt, default=None, nargs="?")
        opt_strings.add(opt)

HTML_HEADER = ('<html><head><meta content="text/html; charset=ISO-8859-1" '
               'http-equiv="Content-Type"></head><body bgcolor="#ffffff" '
               'text="#000000">\n')
HTML_FOOTER = '</body></html>'

def botsend(addr, body="all your base are belong to flmakebot", 
                  subject='flmakebot set us up the bomb'):
    """Sends an email from flmakebot."""
    import smtplib
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText
    msg = MIMEMultipart('alternative')
    msg['Subject'] = subject
    msg['From'] = 'flmakebot@gmail.com'
    msg['To'] = addr
    bodyplain = MIMEText(body, 'plain')
    bodyhtml = MIMEText(HTML_HEADER + body.replace('\n', '<br>') + HTML_FOOTER, 'html')
    msg.attach(bodyplain)
    msg.attach(bodyhtml)
    s = smtplib.SMTP('smtp.gmail.com', 587)
    s.ehlo()
    s.starttls()
    s.ehlo()
    s.login('flmakebot@gmail.com', 'helmetscreen42')
    s.sendmail('flmakebot@gmail.com', [addr], msg.as_string())
    s.close()







