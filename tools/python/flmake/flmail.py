import os
import glob
import sys
import datetime
import subprocess

# relative imports needed
from .utils import botsend
from .. import _metadata

USAGE = ("Sends an email from flmakebot. Defaults are\n"
         "provided as well as standard messages for\n"
         "execution starting and stopping.\n\n"
         "usage: flmake email [-h] [--start] [--stop] [-b BODY] [-s SUBJECT] [-t ADDR]")

def _start_msg():
    user = os.getenv('USER', "<unkown>")
    host = os.uname()[1]
    cdir = os.getcwd()
    curt = str(datetime.datetime.now())
    subject = "Started FLASH run {0}@{1}:{2} at {3}"
    subject = subject.format(user, host, cdir, curt)
    body  = "Run started at {0}\n".format(curt)
    body += "Run user was {0}\n".format(user)
    body += "Run on {0}\n".format(host)
    body += "Run located in {0}\n".format(cdir)
    if os.path.isfile('flash.par'):
        body += '\nRuntime Parameters\n\n<div><font face="courier new, monospace">'
        with open('flash.par') as f:
            body += f.read()
        body += '</font></div>'
    return body, subject

def _stop_msg():
    user = os.getenv('USER', "<unkown>")
    host = os.uname()[1]
    cdir = os.getcwd()
    curt = str(datetime.datetime.now())
    subject = "Stopped FLASH run {0}@{1}:{2} at {3}"
    subject = subject.format(user, host, cdir, curt)
    body  = "Run stopped at {0}\n".format(curt)
    body += "Run user was {0}\n".format(user)
    body += "Run on {0}\n".format(host)
    body += "Run located in {0}\n".format(cdir)
    logs = sorted([(len(l), l) for l in glob.glob('*.log')])
    if 0 < len(logs):
        log = logs[0][1]
        body += "\nLast 100 lines of {0}\n\n".format(log)
        body += '<div><font face="courier new, monospace">'
        with open(log) as f:
            lines = f.readlines()[-100:]
        body += "".join(lines)
        body += '</font></div>'
    return body, subject

def main(ns, rc):
    """sends an email from flmakebot"""
    # handle when we are in a run directory
    if 0 == len(rc) and os.path.isfile('../flashrc.py'):
        execfile('../flashrc.py', rc, rc)
    default_subject = 'flmakebot set us up the bomb'
    default_body = "all your base are belong to flmakebot"

    addrs = [ns.addr, rc.get('email', None), _metadata.get('email', None)]
    addr = [a for a in addrs if a is not None]
    if 0 == len(addr):
        print "No email address found!"
        return 1
    addr = addr[0]

    start_body, start_subject = _start_msg() if ns.start else (None, None) 
    stop_body, stop_subject = _stop_msg() if ns.stop else (None, None) 

    bodies = [ns.body, stop_body, start_body, default_body]
    body = [s for s in bodies if s is not None][0]

    subjects = [ns.subject, stop_subject, start_subject, default_subject]
    subject = [s for s in subjects if s is not None][0]

    botsend(addr, body, subject)
    return 0
