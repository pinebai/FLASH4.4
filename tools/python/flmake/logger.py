import os
import uuid
import logging

logging.basicConfig(
    filename = 'flash.log',
    format = "%(created)f,%(cmd)s,%(user)s,%(logid)s,%(runid)s,%(rundir)s,%(msg)r",
    level=logging.INFO,
    )


def info(msg, cmd=None, runid=None, rundir=None):
    """Log to flmake with INFO status."""
    cmd = "<no-cmd>" if cmd is None else cmd
    user = os.getenv('USER', "<no-user>")
    logid = uuid.uuid4()
    runid = "<no-id>" if runid is None else runid
    rundir = "<no-dir>" if rundir is None else rundir
    logging.info(msg, extra={'cmd': cmd, 'user': user, 'logid': logid, 
                             'runid': runid, 'rundir': rundir})
