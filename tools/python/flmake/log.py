import os
from datetime import datetime

USAGE = ("Displays the history of all (or N) previously\n"
         "executed flmake commands and their metadata.\n\n"
         "usage: flmake log [-n <N>]")

def _parse_row(row):
    t, cmd, user, logid, runid, d, msg = row.split(',', 6)
    return float(t), cmd, user, logid, runid, d, msg[1:-1]


def main(ns, rc):
    """Displays the history of flmake commands."""
    if not os.path.exists('flash.log'):
        return

    with open('flash.log') as f:
        loglines = f.readlines()[::-1]

    logtemplate = ("Run id: {runid}\n"
                   "Run dir: {rundir}\n"
                   "Command: {cmd}\n"
                   "User: {user}\n"
                   "Date: {dt}\n"
                   "Log id: {logid}\n\n"
                   "    {msg}\n\n"
                   )

    logstr = ""
    for row in loglines[:ns.n]:
        t, cmd, user, logid, runid, d, msg = _parse_row(row[:-1])
        dt = datetime.fromtimestamp(t).strftime("%c")
        kwlog = {'dt': dt, 'cmd': cmd, 'user': user, 'logid': logid, 
                 'runid': runid, 'rundir': d, 'msg': msg}
        logstr += logtemplate.format(**kwlog)
    logstr = logstr[:-1]

    print logstr
