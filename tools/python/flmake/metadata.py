import os
import subprocess

USAGE = ("Displays or edits the metadata file.\n\n"
         "usage: flmake metadata [-e/--edit]")
         

def main(ns, rc):
    """Displays or edits the metadata file."""
    metadata_path = os.path.split(os.path.split(__file__)[0])[0] + "/metadata.json"

    if ns.edit:
        editor = os.getenv("EDITOR")
        subprocess.call([editor, metadata_path])
    else:
        with open(metadata_path, 'r') as f:
            print f.read()
