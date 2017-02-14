import os
import json

# load cached metadata
with open(os.path.split(__file__)[0] + "/metadata.json", 'r') as f:
    _metadata = json.load(f)

# set version string
__version__ = _metadata['version']

# set flash source code home directory
if os.getenv('FLASH_SRC_DIR'):
    FLASH_SRC_DIR = os.environ['FLASH_SRC_DIR']
    _metadata['FLASH_SRC_DIR'] = FLASH_SRC_DIR
else:
    FLASH_SRC_DIR = _metadata['FLASH_SRC_DIR']

# set flash clean source code home directory
if os.getenv('FLASH_CLEAN_SRC_DIR'):
    FLASH_CLEAN_SRC_DIR = os.environ['FLASH_CLEAN_SRC_DIR']
    _metadata['FLASH_CLEAN_SRC_DIR'] = FLASH_CLEAN_SRC_DIR
else:
    FLASH_CLEAN_SRC_DIR = _metadata['FLASH_CLEAN_SRC_DIR']

# set mpirun command
if os.getenv('MPIRUN_CMD'):
    MPIRUN_CMD = os.environ['MPIRUN_CMD']
    _metadata['MPIRUN_CMD'] = MPIRUN_CMD
else:
    MPIRUN_CMD = _metadata['MPIRUN_CMD']
