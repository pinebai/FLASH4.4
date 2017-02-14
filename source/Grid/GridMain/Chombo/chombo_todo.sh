#!/bin/bash
find . | grep -Fv .svn | xargs grep DEV | grep -Fv $0
