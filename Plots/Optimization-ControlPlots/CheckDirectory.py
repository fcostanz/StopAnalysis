import sys, os
import commands

def CheckDirectory(path):
    if not os.path.exists(path):
        os.makedirs(path)
