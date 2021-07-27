import fileinput
import subprocess
import sys

def findReplace(filename, toReplace, replacement):
    for line in fileinput.input(filename, inplace=True):
        print line.replace(toReplace, replacement),

def getCompressionLine(filename, suffix):
    cmd = "cat " + filename + " | " + "grep " + suffix
    compressionPath = subprocess.check_output(cmd, shell=True)
    return compressionPath

################################################################################
argc = len(sys.argv)
if argc != 3:
    print "You must pass two command line arguments: a filename, and a compression ratio!"
    sys.exit(-1)

# Get the arguments list 
cmdargs = str(sys.argv)

# By convention, sys.argv[0] stores the name of this script, so we skip to index 1
filename = sys.argv[1]
compressionRatio = sys.argv[2]

# Find and replace the compression path
suffix = "to1/\n"
replacement = "compression path = ./data/reduced.stam.64/" + str(compressionRatio) + suffix
toReplace = getCompressionLine(filename, suffix)

findReplace(filename, toReplace, replacement)

# Find and replace the movie path
suffix = "to1.mov\n"
replacement = "preview movie = ./movies/compressed." + str(compressionRatio) + suffix
toReplace = getCompressionLine(filename, suffix)

findReplace(filename, toReplace,replacement)
