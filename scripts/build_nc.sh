#!/bin/bash
# ============================================================= #
#				Points2Grid Build Script						#
#	This is script will be used to track compilation steps for  #
#	Parallel points2grid on ROGER. This is mainly useful when   #
#	a dependency has broken/moved and we need to see recompile  #
#	quickly.													#
#
#	by: Nathan Casler											#
# ============================================================= #

module load gdal2-stack
module load mpich
module load boost
module load curl

BASE=$(dirname "$(dirname "$(readlink -f "$0")" )")

cd "$BASE"

if [ -f "Makefile" ]; then
	echo "Makefile found"
else
	echo "Failed to find Makefile"
	exit 1
fi

make clean && make





