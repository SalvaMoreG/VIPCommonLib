
source ~/bin/PATHS/setALL_CHLOE01_PATHS_G6.sh

export Ubuntu16ROOTDIR=/CHLOE/VIP/MY_GAMOS6_VIP09/external/root/6.18.00/root/
export Ubuntu16ROOTINCDIR=${Ubuntu16ROOTDIR}/include/

export Ubuntu16ROOTLIBS=$(root-config --libs)

export TEST=$( echo ${Ubuntu16ROOTLIBS} | awk -v rd=${Ubuntu16ROOTDIR} '{ gsub("rootbuild",rd); print $0 }' )

export Ubuntu16ROOTLIBS=$TEST

