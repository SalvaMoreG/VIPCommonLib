 
source ~/bin/PATHS/setALL_CHLOE01_PATHS_G6.sh

export Ubuntu18ROOTDIR=/CHLOE/VIP/MY_GAMOS6_VIP06/external/root/6.12.06/root
export Ubuntu18ROOTINCDIR=${Ubuntu18ROOTDIR}/include/

export Ubuntu18ROOTLIBS=$(root-config --libs)

export TEST=$( echo ${Ubuntu18ROOTLIBS} | awk -v rd=${Ubuntu18ROOTDIR} '{ gsub("rootbuild",rd); print $0 }' )
#	echo TEST: ${TEST}

export Ubuntu18ROOTLIBS=$TEST

echo Ubuntu18ROOTLIBS: ${Ubuntu18ROOTLIBS}
echo Ubuntu18ROOTINCDIR: ${Ubuntu18ROOTINCDIR}


