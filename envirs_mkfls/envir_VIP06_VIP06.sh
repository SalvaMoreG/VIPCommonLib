 
source ~/bin/setALL_PATHS.sh

export VIP06ROOTDIR=/export/VIP/MY_GAMOS6_VIP06/external/root/6.12.06/root
export VIP06ROOTINCDIR=${VIP06ROOTDIR}/include/
export VIP06ROOTLIBDIR=${VIP06ROOTDIR}/lib/

export VIP06ROOTLIBS=$(root-config --libs)
#	echo ROOT dirs: ${VIP06ROOTINCDIR} " " ${VIP06ROOTLIBDIR} 

export TEST=$( echo ${VIP06ROOTLIBS} | awk -v rd=${VIP06ROOTDIR} '{ gsub("rootbuild",rd); print $0 }' )
#	echo TEST: ${TEST}

export VIP06ROOTLIBS=$TEST

echo VIP06ROOTLIBS: ${VIP06ROOTLIBS}


