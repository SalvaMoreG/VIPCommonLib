 
source ~/bin/setALL_PATHS.sh

export VIP39ROOTDIR=/export/VIP/MY_GAMOS6_VIP06/external/root/6.12.06/root
export VIP39ROOTINCDIR=${VIP39ROOTDIR}/include/
export VIP39ROOTLIBDIR=${VIP39ROOTDIR}/lib/

export VIP39ROOTLIBS=$(root-config --libs)
#	echo ROOT dirs: ${VIP39ROOTINCDIR} " " ${VIP39ROOTLIBDIR} 

export TEST=$( echo ${VIP39ROOTLIBS} | awk -v rd=${VIP39ROOTDIR} '{ gsub("rootbuild",rd); print $0 }' )
#	echo TEST: ${TEST}

export VIP39ROOTLIBS=$TEST

echo VIP39ROOTLIBS: ${VIP39ROOTLIBS}
echo VIP39ROOTINCDIR: ${VIP39ROOTINCDIR}


