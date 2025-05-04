 
source ~/bin/setALL_PATHS.sh

export VIP12ROOTDIR=/export/VIP/MY_GAMOS6_VIP06/external/root/6.12.06/root
export VIP12ROOTINCDIR=${VIP12ROOTDIR}/include/
export VIP12ROOTLIBDIR=${VIP12ROOTDIR}/lib/

export VIP12ROOTLIBS=$(root-config --libs)
#	echo ROOT dirs: ${VIP12ROOTINCDIR} " " ${VIP12ROOTLIBDIR} 

export TEST=$( echo ${VIP12ROOTLIBS} | awk -v rd=${VIP12ROOTDIR} '{ gsub("rootbuild",rd); print $0 }' )
#	echo TEST: ${TEST}

export VIP12ROOTLIBS=$TEST

echo VIP12ROOTLIBS: ${VIP12ROOTLIBS}
echo VIP12ROOTINCDIR: ${VIP12ROOTINCDIR}


