source ~/bin/setALL_PATHS.sh

export VIP07U22ROOTDIR=/export/VIP/MY_GAMOS6_U22/external/root/6.26.06/root/
export VIP07U22ROOTINCDIR=${VIP07U22ROOTDIR}/include/
export VIP07U22ROOTLIBDIR=${VIP07U22ROOTDIR}/lib/

export VIP07U22ROOTLIBS=$(root-config --libs)

echo ROOT lib dir: ${VIP07U22ROOTLIBDIR} 
echo ROOT libs: ${VIP07U22ROOTLIBS} 
echo ROOT inc dir: ${VIP07U22ROOTINCDIR} 

export TEST=$( echo ${VIP07U22ROOTLIBS} | awk -v rd=${VIP07U22ROOTDIR} '{ gsub("rootbuild",rd); print $0 }' )
echo TEST: ${TEST}

export VIP07U22ROOTLIBS=$TEST

echo VIP07U22ROOTLIBS: ${VIP07U22ROOTLIBS}






