
source ~/bin/setALL_PATHS_G6.sh

export VIP09ROOTDIR=/export/VIP/MY_GAMOS6_VIP09/external/root/6.18.00/root/
export VIP09ROOTINCDIR=${VIP09ROOTDIR}/include/
export VIP09ROOTLIBDIR=${VIP09ROOTDIR}/lib/

export VIP09ROOTLIBS=$(root-config --libs)

export TEST=$( echo ${VIP09ROOTLIBS} | awk -v rd=${VIP09ROOTDIR} '{ gsub("rootbuild",rd); print $0 }' )

export VIP09ROOTLIBS=$TEST

