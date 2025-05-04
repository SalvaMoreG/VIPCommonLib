
source ~/bin/setALL_PATHS.sh

export VIPHOMEROOTDIR=/home/mkolstein/Programs/GAMOS/external/root/6.18.00/root/
export VIPHOMEROOTINCDIR=${VIPHOMEROOTDIR}/include/
export VIPHOMEROOTLIBDIR=${VIPHOMEROOTDIR}/lib/

export VIPHOMEROOTLIBS=$(root-config --libs)

echo ROOT lib dir: ${VIPHOMEROOTLIBDIR} 
echo ROOT libs: ${VIPHOMEROOTLIBS} 
echo ROOT inc dir: ${VIPHOMEROOTINCDIR} 

export TEST=$( echo ${VIPHOMEROOTLIBS} | awk -v rd=${VIPHOMEROOTDIR} '{ gsub("rootbuild",rd); print $0 }' )
echo TEST: ${TEST}

export VIPHOMEROOTLIBS=$TEST

echo VIPHOMEROOTLIBS: ${VIPHOMEROOTLIBS}






