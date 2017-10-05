# CMinPack 
include_directories("${PROJECT_SOURCE_DIR}/cminpack/include")
set (cminpack_srcs
  cminpack/include/cminpackP.h
  cminpack/src/chkder.c
  cminpack/src/enorm.c
  cminpack/src/hybrd1.c
  cminpack/src/hybrj.c
  cminpack/src/lmdif1.c
  cminpack/src/lmstr1.c
  cminpack/src/qrfac.c
  cminpack/src/r1updt.c
  cminpack/src/dogleg.c
  cminpack/src/fdjac1.c
  cminpack/src/hybrd.c
  cminpack/src/lmder1.c
  cminpack/src/lmdif.c
  cminpack/src/lmstr.c
  cminpack/src/qrsolv.c
  cminpack/src/rwupdt.c
  cminpack/src/dpmpar.c
  cminpack/src/fdjac2.c
  cminpack/src/hybrj1.c
  cminpack/src/lmder.c
  cminpack/src/lmpar.c
  cminpack/src/qform.c
  cminpack/src/r1mpyq.c
  cminpack/src/covar.c
  cminpack/src/covar1.c
  cminpack/src/chkder_.c 
  cminpack/src/enorm_.c  
  cminpack/src/hybrd1_.c 
  cminpack/src/hybrj_.c  
  cminpack/src/lmdif1_.c 
  cminpack/src/lmstr1_.c 
  cminpack/src/qrfac_.c  
  cminpack/src/r1updt_.c
  cminpack/src/dogleg_.c 
  cminpack/src/fdjac1_.c 
  cminpack/src/hybrd_.c  
  cminpack/src/lmder1_.c 
  cminpack/src/lmdif_.c  
  cminpack/src/lmstr_.c  
  cminpack/src/qrsolv_.c 
  cminpack/src/rwupdt_.c
  cminpack/src/dpmpar_.c 
  cminpack/src/fdjac2_.c 
  cminpack/src/hybrj1_.c 
  cminpack/src/lmder_.c  
  cminpack/src/lmpar_.c  
  cminpack/src/qform_.c  
  cminpack/src/r1mpyq_.c 
  cminpack/src/covar_.c
  )
set (cminpack_hdrs
  cminpack/include/cminpack.h 
  cminpack/include/minpack.h)
add_definitions(-DUSE_CBLAS)
add_definitions(-DUSE_LAPACK)
if (MSVC)
  add_definitions(-DCMINPACK_NO_DLL)
  source_group("CMinPack Sources" FILES ${cminpack_srcs})
  source_group("CMinPack Headers" FILES ${cminpack_hdrs})
endif(MSVC)