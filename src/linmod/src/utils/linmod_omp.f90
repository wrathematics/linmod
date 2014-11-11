! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module linmod_omp
  !$ use :: omp_lib
  implicit none
  
  
  integer, public, parameter :: linmod_omp_minsize = 2500
  !$ logical, public, parameter :: linmod_omp_hassimd = (openmp_version >= 201307)
  
end module
