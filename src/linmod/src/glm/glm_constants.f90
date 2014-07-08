! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_constants
  implicit none
  
  ! Family parameters
  integer, public, parameter :: glm_family_unsupported = -1
  
  integer, public, parameter :: glm_family_gaussian = 1
  integer, public, parameter :: glm_family_binomial = 2
  integer, public, parameter :: glm_family_poisson = 3
  integer, public, parameter :: glm_family_gamma = 4
  
  ! Link function parameters
  integer, public, parameter :: glm_link_unsupported = -2
  
  integer, public, parameter :: glm_link_cloglog = 1
  integer, public, parameter :: glm_link_identity = 2
  integer, public, parameter :: glm_link_inverse = 3
  integer, public, parameter :: glm_link_log = 4
  integer, public, parameter :: glm_link_logit = 5
  integer, public, parameter :: glm_link_sqrt = 6
  
  ! Intercept parameters
  integer, public, parameter :: glm_intercept_no = 0
  integer, public, parameter :: glm_intercept_yes = 1
  
  ! Stoprule
  integer, public, parameter :: glm_stoprule_maxiter = 1
  integer, public, parameter :: glm_stoprule_coefs = 2
  integer, public, parameter :: glm_stoprule_deviance = 3
  
  ! Convergence
  integer, public, parameter :: glm_convergence_noconvergence = -1
  integer, public, parameter :: glm_convergence_converged = 1
  integer, public, parameter :: glm_convergence_infparams = 2
  integer, public, parameter :: glm_convergence_nochange = 3
  
end module
