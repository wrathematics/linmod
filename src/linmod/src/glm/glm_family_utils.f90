!!! TODO rename to glm_internal ?
!!! TODO include link definitions

module glm_family_utils
  implicit none
  
  integer, parameter :: glm_family_unsupported = -1
  
  integer, parameter :: glm_family_gaussian = 1
  integer, parameter :: glm_family_binomial = 2
  integer, parameter :: glm_family_poisson = 3
  integer, parameter :: glm_family_gamma = 4
  
end module

