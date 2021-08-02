module wam_ifp_class

      type farr_t
        real, allocatable, dimension(:) :: f107
        real, allocatable, dimension(:) :: kp
        real, allocatable, dimension(:) :: kpa
        real, allocatable, dimension(:) :: f107d
        real, allocatable, dimension(:) :: nhp
        real, allocatable, dimension(:) :: nhpi
        real, allocatable, dimension(:) :: shp
        real, allocatable, dimension(:) :: shpi
        real, allocatable, dimension(:) :: swbz
        real, allocatable, dimension(:) :: swvel
        real, allocatable, dimension(:) :: swbt
        real, allocatable, dimension(:) :: swang
        real, allocatable, dimension(:) :: swden
      end type farr_t

      type forcing_t
        real :: f107
        real :: f107d
        real :: kp
        real :: kpa
        real :: nhp
        real :: nhpi
        real :: shp
        real :: shpi
        real :: swbz
        real :: swvel
        real :: swbt
        real :: swang
        real :: swden
      end type forcing_t

      type param_t
        integer :: ifp_interval
        integer :: kdt_start
        integer :: skip
        integer :: ifp_realtime_interval
      end type param_t

end module wam_ifp_class
