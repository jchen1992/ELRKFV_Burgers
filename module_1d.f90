 module module_1d
    ! data structure
    
    !*******************************
    !  Background cell X(i)
    type, public :: merged_cell
        sequence
        real :: coor_left,coor_right
        real :: cell_itg
        real :: left_char,right_char
        real :: left_flux,right_flux
    end type

    type(merged_cell),allocatable,target,public :: merged_x(:)
    !*******************************
     !merge info
    type, public :: mgcell_info
        sequence
        integer :: cell_id
        integer :: numl,numr
        real :: xt
        real :: value(1:7)
    end type

    type(mgcell_info),allocatable,target,public :: mgcell_info_x(:)
    !*******************************
      ! Upstream cell X_star(i)
    type, public :: downstream_cell  ! downstream points (not cell, bad variable name)
        sequence
        real :: coor
        integer :: id
        real :: char
        real :: flux0, flux1,flux2
    end type

    type(downstream_cell),allocatable,target,public :: ds_x_star(:)

    !*******************************
    type, public :: element1d_downstream
        sequence
        type(downstream_cell) :: point_left,point_right
        real :: up_itg
        real :: int_len
    end type

    type(element1d_downstream),allocatable,target,public :: element_x_star(:)

!*****************************************
    type, public :: F_tail
    sequence
     
    real :: flux_left,flux_right
    end type
    
     type(F_tail),allocatable,target,public :: F_tail_x(:,:)
 !*************************************************   
     
    end module module_1d