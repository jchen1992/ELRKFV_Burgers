!**********************************************
! New ELFV for shocks
! Burgers' equation u_t+1/2 *(u^2)_x=0;
!*********************************************
    
    program SLRK_FV_1d
    use module_1d
    
    implicit none
   integer, parameter :: nl =1, N_x = 200, nd =20, nxd = nl*N_x+nd+1
    integer :: i,j,k,l,knt
    integer :: nx,nt
    real :: dx
    real :: pi,eps
    real :: xg(6),wg(6)
    real :: xleft,xright
    real :: dt,tt,T,cfl,tvd, umax, umin,udif
    real,allocatable :: X(:)
    real,allocatable :: Xmid(:)

    real,allocatable :: u(:,:)
    real :: exact_ave(1-nd:nxd-1)
    real :: error1(1:nl),error2(1:nl),error3(1:nl)
    real :: ord(1:nl)
    integer :: merged_cell_total,iexample

    real,allocatable :: us_n(:)
    integer :: i_stage
    real :: er1,er2,er3
    real :: sum_out, mass,mass_temp,mass_max

    real,allocatable :: speed(:)

 !   integer :: iexample
    
    pi = 4.D0 * datan(1.D0)
    eps = 1.e-8
    xg(1)=-0.466234757101576013906150777246997304567e0
    xg(2)=-0.330604693233132256830699797509952673503e0
    xg(3)=-0.119309593041598454315250860840355967709e0
    xg(4)=-xg(3)
    xg(5)=-xg(2)
    xg(6)=-xg(1)
    wg(1)=1.71324492379170345040296142172733e-1/2e0
    wg(2)=3.60761573048138607569833513837716e-1/2e0
    wg(3)=4.67913934572691047389870343989551e-1/2e0
    wg(4)=wg(3)
    wg(5)=wg(2)
    wg(6)=wg(1)
    
        umax = 1.
        umin = -1.
        udif = umax-umin
    
    call setup
    
    do cfl = 0.1,40.,0.03
    do l=1,nl
    nx = l * N_x
    !*******************************
      !allocate_variables
        allocate( merged_x(-nd:nx+nd+1) )
        allocate( speed(-nd:nx+nd+1) )
        allocate( X(-nd:nx+nd+1) )
        allocate( Xmid(-nd:nx+nd) )
        allocate( element_x_star(-nd:nx+nd+1) )
        allocate( u(-nd:nx+nd+1,0:1) )
        allocate(F_tail_x(1:nx+1,0:0))
        allocate( us_n(1:nx) )
        !*******************************
        
    call init
 !        call CPU_time(t_start)
  tvd = 0.
        do i = 1, nx-1
         tvd= tvd+abs(u(i,0)-u(i+1,0))
        enddo
        Write(1000,*) 'Constant C=', cfl

    		 tt = 0.D0
		 do knt=1,nt
            if(tt+dt.gt.T)then
		    dt = T - tt
		    endif 
         call evolve
         tt = tt + dt  

       tvd = 0.
        do i = 1, nx-1
         tvd = tvd+abs(u(i,0)-u(i+1,0))
        enddo
        Write(1000,*) 'step:',knt, tvd 
        
         write(*,*) 'knt=',knt
        enddo

  !       write(*,*) 'order start'
		call order(l)
  !      write(*,*) 'order end'
!		call CPU_time(t_stop)
!		write(*,*) 'N=', l*N_x , 'cpu time is', t_stop - t_start
!************************************
     !deallocate_variables
        deallocate( merged_x )
        deallocate( speed )
        deallocate( X )
        deallocate( Xmid )
        deallocate( element_x_star )
        deallocate( u )
        deallocate( F_tail_x )
        deallocate( us_n )
    !end deallocate_variables
!*******************************
    enddo
    !call err_ord
    write(7000,*) cfl/udif*umax, error1(1)
    enddo
    
     contains
     
!*******************************************
     subroutine setup
     implicit none

		!write(*,*)'T, e.g. T = 2.D0*pi/5.D0 for test problem'
		!read(*,*) T
         T = 0.5
		!write(*,*)'cfl'
		!read(*,*) cfl
!        cfl = 3.0
 !       write(*,*) 'cfl=',cfl
        !write(*,*) '1.rarfaction, 2.shock, 3. sin, 4. Buckley-Leverett, 5. Non-Convex, 6. linear adv(sin), 7. variable coeff  '
        !read(*,*) iexample
        iexample = 3
        
        if (iexample .eq. 3)then
  !      xleft = pi/2.
		!xright = 5.*pi/2.
        xleft = 0.
		xright = 2.*pi
        elseif(iexample .eq. 4)then
        xleft = 0.
		xright = 1.
        elseif(iexample .eq. 5)then
        xleft = 0.
		xright = 2.   
        else
        xleft = -pi
		xright = pi
        endif
        
        end subroutine setup
!********************************************
        
        subroutine init
		implicit none

        real(8) :: alpha, z, delta, a, beta
        real :: utemp
        integer :: lx
        alpha = 10.D0
        delta = 0.005D0
        beta = log(2.D0)/(36.D0* delta**2.D0)
        a = 0.5D0
        z = -0.7D0
        
        mass = 0.
         i_stage = 0
         merged_cell_total = nx
         
		dx = (xright-xleft)/nx
 		if(iexample .eq. 2)then
        
        dt = cfl * dx/udif
        else
        dt = cfl * dx/udif
        endif
        nt = ceiling(T/dt)
        write(*,*) 'nt=',nt, 'dt/dx=', dt/dx
		do i = -nd, nx+nd+1
! interval end pts
        X(i) = xleft + (i-1)*dx
        enddo
        
        do i =-nd, nx+nd
! mid pts
		Xmid(i) = xleft + i*dx -dx/2.D0
		enddo
        
    do i = 1, nx
          utemp = 0.
        do lx = 1, 6
            utemp = utemp + exact(Xmid(i)+xg(lx)*dx,0.)*wg(lx)
        enddo
        u(i,0) = utemp

        mass = mass + u(i,0)  
    enddo 
    
    call bc(i_stage)
    
         do i = -nd, nx+nd
           merged_x(i)%coor_left = xleft + (i-1) * dx
           merged_x(i)%coor_right = xleft + i * dx
           merged_x(i)%cell_itg = u(i,0)*dx     

        enddo
        
		end subroutine init
!*******************************  
 subroutine mass_check(mass_in)
  implicit none
  real,intent(in) :: mass_in(1:merged_cell_total)
  integer :: i
  
  mass_temp = 0.
  do i = 1, merged_cell_total
      mass_temp = mass_temp+ mass_in(i)
  
  enddo
        
   write(*,*) 'mass error=', abs(mass*dx-mass_temp)
  
    end subroutine  mass_check   
!***********************************
    real function exact(x,t)
    implicit none
    real,intent(in) :: x,t

    if(iexample == 1)then
   if (t .eq. 0.)then 
    if(x .le. 0. )then
    exact = -1.
    else
    exact = 1.
    endif
    else
    
    if(x .lt. -t )then
    exact = -1.
    elseif(x .gt. t )then
    exact = 1.
    else
    exact = x/t
    endif
    endif
    
    elseif (iexample == 2)then
    if (t .eq. 0.)then
    if(x .le. 0. )then
    exact = 2.
    else
    exact = -2.
    endif
    else
    
    if (x .le. 0.*t ) then
    exact = 2.
    elseif (x .gt. 0.*t ) then
    exact = -2.
    endif
    endif
    
    elseif (iexample == 3 .or. iexample == 6)then
      exact = sin( x-t )
    elseif (iexample == 4)then
        if(x .ge. 0. .and. x.le. 0.05 )then
            exact = 1.-20.*x
        elseif(x.ge.0.25 .and. x.le. 0.4) then
            exact = 0.5    
        else
            exact =0.
        endif
    elseif (iexample == 5)then
        if(x .ge. 0. .and. x.le. 0.05 )then
            exact = 1.-20.*x
        elseif(x.ge.0.25 .and. x.le. 0.4) then
            exact = 0.5    
        else
            exact =0.
        endif
    elseif (iexample .eq. 7)then
        exact = 1.
    endif


    end function exact
    !*******************************
     subroutine bc(i_stage)
		 implicit none
        integer,intent(in)  :: i_stage
     
    if (iexample .eq. 3 .or. iexample == 6 .or. iexample .eq. 7) then
		 do i=0, nd
		 u(-i,i_stage) = u(merged_cell_total-i,i_stage)
		 u(merged_cell_total+1+i,i_stage) = u(i+1,i_stage)
		 enddo
    elseif (iexample .eq. 4) then
        do i=0, nd
		 u(-i,i_stage) = 1.  
		 u(merged_cell_total+1+i,i_stage) = 0.
        enddo
    elseif (iexample .eq. 5) then
        do i=0, nd
		 u(-i,i_stage) = 1.  
		 u(merged_cell_total+1+i,i_stage) = 0.
        enddo
    else    
        do i=0, nd
		 u(-i,i_stage) = u(-i+1,i_stage)
		 u(merged_cell_total+1+i,i_stage) = u(merged_cell_total,i_stage)
        enddo
    endif    
		 end subroutine bc
!***********************************
         subroutine bc_merged
         implicit none
     if (iexample .eq. 3 .or. iexample == 6 .or. iexample .eq. 7) then    
        do i=0, nd
		   merged_x(-i)%cell_itg = merged_x(merged_cell_total-i)%cell_itg
		   merged_x(merged_cell_total+1+i)%cell_itg = merged_x(i+1)%cell_itg
        enddo
    else
        do i=0, nd
		   merged_x(-i)%cell_itg = merged_x(-i+1)%cell_itg
           merged_x(merged_cell_total+1+i)%cell_itg = merged_x(merged_cell_total)%cell_itg       
        enddo
    endif
    
         end subroutine bc_merged
!************************************
    subroutine bc_element_x_star
    implicit none
    integer :: i

    if (iexample .eq.3 .or. iexample == 6 .or. iexample .eq. 7) then
    do i=0, nd
        element_x_star(-i)%up_itg = element_x_star(merged_cell_total-i)%up_itg
        element_x_star(merged_cell_total+1+i)%up_itg = element_x_star(i+1)%up_itg
    enddo
    else
    do i=0, nd
        element_x_star(-i)%up_itg = element_x_star(-i+1)%up_itg
        element_x_star(merged_cell_total+1+i)%up_itg = element_x_star(merged_cell_total)%up_itg
    enddo
    endif
    
    end subroutine bc_element_x_star   
         
!**************************************
    subroutine evolve
		 implicit none

         type(element1d_downstream), pointer :: p(:)
         integer :: i

         merged_cell_total = nx
         
         i_stage = 0
         call bc(i_stage)
       
        do i = -nd, nx+nd
           merged_x(i)%coor_left = xleft + (i-1) * dx
           merged_x(i)%coor_right = xleft + i * dx
           merged_x(i)%cell_itg = u(i,0)*dx     
           
           element_x_star(i)%up_itg = merged_x(i)%cell_itg 
        enddo

         call char_comput
         call downstream_locate(0.)
    ! Step 1 ----------------------------------------------------------------------------------------    
         do i = 1, nx
          p => element_x_star(i-3:i+3)
         
		 call flux_comput(p,F_tail_x(i,0)%flux_left,F_tail_x(i,0)%flux_right)  !reconstruction of F_tail at each approx chara line

          enddo 
          
        call downstream_locate(dt) 
        call merging
        call bc_merged
        call bc_element_x_star
        
         do i = 1, merged_cell_total
      
         us_n(i) = merged_x(i)%cell_itg

         element_x_star(i)%up_itg = us_n(i)-dt*(element_x_star(i)%point_right%flux0 - element_x_star(i)%point_left%flux0)
        
        enddo        
!reaverage -------------------------------
        
        call bc_element_x_star
        
         call reproject(u(1:nx,1))
         u(1:nx,0) = u(1:nx,1)
                
		 end subroutine evolve
!**************************************************
     subroutine merging
    implicit none
    integer :: merge_rec(1:nx), i, j, k, m
    real :: cross_t(1:nx)
    integer :: merge_num,cross_num
    type(merged_cell),pointer :: p_merged
    type(mgcell_info), pointer ::  p_cell, p2
    integer :: left_id, right_id

        merged_cell_total = nx
        merge_num = 2

        do i = 1,nx
           merged_x(i)%left_char = speed(i)
           element_x_star(i)%point_left%char = speed(i)
           merged_x(i)%right_char = speed(i+1)
           element_x_star(i)%point_right%char = speed(i+1)
           merged_x(i)%left_flux = F_tail_x(i,0)%flux_left
           element_x_star(i)%point_left%flux0 = F_tail_x(i,0)%flux_left
           merged_x(i)%right_flux = F_tail_x(i,0)%flux_right
           element_x_star(i)%point_right%flux0 = F_tail_x(i,0)%flux_right
        enddo
        call bc_element_x_star
        call bc_merged	

!100     continue
            cross_t(:) = T*2
            merge_rec(:) = -2*nd
            cross_num = 0
            !p_cell = -2*nd

            cross_num=0
        do j=1,merged_cell_total
        if (element_x_star(j)%point_left%coor .ge. element_x_star(j)%point_right%coor) then
            cross_num = cross_num +1
            merge_rec(cross_num) = j
            p_merged => merged_x(j) 
            call cross_time(p_merged, cross_t(cross_num))
            !write(*,*) 'the', cross_num,  '-th cross is at', j, 'time=', cross_t(cross_num)
            write(*,*) knt, 'cross num=', cross_num
        endif
        enddo
        if(cross_num .eq. 0)then
        !write(*,*) 'No merging'    
        return
        endif
        
        allocate(mgcell_info_x(1:cross_num))
        
        do i = 1, cross_num
        mgcell_info_x(i)%cell_id = merge_rec(i)
        mgcell_info_x(i)%xt = cross_t(i)
        mgcell_info_x(i)%value(1:7) = u(mgcell_info_x(i)%cell_id-3:mgcell_info_x(i)%cell_id+3,0)
        p2 => mgcell_info_x(i)
        call mgnum_criteria(p2)
        Write(2008,*) cfl, knt, mgcell_info_x(i)%cell_id,mgcell_info_x(i)%numl,mgcell_info_x(i)%numr
        enddo
      
        i=1
        do while (i .le. cross_num )
        p_cell => mgcell_info_x(i)
        if ( i .lt. cross_num ) then
        if ( (mgcell_info_x(i+1)%cell_id - p_cell%cell_id) .le. 1 ) then
            i = i+1   
        endif
        endif
            
            k=0
            left_id = p_cell%cell_id - p_cell%numl 
            right_id = p_cell%cell_id + p_cell%numr 
            k = i+1
  
            do while(k .le. cross_num)
            if ((mgcell_info_x(k)%cell_id-p_cell%cell_id) .le. (mgcell_info_x(k)%numl + p_cell%numr))then
            p_cell => mgcell_info_x(k)
            i = i+1
                if(k .lt. cross_num) then
                if((mgcell_info_x(k+1)%cell_id-p_cell%cell_id) .le. 1)then         
                i=i+1
                k=k+1
                endif
                endif
                
            right_id = p_cell%cell_id + p_cell%numr     
            k=k+1
            i=i+1
            else
            exit
            endif
            enddo
            i=i+1
        !write(*,*)  p_cell%cell_id, left_id, right_id
        if (cross_num .gt. 0 ) then
            merged_x(left_id)%coor_right = merged_x(right_id)%coor_right
            element_x_star(left_id)%point_right = element_x_star(right_id)%point_right
            merged_x(left_id)%right_char = merged_x(right_id)%right_char  ! characteristic of right end pt
            merged_x(left_id)%right_flux = merged_x(right_id)%right_flux
            merged_x(left_id)%cell_itg = sum(merged_x(left_id : right_id)%cell_itg)
            element_x_star(left_id)%int_len = element_x_star(right_id)%point_right%coor &
            & -element_x_star(left_id)%point_left%coor
            merged_cell_total = nx - (right_id-left_id)
        write(*,*) 'merged_cell_total=',merged_cell_total

        do j = left_id+1, merged_cell_total+nd 
           merged_x(j) = merged_x(j+right_id-left_id)
           element_x_star(j) = element_x_star(j+right_id-left_id)
        enddo 
        
        do m = i,cross_num
            mgcell_info_x(m)%cell_id = mgcell_info_x(m)%cell_id - (right_id-left_id)
        enddo
        endif
    enddo 

     deallocate(mgcell_info_x)
    !if (cross_num .gt. 0)then
    !goto 100
    !endif
    
end subroutine merging
!****************************************
  subroutine mgnum_criteria(pc)
    implicit none
    type(mgcell_info),pointer :: pc
    !integer,intent(out) :: num_out
    
    if((pc%value(6)+pc%value(7)) .lt. (umax+3.*umin)/2. .and. (pc%value(3)+pc%value(4) &
    & +pc%value(5)) .ge. (7.*umax+5.*umin)/4.)then
    pc%numl = 2
    pc%numr = 3        !  i-2,....,i+3
    elseif((pc%value(1)+pc%value(2)) .gt. (3.*umax+umin)/2. .and. (pc%value(3)+pc%value(4) &
    & +pc%value(5)) .le. (5.*umax+7.*umin)/4.)then
    pc%numl = 3
    pc%numr = 2       ! i-3,...,i+2
    else
    pc%numl = 2
    pc%numr = 2          ! i-2,...,i+2
    endif
    end subroutine mgnum_criteria
!*************************************************************************
        subroutine cross_time(pc,i_time)
        implicit none
        type(merged_cell),pointer :: pc
        real,intent(out) :: i_time
        
        i_time=(pc%coor_right-pc%coor_left)/(pc%left_char-pc%right_char)
        
       end subroutine 
!***************************************
         subroutine char_comput
		 implicit none

		do i = 1-nd,nx+nd+1

           call roe_speed( u(i-1,0),u(i,0),X(i),speed(i) )

        enddo
          
          do i = 1-nd,nx+nd
           element_x_star(i)%point_left%char = speed(i)
           element_x_star(i)%point_right%char = speed(i+1)
          enddo

		 end subroutine char_comput    
!************************************************
     subroutine downstream_locate(dt_local)    
      implicit none
      real,intent(in) :: dt_local
      integer :: i
      do i = 1-nd,merged_cell_total+nd
           
        element_x_star(i)%point_left%coor = merged_x(i)%coor_left+element_x_star(i)%point_left%char*dt_local
        element_x_star(i)%point_right%coor = merged_x(i)%coor_right+element_x_star(i)%point_right%char*dt_local
           
        element_x_star(i)%int_len = element_x_star(i)%point_right%coor-element_x_star(i)%point_left%coor
            
        element_x_star(i)%point_left%id = ceiling( (element_x_star(i)%point_left%coor-xleft)/dx )
        element_x_star(i)%point_right%id = ceiling( (element_x_star(i)%point_right%coor-xleft)/dx )
        
     enddo
     
     end subroutine downstream_locate
!***************************************************
    subroutine roe_speed( uleft,uright,x,speed )
    implicit none
    real,intent(in) :: uleft,uright,x
    real,intent(out) :: speed

    if( abs(uright - uleft)<eps )then
        speed = fp(uleft*0.5+uright*0.5,x)
    else
        speed = ( f(uright,x) - f(uleft,x) )/(uright-uleft)
    endif

    end subroutine  roe_speed
!***********************************
    
    real function fp(u,x)
    implicit none
    real,intent(in) :: u,x

    if(iexample.eq.4)then
        fp = (2.*u*(1-u))/(2.*u**2.-2.*u+1)**2.
    elseif(iexample.eq.5)then
    if(u .le. 0.01)then
        fp = 10.
    else
        fp = 1./(sqrt(u)*2.)        
    endif
    elseif(iexample.eq.6)then
        fp = 1.
    elseif(iexample.eq.7)then
        fp = sin(x)
    else
        fp = u
    endif
    
    end function fp
!*******************************

    real function f(u,x)
    implicit none
    real,intent(in) :: u,x

    if(iexample.eq.4)then
        f = u**2./(u**2.+(1.-u)**2.)
    
    elseif(iexample.eq.5)then
        if(u .le. 0.01)then
            f = 10.*u
        else
            f = sqrt(u)        
        endif
    elseif(iexample.eq.6)then
         f = u
    elseif(iexample.eq.7)then
         f = sin(x)*u
    else
        f = 0.5*u**2.
    endif
    
    end function f
!************************************
    
!***********************************************************
    subroutine flux_comput(p,fleft_out,fright_out)
    implicit none
    real,intent(out) :: fleft_out,fright_out
    type(element1d_downstream), pointer :: p(:)
    integer :: i0
    real ::  a3m,b3m,c3m
    real :: a3p,b3p,c3p
    real :: fleft_minus,fleft_plus
    real :: fright_minus,fright_plus
    real :: xi,lambda
    real :: fleft, fright
    
    
!F_{j-1/2}  
    
    i0 = 3  
    call  weno3(p,i0,0.5,fleft_minus)

    i0 = 4
    call  weno3(p,i0,-0.5,fleft_plus)
    
    call lambda_local( fleft_minus,fleft_plus, p(4)%point_left%char, p(4)%point_left%coor,lambda  )
    
    fleft = 0.5*(f(fleft_minus,p(4)%point_left%coor) - p(4)%point_left%char *fleft_minus &
        & + f(fleft_plus,p(4)%point_left%coor)-p(4)%point_left%char *fleft_plus) - 0.5*lambda*(fleft_plus-fleft_minus )
    
    ! F_{j+1/2}
    i0 = 4
        call  weno3(p,i0,0.5,fright_minus)

    i0 = 5
    call  weno3(p,i0,-0.5,fright_plus)
    
    call lambda_local( fright_minus,fright_plus, p(4)%point_right%char, p(4)%point_right%coor,lambda)

    fright = 0.5*(f(fright_minus,p(4)%point_right%coor)-p(4)%point_right%char *fright_minus &
         & + f(fright_plus,p(4)%point_right%coor)-p(4)%point_right%char*fright_plus)-0.5*lambda*(fright_plus-fright_minus)
    
    fleft_out = fleft
    fright_out = fright
    
    
    end subroutine flux_comput
!************************************************************************
    subroutine weno3(p,i0,xi,f_out)
    implicit none
    real,intent(in) :: xi
    integer,intent(in) :: i0
    real,intent(out) :: f_out
    type(element1d_downstream), pointer :: p(:)
    real :: avg_l2,avg_l1,avg_c,avg_r1,avg_r2
    real :: a3,b3,c3
    real :: ptl2_mid,ptl1_mid, ptr1_mid,ptr2_mid,ptc_mid
    real ::  avg_m, avg_l, avg_r
    real :: lenl1,lenl2,len3, lenr1,lenr2
    
    lenl1 = p(i0-1)%point_right%coor - p(i0-1)%point_left%coor
    lenr1 = p(i0+1)%point_right%coor - p(i0+1)%point_left%coor
    len3 = p(i0)%point_right%coor - p(i0)%point_left%coor
    lenl2 = p(i0-2)%point_right%coor - p(i0-2)%point_left%coor
    lenr2 = p(i0+2)%point_right%coor - p(i0+2)%point_left%coor
    
     avg_l2 = p(i0-2)%up_itg/lenl2
     avg_l1 = p(i0-1)%up_itg/lenl1
     avg_c = p(i0)%up_itg/len3
     avg_r1 = p(i0+1)%up_itg/lenr1
     avg_r2 = p(i0+2)%up_itg/lenr2
     
    ptl2_mid = (p(i0-2)%point_right%coor + p(i0-2)%point_left%coor)/2.
    ptl1_mid = (p(i0-1)%point_right%coor + p(i0-1)%point_left%coor)/2.
    ptc_mid = (p(i0)%point_right%coor + p(i0)%point_left%coor)/2.
    ptr1_mid = (p(i0+1)%point_right%coor + p(i0+1)%point_left%coor)/2.
    ptr2_mid = (p(i0+2)%point_right%coor + p(i0+2)%point_left%coor)/2.
   
!ui-2, ui-1, ui ************************************   
    if(abs((avg_c-avg_l1)/(ptc_mid-ptl1_mid)) .le. abs((avg_c-avg_r1)/(ptc_mid-ptr1_mid))&
    & .and. abs(((avg_l2-avg_l1)/(ptl2_mid-ptl1_mid)-(avg_l1-avg_c)/(ptl1_mid-ptc_mid))/(ptl2_mid-ptc_mid)) &
    & .le. abs(((avg_l1-avg_c)/(ptl1_mid-ptc_mid)-(avg_c-avg_r1)/(ptc_mid-ptr1_mid))/(ptl1_mid-ptr1_mid)))then
!ui-2, ui-1, ui
    avg_l = p(i0-2)%up_itg/len3
    avg_m = p(i0-1)%up_itg/len3   
    avg_r = p(i0)%up_itg/len3  
    
    a3 = (3.*lenl1*len3**3.*(lenl1+len3)*avg_l-3.*lenl2*len3**3.*(lenl2+2.*lenl1+len3)*avg_m &
    & +3.*lenl2*lenl1*(lenl2+lenl1)*len3**2.*avg_r)/(lenl2*lenl1*(lenl2+lenl1)*(lenl1+len3)*(lenl2+lenl1+len3))
    
    b3 = (lenl1*len3**2.*(2.*lenl1**2.+3.*lenl1*len3+len3**2.)*avg_l-lenl2*len3**2.*(2.*lenl2**2.+6.*lenl1**2.+6.*lenl1*len3 &
    &+len3**2.+3.*lenl2*(2.*lenl1+len3))*avg_m+lenl2*lenl1*(lenl2+lenl1)*len3*(2.*lenl2+4.*lenl1+3.*len3)*avg_r) &
    & /(lenl2*lenl1*(lenl2+lenl1)*(lenl1+len3)*(lenl2+lenl1+len3))
    
    c3 = (-1.*lenl1*len3**3.*(lenl1+len3)*avg_l+lenl2*len3**3.*(lenl2+2.*lenl1+len3)*avg_m &
    & +lenl2*lenl1*(4.*lenl2**2.*(lenl1+len3)+lenl1*(4.*lenl1**2.+8.*lenl1*len3+3.*len3**2.) &
    & +lenl2*(8.*lenl1**2.+12.*lenl1*len3+3.*len3**2.))*avg_r) &
    & /(4.*(lenl2*lenl1*(lenl2+lenl1)*(lenl1+len3)*(lenl2+lenl1+len3)))
    
   f_out = a3*xi**2. + b3*xi + c3

!ui, ui+1, ui+2 *****************************************************   
    elseif(abs((avg_c-avg_l1)/(ptc_mid-ptl1_mid)) .ge. abs((avg_c-avg_r1)/(ptc_mid-ptr1_mid))&
    & .and. abs(((avg_l1-avg_c)/(ptl1_mid-ptc_mid)-(avg_r1-avg_c)/(ptr1_mid-ptc_mid))/(ptl1_mid-ptr1_mid)) &
    & .ge. abs(((avg_c-avg_r1)/(ptc_mid-ptr1_mid)-(avg_r1-avg_r2)/(ptr1_mid-ptr2_mid))/(ptc_mid-ptr2_mid)))then
  
    avg_l = p(i0)%up_itg/len3
    avg_m = p(i0+1)%up_itg/len3   
    avg_r = p(i0+2)%up_itg/len3 
    
    a3 = (3.*len3**2.*lenr1*lenr2*(lenr1+lenr2)*avg_l-3.*len3**3.*lenr2*(len3+2.*lenr1+lenr2)*avg_m &
    & +3.*len3**3.*lenr1*(len3+lenr1)*avg_r)/(lenr1*(len3+lenr1)*lenr2*(lenr1+lenr2)*(len3+lenr1+lenr2))
    
    b3 = (-1.*len3*lenr1*lenr2*(lenr1 + lenr2)*(3.*len3+4.*lenr1+2.*lenr2)*avg_l &
    & +len3**2.*lenr2*(len3**2.+6.*len3*lenr1+6.*lenr1**2.+3.*len3*lenr2+6.*lenr1*lenr2 & 
    & +2.*lenr2**2.)*avg_m-len3**2.*lenr1*(len3**2.+3.*len3*lenr1+2.*lenr1**2.)*avg_r) &
    & /(lenr1*(len3+lenr1)*lenr2*(lenr1+lenr2)*(len3+lenr1+lenr2))
    
    c3 = (lenr1*lenr2*(lenr1+lenr2)*(3.*len3**2.+4.*lenr1*(lenr1+lenr2)+4.*len3*(2.*lenr1+lenr2))*avg_l &
    & +len3**3.*lenr2*(len3+2.*lenr1+lenr2)*avg_m-len3**3.*lenr1*(len3+lenr1)*avg_r) &
    & /(lenr1*(len3+lenr1)*lenr2*(lenr1+lenr2)*(len3+lenr1+lenr2)*4.)
    
    f_out = a3*xi**2. + b3*xi + c3

    else
!ui-1,ui,ui+1**************************************
    avg_l = p(i0-1)%up_itg/len3
    avg_m = p(i0)%up_itg/len3   
    avg_r = p(i0+1)%up_itg/len3  
    
    
   a3 = (3.*len3**3.*lenr1*(len3+lenr1)*avg_l &
        & -3.*lenl1*len3**2.*lenr1*(lenl1+2.*len3 &
        & +lenr1)*avg_m &
        & +3.*lenl1*len3**3.*(lenl1+len3)*avg_r) &
        & /(lenl1*(lenl1+len3)*lenr1*(len3+lenr1) &
        & *(lenl1+len3+lenr1))

    b3 = (-1.*len3**2.*lenr1*(len3**2.+3.*len3*lenr1+2.*lenr1**2.) &
        & *avg_l-lenl1*len3*(lenl1-lenr1)*lenr1*(2.*lenl1 &
        & +3.*len3+2.*lenr1)*avg_m +lenl1*len3**2.*(2.*lenl1**2.&
        & +3.*lenl1*len3+len3**2.)*avg_r)/(lenl1*(lenl1 &
        & +len3)*lenr1*(len3+lenr1)*(lenl1+len3+lenr1))

    c3 = (-1.*len3**3.*lenr1*(len3+lenr1)*avg_l &
        & +lenl1*lenr1*(4.*lenl1**2.*(len3+lenr1)+lenl1 &
        & *(3.*len3+2.*lenr1)**2.+len3*(6.*len3**2.+9.*len3*lenr1 &
        & +4.*lenr1**2.))*avg_m-lenl1*len3**3.*(lenl1 &
        & +len3)*avg_r)/(4.*lenl1*(lenl1+len3) &
        & *lenr1*(len3+lenr1)*(lenl1+len3+lenr1))
    f_out = a3*xi**2. + b3*xi + c3
    endif
    
    end subroutine weno3
    
!************************************************************************************
 subroutine weno3_downstream_integral(p,i0,pt_left,pt_right,f_out)
    implicit none
    real,intent(in) :: pt_left,pt_right
    integer,intent(in) :: i0
    real,intent(out) :: f_out
    type(element1d_downstream), pointer :: p(:)
    real :: avg_l2,avg_l1,avg_c,avg_r1,avg_r2
    real :: a3,b3,c3
    real :: ptl2_mid,ptl1_mid, ptr1_mid,ptr2_mid,ptc_mid
    real ::  avg_m, avg_l, avg_r
    real :: lenl1,lenl2,len3, lenr1,lenr2
    real ::  kesi_left,kesi_right
     
    lenl1 = p(i0-1)%point_right%coor - p(i0-1)%point_left%coor
    lenr1 = p(i0+1)%point_right%coor - p(i0+1)%point_left%coor
    len3 = p(i0)%point_right%coor - p(i0)%point_left%coor
    lenl2 = p(i0-2)%point_right%coor - p(i0-2)%point_left%coor
    lenr2 = p(i0+2)%point_right%coor - p(i0+2)%point_left%coor
    
     avg_l2 = p(i0-2)%up_itg/lenl2
     avg_l1 = p(i0-1)%up_itg/lenl1
     avg_c = p(i0)%up_itg/len3
     avg_r1 = p(i0+1)%up_itg/lenr1
     avg_r2 = p(i0+2)%up_itg/lenr2
     
    ptl2_mid = (p(i0-2)%point_right%coor + p(i0-2)%point_left%coor)/2.
    ptl1_mid = (p(i0-1)%point_right%coor + p(i0-1)%point_left%coor)/2.
    ptc_mid = (p(i0)%point_right%coor + p(i0)%point_left%coor)/2.
    ptr1_mid = (p(i0+1)%point_right%coor + p(i0+1)%point_left%coor)/2.
    ptr2_mid = (p(i0+2)%point_right%coor + p(i0+2)%point_left%coor)/2.
    
    kesi_left = (pt_left - ptc_mid)/len3
    kesi_right = (pt_right - ptc_mid)/len3
    
!ui-2, ui-1, ui ************************************   
    if(abs((avg_c-avg_l1)/(ptc_mid-ptl1_mid)) .le. abs((avg_c-avg_r1)/(ptc_mid-ptr1_mid))&
    & .and. abs(((avg_l2-avg_l1)/(ptl2_mid-ptl1_mid)-(avg_l1-avg_c)/(ptl1_mid-ptc_mid))/(ptl2_mid-ptc_mid)) &
    & .le. abs(((avg_l1-avg_c)/(ptl1_mid-ptc_mid)-(avg_c-avg_r1)/(ptc_mid-ptr1_mid))/(ptl1_mid-ptr1_mid)))then
!ui-2, ui-1, ui
    avg_l = p(i0-2)%up_itg/len3
    avg_m = p(i0-1)%up_itg/len3   
    avg_r = p(i0)%up_itg/len3  
    
    a3 = (3.*lenl1*len3**3.*(lenl1+len3)*avg_l-3.*lenl2*len3**3.*(lenl2+2.*lenl1+len3)*avg_m &
    & +3.*lenl2*lenl1*(lenl2+lenl1)*len3**2.*avg_r)/(lenl2*lenl1*(lenl2+lenl1)*(lenl1+len3)*(lenl2+lenl1+len3))
    
    b3 = (lenl1*len3**2.*(2.*lenl1**2.+3.*lenl1*len3+len3**2.)*avg_l-lenl2*len3**2.*(2.*lenl2**2.+6.*lenl1**2.+6.*lenl1*len3 &
    &+len3**2.+3.*lenl2*(2.*lenl1+len3))*avg_m+lenl2*lenl1*(lenl2+lenl1)*len3*(2.*lenl2+4.*lenl1+3.*len3)*avg_r) &
    & /(lenl2*lenl1*(lenl2+lenl1)*(lenl1+len3)*(lenl2+lenl1+len3))
    
    c3 = (-1.*lenl1*len3**3.*(lenl1+len3)*avg_l+lenl2*len3**3.*(lenl2+2.*lenl1+len3)*avg_m &
    & +lenl2*lenl1*(4.*lenl2**2.*(lenl1+len3)+lenl1*(4.*lenl1**2.+8.*lenl1*len3+3.*len3**2.) &
    & +lenl2*(8.*lenl1**2.+12.*lenl1*len3+3.*len3**2.))*avg_r) &
    & /(4.*(lenl2*lenl1*(lenl2+lenl1)*(lenl1+len3)*(lenl2+lenl1+len3)))

   f_out = a3*(kesi_right**3.-kesi_left**3.)/3. + 0.5*b3*(kesi_right**2.-kesi_left**2.) + c3*(kesi_right-kesi_left) 
!ui, ui+1, ui+2 *****************************************************   
    elseif(abs((avg_c-avg_l1)/(ptc_mid-ptl1_mid)) .ge. abs((avg_c-avg_r1)/(ptc_mid-ptr1_mid))&
    & .and. abs(((avg_l1-avg_c)/(ptl1_mid-ptc_mid)-(avg_r1-avg_c)/(ptr1_mid-ptc_mid))/(ptl1_mid-ptr1_mid)) &
    & .ge. abs(((avg_c-avg_r1)/(ptc_mid-ptr1_mid)-(avg_r1-avg_r2)/(ptr1_mid-ptr2_mid))/(ptc_mid-ptr2_mid)))then
  
    avg_l = p(i0)%up_itg/len3
    avg_m = p(i0+1)%up_itg/len3   
    avg_r = p(i0+2)%up_itg/len3 
    
    a3 = (3.*len3**2.*lenr1*lenr2*(lenr1+lenr2)*avg_l-3.*len3**3.*lenr2*(len3+2.*lenr1+lenr2)*avg_m &
    & +3.*len3**3.*lenr1*(len3+lenr1)*avg_r)/(lenr1*(len3+lenr1)*lenr2*(lenr1+lenr2)*(len3+lenr1+lenr2))
    
    b3 = (-1.*len3*lenr1*lenr2*(lenr1 + lenr2)*(3.*len3+4.*lenr1+2.*lenr2)*avg_l &
    & +len3**2.*lenr2*(len3**2.+6.*len3*lenr1+6.*lenr1**2.+3.*len3*lenr2+6.*lenr1*lenr2 & 
    & +2.*lenr2**2.)*avg_m-len3**2.*lenr1*(len3**2.+3.*len3*lenr1+2.*lenr1**2.)*avg_r) &
    & /(lenr1*(len3+lenr1)*lenr2*(lenr1+lenr2)*(len3+lenr1+lenr2))
    
    c3 = (lenr1*lenr2*(lenr1+lenr2)*(3.*len3**2.+4.*lenr1*(lenr1+lenr2)+4.*len3*(2.*lenr1+lenr2))*avg_l &
    & +len3**3.*lenr2*(len3+2.*lenr1+lenr2)*avg_m-len3**3.*lenr1*(len3+lenr1)*avg_r) &
    & /(lenr1*(len3+lenr1)*lenr2*(lenr1+lenr2)*(len3+lenr1+lenr2)*4.)
    
    f_out = a3*(kesi_right**3.-kesi_left**3.)/3. + 0.5*b3*(kesi_right**2.-kesi_left**2.) + c3*(kesi_right-kesi_left) 
    else
!ui-1,ui,ui+1**************************************
    avg_l = p(i0-1)%up_itg/len3
    avg_m = p(i0)%up_itg/len3   
    avg_r = p(i0+1)%up_itg/len3  
    
    
   a3 = (3.*len3**3.*lenr1*(len3+lenr1)*avg_l &
        & -3.*lenl1*len3**2.*lenr1*(lenl1+2.*len3 &
        & +lenr1)*avg_m &
        & +3.*lenl1*len3**3.*(lenl1+len3)*avg_r) &
        & /(lenl1*(lenl1+len3)*lenr1*(len3+lenr1) &
        & *(lenl1+len3+lenr1))

    b3 = (-1.*len3**2.*lenr1*(len3**2.+3.*len3*lenr1+2.*lenr1**2.) &
        & *avg_l-lenl1*len3*(lenl1-lenr1)*lenr1*(2.*lenl1 &
        & +3.*len3+2.*lenr1)*avg_m +lenl1*len3**2.*(2.*lenl1**2.&
        & +3.*lenl1*len3+len3**2.)*avg_r)/(lenl1*(lenl1 &
        & +len3)*lenr1*(len3+lenr1)*(lenl1+len3+lenr1))

    c3 = (-1.*len3**3.*lenr1*(len3+lenr1)*avg_l &
        & +lenl1*lenr1*(4.*lenl1**2.*(len3+lenr1)+lenl1 &
        & *(3.*len3+2.*lenr1)**2.+len3*(6.*len3**2.+9.*len3*lenr1 &
        & +4.*lenr1**2.))*avg_m-lenl1*len3**3.*(lenl1 &
        & +len3)*avg_r)/(4.*lenl1*(lenl1+len3) &
        & *lenr1*(len3+lenr1)*(lenl1+len3+lenr1))
    f_out = a3*(kesi_right**3.-kesi_left**3.)/3.+0.5*b3*(kesi_right**2.-kesi_left**2.)+c3*(kesi_right-kesi_left)

    endif
   
    end subroutine weno3_downstream_integral
!************************************************************************************

!*************************************************************************************
    subroutine lambda_local(minus,plus,gamma,endpt,lambda_out)
    implicit none
    real,intent(in) :: minus,plus,gamma,endpt
    real,intent(out) :: lambda_out
    real :: local
    integer :: i

     local = 0.
        local = max(abs( fp(minus,endpt)-gamma), abs(fp(plus,endpt)-gamma) )

    lambda_out = local

    end subroutine lambda_local
    
!**************************************************************************
    subroutine reproject(result_out)
    implicit none
    real,intent(out) :: result_out(1:nx)
    type(element1d_downstream), pointer :: p1,p2
    type(element1d_downstream), pointer :: p(:)
    integer :: i,j,j_note
    real :: sum, int_left, int_right,int_cell
    integer :: i0
    
    do i =1,nx
    result_out(i) = 0.
    sum = 0.
    j_note = 0
    int_left = 0.
    int_right = 0.
    int_cell = 0.
    
    do j = 2-nd, merged_cell_total+nd
       p1 => element_x_star(j)
       p2 => element_x_star(j-1)
       if (element_x_star(j)%point_left%id .eq. i) then      
       
       if (j_note .eq. 0) then
       j_note = j_note + 1 
       
       p => element_x_star(j-3:j+1)
       i0 = 3
       call weno3_downstream_integral(p,i0,X(p1%point_left%id),p1%point_left%coor,int_left)
       sum = sum + int_left*(p2%int_len)
        
       endif
       
       if (element_x_star(j)%point_right%id .eq. i) then
       
       sum =  sum + p1%up_itg
        
       elseif( element_x_star(j)%point_right%id .gt. i ) then
       
       p => element_x_star(j-2:j+2) 
       i0 = 3
      call weno3_downstream_integral(p,i0,p1%point_left%coor,X(p1%point_left%id+1),int_right)     
       sum = sum + int_right* p1%int_len
               
       endif
       elseif (element_x_star(j)%point_left%id .lt. i .and. element_x_star(j)%point_right%id .gt. i) then
  
         p => element_x_star(j-2:j+2)
       i0 = 3
       call weno3_downstream_integral(p,i0,X(i),X(i+1),int_cell)
       
       sum = sum + int_cell*p1%int_len
        
    endif
    enddo
    result_out(i) = sum/dx
    enddo
    
    end subroutine reproject 
!***************************************************************************
    subroutine order(l)
    implicit none
    
    integer :: kk0,ig
    real :: xrg
    real :: pt1, pt2, pt3, gw1, gw2, gw3
    integer :: l
    integer :: i
    
    !do i=1, nx
    !if(l.eq.nl)then
    !write(1003,*) Xmid(i), u(i,0)
    !endif
    !enddo
    !close(1003)

        pt1 = 0.
        pt2 = sqrt(3./5.)
        pt3 = -pt2
        gw1 = 8./9.
        gw2 = 5./9.
        gw3 = gw2
        
    error1(l)=0.
    error2(l)=0.
    error3(l)=0.
    do kk0=1,nx
    exact_ave(kk0) = 0.
        do ig = 1 , 6
            xrg = Xmid(kk0) + dx * xg(ig)
            !if (iexample .eq. 3 .and. T .lt. 1.) then
            if (iexample .eq. 3) then
                exact_ave(kk0) = exact_ave(kk0) + burgex( xrg,T ) * wg(ig)
            elseif (iexample .eq. 1 .or. iexample .eq. 2 .or. iexample .eq. 6)then
                exact_ave(kk0) = exact_ave(kk0) + exact( xrg,T ) * wg(ig)
            elseif(iexample .eq. 7)then
                exact_ave(kk0) = exact_ave(kk0) + varex( xrg,T ) * wg(ig)
            else
            return
            endif
        enddo
  enddo
        do i=1, nx
        if(l.eq.nl)then
        write(2000, *) Xmid(i), exact_ave(i)
        endif
        enddo
        close(2000)
         
        if (iexample .eq. 3 .and. T .ge. 1) then
        do i=1,nx
            if(Xmid(i) .ge. pi-0.1 .and. Xmid(i) .le. pi+0.1) then
                exact_ave(i) = 0.
                u(i,i_stage) = 0.
            endif 
        enddo
        endif
        do kk0=1,nx
        error1(l)=error1(l)+abs( exact_ave(kk0)-u(kk0,i_stage)   )
        error3(l)=max(error3(l),abs(exact_ave(kk0)-u(kk0,i_stage)  ))
        error2(l)=error2(l)+( exact_ave(kk0)-u(kk0,i_stage)  )**2
    enddo
    error1(l)=error1(l)/nx
    error2(l)=sqrt(error2(l)/nx)
    
        if(l.eq.nl)then
        write(8000, *) cfl/udif*umax,error1(l)
        endif

    end subroutine order
    !*******************************
    real function varex(x,t)
    implicit none
    real,intent(in) :: x,t

    varex = sin( 2.*atan(exp(-t)*tan(x/2.) ) )/sin(x)

    end function varex
    !*******************************
    real function burgex(xx,tt)
    implicit none
    real, intent(in) :: xx,tt
    integer :: nn,ii
    real :: q0,q1,qtemp
    real :: xxx(500),qq(500),xux(500)
    real :: ptl,ptr, ptc
    nn = 400
    
    ptl = -2.
    ptr = 2.
    do while(abs(ptl-ptr) .gt. 0.2 )
    ptc = (ptl+ptr)/2.
    if((ptl-sin(xx-ptl*tt)) * (ptc-sin(xx-ptc*tt)) .lt. 0. )then
       ptr=ptc
    elseif((ptr-sin(xx-ptr*tt)) * (ptc-sin(xx-ptc*tt)) .lt. 0. )then
       ptl=ptc
    elseif((ptc-sin(xx-ptc*tt)) .eq. 0. )then
    burgex = ptc
    return
    else
    Write(*,*) 'No suitable condition at', xx
    endif
    enddo
    q0 = ptc 
    
    ! q1 = q0 - f(q0)/f'(q0)
    q1 = 10000.
    !q0 = 0.
    qtemp = q0
    
    do while(abs(q0-q1)>eps )
   !  do while(abs(q1-sin(xx-q1*tt)) > eps )
        q0 = qtemp
        q1 = q0 - ( q0-sin(xx-q0*tt) )/ ( 1. + tt*cos( xx-q0*tt ))
        qtemp = q1
    enddo
    burgex = q1

    end function burgex
    !*******************************
    subroutine err_ord
    implicit none
    integer :: l  
    
         l = 1
		 write(*,*) l, error1(l)
		 do l=2, nl
		 ord(l) = -log(error1(l)/error1(l-1))/log((l+0.)/(l-1.))
		 write(*,*) l, error1(l), ord(l)
		 enddo
		 write(*,*) 'T=',T,'tt=',tt,'nt=',nt,'dt=',dt
         pause

    end subroutine err_ord
    !******************************
    end program SLRK_FV_1d   

   