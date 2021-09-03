!module srfc
  !
  ! ----------------------------------------------------------------------
  ! Subroutine srfcscls:  returns scale values based on Businger/Dye
  ! similarity functions.
  !
  subroutine srfcscls(n2,n3,z,u,dth,drt,ustar,tstar,rstar)

    !use mpi_interface, only : myid

    implicit none

    real   :: z0 = 0.001   ! momentum roughness height
    real   :: th00 = 350.1       ! basic state temperature

    character (len=7) :: runtype = 'HISTORY'
    real, parameter :: vonk   = 0.40
    real, parameter :: g      = 9.8
    real, parameter :: R      = 287.04
    real, parameter :: Rm     = 461.5
    real, parameter :: ep2    = Rm/R - 1.

    !real, parameter     :: ah   =  4.7   ! stability function parameter
    !real, parameter     :: bh   = 16.0   !   "          "         "
    !real, parameter     :: am   =  4.7   !   "          "         "
    !real, parameter     :: bm   = 16.0   !   "          "         "
    !real, parameter     :: ah   =  7.8   ! stability function parameter
    !real, parameter     :: bh   = 12.0   !   "          "         "
    !real, parameter     :: am   =  4.8   !   "          "         "
    !real, parameter     :: bm   = 19.3   !   "          "         "
    real, parameter     :: pr   =  0.74  ! prandlt number
    real, parameter     :: eps  = 1.e-10 ! non-zero, small number

    integer, intent(in) :: n2,n3   ! span of indicies covering plane
    real, intent(in)    :: z             ! height where u & T locate
    real, intent(in)    :: u(n2,n3)      ! velocities at z
    real, intent(in)    :: dth(n2,n3)    ! theta (th(z) - th(z0))
    real, intent(inout) :: ustar(n2,n3)  ! scale velocity
    real, intent(inout) :: tstar(n2,n3)  ! scale temperature
    real, intent(in)     :: drt(n2,n3)    ! qt(z) - qt(z0)
    real, intent(inout) :: rstar(n2,n3)  ! scale value of qt

    logical, save :: first_call=.True.
    integer       :: i,j,iterate,iter
    real          :: lnz, klnz, betg
    real          :: zeta, lmo, dtv
    real          :: thv,Rib,Lstart,Lend,Lold,fx,fxdif,Ldif,zeff
    logical       :: exititer

    real, dimension(n2,n3) :: obl    ! Obukhov Length

    real, external :: psim
    real, external :: psih

    lnz   = log(z/z0)
    klnz  = vonk/lnz
    betg  = th00/g


    do j=1,n3
       do i=1,n2
          dtv = dth(i,j) + ep2*th00*drt(i,j)
          !
          ! Neutral case
          ! 
          if (dtv == 0.) then
            ustar(i,j) =  vonk*u(i,j)/lnz
            tstar(i,j) =  vonk*dtv/(pr*lnz)
            lmo        = -1.e10

          !
          ! start iterations from values at previous tstep, 
          ! unless the sign has changed or if it is the first call, then 
          ! use neutral values.
          !
          else

             if ((runtype=='INITIAL' .and. first_call) .or.( tstar(i,j)*dtv <= 0.)) then
               tstar(i,j) =  vonk*dtv/(pr*lnz)
               ustar(i,j) =  vonk*u(i,j)/lnz
               lmo        = -1.e10
             end if

             if(ustar(i,j) == 0) ustar(i,j) = 0.1

             Lold  = 1e9
             Ldif  = 1e9
             iter  = 0
             exititer = .false.
             !do iterate = 1,100
             do while(abs(Ldif)>0.1)
               lmo        = betg*ustar(i,j)**2/(vonk*tstar(i,j))
               Ldif       = lmo - Lold
               Lold       = lmo

               if ((dtv < 0) .and. (lmo > -0.001)) lmo = -0.001999
               if ((dtv > 0) .and. (lmo < +0.001)) lmo = +0.001777

               ! BvS : Following ECMWF, limit z/L for very stable conditions
               if(z/lmo > 5.) then
                 zeff = lmo * 5.
                 exititer = .true.
               else
                 zeff = z
               end if

               zeta       = zeff/lmo
               ustar(i,j) = u(i,j)*vonk/(log(zeff/z0) - psim(zeta))
               if(ustar(i,j)<0.) ustar(i,j) = 0.1
               tstar(i,j) = (dtv*vonk/pr)/(log(zeff/z0) - psih(zeta))

               if(exititer) then
                 lmo        = zeff/5. 
                 exit
               end if

               iter = iter + 1

               ! Limit L for day/night transitions
               if(lmo > 1e6)  lmo = 1e6
               if(lmo < -1e6) lmo = -1e6 

               if(iter>10000) then
                ! print*,'Obukh. length not converged, myid=',myid,'i,j=',i,j
                 stop
               end if
             end do
          end if
 
          obl(i,j) = lmo
          rstar(i,j) = tstar(i,j)*drt(i,j)/(dtv + eps)
          tstar(i,j) = tstar(i,j)*dth(i,j)/(dtv + eps)
       end do
    end do

    first_call = .False.

    return 
  end subroutine srfcscls


  !
  ! ----------------------------------------------------------
  ! Malte: Integrated stability function for momentum
  !
  function psim(zeta)
    implicit none

    real             :: psim
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      x     = (1. - 15. * zeta) ** (0.25)
      !psim = 3.14159265/2. - 2. *atan(x) + log((1.+x)** 2. * (1. + x**2.)/ 8.)
      psim  = 3.14159265/2. - atan(x) + 2.*log((1+x)/2.) + log((1+x*x)/2.)
    else
      psim  = - 4.7 * zeta
    end if

    return 
  end function psim

  !
  ! ----------------------------------------------------------------------
  ! Malte: Integrated stability function for heat
  !
  function psih(zeta)

    implicit none

    real             :: psih
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      x     = (1. - 15. * zeta) ** (0.25)
      psih  = 2. * log( (1. + x ** 2.) / 2. )
    else
      psih  = - 4.7 * zeta
    end if

    return 
  end function psih


