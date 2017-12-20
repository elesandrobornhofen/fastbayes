module misc
!
! DESCRIPTION:
!   collection of miscellaneous functions
!
! REVISIONS:
!   09/06/09 wittenburg
!
  use kindnumber
  use cdf_normal_mod

  implicit none

contains

  function exp_trunc(mean,var,point,lower)
  !expectation of lower (lower=TRUE) truncated normal distribution at truncation point 'point'

    real(rk8),intent(in) :: point
    real(rk8),intent(in) :: mean,var
    logical,intent(in) :: lower
    real(rk8) :: sd,denom,exp_trunc

    sd=sqrt(var)
    lowertrunc: if(lower == .TRUE.) then
       denom=-cum_normal(point,mean,sd)
    else
       denom=1-cum_normal(point,mean,sd)
    end if lowertrunc
    exp_trunc=sd/denom*normalpdf((point-mean)/sd,0d0,1d0)+mean

  end function exp_trunc


  function secMoment(mu,sd,lowerbound,upperbound)
!!! return the second moment of a truncated normal distribution
!!!  (Phoebus J. Dhrymes, 2005)

    real(rk8), intent(in) :: mu, sd, lowerbound, upperbound
    real(rk8) :: secMoment, lower_std, upper_std, term(3),dnormU&
         &,dnormL,pnormU,pnormL

    lower_std=(lowerbound-mu)/sd
    upper_std=(upperbound-mu)/sd

    dnormU=normalpdf(upper_std,0d0,1d0)
    dnormL=normalpdf(lower_std,0d0,1d0)
    pnormU=cum_normal(upper_std,0d0,1d0)
    pnormL=cum_normal(lower_std,0d0,1d0)

    term(1)=mu**2
    term(2)=2*mu*sd*(dnormU-dnormL)/(pnormU-pnormL)
    term(3)=sd**2*(1-(upper_std*dnormU-lower_std*dnormL)/(pnormU&
         &-pnormL))

    secMoment=sum(term)

  end function secMoment


  function normalpdf(x,mean,sd)
  !probability density of normal distribution

    real(rk8),intent(in) :: x
    real(rk8),intent(in) :: mean,sd
    real(rk8) :: normalpdf,var
    real(rk8),parameter :: pi=acos(0.0)*2.

    var=sd**2
    normalpdf=1./sqrt(2.*pi*var)*exp(-(x-mean)**2/2./var)

  end function normalpdf


  function prop_exp_trunc(mean,var,point,lower)
  !expectation of truncated normal distribution at truncation point
    ! 'point' neglecting the factor of proportionality

    integer(ik2) :: factor
    real(rk8),intent(in) :: point
    real(rk8),intent(in) :: mean,var
    logical,intent(in) :: lower
    real(rk8) :: sd,denom,prop_exp_trunc

    sd=sqrt(var)
    if(lower==.true.) then
       factor=-1
       denom=cum_normal(point,mean,sd)
    else
       factor=1
       denom=1-cum_normal(point,mean,sd)
    end if
    prop_exp_trunc=factor*sd*normalpdf((point-mean)/sd,0d0,1d0)+mean*denom

  end function prop_exp_trunc


  subroutine makevector(XTX,spalten,Element,m)
  !Umspeicherung von XTX in Array-Schreibweise für Subroutine DKMWHF

    integer(ik4),intent(in) :: spalten,m
    integer(ik4) :: i,j,h
    real(rk8),dimension(spalten,spalten),intent(in) :: XTX
    real(rk8),dimension(m),intent(out) :: Element

    do i = 1,spalten
       do j = 1,spalten
          if(i .ge. j) then
             h=i+(2*spalten-j)*(j-1)/2
             Element(h)=XTX(i,j)
          end if
       end do
    end do

  end subroutine makevector


  subroutine makematrix(Element,spalten,XTX,m)
  !Speichern der Inversen von XTX in XTX

    integer(ik4),intent(in) :: spalten,m
    integer(ik4) :: i,j,h
    real(rk8),dimension(spalten,spalten),intent(out) :: XTX
    real(rk8),dimension(m),intent(in) :: Element

    do i = 1,spalten
       do j = 1,spalten
          if(i .ge. j) then
             h=i+(2*spalten-j)*(j-1)/2
             XTX(i,j)=Element(h)
             XTX(j,i)=XTX(i,j)
          end if
       end do
    end do

  end subroutine makematrix


  function findIndex(row,column,nLoc,solveAll)
!!! calculate the index for the vector of estimated genetic
!!!  effects corresponding to a given pair of loci
!!! if only additive effects are involved -> solveAll=.false.

    integer(ik4),intent(in) :: row,column,nLoc
    logical,intent(in) :: solveAll
    integer(ik4) :: findIndex(4)

    findIndex=0

    if((row<column).and.(column<=nLoc).and.(.not.solveAll)) then
       !only additive effects are modelled
       findIndex(1)=column-1+(2*(nLoc-1)-row)*(row-1)/2
       findIndex(1)=findIndex(1)+nLoc !nLoc main effects are stored first
    else if ((row<column).and.(column<=nLoc).and.solveAll) then
       !first triangle (additive x additive)
       findIndex(1)=column-1+(2*(nLoc-1)-row)*(row-1)/2+(row-1)*(nLoc-1)
       !upper triangle in upper right corner (additive x dominant)
       findIndex(2)=column-1+(2*(nLoc-1)-row)*(row-1)/2+row*(nLoc-1)
       !lower triangle in upper right corner (dominant x additive)
       findIndex(3)=row+nLoc-1+(2*(2*nLoc-1)-column)*(column-1)/2&
            &-(column-1)
       !second triangle (dominant x dominant)
       findIndex(4)=column+nLoc-1+(2*(2*nLoc-1)-(row+nLoc))*(row+nLoc&
            &-1)/2-nLoc
       findIndex=findIndex+2*nLoc !2*nLoc main effects are stored first
    else
       print*,'CHECK:',row,'<',column,'or',column,'<=',nLoc
    end if

  end function findIndex

end module misc
