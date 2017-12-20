module estimate
!
! DESCRIPTION:
!   functions on estimating variance components and proportion of non-zero genetic effects (i.e. gamma)
! REVISIONS:
!  06/07/09 wittenburg
!  12/05/14 wittenburg; add gamma estimation based on single-locus model
!
  use kindnumber
  use sharedTypes
  use misc

  implicit none

contains

  function singleGammaEst(bigY,lambda,sigma2,a)
!!! estimate posterior expectation E(gamma|y) based on 1-locus model

    real(rk8),intent(in) :: bigY,lambda,sigma2,a
    real(rk8) :: singleGammaEst,bigYp,bigYm,lambdaY,term(3),denom,num,skal(3)
    real(rk8),parameter :: maxValue=50.

    term=0.
    skal=0.

    bigYp=bigY+lambda*sigma2
    bigYm=bigY-lambda*sigma2
    lambdaY=lambda*bigY
    if(lambdaY>maxValue) lambdaY=maxValue !numerical stabilisation
    if(lambdaY<-maxValue) lambdaY=-maxValue

    skal(1)=0.5*lambda*exp(0.5*lambda**2*sigma2+lambdaY)
    term(1)=skal(1)*cum_normal(0d0,bigYp,sqrt(sigma2))
    skal(2)=(a+1)/2
    term(2)=skal(2)*normalpdf(bigY,0d0,sqrt(sigma2))
    skal(3)=0.5*lambda*exp(0.5*lambda**2*sigma2-lambdaY)
    term(3)=skal(3)*(1-cum_normal(0d0,bigYm,sqrt(sigma2)))
    num=2/(a+3)*sum(term)

    skal(2)=a+1
    term(2)=skal(2)*normalpdf(bigY,0d0,sqrt(sigma2))
    denom=sum(term)

    singleGammaEst=num/denom
    if(.not.abs(singleGammaEst)>0) singleGammaEst=0.

  end function singleGammaEst


  function sigmaEst(effects,col,freq,p)
!!! estimate additive and dominance variance components
!!! depend on reparametrisation method chosen earlier

    integer(ikxl),intent(in) :: col
    integer(ik4),intent(in) :: p
    real(rk8),dimension(col),intent(in) :: effects
    real(rk8),dimension(p),intent(in) :: freq
    integer(ik2) :: eins(p)
    integer(ik4) :: j,k,find(4)
    real(rk8) :: sigmaEst(6),varAdd,varDom,VarAA,VarAD,VarDA,VarDD&
         &,hetero(p),AddLoc1(3),AddLoc2(3),DomLoc1(3),DomLoc2(3)&
         &,geno(9),nLoc,meanLoc1,meanLoc2,varLoc1,varLoc2,allelsub(p)

    where((freq(:)<=maf).or.(1-freq(:)<=maf))
       eins(:)=0
    elsewhere
       eins(:)=1
    end where

    nLoc=sum(real(eins,rk8))

    VarDom=0.; VarAA=0.; VarAD=0.; VarDD=0.
    hetero(:)=2.*freq(:)*(1.-freq(:))*eins(:)

    select case (reparamMethod)
    case(1,2,3)

!!! since columns of design matrices are standardised, the
!!!  calculation of variance components simplifies

       !total additive genetic variance
       VarAdd=sum(effects(1:p)**2)
      
       !total dominance variance
       if(solveDom) then
          varDom=sum(effects((p+1):(2*p))**2)
          if(reparamMethod==1) VarDom=VarDom/4.
       end if

       !epistatic variance components
       ifEpi: if(solveEpi) then 
          loopLoc1: do j=1,p-1

             loopLoc2: do k=j+1,p
                find=findIndex(j,k,p,solveDom)

                !additive x additive interaction
                VarAA=VarAA+effects(find(1))**2

                DomInterctions: if(solveDom) then
                   !additive x dominance interaction
                   VarAD=VarAD+effects(find(2))**2

                   !dominance x additive interaction
                   VarDA=VarDA+effects(find(3))**2

                   !dominance x dominance interaction
                   VarDD=VarDD+effects(find(4))**2
                end if DomInterctions

             end do loopLoc2
          end do loopLoc1
       end if ifEpi

!!! in case of no reparametrisation
    case default 
       ! additive genetic and dominance variance
       if(solveDom) then
          allelsub=effects(1:p)+effects((p+1):(2*p))*(1-2*freq(:))
          varAdd=sum(hetero(:)*allelsub(:)**2)
          varDom=sum(hetero(:)**2*effects((p+1):(2*p))**2)
       else          
          varAdd=sum(hetero(:)*effects(1:p)**2)
       end if

    end select

    sigmaEst=(/varAdd,varDom,VarAA,VarAD,VarDA,VarDD/)

  end function sigmaEst


end module estimate
