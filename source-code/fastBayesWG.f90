module fastBayesWeight
!
! DESCRIPTION:
!   implementation of fBayesB (Meuwissen et al., 2009, GSE)
!
! REVISIONS:
!   06/07/09 wittenburg
!   04/01/10 wittenburg; fix hyperparameter lambda
!   15/02/10 wittenburg; normalisation for calculation of epistatic effects
!   29/06/10 wittenburg; include weights per SNP for estimating
  !    main genetic effects; ToDo -> weights for epistasis???
!   13/04/11 wittenburg; solve nuisance effects via BLUE
!   18/11/13 wittenburg; include credibility interval for testing genetic effects -> Bayes factor deleted: no useful results
!   21/01/14 wittenburg; weight=0 -> exclude SNP
!   07/05/14 wittenburg; estimate gamma; delete gPostVar (its not validated so far)
!   26/01/18 wittenburg; calculation of bayesfactor 
!
  use kindnumber
  use sharedTypes
  use misc
  use cdf_normal_mod
  use inout
  use inverse
  use estimate

  implicit none

  real(rk8),parameter :: epsilon=1d-8

contains

  subroutine ice(y,X,D,n,p,genefreq,weight,desNuis,q)
!!! ICE algorithm of Meuwissen et al., 2009

!!! notations:
!!! ycorr=y-X*gHat-W*bHat

    integer(ik2) :: io,system,status,last
    integer(ik4),intent(in) :: n,p,q
    integer(ik0),allocatable :: sign(:)
    integer(ik4) :: i,j,k,l,iter,dim,nonNull(6),temp(2),datacol
    integer(ikxl) :: counter,col,nEpi
    real(rk8),intent(in) :: y(n),X(n,p),D(n,p),desNuis(n,q)&
         &,genefreq(p),weight(p)
    real(rk8),allocatable :: gtemp(:),diff(:),ycorr(:),Wb(:)&
         &,g(:),z(:),solution(:,:),t(:,:),gammaProposal(:),moe(:),BF(:)
    real(rk8) :: lambda,gamma,sigma2,bigY,alphaL,alphaU&
         &,hetero(p),euklid,sigmaG2,sigmaD2,sigmaEpi2(1:4)&
         &,div,mu,sigmaE2,varcomp(6),precision=1d-8&
         &,squaredWeight,XtXXt(q,n),lower,upper,frac,aMain,aEpi
    character(len=100) :: shell,text
    real(rk8),parameter :: threshold=1d-4 ! for storing epistatic effects

    alphaL=errorlevel/2
    alphaU=1-alphaL

    if(solveEpi) then
       precision=1d-6
       nEpi=int(p*(p-1)/2,ikxl)
    end if
    
    !number of genetic effects
    if(solveDom) then
       dim=2*p
       if(solveEpi) col=int(dim*(dim+1)/2-p,ik4) !-p because no add x dom interaction at the same locus
    else
       dim=p
       if(solveEpi) col=int(dim*(dim+1)/2,ik4)
    end if
    if(.not.solveEpi) col=dim

    allocate(g(col),gtemp(col),diff(col),ycorr(n),z(n),t(col,2),&
         &sign(col),gammaProposal(col),moe(col),BF(col))

    !starting values
    g=0.; sign=0; t(:,1)=-999.; t(:,2)=999.; moe=1.; BF=1.
    gtemp=0.; g=gtemp
    ycorr=y

    sigmaE2=1.
    do j=1,dim
       gammaproposal(j)=nonzeroMain
    end do
    if(solveEpi) then
       do j=(dim+1),col
          gammaproposal(j)=nonzeroEpi
       end do
    end if

    aMain=1/nonZeroMain-2
    aEpi=1/nonZeroEpi-2

    write(*,'(1x,A,2F15.8)') 'START with: gamma', nonZeroMain, nonZeroEpi
    if(.not.GammaFix)  write(*,'(1x,A,2F15.4)') 'hyper-parameter for gamma (a):', aMain, aEpi

    workfile: if(solveNuisance) then
       allocate(Wb(n)); Wb=0.

       shell='cp '//trim(datafile)//' '//trim(corrfile)
       io=system(shell)
       if(io /= 0) write(*,*) 'ERROR copying data file'

       temp=finddim(corrfile,skiplines)
       datacol=temp(2)
       if(datacol<posobs) write(*,*) 'Warning: wrong column number for trait vector.'

    end if workfile

    last=0
    hetero(:)=2*genefreq(:)*(1-genefreq(:))

    iterationen: do iter=1,maxit
       write(*,'(1x,A,I6)')'ITERATION',iter

       !estimate nuisance parameters
       nuisance: if (solveNuisance) then

          forall(i=1:n)
             ycorr(i)=ycorr(i)+Wb(i)
          end forall
   
          if(iter==1) XtXXt=coefMatrix(desNuis,n,q)
          call runBLUE(nuisanceFile,ycorr,desNuis,XtXXt,n,q,Wb)

          forall(i=1:n)
             ycorr(i)=ycorr(i)-Wb(i)
          end forall

       else
          if(iter==1) then
             mu=sum(y)/n
             ycorr(:)=ycorr(:)-mu
          end if

       end if nuisance

       !estimate g_j
       add_et_dom: do j=1,dim
          if(j.le.p) then !additive effects
             if((genefreq(j)<=maf).or.(genefreq(j)>=1-maf)) cycle   
             if(weight(j)<1d-5) cycle
             z(:)=X(:,j)
             squaredWeight=weight(j)
          else !dominance effects
             if((genefreq(j-p)<=maf).or.(genefreq(j-p)>=1-maf)) cycle
             if(weight(j-p)<1d-5) cycle
             z(:)=D(:,j-p)
             squaredWeight=weight(j-p)
          end if
          gamma=gammaproposal(j)
          lambda=sqrt(p*gamma*2)

          !estimate sigma2; e=ycorr
          sigmaE2=dot_product(ycorr,ycorr)/(n-q)
          div=dot_product(z(:),z(:))*squaredWeight
          sigma2=sigmaE2/div

          !determine y_{-j}=y-W*b-X*g[j=0]
          ycorr=ycorr+z(:)*gtemp(j)
          bigY=dot_product(z(:),ycorr)*squaredWeight/div

          !update SNP effect
          gtemp(j)=updateSNP(bigY,lambda,gamma,sigma2)
          ycorr=ycorr-z(:)*gtemp(j)

          !posterior expectation of gamma based on single locus model
          if(.not.GammaFix) gammaProposal(j)=singleGammaEst(bigY,lambda,sigma2,aMain)
       
          !last round: determine 5%-credible interval of g's
           if(last==1) then
              if(j==1) call timestamp()
              if(stepnum>0) then !if credible interval is desired
                 t(j,1)=quantile(probg,bigY,lambda,gamma,sigma2,alphaL,lower,upper,0)
                 t(j,2)=quantile(probg,bigY,lambda,gamma,sigma2,alphaU,lower,upper,0)
              end if
              ! measure of evidence (Pereira & Stern, 1999), but assuming symmetric density p(g|y)
              moe(j)=measureOfEvidence(0d0,gtemp(j),bigY,lambda,gamma,sigma2,lower,upper)
              if(moe(j)<=1-alphaU+alphaL) sign(j)=1
              BF(j)=bayesfactor(z,n,ycorr,gtemp(j),sigmaE2)
           end if

       end do add_et_dom

       counter=dim

       epistatsis: if(solveEpi) then

          do k=1,dim-1
             do l=k+1,dim

                if(l==k+p) cycle
                counter=counter+1
                gamma=gammaproposal(counter)
                lambda=sqrt(p*(p-1)/2*gamma*2)

                choose_couple: if((k.le.p).and.(l.le.p)) then
                   if((genefreq(k)<=maf).or.(genefreq(k)>=1-maf))&
                        & cycle                   
                   if((genefreq(l)<=maf).or.(genefreq(l)>=1-maf))&
                        & cycle
                   if(weight(k)<1d-5) cycle
                   if(weight(l)<1d-5) cycle
                   z(:)=X(:,k)*X(:,l) !1: add x add                   
                   squaredWeight=max(weight(k),weight(l))
                else if((k.le.p).and.(l.gt.p)) then
                   if((genefreq(k)<=maf).or.(genefreq(k)>=1-maf))&
                        & cycle                   
                   if((genefreq(l-p)<=maf).or.(genefreq(l-p)>=1-maf))&
                        & cycle
                   if(weight(k)<1d-5) cycle
                   if(weight(l-p)<1d-5) cycle
                   z(:)=X(:,k)*D(:,l-p) !2: add  x dom
                   squaredWeight=max(weight(k),weight(l-p))
                else
                   if((genefreq(k-p)<=maf).or.(genefreq(k-p)>=1-maf))&
                        & cycle
                   if((genefreq(l-p)<=maf).or.(genefreq(l-p)>=1-maf))&
                        & cycle
                   if(weight(k-p)<1d-5) cycle
                   if(weight(l-p)<1d-5) cycle
                   z(:)=D(:,k-p)*D(:,l-p) !3: dom x dom
                   squaredWeight=max(weight(k-p),weight(l-p))
                end if choose_couple

                !estimate sigma2; e=ycorr
                sigmaE2=dot_product(ycorr,ycorr)/(n-q)
                div=dot_product(z(:),z(:))*squaredWeight
                sigma2=sigmaE2/div

                !determine y_{-j}=y-W*b-X*g[j=0]
                ycorr=ycorr+z(:)*gtemp(counter)
                bigY=dot_product(z(:),ycorr)*squaredWeight/div

                !update SNP effect
                gtemp(counter)=updateSNP(bigY,lambda,gamma,sigma2)
                ycorr=ycorr-z(:)*gtemp(counter)

                if(.not.GammaFix) gammaProposal(counter)=singleGammaEst(bigY,lambda,sigma2,aEpi)

                !last round
                if(last==1) then
                   if(stepnum>0) then
                      t(counter,1)=quantile(probg,bigY,lambda,gamma,sigma2,alphaL,lower,upper,0)
                      t(counter,2)=quantile(probg,bigY,lambda,gamma,sigma2,alphaU,lower,upper,0)
                   end if
                   moe(counter)=measureOfEvidence(0d0,gtemp(counter),bigY,lambda,gamma,sigma2,lower,upper)
                   if(moe(counter)<=1-alphaU+alphaL) sign(counter)=1
                   BF(j)=bayesfactor(z,n,ycorr,gtemp(counter),sigmaE2)
                end if

             end do
          end do

       end if epistatsis

       diff=gtemp-g; euklid=dot_product(diff,diff)/dot_product(gtemp&
            &,gtemp)
       g=gtemp

       !estimate genetic variance components
       varcomp=sigmaEst(g,col,genefreq,p)
      
       !some output variables
       write(*,'(1x,6A15)') &
            'sigmaA2', 'sigmaD2','sigmaAA2','sigmaAD2','sigmaDA2','sigmaDD2'
       write(*,'(1x,6F15.8)') varcomp
       write(*,*) 'norm of residuals',euklid

       if(.not.euklid.gt.0) then
          write(*,*) 'Unexpected ERROR occured.'
          stop
       end if

       if ((euklid<precision).and.(last==1)) exit
       if ((euklid<precision).and.(iter>=9)) then
          last=1
          lower=minval(g)*10.
          upper=maxval(g)*10.
          write(*,'(1x,A,1x,2F10.6)') 'lower and upper bound for credibility interval',lower,upper
       end if
         
       backup: if(mod(iter,100)==0) then
          write(*,*) '%%%%% It is time for a backup. %%%%%'
          open(unit=9,file='backup.tmp',action='write',status='replace'&
               &,iostat=status)
          if(status==0) then
             do j=1,col
                if(abs(g(j))>threshold) write(9,'(1X,I10,F18.10)') j, g(j)
             end do
          end if
          close(9)
       end if backup

    end do iterationen

    !write yhat and residuals
    allocate(solution(n,2))
    forall(i=1:n) 
       solution(i,1)=y(i)-ycorr(i)
       solution(i,2)=ycorr(i)
    end forall
    call writeMatrixCol(solution,n,2,resfile)

    !write genetic effects
    open(unit=3,file=outfile,status='replace',action='write',iostat&
         &=status)
    file_io: if(status==0) then
       write(3,'(1X,9A12)') 'loc_1','loc_2','gPostExp','gamma','measOfEvid','sign',&
            &'cred_t1','cred_t2','bayesfactor'
       write_main: do j=1,dim
          write(3,1011) j,j,g(j),gammaproposal(j),moe(j),sign(j),t(j,1),t(j,2),BF(j)
       end do write_main

1011   format(1X,2I10,3F18.10,1x,I2,1x,3F18.10)

       epi_out: if(solveEpi) then
       	  write(*,1013) 'Estimates of epistatic effects being greater than ',&
       	  & threshold,' are stored.'
1013      format(1X,A,F10.6,A)
          counter=dim
          do k=1,dim-1
             do l=k+1,dim !!! add x add, add x dom (=dom x add -- if l>p>k), dom x dom
                if(l==k+p) cycle
                counter=counter+1
                if(abs(g(counter))>threshold) write(3,1011) k,l,g(counter)&
                     &,gammaproposal(counter),moe(counter),sign(counter)&
                     &,t(counter,1),t(counter,2),BF(counter)
             end do
          end do
       end if epi_out

       write(*,*) 'gHat is written to ', trim(outfile)
    else
       write(*,*) 'ERROR opening ',outfile
    end if file_io
    close(3)


    write(*,*)
    write(*,*) '--- RESULTS ---->'
   
    write(*,'(1X,A,F10.6)') 'Proportion of significant genetic effects:'
    frac=sum(real(sign(1:p),rk8))/p
    write(*,'(A20,F10.6)') 'additive',frac
    if(solveDom) then
       frac=sum(real(sign((p+1):dim),rk8))/p
       write(*,'(A20,F10.6)') 'dominance',frac
    end if
    if(solveEpi) then
       frac=sum(real(sign((dim+1):counter),rk8))/nEpi/4
       write(*,'(A20,F10.6)') 'epistatic',frac
    end if
    write(*,*)


    global_mu: if(.not.solveNuisance) then
       write(*,'(1X,A25,F15.6)') 'global mean',mu
       write(*,*)
    end if global_mu

    write(*,'(1x,A)') '*** Genetic variance components are calculated as plug-in estimates under LE ***'
    write(*,'(1x,A)') '*** More realistic Var(genetic values) should be calculated elsewhere ***'
    select case (reparamMethod) 
    case (1,2,3)
       write(*,*) 
    case default
       write(*,*) '*** Epistatic variance components were not calculated! ***'
       write(*,*) 
    end select

    write(*,'(1x,A25,F15.6)') 'total genetic variance',sum(varcomp)
    write(*,'(1X,A25,F15.6)') 'residual variance',sigmaE2
    write(*,*)

    finish: if(iter.ge.maxit) then
       write(*,*) 'Warning: ice algorithm could not converge within'&
            &,iter-1,'steps.'
    else
       write(*,*) 'Finished: ice algorithm converged within',iter,'ste&
            &ps.'
    end if finish


    deallocate(g,gtemp,ycorr,diff,z,solution,t,sign,gammaProposal,moe,BF)
    if(solveNuisance) then
       deallocate(Wb)
       call deleteOld(corrfile)
    end if

  end subroutine ice


  function updateSNP(bigY,lambda,gamma,sigma2)
!!! calculate estimate of SNP effect as mean of its posterior distribution

    real(rk8),intent(in) :: bigY,lambda,gamma,sigma2
    real(rk8) :: bigYp,bigYm,thetaL,thetaU,term(4),num,denom,lambdaY&
         &,updateSNP
    real(rk8), parameter :: maxValue=50.

    term=0.

    bigYp=bigY+lambda*sigma2
    bigYm=bigY-lambda*sigma2
    lambdaY=lambda*bigY
    if(lambdaY>maxValue) lambdaY=maxValue !numerical stabilisation
    if(lambdaY<-maxValue) lambdaY=-maxValue
    thetaL=prop_exp_trunc(bigYp,sigma2,0d0,.true.)
    thetaU=prop_exp_trunc(bigYm,sigma2,0d0,.false.)
    term(1)=exp(-lambdaY)
    term(2)=1-cum_normal(0d0,bigYm,sqrt(sigma2))
    term(3)=1/term(1)
    term(4)=cum_normal(0d0,bigYp,sqrt(sigma2))
    num=term(1)*thetaU+term(3)*thetaL
    denom=term(1)*term(2)+term(3)*term(4)+2.*(1-gamma)/&
         &(gamma*lambda)*exp(-0.5*lambda**2*sigma2)*&
         &normalpdf(bigY,0d0,sqrt(sigma2))
    updateSNP=num/denom

    if(.not.abs(updateSNP)>0) updateSNP=0.
   
  end function updateSNP


  function postVar(bigY,lambda,gamma,sigma2)
!!! if required: postVar calculates the posterior variance of a single SNP effect
!!! is not validated yet

    real(rk8),intent(in) :: bigY,lambda,gamma,sigma2
    real(rk8) :: postVar,bigYp,bigYm,denom,num,expect,term(6)

    expect=updateSNP(bigY,lambda,gamma,sigma2)
    term=0.

    bigYp=bigY+lambda*sigma2
    bigYm=bigY-lambda*sigma2
    term(1)=exp(-lambda*bigY)
    term(2)=1-cum_normal(0d0,bigYm,sqrt(sigma2))
    term(3)=secMoment(bigYm,sqrt(sigma2),0d0,100000d0)! upper
    ! truncated normal
    term(4)=cum_normal(0d0,bigYp,sqrt(sigma2))
    term(5)=secMoment(bigYp,sqrt(sigma2),-100000d0,0d0)
    num=term(1)*term(2)*term(3)+1/term(1)*term(4)*term(5)
    term(6)=(1-gamma)*2./gamma/lambda*exp(-0.5*lambda**2*sigma2)
    denom=term(1)*term(2)+1/term(1)*term(4)+term(6)*normalpdf(bigY&
         &,0d0,sqrt(sigma2))
    postVar=num/denom-expect**2

  end function postVar


  function Eg_y(y,lambda,Pqtl,ve)
  !!! function from Meuwissen's gwblupcgm.f90

    real(rk8),intent(in) :: y,lambda,Pqtl,ve
    real(rk8) :: noemer,lambday,Eg_y,term,teller,phi,phi2,si,si2,phis,phis2

    phi=cum_normal((y-lambda*ve)/sqrt(ve));
    phi2=cum_normal(-(y+lambda*ve)/sqrt(ve));
    phis=exp(-0.5*(y-lambda*ve)**2/ve)/sqrt(2*3.1415927*ve);
    phis2=exp(-0.5*(y+lambda*ve)**2/ve)/sqrt(2*3.1415927*ve);
    if(phi>0.0)then
       si=phis/phi
    else
       si=10.0 !high intensity; irrelevant what is filled in here, since multiplied by PHI
    end if
    if(phi2>0.0)then
       si2=phis2/phi2
    else
       si2=10.0 !high intensity; irrelevant what is filled in here, since multiplied by PHI
    end if

    lambday=lambda*y  !perhaps some numerical stabilisation needed some
    if(abs(lambday)<1.e-20)lambday=0.0
    if(lambday>50.)lambday=50.0
    if(lambday<-50.)lambday=-50.0

    term=Pqtl*lambda/2.
    teller=exp(-lambday)*phi*(y-lambda*ve+sqrt(ve)*si);
    teller=teller+exp(lambday)*phi2*(y+lambda*ve-sqrt(ve)*si2);
    teller=teller*term
    noemer=term*(exp(-lambday)*phi+exp(lambday)*phi2)  &
         +exp(-0.5*lambda**2*ve)*(1.-pqtl)*exp(-0.5*y**2/ve)/sqrt(2*3.1415927*ve);
    Eg_y=teller/noemer;

  end function Eg_y

  
  subroutine deleteOld(workfile,extension)

    character(len=100),intent(in) :: workfile
    character(len=5),intent(in),optional :: extension
    character(len=100) :: shell,delfile
    integer(ik2) :: io,system
    logical :: ex

    if(present(extension)) then
       delfile=trim(workfile)//trim(extension)
    else
       delfile=trim(workfile)
    end if

    inquire(file=delfile,exist=ex)

    if(ex) then
       shell='rm '//trim(delfile)
       io=system(shell)
       if(io/=0) write(*,*) 'ERROR: deleting ',trim(delfile)
    else
       write(*,*) 'Warning: file not found: ',trim(delfile)
    end if

  end subroutine deleteOld


  subroutine runBLUE(workfile,y,W,XXX,rows,cols,Wb)

    integer(ik4),intent(in) :: rows,cols
    real(rk8),intent(in) :: y(rows),W(rows,cols),XXX(cols,rows)
    real(rk8) :: effects(cols)
    real(rk8),intent(out) :: Wb(rows)
    character(len=100),intent(in) :: workfile

    !least square estimate
    effects=matmul(XXX,y)
    Wb=matmul(W,effects)

    call writeVektor(effects,cols,workfile,.false.)

  end subroutine runBLUE


  function coefMatrix(W,rows,cols)

    integer(ik4),intent(in) :: rows,cols
    integer(ik4) :: io,temp(2),lang,breit,nrang
    integer(ik4),allocatable :: iflag(:)
    real(rk8),intent(in) :: W(rows,cols)
    real(rk8) :: det,zero=1d-10
    real(rk8),allocatable :: XtX(:,:),matVec(:),Vvec(:),Wvec(:)
    real(rk8) :: coefMatrix(cols,rows)

    allocate(XtX(cols,cols))

    XtX=matmul(transpose(W),W)
    temp=shape(XtX)
    breit=temp(1)
    lang=int(breit*(breit+1)/2)

    ! coefMatrix = (X'X)^{-1}X'
    allocate(matVec(lang))
    call makevector(XtX,breit,matVec,lang)
    allocate(Vvec(breit),Wvec(breit),iflag(breit))
    call DKMWHF(matVec,Vvec,Wvec,det,zero,iflag,nrang,breit,lang,0)
    call makematrix(matVec,breit,XtX,lang)

    coefMatrix=matmul(XtX,transpose(W))

    deallocate(XtX,matVec,Vvec,Wvec,iflag)

  end function coefMatrix


  function probg(x,bigY,lambda,gamma,sigma2)
!!! calculate posterior probability of SNP effect p(g|Y)

    real(rk8),intent(in) :: x,bigY,lambda,gamma,sigma2
    real(rk8) :: probg,bigYp,bigYm,lambdaY,term(3),denom,num,skal(3)
    real(rk8),parameter :: maxValue=50.

    term=0.
    skal=0.

    bigYp=bigY+lambda*sigma2
    bigYm=bigY-lambda*sigma2
    lambdaY=lambda*bigY
    if(lambdaY>maxValue) lambdaY=maxValue !numerical stabilisation
    if(lambdaY<-maxValue) lambdaY=-maxValue

    skal(1)=0.5*gamma*lambda*exp(0.5*lambda**2*sigma2+lambdaY)
    term(1)=skal(1)*cum_normal(0d0,bigYp,sqrt(sigma2))
    term(2)=(1-gamma)*normalpdf(bigY,0d0,sqrt(sigma2))
    skal(3)=0.5*gamma*lambda*exp(0.5*lambda**2*sigma2-lambdaY)
    term(3)=skal(3)*(1-cum_normal(0d0,bigYm,sqrt(sigma2)))
    denom=sum(term)

    if(x<=0) then
       num=skal(1)*cum_normal(x,bigYp,sqrt(sigma2))
    else
       num=term(1)+term(2)+skal(3)*(cum_normal(x,bigYm,sqrt(sigma2))&
            &-cum_normal(0d0,bigYm,sqrt(sigma2)))
    end if

    probg=num/denom
    if(.not.abs(probg)>0) probg=0.

  end function probg


  function quantile(func,bigY,lambda,gamma,sigma2,alpha,lower,upper,iprint)
    !sequencial search for quantile

    real(rk8),intent(in) :: alpha,lower,upper,bigY,lambda,gamma,sigma2
    integer(ik2),intent(in) :: iprint
    integer(ik4) :: i
    real(rk8) :: step,quantile,yval,xval,minx,miny
    real(rk8),external :: func

    step=(upper-lower)/stepnum

    do i=1,stepnum
       xval=lower+(i-1)*step
       yval=func(xval,bigY,lambda,gamma,sigma2)
       if(yval.gt.alpha) exit
       miny=yval
       minx=xval
    end do

    if(iprint>0) then
       print*,'stepsize',step
       print*,'quantile found after',i,'cycles'
       print*,'quantile at',minx,miny
    end if
    quantile=minx

  end function quantile


  function measureOfEvidence(x,gEst,bigY,lambda,gamma,sigma2,lower,upper)
!!! Pereira & Stern, Entropy, 1999
!!! x=theta0 value under the null hypothesis
    
    real(rk8),intent(in) :: x,gEst,bigY,lambda,gamma,sigma2,lower,upper
    integer(ik2) :: i,nstep=2000
    real(rk8) :: bigYp,bigYm,lambdaY,measureOfEvidence,phi,step,xval,yval,kappa,x0
    real(rk8),parameter :: maxValue=50.

    bigYp=bigY+lambda*sigma2
    bigYm=bigY-lambda*sigma2
    lambdaY=lambda*bigY
    if(lambdaY>maxValue) lambdaY=maxValue !numerical stabilisation
    if(lambdaY<-maxValue) lambdaY=-maxValue

    ! posterior density p(g|Y)
    phi=density(x)
      
    step=(upper-lower)/nstep

    do i=1,nstep
       if(gEst>0) then
          xval=lower+(i-1)*step
          yval=density(xval)
          if((yval.le.phi).and.(xval.gt.x).and.(xval.ge.gEst)) then
             x0=x+epsilon ! calculate mass of *OPEN* interval
             kappa=probg(xval,bigY,lambda,gamma,sigma2)-probg(x0,bigY,lambda,gamma,sigma2) 
             exit
          end if
       else
          xval=upper-(i-1)*step
          yval=density(xval)
          if((yval.le.phi).and.(xval.lt.x).and.(xval.le.gEst)) then
             xval=xval+epsilon ! calculate mass of *OPEN* interval
             kappa=probg(x,bigY,lambda,gamma,sigma2)-probg(xval,bigY,lambda,gamma,sigma2) 
             exit
          end if
       end if
    end do

    measureOfEvidence=1-kappa !bayesian analogue of p-value

  contains

    function density(x)
      ! posterior density p(g|Y)

      real(rk8),intent(in) :: x
      real(rk8) :: density,num,denom,skal(3),term(3)

      bigYp=bigY+lambda*sigma2
      bigYm=bigY-lambda*sigma2
      lambdaY=lambda*bigY
      if(lambdaY>maxValue) lambdaY=maxValue !numerical stabilisation
      if(lambdaY<-maxValue) lambdaY=-maxValue

      skal(1)=0.5*gamma*lambda*exp(0.5*lambda**2*sigma2+lambdaY)
      term(1)=skal(1)*cum_normal(0d0,bigYp,sqrt(sigma2))
      term(2)=(1-gamma)*normalpdf(bigY,0d0,sqrt(sigma2))
      skal(3)=0.5*gamma*lambda*exp(0.5*lambda**2*sigma2-lambdaY)
      term(3)=skal(3)*(1-cum_normal(0d0,bigYm,sqrt(sigma2)))
      denom=sum(term)

      if(x<0d0) then
         num=skal(1)*normalpdf(x,bigYp,sqrt(sigma2))
      else if (x.eq.0d0) then
         num=(1-gamma)*normalpdf(bigY,0d0,sqrt(sigma2))
      else
         num=skal(3)*normalpdf(x,bigYm,sqrt(sigma2))
      end if
      density=num/denom

    end function density

  end function measureOfEvidence

 
  function bayesfactor(z,n,ycorr,effect,sigmaE2)
!!! calculate the bayes factor for a single effect at a time
    
    integer(ik4),intent(in) :: n
    real(rk8),intent(in) :: z(n),ycorr(n),sigmaE2,effect
    real(rk8) :: bayesfactor,ytemp(n),gtemp(n)

    gtemp=z(:)*effect
    ytemp=ycorr+gtemp
    bayesfactor=exp((2*dot_product(ytemp,gtemp)-dot_product(gtemp,gtemp))/2/sigmaE2)
    if(bayesfactor>1d+6) bayesfactor=99999

  end function bayesfactor


end module fastBayesWeight
