module nuisance

  use kindnumber
  use inout

  implicit none

contains

  function designNuisance(n,datafile,skiplines)

    integer(ik2) :: io,regDegree
    integer(ik4),intent(in) :: n
    integer(ik4),parameter :: maxEff=10
    integer(ik4) :: counter,k,i,l,nEff,& ! number of other/nuisance effects
         &nLevA,nLevI,levelsI(n) !number of levels and levels of a factor
    integer(ik4) :: columnDat(maxEff),& !which columns for other/nuisance effects
         &setZero(maxEff),& !which level is set to zero for each nuisance effect
         &factorI(n),index
    real(rk8),allocatable :: WX(:,:)
    real(rk8) :: factorR(n)
    real(rk8),dimension(:,:),pointer :: designNuisance
    character(len=2) :: effectType(maxEff),regr !type of effect of other
    !/nuisance effects (I/A Integer/Alpha factor, R regressor)
    character(len=20) :: factorA(n),levelsA(n)
    character(len=100) :: levelFile='coding_of_level.txt'

    integer(ik4),intent(in) :: skiplines
    integer(ik4) :: dfNuisance
    character(len=100),intent(in) :: datafile
   
    namelist / parameterList / nEff,columnDat,effectType,setZero

    !read in list of parameters
    open(unit=1,file='nuisance.par',status='old',action='read',iostat=io)
    readparameter: if(io == 0) then
       read(1,nml=parameterList,iostat=io)
       write(*,*) 'Parameters are read from nuisance.par'
    else
       write(*,*) 'ERROR reading parameter list'
       stop
    end if readparameter
    close(1)

    print*,'Maximum nuisance effects possible:',maxEff
    if(maxEff<nEff) then
       print*,'If more effects required, change >maxEff< and re-compile!'
       stop
    end if

    allocate(WX(n,n),stat=io)
    if(io/=0) then
       print*,'ERROR: unsufficient memory for nuisance design matrix.'
       stop
    end if

    open(unit=2,file=levelfile,status='unknown',action='write',iostat&
         &=io)
    write(2,'(1x,A8,2A20)') 'factor','level_numeric','level_original'
   
    WX=0.; WX(:,1)=1. !global mean
    write(2,1012) 0,1,'mu'

    counter=1; dfNuisance=0
    do k=1,nEff
       regr=effectType(k)
       select case (regr(1:1))
       case('A')
          call readColumnAlpha(datafile,n,columnDat(k),skiplines&
               &,factorA)
          call countLevelsA(factorA,n,nlevA,levelsA)
          write(*,1011),'Factor Nr',k,'is alphanumeric,',nlevA,'levels'
          !make columns in design matrix
          do i=1,n
             index=findLevA(factorA(i),levelsA,nlevA)
             WX(i,counter+index)=1.
          end do
          if(setZero(k)>0) then
             WX(:,counter+setZero(k))=0.
             dfNuisance=dfNuisance-1
          end if
          counter=counter+nLevA
          writeLevelsA: do i=1,nlevA
             write(2,1012) k,i,levelsA(i)
          end do writeLevelsA
       case('I')
          call readColumnInt(datafile,n,columnDat(k),skiplines&
               &,factorI)
          call countLevels(factorI,n,nlevI,levelsI)
          write(*,1011),'Factor Nr',k,'is integer,',nlevI,'levels'
          !make columns in design matrix
          do i=1,n
             index=findLevI(factorI(i),levelsI,nlevI)
             WX(i,counter+index)=1.
          end do
          if(setZero(k)>0) then
             WX(:,counter+setZero(k))=0.
             dfNuisance=dfNuisance-1
          end if
          counter=counter+nLevI
          writeLevelsI: do i=1,nlevI
             write(2,1013) k,i,levelsI(i)
          end do writeLevelsI
       case('R')
          regDegree=iachar(regr(2:2))-48 !ascii
          call readColumn(datafile,n,columnDat(k),skiplines,factorR)
          if((regDegree>9).or.(regDegree<=0)) then
             regDegree=1
          end if
          write(*,1011),'Factor Nr',k,'is regressor,',regDegree,'degrees'
          do l=1,regDegree
             WX(:,counter+1)=factorR(:)**l
             counter=counter+1
             write(2,1014) k,l,'x**',l
          end do
       case default
          print*,'Factor Nr',k,'has unknown type'
       end select
    end do
1011 format(1x,A,1x,I6,1x,A,1x,I6,1x,A)
1012 format(1x,2I6,1x,A20)
1013 format(1x,2I6,1x,I6)
1014 format(1x,2I6,1x,A3,I1)
    close(2)
    dfNuisance=dfNuisance+counter
    print*,'Degrees of freedom of fixed effects:',dfNuisance

    allocate(designNuisance(n,counter),stat=io)
    if(io/=0) then
       print*,'ERROR: unsufficient memory for nuisance design matrix.'
       stop
    end if
    designNuisance=WX(1:n,1:counter)
    deallocate(WX)

  end function designNuisance


  subroutine countLevels(factor,q,nlevels,levels)
!!! count levels of integer-valued factor

    integer(ik4),intent(in) :: q,factor(q)
    integer(ik4),intent(out) :: nlevels,levels(q)
    integer(ik4) :: i,j,count

    nlevels=1
    levels(1)=factor(1)
    if(q>1) then
       observedFactor: do j=2,q
          count=0
          uniqueLevel: do i=1,nlevels
             if(factor(j)/=levels(i)) count=count+1
          end do uniqueLevel
          if(count==nlevels) then             
             nlevels=nlevels+1
             levels(nlevels)=factor(j)
          end if
       end do observedFactor
    end if
    
  end subroutine countLevels


  subroutine countLevelsA(factor,q,nlevels,levels)
!!! count levels of alphanumeric factor

    integer(ik4),intent(in) :: q
    character(len=20),intent(in) :: factor(q)
    integer(ik4),intent(out) :: nlevels
    character(len=20),intent(out) :: levels(q)
    integer(ik4) :: i,j,count

    nlevels=1
    levels(1)=factor(1)
    if(q>1) then
       observedFactor: do j=2,q
          count=0
          uniqueLevel: do i=1,nlevels
             if(factor(j)/=levels(i)) count=count+1
          end do uniqueLevel
          if(count==nlevels) then             
             nlevels=nlevels+1
             levels(nlevels)=factor(j)
          end if
       end do observedFactor
    end if
    
  end subroutine countLevelsA


  function findLevI(factorValue,levels,nlev)

    integer(ik4),intent(in) :: nLev,levels(nlev),factorValue
    integer(ik4) :: findLevI,l

    findLevI=0
    do l=1,nLev
       if(factorValue==levels(l)) then
          findLevI=l
          exit
       end if
    end do

    if(findLevI==0) print*,'Warning: factor level',factorValue,&
         &'was not found'

  end function findLevI


  function findLevA(factorValue,levels,nlev)

    integer(ik4),intent(in) :: nLev
    character(len=20),intent(in) :: levels(nlev),factorValue
    integer(ik4) :: findLevA,l

    findLevA=0
    do l=1,nLev
       if(factorValue==levels(l)) then
          findLevA=l
          exit
       end if
    end do

    if(findLevA==0) print*,'Warning: factor level',factorValue,&
         &'was not found'

  end function findLevA


end module nuisance
