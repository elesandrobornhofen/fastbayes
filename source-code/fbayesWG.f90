program fbayesWeight
!
! DESCRIPTION:
!  estimation of SNP effects using fBayesB, 02/06/09 wittenburg
!
! REVISIONS:
!   05/10/09 wittenburg; include epistatic effects
!   26/04/10 wittenburg; include weighting if SNP positions
!   12/04/11 wittenburg; include estimation of nuisance effects via BLUE
!   18/11/13 wittenburg; improve request for optional weights
!   14/04/14 wittenburg; allow for "no parametrisation" of additive effects
!
  use kindnumber
  use sharedTypes
  use inout
  use fastBayesWeight
  use nuisance
 
  implicit none

  integer(ik2) :: io,io2,system ! control parameter
  integer(ik4) :: k,temp(2),&
       &snp,& ! number of SNP 
       &nRecords,& ! number of individuals with phenotype
       &snpobs,& ! number of individuals with genotype for additive
       !  effects
       &mainEffects,& !number of genetic effects
       &nuisanceEffects !number of nuisance effects
  real(rk8),allocatable ::weight(:),& !vector of weights
       &Genoadd(:,:),Genodom(:,:) ! SNP genotypes
  real(rk8),allocatable :: obs(:),& ! vector of observations
       &genefreq(:),genoProb(:,:) ! allele and genotype frequencies
  real(rk8),dimension(:,:),pointer :: W !design matrix for nuisance effects
  logical :: withWeight
  character(len=100) :: shell,outPrefix

  namelist / input / snpfile,datafile,skiplines,posobs&
       &,weightFile,outPrefix,maxit,stepnum&
       &,reparamMethod,solveNuisance&
       &,solveDom,solveEpi,maf,nonZeroMain,nonZeroEpi&
       &,gammaFix,lambdaFix,writeFreq,errorlevel

  call timeStamp()

  !read in list of global parameters
  open(unit=1, file='inputWG.par',status='old',action='read',iostat=io)
  readparameter: if (io == 0) then
     read(1, nml = input, iostat = io)
     if (io == 0) then
        print '(1XA)', 'Read parameters from inputWG.par'
     else
        print '(1XA)', 'ERROR reading from inputWG.par'
        stop
     end if
     close(1)
  else
     write(*,*) 'ERROR: no input file given.'
  end if readparameter

  nuisanceFile=trim(outprefix)//'.nui'
  outfile=trim(outprefix)//'.eff'
  resFile=trim(outprefix)//'.res'
  logFile=trim(outprefix)//'.log'

  open(unit=99,file=logFile,action='write',status='replace')
  write(99,*) '********************************************'

  !read in phenotypic values
  temp=finddim(datafile,skiplines)
  nRecords=temp(1)
  allocate(obs(nRecords))
  call readColumn(datafile,nRecords,posobs,skiplines,obs)
 
  !read in SNP genotypes (additive)
  temp=finddim(snpfile,0)
  snpobs=temp(1)
  snp=temp(2)
  if(snpobs .lt. nRecords) then
     write(*,*) 'ERROR: too less genotypes'
     stop
  end if
  allocate(Genoadd(nRecords,snp)); Genoadd=0.
  call readin(snpfile,nRecords,snp,0,Genoadd) 
  write(*,*)

  !read in SNP weights
  allocate(weight(snp))
  inquire(file=weightFile,exist=withWeight)
  if(withWeight) then
     temp=finddim(weightFile,0)
     print*,'SNPs with weight<1e-5 are excluded.'
     if(temp(1) .ne. snp) then
        print*,'ERROR in dimension of weights',temp(1)
        stop
     else
        call readin(weightFile,snp,1,0,weight) 
     end if
  else
     write(*,*) 'Warning: fix weight=1'
     weight(1:snp)=1.
  end if
  write(*,*)

  !reparametrise SNP genotypes
  allocate(genefreq(snp),genoProb(snp,3))
  allocate(Genodom(nRecords,snp)); Genodom=0. 
  dominance: if (.not.solveDom) then
     select case (reparamMethod)
     case (1,2,3)
        call reparam_add(Genoadd,nRecords,snp,genefreq,maf,.false.)
        call standardiseGenoMat(Genoadd,Genodom,nRecords,snp,genefreq)
     case default
        call reparam_nothing(Genoadd,Genodom,nRecords,snp&
             &,genefreq,maf,.false.)
        genefreq=0.5 ! *experimental use for regression analysis*
        write(*,1011) 'NOTHING'
     end select
  else
     select case (reparamMethod)
     case (1)
        call reparam_contrast1(Genoadd,Genodom,nRecords,snp&
             &,genefreq,maf,.false.)
        write(*,1011) 'Cockerham, Genetics, 1954'
        call standardiseGenoMat(Genoadd,Genodom,nRecords,snp,genefreq)
     case(2)
        call reparam_contrast2(Genoadd,Genodom,nRecords,snp&
             &,genefreq,maf,genoProb,.false.)
        write(*,1011) 'Alvarez-Castro & Carlborg, Genetics, 2007'
        call standardiseGenoMat(Genoadd,Genodom,nRecords,snp,genefreq)
     case(3)
        call reparam_contrast3(Genoadd,Genodom,nRecords,snp&
             &,genefreq,maf,.false.)
        write(*,1011) 'Zeng et al., Genetics, 2005'
        call standardiseGenoMat(Genoadd,Genodom,nRecords,snp,genefreq)
     case default
        call reparam_nothing(Genoadd,Genodom,nRecords,snp&
             &,genefreq,maf,.false.)
        write(*,1011) 'NOTHING'
     end select
     if(reparamMethod/=2) genoProb=0. !not used yet
  end if dominance
1011 format(1X,'Orthogonalisation of genetic effects according to *&
          &**',A,'***')

  writeAlleleFreq: if(writeFreq) then
     write(*,*) 'Writing allele frequencies...'
     open(5,file='alleleFreq.txt',iostat=io)
     if(reparamMethod==2) then
        write(*,*) 'Writing genotype probabilities...'
        open(91,file='genoProb.txt',iostat=io2)
     end if

     if(io==0) then
        do k=1,snp
           write(5,'(1X,I6,F10.6)') k,genefreq(k)
           if((reparamMethod==2).and.(io2==0)) &
                &write(91,'(1x,I6,3F15.8)') k,genoProb(k,:)
        end do
     else
        write(*,*) 'ERROR writing frequencies.'
     end if
     close(5)
     if(reparamMethod==2) close(91)
     write(*,*)
  end if writeAlleleFreq

  call timeStamp()

  !read in other/nuisance factors
  if(solveNuisance) then
     print*,'Involve nuisance effects...'
     W=>designNuisance(nRecords,datafile,skiplines)
     temp=shape(W)
     nuisanceEffects=temp(2)
     print*,'Total number of nuisance effects',nuisanceEffects
     print*
  else
     nuisanceEffects=1
     allocate(W(nRecords,nuisanceEffects)); W=1.
  end if
 
  !estimate genotypic value
  write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  write(*,*)
  write(*,*) '  Welcome to fBayesB (Meuwissen et al., 2009, GSE)  '
  write(*,*)
  write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'


  call ice(obs,Genoadd,Genodom,nRecords,snp,genefreq,weight,W&
       &,nuisanceEffects)
  call timeStamp() 

 
  deallocate(obs,Genoadd,Genodom,genefreq,genoProb)
 
  close(99)

end program fbayesWeight


