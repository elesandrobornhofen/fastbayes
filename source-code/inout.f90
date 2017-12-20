module inout
!
! DESCRIPTION:
!  read in/write out numerical data
!  reparametrize columns of a matrix
!
! REVISIONS:
!   02/06/09 wittenburg
!   28/03/10 wittenburg; writeMatrix; readinAlpha etc.
!   07/05/14 wittenburg; set output criterion "var<1d-6" in reparam_conrast2 to have some information in logfile
!
  use kindnumber
  implicit none

contains

  function finddim(filename,skip)
  !count number of rows in filename

    integer(ik2) :: status
    integer(ik4) :: i,rows,columns,finddim(2),ilen
    integer(ik4),intent(in) :: skip
    character(len=100),intent(in) :: filename
    character(len=500000) :: inputLine

    rows=0
    open(unit=1,file=filename,status='old',action='read',iostat&
         &=status)

    openif: if(status ==0) then
       skiplines: if(skip>0) then
          do i=1,skip
             read(1,*,iostat=status)
             if(status/=0) exit
          end do
       end if skiplines

       readloop: do
          read(1,'(A)',iostat=status) inputLine
          if (status /=0) exit
          rows=rows+1
          numColumns: if(rows==1) then
             ilen=len_trim(inputLine)
             if (ilen >0) columns=countWords(inputLine,ilen)
          end if numColumns
       end do readloop

       readif: if (status >0) then
          write(*,*) 'ERROR occured in line',rows+1
       else
          write(*,1020) trim(filename),rows,columns
1020      format(1X,'Reading ',A,' was successful. Dimension: ',2I8)
       end if readif

    else openif
       write(*,*) 'ERROR opening ',trim(filename)
    end if openif

    close(unit=1)
    finddim=(/rows,columns/)

  end function finddim


  subroutine readin(filename,rows,columns,skip,output)
  !read in data array of type real(rk8)

    integer(ik2) :: status
    integer(ik4) :: i
    integer(ik4),intent(in) :: rows,columns,skip
    real(rk8),dimension(rows,columns),intent(out) :: output
    character(len=100),intent(in) :: filename

    open(unit=1,file=filename,status='old',action='read',iostat&
         &=status)
    openif: if(status ==0) then
       skiplines: if(skip>0) then
          do i=1,skip
             read(1,*,iostat=status)
          end do
       end if skiplines
       do i=1,rows
          if (status /=0) exit
          read(1,*,iostat=status) output(i,1:columns)
       end do
    else openif
       write(*,*) 'ERROR opening ',trim(filename)
    end if openif

    close(unit=1)

  end subroutine readin


  subroutine readinAlpha(filename,rows,columns,skip,output)
  !read in data array of type character

    integer(ik2) :: status
    integer(ik4) :: i
    integer(ik4),intent(in) :: rows,columns,skip
    character(len=20),dimension(rows,columns),intent(out) :: output
    character(len=100),intent(in) :: filename

    open(unit=1,file=filename,status='old',action='read',iostat&
         &=status)
    openif: if(status ==0) then
       skiplines: if(skip>0) then
          do i=1,skip
             read(1,*,iostat=status)
          end do
       end if skiplines
       do i=1,rows
          if (status /=0) exit
          read(1,*,iostat=status) output(i,1:columns)
       end do
    else openif
       write(*,*) 'ERROR opening ',trim(filename)
    end if openif

    close(unit=1)

  end subroutine readinAlpha


  subroutine readinInt(filename,rows,columns,skip,output)
  !read in data array of type integer(ik0)

    integer(ik2) :: status
    integer(ik4) :: i
    integer(ik4),intent(in) :: rows,columns,skip
    integer(ik0),dimension(rows,columns),intent(out) :: output
    character(len=100),intent(in) :: filename

    open(unit=1,file=filename,status='old',action='read',iostat&
         &=status)

    openif: if(status ==0) then
       skiplines: if(skip>0) then
          do i=1,skip
             read(1,*,iostat=status)
          end do
       end if skiplines
       do i=1,rows
          if(status/=0) exit
          read(1,*,iostat=status) output(i,1:columns)
       end do
    else openif
       write(*,*) 'ERROR opening ',trim(filename)
    end if openif

    close(unit=1)

  end subroutine readinInt


  subroutine readColumn(filename,rows,column,skip,output)
  !read in column of type real(rk8)

    integer(ik2) :: status
    integer(ik4) :: i
    integer(ik4),intent(in) :: rows,column,skip
    real(rk8),dimension(rows),intent(out) :: output
    character(len=100),intent(in) :: filename
    character(len=100) :: dummy(column-1)

    open(unit=1,file=filename,status='old',action='read',iostat&
         &=status)
    openif: if(status ==0) then
       skiplines: if(skip>0) then
          do i=1,skip
             read(1,*,iostat=status)
          end do
       end if skiplines
       do i=1,rows
          if (status /=0) exit
          read(1,*,iostat=status) dummy(:), output(i)
       end do
    else openif
       write(*,*) 'ERROR opening ',trim(filename)
    end if openif

    close(unit=1)

  end subroutine readColumn


  subroutine readColumnInt(filename,rows,column,skip,output)
  !read in integer column

    integer(ik2) :: status
    integer(ik4) :: i
    integer(ik4),intent(in) :: rows,column,skip
    integer(ik4),dimension(rows),intent(out) :: output
    character(len=100),intent(in) :: filename
    character(len=100) :: dummy(column-1)

    open(unit=1,file=filename,status='old',action='read',iostat&
         &=status)
    openif: if(status ==0) then
       skiplines: if(skip>0) then
          do i=1,skip
             read(1,*,iostat=status)
          end do
       end if skiplines
       do i=1,rows
          if (status /=0) exit
          read(1,*,iostat=status) dummy(:), output(i)
       end do
    else openif
       write(*,*) 'ERROR opening ',trim(filename)
    end if openif

    close(unit=1)

  end subroutine readColumnInt


  subroutine readColumnAlpha(filename,rows,column,skip,output)
  !read in  alphanumeric column

    integer(ik2) :: status
    integer(ik4) :: i
    integer(ik4),intent(in) :: rows,column,skip
    character(len=20),dimension(rows),intent(out) :: output
    character(len=100),intent(in) :: filename
    character(len=100) :: dummy(column-1)

    open(unit=1,file=filename,status='old',action='read',iostat&
         &=status)
    openif: if(status ==0) then
       skiplines: if(skip>0) then
          do i=1,skip
             read(1,*,iostat=status)
          end do
       end if skiplines
       do i=1,rows
          if (status /=0) exit
          read(1,*,iostat=status) dummy(:), output(i)
       end do
    else openif
       write(*,*) 'ERROR opening ',trim(filename)
    end if openif

    close(unit=1)

  end subroutine readColumnAlpha


  subroutine reparam_add(matrix,rows,columns,freq,maf,exFreq)
  ! center columns of matrix for additive genetic effect -> mean(X)=0

    integer(ik4) :: i,k
    integer(ik4),intent(in) :: rows,columns
    real(rk8),dimension(rows,columns),intent(inout) :: matrix
    real(rk8),dimension(columns),intent(inout) :: freq
    real(rk8),intent(in) :: maf
    logical,intent(in) :: exFreq

    write(*,*) '... reparametrization of X is running ...'

    do k=1,columns
       if(.not.exFreq) freq(k)=sum(matrix(:,k))/(2*rows) !count '#2' allele

       reparametrise:  if((freq(k)>=maf).and.(freq(k)<=1-maf)) then
          do i=1,rows
             matrix(i,k)=(matrix(i,k)-2*freq(k))!/sqrt(var)
          end do
       else
          write(99,1021) k
1021      format(1X,'Warning: no SNP variation (freq < MAF) in column ',I6)
       end if reparametrise
    end do

  end subroutine reparam_add


  subroutine reparam_contrast1(matadd,matdom,rows,columns,freq,maf,exFreq)
  ! reparametrise columns of matrix for additive and dominance effect
    ! according to Li "First course in population genetics" (equals
    ! the Cockerham model)

    integer(ik2) :: geno
    integer(ik4) :: i,k
    integer(ik4),intent(in) :: rows,columns
    real(rk8),dimension(rows,columns),intent(inout) :: matadd,matdom
    real(rk8),dimension(columns),intent(inout) :: freq
    real(rk8),intent(in) :: maf
    logical,intent(in) :: exFreq

    write(*,*) '... reparametrization of X and D is running ...'

    do k=1,columns
       if(.not.exFreq) freq(k)=sum(matadd(:,k))/(2*rows) !count '#2' allele

       reparametrise: if((freq(k)>=maf).and.freq(k)<=(1-maf)) then
          do i=1,rows
             geno=nint(matadd(i,k))
             select case (geno)
             case(0)
                matdom(i,k)=freq(k)**2
             case(1)
                matdom(i,k)=-freq(k)*(1-freq(k))
             case(2)
                matdom(i,k)=(1-freq(k))**2
             end select
             matadd(i,k)=matadd(i,k)-2*freq(k)
          end do
       else
          write(99,1022) k
1022      format(1X,'Warning: no SNP variation (freq < MAF) in column ',I6)
       end if reparametrise

    end do

  end subroutine reparam_contrast1


  subroutine reparam_contrast2(matadd,matdom,rows,columns,freq,maf,genoProb,exFreq)
  ! reparametrise columns of matrix for additive and dominance effect
    ! according to Alvarez-Castro & Carlborg, 2007

    integer(ik2) :: io
    integer(ik4),intent(in) :: rows,columns
    integer(ik4) :: i,k,count(rows),geno(rows)
    real(rk8),dimension(rows,columns),intent(inout) :: matadd,matdom
    real(rk8),dimension(columns),intent(inout) :: freq
    real(rk8),dimension(columns,3),intent(inout) :: genoProb
    real(rk8),intent(in) :: maf
    real(rk8) :: var,p11,p12,p22
    real(rk8),parameter :: minValue=1d-8
    logical,intent(in) :: exFreq
   
    write(*,*) '... reparametrization of X and D is running ...'   
   
    do k=1,columns
       geno(:)=nint(matadd(:,k))

       if(.not.exFreq) then
          freq(k)=sum(matadd(:,k))/(2*rows) !count '#2' allele
       
          !genotype probabilities
          where(geno(:)==0)
             count(:)=1
          elsewhere
             count(:)=0
          end where
          p11=sum(count)/real(rows,rk8)

          where(geno(:)==2)
             count(:)=1
          elsewhere
             count(:)=0
          end where
          p22=sum(count)/real(rows,rk8)
          p12=1-p11-p22       
         
          if(p11<=minValue) p11=minValue      ! numerical stabilisation 
          if(p22<=minValue) p22=minValue
          if(p12<=minValue) p12=minValue
          genoProb(k,1:3)=(/p11,p12,p22/)
       else
          p11=genoProb(k,1)
          p12=genoProb(k,2)
          p22=genoProb(k,3)
       end if
       var=p11+p22-(p11-p22)**2

       if((freq(k)>=maf).and.freq(k)<=(1-maf)) then
          do i=1,rows
             select case (geno(i))
             case(0)
                matadd(i,k)=-p12-2*p22
                matdom(i,k)=-2*p12*p22/var
             case(1)
                matadd(i,k)=1-p12-2*p22
                matdom(i,k)=4*p11*p22/var
             case(2)
                matadd(i,k)=2-p12-2*p22
                matdom(i,k)=-2*p11*p12/var
             end select
          end do
          if(var<1d-6) write(99,1021) k ! shouldn't happen actually
       else
          write(99,1022) k
          cycle
       end if
    end do

1021 format(1X,'ERROR: critical reparametrisation in column ',I6,&
            &' but maf is ok')
1022 format(1X,'Warning: no SNP variation (freq < MAF) in column ',I6)

  end subroutine reparam_contrast2


  subroutine reparam_contrast3(matadd,matdom,rows,columns,freq,maf,exFreq)
  ! reparametrise columns of matrix for additive and dominance effect
    ! according to Zeng et al. (2005)

    integer(ik2) :: geno
    integer(ik4) :: i,k
    integer(ik4),intent(in) :: rows,columns
    real(rk8),dimension(rows,columns),intent(inout) :: matadd,matdom
    real(rk8),dimension(columns),intent(inout) :: freq
    real(rk8),intent(in) :: maf
    logical,intent(in) :: exFreq
 
    write(*,*) '... reparametrization of X and D is running ...'

    do k=1,columns
       if(.not.exFreq) freq(k)=sum(matadd(:,k))/(2*rows) !count '#2' allele
      
       reparametrise: if((freq(k)>=maf).and.freq(k)<=(1-maf)) then
          do i=1,rows
             geno=nint(matadd(i,k))
             select case (geno)
             case(0)
                matdom(i,k)=-2*freq(k)**2
             case(1)
                matdom(i,k)=2*freq(k)*(1-freq(k))
             case(2)
                matdom(i,k)=-2*(1-freq(k))**2
             end select
             matadd(i,k)=matadd(i,k)-2*freq(k)
          end do
       else
          write(99,1022) k
1022      format(1X,'Warning: no SNP variation (freq < MAF) in column ',I6)
       end if reparametrise

    end do

  end subroutine reparam_contrast3


  subroutine reparam_nothing(matadd,matdom,rows,columns,freq,maf,exFreq)
  ! no reparametrisation but recording gene frequencies

    integer(ik2) :: geno
    integer(ik4) :: i,k
    integer(ik4),intent(in) :: rows,columns
    real(rk8),dimension(rows,columns),intent(inout) :: matadd,matdom
    real(rk8),dimension(columns),intent(inout) :: freq
    real(rk8),intent(in) :: maf
    logical,intent(in) :: exFreq

    write(*,*) '... counting of gene frequencies is running ...'

    do k=1,columns
       if(.not.exFreq) freq(k)=sum(matadd(:,k))/(2*rows) !count '#2' allele

       reparametrise: if((freq(k)>=maf).and.freq(k)<=(1-maf)) then
          do i=1,rows
             geno=nint(matadd(i,k))
             select case (geno)
             case(1)
                matdom(i,k)=1
             case default
                matdom(i,k)=0
             end select
          end do
       else
          write(99,1022) k
1022      format(1X,'Warning: no SNP variation (freq < MAF) in column ',I6)
       end if reparametrise

    end do

  end subroutine reparam_nothing


  subroutine standardiseGenoMat(X,D,n,p,genefreq)

    integer(ik4),intent(in) :: n,p
    integer(ik4) :: i,j
    real(rk8),intent(in) :: genefreq(p)
    real(rk8),intent(inout) :: X(n,p),D(n,p)
    real(rk8) :: hetero

    do j=1,p
       hetero=2*genefreq(j)*(1-genefreq(j))
       if(hetero<1d-3) cycle
       do i=1,n
          X(i,j)=X(i,j)/sqrt(hetero)
          D(i,j)=D(i,j)/hetero
       end do
    end do

  end subroutine standardiseGenoMat


  function filterMAF(genefreq,p,maf)

    integer(ik4),intent(in) :: p
    real(rk8),intent(in) :: genefreq(p),maf
    integer(ik0) :: filterMAF(p)

    where((genefreq(:)<=maf).or.(genefreq(:)>=1-maf))
       filterMAF(:)=1
    elsewhere
       filterMAF(:)=0
    end where

  end function filterMAF


  subroutine jointGenoFreq(genofile,rows,columns,freq)

    character(len=100),intent(in) :: genofile
    integer(ik4),intent(in) :: rows, columns
    real(rk8), intent(out) :: freq(columns,columns,9)
    integer(ik0), allocatable :: matadd(:,:)
    integer(ik4) :: j,k,geno1,geno2,counter,f(rows)

    allocate(matadd(rows,columns))
    call readinInt(genofile,rows,columns,0,matadd)
    freq=0

    loc1: do j=1,columns-1
       if(mod(j,500)==0) write(*,*) '... calculation at locus #',j
       loc2: do k=j+1,columns
          counter=0
          do geno2=0,2
             do geno1=0,2
                counter=counter+1
                where((matadd(:,j)==geno1).and.(matadd(:,k)==geno2))
                   f(:)=1
                elsewhere
                   f(:)=0
                end where
                freq(j,k,counter)=sum(f)/real(rows,rk8)
             end do
          end do
       end do loc2
    end do loc1

    deallocate(matadd)

  end subroutine jointGenoFreq


  subroutine writeout(filename,rows,columns,skip,neu,pos)
  !substitute column no. 'pos' by 'neu' in filename

    integer(ik2) :: status
    integer(ik4) :: i
    integer(ik4),intent(in) :: rows,columns,pos,skip
    real(rk8),dimension(rows),intent(in) :: neu
    real(rk8),allocatable :: tempmatrix(:,:)
    character(len=100),intent(in) :: filename
    character(len=8) :: kopf(skip,columns)
    character(len=100) :: styleA,styleF

    allocate(tempmatrix(rows,columns))
    call readin(filename,rows,columns,skip,tempmatrix)
    tempmatrix(:,pos)=neu

    skiplines: if(skip>0) then
       open(unit=2,file=filename,status='old',action='read',iostat&
            &=status)
       if(status==0) then
          do i=1,skip
             read(2,*,iostat=status) kopf(skip,1:columns)
          end do
       else
          write(*,*) 'ERROR opening ',trim(filename)
       end if
       close(2)
    end if skiplines

    styleA='(1X'
    styleF=styleA
    formatA: do i=1,columns
       styleA=trim(styleA)//',A10'
       styleF=trim(styleF)//',F20.10'
    end do formatA
    styleA=trim(styleA)//')'
    styleF=trim(styleF)//')'

    open(unit=2,file=filename,status='replace',action='write',iostat&
         &=status)
    if (status == 0) then
       if(skip>0) then
          do i=1,skip
             write(2,styleA) kopf(skip,1:columns)
          end do
       end if
       do i=1,rows
          write(2,styleF) tempmatrix(i,1:columns)
       end do
    else
       write(*,*) 'ERROR opening ',trim(filename)
    end if
    close(2)

    deallocate(tempmatrix)

  end subroutine writeout


  subroutine writeoutAlpha(filename,rows,columns,skip,neu,pos)
  !substitute column no. 'pos' by 'neu' in filename

    integer(ik2) :: status
    integer(ik4) :: i
    integer(ik4),intent(in) :: rows,columns,pos,skip
    real(rk8),dimension(rows),intent(in) :: neu
    character(len=100),intent(in) :: filename
    character(len=20),allocatable :: tempmatrix(:,:)   
    character(len=20) :: kopf(skip,columns)
    character(len=100) :: styleA,styleF

    allocate(tempmatrix(rows,columns))
    call readinAlpha(filename,rows,columns,skip,tempmatrix)

    skiplines: if(skip>0) then
       open(unit=2,file=filename,status='old',action='read',iostat&
            &=status)
       if(status==0) then
          do i=1,skip
             read(2,*,iostat=status) kopf(skip,1:columns)
          end do
       else
          write(*,*) 'ERROR opening ',trim(filename)
       end if
       close(2)
    end if skiplines

    styleA='(1X'
    styleF=styleA
    formatA: do i=1,pos-1
       styleA=trim(styleA)//',A21'
    end do formatA
    styleF=trim(styleA)//',F20.10,2x'
    styleA=trim(styleA)//',A21'
    formatF: do i=pos,columns
       styleA=trim(styleA)//',A21'
       styleF=trim(styleF)//',A21'
    end do formatF
    styleA=trim(styleA)//')'
    styleF=trim(styleF)//')'

    open(unit=2,file=filename,status='replace',action='write',iostat&
         &=status)
    if (status == 0) then
       if(skip>0) then
          do i=1,skip
             write(2,styleA) kopf(skip,1:columns)
          end do
       end if
       do i=1,rows
          write(2,styleF) tempmatrix(i,1:(pos-1)),neu(i),&
               &tempmatrix(i,(pos+1):columns)
       end do
    else
       write(*,*) 'ERROR opening ',trim(filename)
    end if
    close(2)

    deallocate(tempmatrix)

  end subroutine writeoutAlpha


  subroutine readMat(matfile,rows,Mat)
!!! read matrix which was written to a file in 3 columns: index_row
!!!  index_column value

    character(len=100),intent(in) :: matfile
    integer(ik4),intent(in) :: rows
    integer(ik4) :: i,j,k
    integer(ik2) :: status
    real(rk8),intent(out) :: Mat(rows,rows)

    Mat=0.

    open(unit=1,file=matfile,status='old',action='read',iostat&
         &=status)
    openif: if(status ==0) then
       do
          if (status /=0) exit
          read(1,*,iostat=status) i,j,Mat(i,j)
          Mat(j,i)=Mat(i,j)
       end do
    else openif
       write(*,*) 'ERROR opening ',trim(matfile)
    end if openif

    close(unit=1)

  end subroutine readMat


  subroutine writeMat(matrix,nRow,matFile,cut)
!!! write out *symmetric* matrix in 3 columns: index_row index_column value

    integer(ik4),intent(in) :: nRow
    integer(ik4) :: io, i,j
    real(rk8),intent(in) :: matrix(nRow,nRow)
    real(rk8) :: limit=1d-8
    character(len=100),intent(in) :: matFile
    logical,intent(in) :: cut

    open(unit=9,file=matFile,action='write',status='unknown',iostat=io)
    if(io==0) then
       write(*,*) '... write matrix to ',trim(matFile)
       do i=1,nRow
          do j=1,i
             if(cut.and.(abs(matrix(i,j))<limit)) cycle
             write(9,'(1x,2I8,F20.6)') i,j,matrix(i,j)
          end do
       end do
    else
       write(*,*) 'ERROR opening ',trim(matFile)
    end if
    close(9)
    
  end subroutine writeMat


  subroutine writeMatrix(matrix,nRow,nColumn,matFile,cut)
!!! write out matrix in 3 columns: index_row index_column value

    integer(ik4),intent(in) :: nRow,nColumn
    integer(ik4) :: io, i,j
    real(rk8),intent(in) :: matrix(nRow,nColumn)
    real(rk8) :: limit=1d-8
    character(len=100),intent(in) :: matFile
    logical,intent(in) :: cut

    open(unit=9,file=matFile,action='write',status='unknown',iostat=io)
    if(io==0) then
       write(*,*) '... write matrix to ',trim(matFile)
       do i=1,nRow
          do j=1,nColumn
             if(cut.and.(abs(matrix(i,j))<limit)) cycle
             write(9,'(1x,2I8,F20.6)') i,j,matrix(i,j)
          end do
       end do
    else
       write(*,*) 'ERROR opening ',trim(matFile)
    end if
    close(9)
    
  end subroutine writeMatrix


  subroutine writeMatrixCol(matrix,nRow,nColumn,matFile)
!!! write out matrix 

    integer(ik4),intent(in) :: nRow,nColumn
    integer(ik4) :: io, i,j
    real(rk8),intent(in) :: matrix(nRow,nColumn)
    real(rk8) :: limit=1d-8
    character(len=100),intent(in) :: matFile
    character(len=100) :: style

    write(style,*) nColumn
    style='(1x,'//trim(adjustl(style))//'F15.6)'

    open(unit=9,file=matFile,action='write',status='unknown',iostat=io)
    if(io==0) then
       write(*,*) '... write matrix to ',trim(matFile)
       do i=1,nRow  
          write(9,style) matrix(i,:)
       end do
    else
       write(*,*) 'ERROR opening ',trim(matFile)
    end if
    close(9)
    
  end subroutine writeMatrixCol


  subroutine writeVektor(vektor,nRow,outFile,cut)
!!! write out vector in 2 columns: index value

    integer(ik4),intent(in) :: nRow
    integer(ik4) :: io, i
    real(rk8),intent(in) :: vektor(nRow)
    real(rk8) :: limit=1d-8
    character(len=100),intent(in) :: outFile
    logical,intent(in) :: cut

    open(unit=9,file=outFile,action='write',status='unknown',iostat=io)
    if(io==0) then
       write(*,*) '... write vektor to ',trim(outFile)
       do i=1,nRow
          if(cut.and.(abs(vektor(i))<limit)) cycle
          write(9,'(1x,I8,F20.6)') i,vektor(i)
       end do
    else
       write(*,*) 'ERROR opening ',trim(outFile)
    end if
    close(9)
    
  end subroutine writeVektor


  function countWords(inputLine, ilen)
!!! from D. Berry & I. Stranden, gsBayesB.f90 (Aug, 2008)

    character(len=500000),intent(in) :: inputLine
    integer(ik4), intent(in) :: ilen
    integer(ik4) :: countWords, i
    logical :: yesSpace

    yesSpace = .TRUE. ! prior to first character there is space
    countWords = 0

    do i = 1,ilen
       if (inputLine(i:i) == ' ') then
          yesSpace = .TRUE. ! word delimiter
       else
          if (yesSpace) countWords = countWords+1 ! this is the first non-space character
          yesSpace = .FALSE. ! we are in a word
       endif
    enddo

  end function countWords


  subroutine reduceMatrix(oldMatrix,rows,columns,del,newMatrix,colNew)
!!! take only every del-th column of oldMatrix

    integer(ik2),intent(in) :: del
    integer(ik4),intent(in) :: columns,rows,colNew
    integer(ik4) :: i,count,dim
    real(rk8),intent(in) :: oldMatrix(rows,columns)
    real(rk8),intent(out) :: newMatrix(rows,colNew)

    dim=floor(real(columns/del,rk8))
    if(dim>colNew) then
       print*,'Only ',colNew,' columns are taken from input matrix.'
    end if

    count=0
    do i=1,columns
       if((mod(i,del)==1).and.(count<colNew)) then
          count=count+1
          newMatrix(:,count)=oldMatrix(:,i)
       end if
    end do

  end subroutine reduceMatrix


  subroutine timeStamp()

    character(len=12) :: date,time

    call date_and_time (date,time)

    print*
    print*,'*** time: ',time(1:2),':',time(3:4),'   date: ',date(7:8)&
         &,'.',date(5:6),'.',date(1:4)
    print*

  end subroutine timeStamp

end module inout
