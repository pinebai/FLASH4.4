!! This utility compares two files written by the function 
!! Grid_dump in FLASH and determines if they are identical. 
!! Grid_dump's interface is Grid_dump(vars,numvars,blockId,gcell), where
!! vars is the list of variables being dumped in a single call, and 
!! numvars is the number of variables in the list. 
!! For this utility to be usable, the call should always be made with 
!! gcell = .false., the utility only compares the interior domain
!! correctly. The variable numvars determines the number on records in 
!! each file. Each variable dumped in a single call to Grid_dump represents a
!! record, and each call to Grid_dump generates a new file.

Program GridDumpCompare

  integer,parameter :: double=SELECTED_REAL_KIND(p=15)
  integer :: xsize,ysize,zsize
  real(kind=double),allocatable,dimension(:,:,:) :: arrOne,arrTwo
  character(len=80):: filename1,filename2,base1,base2
  character(len=4) :: numStr
  integer :: reclen,wordsize,numRecs,i,j,pos1,pos2

  open(1,file="GridDumpCompareInput")
  read(1,*)base1
  read(1,*)base2
  read(1,*)xsize,ysize,zsize
  read(1,*)wordsize
  read(1,*)numRecs
  read(1,*)numFiles
  close(1)
  allocate(arrOne(xsize,ysize,zsize))
  allocate(arrTwo(xsize,ysize,zsize))
  
  reclen = xsize*ysize*zsize*wordsize
  pos1=index(base1,' ')
  pos2=index(base2,' ')
  do j = 0,numfiles-1
     write(numStr,'(i4.4)')j
     filename1=base1(:pos1-1)//"/FL3"//numStr
     filename2=base2(:pos2-1)//"/FL3"//numStr
     print*,base1,base2,pos1,pos2
     print*,filename1,filename2
     open(1,file=filename1,access='direct',form='unformatted',recl=reclen)
     open(2,file=filename2,access='direct',form='unformatted',recl=reclen)

     do i = 1,numRecs
        read(1,rec=i)arrOne
        read(2,rec=i)arrTwo
        print*,'error in ',i,'th record and ',j,'th file', &
             &maxval(arrOne-arrTwo),minVal(arrOne-arrTwo)
     end do
     close(1)
     close(2)
  end do
end Program GridDumpCompare
