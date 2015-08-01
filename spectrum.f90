module spectrum
  !
  use accuracy
  use timer
  !
  implicit none
  !
  private
  public intensity,readinput,verbose
  !
  integer(ik),parameter   :: nfiles_max =1000, max_items = 1000
  !
  integer(ik)   :: GNS=1,npoints=1001,nchar=1,nfiles=1,ipartf=0,verbose=3
  real(rk)      :: temp=298.0,partfunc=0,freql=-small_,freqr= 20000.0,thresh=1.0d-90,halfwidth=1e-2,meanmass=1.0,maxtemp=10000.0
  real(rk)      :: enermax = 1e6
  !
  character(len=cl) :: specttype,enrfilename,intfilename(nfiles_max),proftype="DOPPL",output="output"
  character(4) a_fmt
  character(9) b_fmt
  !
  type selectT
    integer(ik)  :: i
    character(len=cl) :: mask
  end type selectT
  !
  type(selectT) :: upper,lower
  !
  logical :: partfunc_do = .true., filter = .false.
  !
  contains
  !
  subroutine ReadInput
    !
    use  input
    !
    logical :: eof
    character(len=cl) :: w,v
    real(rk)      :: m1,m2
    integer(ik)   :: i
    ! -----------------------------------------------------------
    !
    write(out,"('Read the input')")
    !
    call input_options(echo_lines=.true.,error_flag=1)
    !
    !
    ! default constants
    !
    write(a_fmt,'(I3)') nchar
    write(a_fmt,'("a",a3)') adjustl(a_fmt)
    !
    lower%i = 0
    upper%i = 0
    lower%mask = ""
    upper%mask = ""
    !
    write(my_fmt,'(a,a4,a,a4,a)')  '(1x,2es16.8,1x,f5.1,1x,f12.4," <- ",f5.1,1x,f12.4,5x)'
    !
    do
        call read_line(eof) ; if (eof) exit
        call readu(w)
        select case(w)
        !
        case("STOP","FINISH","END")
          exit
        case("")
          print "(1x)"    !  Echo blank lines
          !
        case ("GNS")
          !
          call readi(gns)
          !
        case ("TEMPERATURE","TEMP")
          !
          call readf(temp)
          !
        case ("PF","QSTAT")
          !
          call readf(partfunc)
          !
          partfunc_do = .false.
          !
        case ("WINDOW","FREQUENCY","FREQUENCIES","WAVENUMBERS","RANGE")
          !
          call readf(freql)
          call readf(freqr)
          !
        case ("NPOINTS","NUMBER-OF-POINTS","NTEMPS")
          !
          call readi(npoints)
          !
          ! must be an odd number
          !
          if (mod(npoints,2)==0) npoints = npoints + 1
          !
        case ("ABSORPTION","EMISSION")
          !
          specttype = trim(w)
          !
        case ("VERBOSE")
          !
          call readi(verbose)
          !
        case ("PARTFUNC","PARTITION-FUNCTION")
          !
          specttype = trim(w)
          !
          call read_line(eof) ; if (eof) exit
          call readu(w)
          !
          do while (trim(w)/="".and.trim(w)/="END")
            !
            select case(w)
            !
            case('TEMPMAX','MAXTEMP','MAX-TEMPERATURE')
              !
              call readf(maxtemp)
              !
              if (maxtemp<small_) call report ("Illegal Max Temperature maxtem=0 "//trim(w),.true.)
              !
            case ("NTEMPS","NPOINTS")
              !
              call readi(npoints)
              !
            case ("MOMENT")
              !
              call readi(ipartf)
              !
            case ("CP")
              !
              ipartf = 3
              !
            case default
              !
              call report ("Unrecognized unit name "//trim(w),.true.)
              !
            end select
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
          enddo
          !
          if (trim(w)/="".and.trim(w)/="END") then
             !
             write (out,"('input: wrong last line in CONTRACTION =',a)") trim(w)
             stop 'input - illigal last line in CONTRACTION'
             !
          endif
          !
        case ("SELECT","FILTER")
          !
          filter  = .true.
          !
          call read_line(eof) ; if (eof) exit
          call readu(w)
          !
          do while (trim(w)/="".and.trim(w)/="END")
            !
            select case(w)
            !
            case('LOWER')
              !
              call readi(lower%i) ; lower%i = lower%i - 4
              call reada(lower%mask)
              !
            case('UPPER')
              !
              call readi(upper%i) ; upper%i = upper%i - 4
              call reada(upper%mask)
              !
            case default
              !
              call report ("Unrecognized unit name "//trim(w),.true.)
              !
            end select
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
          enddo
          !
          if (trim(w)/="".and.trim(w)/="END") then
             !
             write (out,"('input: wrong last line in CONTRACTION =',a)") trim(w)
             stop 'input - illigal last line in CONTRACTION'
             !
          endif
          !
        case ("OUTPUT")
          !
          call reada(output)
          !
        case ("NCHARS","CHARACTERS","LENGH-OF-QN-FIELD","LENGTH-OF-QN")
          !
          call readi(nchar)
          !
          write(a_fmt,'(I3)') nchar
          write(a_fmt,'("a",a3)') adjustl(a_fmt)
          !
          write(my_fmt,'(a,a4,a,a4,a)')  '(1x,2es16.8,1x,f5.1,1x,f12.4," <- ",f5.1,1x,f12.4,5x,',a_fmt,'," <-  ",',a_fmt,'))'
          !
       case ("STATES","STATESFILE","STATEFILE","STATES_FILE","STATES-FILE")
          !
          call reada(enrfilename)
          !
       case ("TRANSITION","TRANSITIONFILES","TRANSITIONS")
         !
         i = 0 
         !
         if (nitems > 1) then 
           !
           nfiles = 1 
           call reada(intfilename(1))
           !
           cycle 
           !
         endif
         !
         call read_line(eof) ; if (eof) exit
         call reada(w)
         v = trim(w)
         call upcase(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           i = i + 1
           !
           if (i>nfiles_max) call report ("Too many files, increase nfiles_max"//trim(w),.true.)
           !
           intfilename(i) = trim(v)
           !
           call read_line(eof) ; if (eof) exit
           call reada(w)
           v = trim(w)
           call upcase(w)
           !
         enddo
         !
         nfiles = i
         !
         if (trim(w)/="".and.trim(w)/="END") then
            !
            call report ("Unrecognized unit name in GRID "//trim(w),.true.)
            !
         endif
          !
       case ("THRESHOLD")
          !
          call readf(thresh)
          !
       case ("ENERMAX")
          !
          call readf(enermax)
          !
       case('GAUSSIAN','GAUSS','DOPPL','DOPPLER','RECT','BOX','BIN','STICKS','STICK')
          !
          proftype = trim(w)
          !
       case ("HWHM","HALFWIDTH")
          !
          call readf(halfwidth)
          !
       case ("MASSES","MASS")
          !
          call readf(m1)
          call readf(m2)
          !
      case default
        !
        call report ("Unrecognized unit name "//trim(w),.true.)
        !
      end select 
      !
    enddo
    !
    write(out,"('...done!'/)")
    !
  end subroutine ReadInput
  !
  subroutine intensity
   !
   use  input
   !
   integer(ik) :: info,ipoint,nlevels,i,itemp,enunit,tunit,sunit,j,ilog,ib,ie,j0,ilevelf,ileveli,indexi,indexf,iline,maxitems,kitem,l
   real(rk)    :: beta,ln2,dtemp,dfreq,temp0,beta0,intband,dpwcoef,x0,tranfreq0,tranfreq,abscoef,dfreq_,xp,xm,de,lor,b,lor2
   real(rk)    :: cmcoef,emcoef,energy,energyf,energyi,elow,jf,ji,acoef,j0rk,alpha
   character(len=cl) :: ioname
   character(wl) :: string_tmp
   !
   real(rk),allocatable :: freq(:),intens(:),jrot(:),pf(:,:),energies(:)
   integer(ik),allocatable :: gtot(:),indecies(:)
   character(len=20),allocatable :: quantum_numbers(:,:)
   !
   character(len=9) :: npoints_fmt  !text variable containing formats for reads/writes
   integer(hik):: Nintens(-60:60)
   !
   logical :: eof
   character(len=cl) :: w
   !
   ln2=log(2.0_rk)
   beta=c2/temp
   cmcoef=1.0_rk/(8.0_rk*pi*vellgt)
   emcoef=1.0_rk*planck*vellgt/(4.0_rk*pi)
   dfreq=(freqr-freql)/real(npoints-1,rk)
   !
   !   half width for Doppler profiling
   !
   dpwcoef=sqrt(2.0*ln2*boltz*avogno)/vellgt
   dpwcoef = dpwcoef*sqrt(temp/meanmass)
   !
   write(npoints_fmt,'(i9)') npoints
   !
   allocate(freq(npoints),stat=info)
   call ArrayStart('frequency',info,size(freq),kind(freq))
   !
   allocate(intens(npoints),stat=info)
   call ArrayStart('intens',info,size(intens),kind(intens))
   !
   forall(ipoint=1:npoints) freq(ipoint)=freql+real(ipoint-1,rk)*dfreq
   !
   !count number of lines (levels) in the Energy files
   !
   write(ioname, '(a)') 'Energy file'
   call IOstart(trim(ioname),enunit)
   open(unit=enunit,file=trim(enrfilename),action='read',status='old')
   i = 0
   iline = 0
   maxitems = 5
   !
   call input_options(echo_lines=.false.,error_flag=1)
   !
   do
      !
      call read_line(eof,enunit) ; if (eof) exit
      !
      call readi(itemp)
      !
      !read(enunit,*,end=10) itemp
      !
      iline = iline + 1
      !
      i = max(i,itemp)
      !
      maxitems = max(nitems,maxitems)
      !
      !cycle
      !
      !10  exit
   end do
   !
   open(unit=enunit,file=trim(enrfilename),action='read',status='old')
   !
   rewind(enunit)
   !
   nlevels = i
   !
   maxitems = maxitems-4
   !
   if (verbose>=3) print*,"nlevel =",nlevels
   !
   allocate(energies(iline),Jrot(iline),gtot(iline),indecies(nlevels),stat=info)
   call ArrayStart('energies',info,size(energies),kind(energies))
   call ArrayStart('Jrot',info,size(Jrot),kind(Jrot))
   call ArrayStart('gtot',info,size(gtot),kind(gtot))
   !
   !
   allocate(quantum_numbers(maxitems,iline),stat=info)
   quantum_numbers = ""
   !
   indecies = 0
   !
   energies = 0
   Jrot = 0
   gtot = 0
   !
   ! In case of specttype = partfunc the partition function will be computed for a series of temperatures. 
   ! npoints stands for the number of temperature steps, and the frequency limits as temperature limits
   !
   if (specttype(1:4)=='PART') then
     !
     dtemp=maxtemp/real(npoints,rk)
     !
     if (ipartf<0.or.ipartf>3) stop "illegal partfunc component, can be only 0,1,2,3"
     !
     if (verbose>=1) print("('0,1,2,3 stand for PF, 1st and 2d moments and Cp:')")
     !
     if (verbose>=1) print('("!",5x,a4,( 1x,'//npoints_fmt//'( 5x,"T=",f8.2,5x) ) )'),'  J ',(i*dtemp,i=1,npoints)
     !
     allocate(pf(0:3,npoints),stat=info)
     call ArrayStart('pf',info,size(pf),kind(pf))
     !
     pf = 0 
     !
   endif
   j0 = 0
   i = 0
   !
   ! start reading the energies 
   !
   do
      !
      call read_line(eof,enunit) ; if (eof) exit
      i = i + 1
      !
      call readi(itemp)
      call readf(energies(i))
      call readi(gtot(i))
      call readf(jrot(i))
      !
      do kitem = 5,nitems
         !
         call reada(quantum_numbers(kitem-4,i))
         !
      enddo
      !
      if (filter.and.i==1) then
        !
        if (lower%i<1.or.lower%i>nitems-4.or.upper%i<1.or.upper%i>nitems-4) then
          !
          write(out,"('wrong filter indices, upper or lower: ',2i)") upper%i,lower%i
          stop 'wrong filter indices, upper or lower'
          !
        endif
        !
      endif 
      !
      indecies(itemp) = i
      !
      j = 2*int(jrot(i))+1
      energy = energies(i)
      !
      if (partfunc_do) then
        if (j/=j0) then 
           !
           if (specttype(1:4)=='PART') then 
             if (verbose>=1.and.npoints<1000) print('("!",f8.1,3(1x,'//npoints_fmt//'es20.8))'),jrot(i),pf(ipartf,1:npoints)
           else 
             if (verbose>=4) print('("|",f8.1,1x,es16.8)'),jrot(i),partfunc
           endif
           !
        endif
        !
        partfunc=partfunc+gtot(i)*exp(-beta*energy)
        !
        j0 = j
        !
       endif
       !
       if (specttype(1:4)=='PART') then
         !
         do itemp = 1,npoints
           !
           if (energy>enermax) cycle
           !
           temp0 = real(itemp,rk)*dtemp
           !
           beta0 = planck*vellgt/(boltz*temp0)
           !
           pf(0,itemp) = pf(0,itemp) + gtot(i)*exp(-beta0*energy)
           pf(1,itemp) = pf(1,itemp) + gtot(i)*exp(-beta0*energy)*(beta0*energy)
           pf(2,itemp) = pf(2,itemp) + gtot(i)*exp(-beta0*energy)*(beta0*energy)**2
           !
           pf(3,itemp) = ( pf(2,itemp)/pf(0,itemp)-(pf(1,itemp)/pf(0,itemp))**2 )
           !
         enddo
         !
       endif
       !
   end do
   close(enunit)
   !
   ! print out the computed part. function objects and finish
   !
   if (specttype(1:4)=='PART') then
     !
     j0rk = real(j0-1)*0.5_rk
     !
     if (verbose>=0) print('("!",f8.1,3(1x,'//npoints_fmt//'es20.8/))'),j0rk,pf(ipartf,1:npoints)
     !
     if (verbose>=0) print('("!0",8x,'//npoints_fmt//'es20.8)'),pf(0,1:npoints)
     if (verbose>=0) print('("!1",8x,'//npoints_fmt//'es20.8)'),pf(1,1:npoints)
     if (verbose>=0) print('("!2",8x,'//npoints_fmt//'es20.8)'),pf(2,1:npoints)
     if (verbose>=0) print('("!3",8x,'//npoints_fmt//'es20.8)'),pf(3,1:npoints)
     !
     if (verbose>=3) print('(/a,1x,es16.8/)'),'! partition function value is',partfunc
     !
     if (verbose>=0) then 
       print('(a)'),"Partition function"
       do itemp = 1,npoints
         temp0 = real(itemp,rk)*dtemp
         print('(f12.2,1x,es20.8)'),temp0,pf(0,itemp)
       enddo
     endif 
     !
     stop
     !
   else
     !
     if (verbose>=4) print('("|",f9.1,1x,es16.8)'),real(j0-1)*0.5,partfunc
     !
   endif
   !
   if (verbose>=2.or.(partfunc_do.and.verbose>0)) print('(1x,a,1x,es16.8/)'),'! partition function value is',partfunc
   !
   !
   ! Intensities
   !
   intens = 0.0
   intband = 0.0
   sunit = 0
   !
   Nintens = 0
   !
   write(ioname, '(a)') 'Transition file'
   call IOstart(trim(ioname),tunit)
   !
   select case (trim(proftype(1:5)))
       !
   case ('GAUSS','DOPPL','LOREN')
       !
       if (verbose>=2) then
          !
          write(out,"(10x,/'Cross-sections using ',a,' profile with HWHM = ',f17.8)") trim(proftype),halfwidth
          write(out,"(10x,'Number of grid points = ',i8)") Npoints
          write(out,"(10x,'Range = ',f18.7,'-',f18.7)") freql,freqr
          write(out,"(10x,'Temperature = ',f18.7)") temp
          write(out,"(10x,'Partition function = ',f17.4)") partfunc
          write(out,"(10x,'Spectrum type = ',a/)") trim(specttype)
          !
       endif
       !
   case ('RECT','BOX')
       !
       if (verbose>=2) then
          !
          write(out,"(10x,'Range = ',f18.7,'-',f18.7)") freql,freqr
          write(out,"(10x,'Temperature = ',f18.7)") temp
          write(out,"(10x,'Partition function = ',f17.4)") partfunc
          write(out,"(10x,'Spectrum type = ',a/)") trim(specttype)
          !
       endif
       !
   case ('STICK')
       !
       write(ioname, '(a)') 'Stick spectrum'
       call IOstart(trim(ioname),sunit)
       open(unit=sunit,file=trim(output)//".stick",action='write',status='replace')
       !
       if (verbose>=2) then
          !
          write(out,"(10x,/'Stick pectra of type ',a,' stronger than ',e18.5)") trim(proftype),thresh
          write(out,"(10x,'Range = ',f18.7,'-',f18.7)") freql,freqr
          write(out,"(10x,'Temperature = ',f18.7)") temp
          write(out,"(10x,'Partition function = ',f17.4)") partfunc
          write(out,"(10x,'Spectrum type = ',a/)") trim(specttype)
          !
       endif
       !
   end select
   !
   do i = 1,nfiles
     !
     open(unit=tunit,file=trim(intfilename(i)),action='read',status='old')
     do
        !   read new line from intensities file
        !
        !read(tunit,*,end=20) ilevelf,ileveli,acoef
        !
        read(tunit,*,end=20) indexf,indexi,acoef
        !
        if (indexf>nlevels.or.indexi>nlevels) cycle
        !
        ilevelf = indecies(indexf)
        ileveli = indecies(indexi)
        !
        if (ilevelf==0.or.ileveli==0) cycle
        !
        energyf = energies(ilevelf)
        energyi = energies(ileveli)
        !
        if (energyf-energyi<0) then 
          write(6,"(4i12,2x,3es16.8)"),ilevelf,ileveli,indexf,indexi,acoef,energyf,energyi
          stop 'wrong order of indices'
        endif 
        !
        if (filter) then
          !
          if (upper%mask/=trim(quantum_numbers(upper%i,ilevelf)).or.lower%mask/=trim(quantum_numbers(lower%i,ileveli))) cycle
          !
        endif
        !
        jf = jrot(ilevelf)
        ji = jrot(ileveli)
        !
        tranfreq = energyf-energyi
        !
        tranfreq0 = tranfreq
        !
        !   half width for Doppler profiling
        if (proftype(1:5)=='DOPPL') then
          halfwidth=dpwcoef*tranfreq
        endif 
        !
        !   if transition frequency is out of selected range
        if (tranfreq0>freqr+10.0*halfwidth.or.tranfreq0<freql-10.0*halfwidth) cycle
        !
        x0 = sqrt(ln2)/halfwidth*dfreq*0.5_rk
        !
        if (tranfreq0<small_) cycle
        !
        select case (trim(specttype))
          !
        case default
          !
          print('(a,2x,a)'),'Illegal key:',trim(specttype)
          stop 'Illegal specttype-key'
          !
        case ('ABSORPTION')
          !
          abscoef=cmcoef*acoef*gtot(ilevelf)*exp(-beta*energyi)*(1.0_rk-exp(-beta*tranfreq))/(tranfreq**2*partfunc)
          !
        case ('emission','EMISSION','EMISS','emiss')
          !
          ! emission coefficient [Ergs/mol/Sr]
          !
          abscoef=emcoef*acoef*gtot(ilevelf)/real(2*ji+1,rk)*real(2*jf+1,rk)*exp(-beta*energyf)*tranfreq/(partfunc)
          !
        end select
        !
        intband = intband + abscoef
        !
        ilog=log10(abscoef)
        !
        ilog = max(min(ubound(Nintens,dim=1),ilog),lbound(Nintens,dim=1))
        !
        Nintens(ilog) = Nintens(ilog)+1
        !
        ! if only stick spectrum needed
        if (proftype(1:5)=='STICK') then
           !
           if (abscoef>thresh)  then 
              !
              ! write to .stick-file
              write(sunit,my_fmt) tranfreq,abscoef
              !
              write(out,my_fmt,advance="no"), &
              tranfreq,abscoef,jrot(ilevelf),energyf, jrot(ileveli),energyi 
              !
              write(out,"(a4)",advance="no"), " <- "
              !
              do kitem = 1,maxitems 
                !
                l = len(trim(quantum_numbers(kitem,ilevelf)))
                !
                b_fmt = "(1x,a3)" ; if (l>3) b_fmt = "(1x,a8)"
                !
                write(out,b_fmt,advance="no"), trim(quantum_numbers(kitem,ilevelf))
                !
              enddo
              !
              write(out,"(a3)",advance="no") " <-"
              !
              do kitem = 1,maxitems 
                !
                l = len(trim(quantum_numbers(kitem,ileveli)))
                !
                b_fmt = "(1x,a3)" ; if (l>3) b_fmt = "(1x,a8)"
                !
                write(out,b_fmt,advance="no"), trim(quantum_numbers(kitem,ileveli))
                !
              enddo
              !
              write(out,"(a1)",advance="yes") " "
              !
           endif 
           !
           cycle
        end if
        !
        ib =  max(nint( ( tranfreq0-halfwidth*10.0-freql)/dfreq )+1,1)
        ie =  min(nint( ( tranfreq0+halfwidth*10.0-freql)/dfreq )+1,npoints)
        !
        select case (trim(proftype(1:5)))
            !
        case ('GAUSS0','DOPPL0')
            !
            alpha = -ln2/halfwidth**2
            !
            !$omp parallel do private(ipoint,dfreq_,de) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               dfreq_=freq(ipoint)-tranfreq0
               !
               de = exp(-alpha*dfreq_)
               !
               intens(ipoint)=intens(ipoint)+abscoef*de
               !
            enddo
            !$omp end parallel do 
            !
        case ('GAUSS','DOPPL')
            !
            !$omp parallel do private(ipoint,dfreq_,xp,xm,de) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               dfreq_=freq(ipoint)-tranfreq0
               !
               xp = sqrt(ln2)/halfwidth*(dfreq_)+x0
               xm = sqrt(ln2)/halfwidth*(dfreq_)-x0
               !
               de = erf(xp)-erf(xm)
               !
               intens(ipoint)=intens(ipoint)+abscoef*0.5_rk/dfreq*de
               !
            enddo
            !$omp end parallel do 
            !
        case ('LOREN')
            !
            ib =  max(nint( ( tranfreq0-halfwidth*200.0-freql)/dfreq )+1,1)
            ie =  min(nint( ( tranfreq0+halfwidth*200.0-freql)/dfreq )+1,npoints)
            !
            !abscoef=abscoef*sqrt(ln2/pi)/halfwidth
            !
            lor = halfwidth/dfreq
            b = 0.25*pi*lor*( 4.0_rk+lor**2**2 )/( 2.0_rk+lor**2 )
            !
            b = 1.0_rk
            !
            lor = 0.5_rk/pi*halfwidth*abscoef*b
            lor2 = 0.25_rk*halfwidth**2
            !
            !$omp parallel do private(ipoint,dfreq_,xp,xm,de) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               dfreq_=freq(ipoint)-tranfreq0
               !
               intens(ipoint)=intens(ipoint)+lor/( dfreq_**2+ lor2 )
               !
            enddo
            !$omp end parallel do 
            !
        case ('RECT','BOX');
          !
          !$omp parallel do private(ipoint,dfreq_) shared(intens) schedule(dynamic)
          do ipoint=ib,ie
             dfreq_=freq(ipoint)-tranfreq0
             !if (abs(dfreq_)>halfwidth*10.0) cycle
             intens(ipoint)=max(intens(ipoint),abscoef)
          enddo
          !$omp end parallel do
          !
        case ('BIN');
          !
          !$omp parallel do private(ipoint,dfreq_) shared(intens) schedule(dynamic)
          do ipoint=ib,ie
             dfreq_=freq(ipoint)-tranfreq0
             !
             intens(ipoint)=intens(ipoint)+abscoef/dfreq
             !
          enddo
          !$omp end parallel do
          !
        end select
        !
        cycle
     20  exit
     enddo
     !
     close(tunit)
     !
   enddo 
   !
   call IOstop(trim(ioname))
   !
   if (proftype(1:5)=='DOPPL'.or.proftype(1:5)=='GAUSS'.or.proftype(1:4)=='RECT'.or.proftype(1:3)=='BIN'.or.proftype(1:3)=='BOX') then 
     !
     write(ioname, '(a)') 'Cross sections or intensities'
     call IOstart(trim(ioname),tunit)
     !
     open(unit=tunit,file=trim(output)//".xsec",action='write',status='replace')
     !
     write(tunit,'(2(1x,es16.8))'),(freq(ipoint),intens(ipoint),ipoint=1,npoints)
     !
     if (verbose>=2) print('("Total intensity  (sum):",es16.8," (int):",es16.8)'), intband,sum(intens)*dfreq
     !
     close(tunit,status='keep')
     !
   else
     !
     if (verbose>=2) print('("Total intensity = ",es16.8)'), intband
     !
   endif
   !
   !call IOStop(trim(ioname))
   !
   if (verbose>=2) then 
      !
      write(ioname, '(a)') 'Intensity distribution...'
      call IOstart(trim(ioname),tunit)
      !
      open(unit=tunit,file=trim(output)//".log",action='write',status='replace')
      !
      if (verbose>=3) print('(/"Intensity distribution...")')
      !
      do i = lbound(Nintens,dim=1),ubound(Nintens,dim=1)
        !
        write(tunit,'(i4,2x,i10)'), i,Nintens(i)
        !
      enddo
      !
      close(tunit,status='keep')
      !
   endif
   !
   if (sunit/=0) close(sunit,status='keep')
   !
   call ArrayStop('frequency')
   call ArrayStop('intens')
   call ArrayStop('energies')
   call ArrayStop('Jrot')
   call ArrayStop('gtot')
   if (specttype(1:4)=='PART') call ArrayStop('pf')
   !
  end subroutine intensity
  !
end module spectrum



