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
  integer(ik),parameter   :: nfiles_max =1000, max_items = 1000, nspecies_max = 10, nquadmax = 101, filtermax = 100
  !
  integer(ik)   :: GNS=1,npoints=1001,nchar=1,nfiles=1,ipartf=0,verbose=3,ioffset = 10,iso=1
  real(rk)      :: temp=298.0,partfunc=-1.0,freql=-small_,freqr= 20000.0,thresh=1.0d-90,halfwidth=1e-2,meanmass=1.0,maxtemp=10000.0
  real(rk)      :: voigt_gamma = 0.05, voigt_n = 0.44, offset = -25.0, pressure = 1.0_rk 
  real(rk)      :: enermax = 1e6, abscoef_thresh = 1.0d-50
  integer(ik)   :: nquad = 20      ! Number of quadrature points 
  !
  character(len=cl) :: specttype="absorption",enrfilename="none",intfilename(nfiles_max),proftype="DOPPL",output="output"
  integer(ik)   :: intJvalue(nfiles_max)
  character(4) a_fmt
  character(9) b_fmt
  !
  type selectT
    integer(ik)  :: i = 0
    character(len=cl) :: mask = ""
  end type selectT
  !
  type speciesT ! broadener
    !
    real(rk)       :: N = 0.5_rk      ! Voigt parameter N
    real(rk)       :: gamma = 0.05_rk ! Voigt parameter gamma
    real(rk)       :: ratio = 0.9_rk  ! Ratio
    real(rk)       :: T0    = 298_rk  ! Reference T, K
    real(rk)       :: P0    = 1.0_rk  ! Reference P, bar
    character(len=cl) :: name         ! Broadener number of quadrature points 
    character(len=cl) :: filename=""  ! File name with QN-dependent broadening parameters 
    character(len=cl) :: model="const"  ! Broadening model
    real(rk),pointer  :: gammaQN(:,:) ! Voigt parameter gamma as a function of QN
    real(rk),pointer  :: nQN(:,:)     ! Voigt parameter n as a function of QN
    !
  end type speciesT
  !
  type(selectT) :: upper(filtermax),lower(filtermax)
  !
  integer(ik)    :: Nspecies = 0, Nfilters = 0
  type(speciesT) :: species(nspecies_max)
  !
  logical :: partfunc_do = .true., filter = .false., histogram = .false., hitran = .false.,  histogramJ = .false.
  !
  contains
  !
  subroutine ReadInput
    !
    use  input
    !
    implicit none
    !
    logical :: eof
    character(len=cl) :: w,v
    real(rk)      :: m1,m2
    integer(ik)   :: i,ifilter
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
        case ("PRESSURE")
          !
          call readf(pressure)
          !
        case ("NQUADRATURE","NQUAD")
          !
          call readi(nquad)
          !
          if (nquad>nquadmax) call report ("nquad is > nquadmax, which then should be icreased ",.true.)
          !
          ! nquad must be odd
          !
          !if (mod(nquad,2)==0) nquad = nquad + 1
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
        case ("ABSORPTION","EMISSION","LIFETIME")
          !
          specttype = trim(w)
          !
       case ("MEM","MEMORY")
         !
         call readf(memory_limit)
          !
        case ("ISO","ISOTOPE")
          !
          call readi(iso)  
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
        case ("HISTOGRAM")
          !
          histogram = .true.
          !
        case ("HISTOGRAM-J")
          !
          histogram = .true.
          histogramJ = .true.
          !
        case ("HITRAN")
          !
          hitran = .true.
          !
        case ("SELECT","FILTER")
          !
          filter  = .true.
          !
          call read_line(eof) ; if (eof) exit
          !call readu(w)
          !
          ifilter = 0 
          !
          call readu(w)
          !
          do while (trim(w)/="".and.trim(w)/="END")
            !
            ifilter = ifilter + 1
            !
            do while (trim(w)/="".and.item<Nitems)
              !
              select case(w)
              !
              case('LOWER')
                !
                call readi(lower(ifilter)%i) ; lower(ifilter)%i = lower(ifilter)%i - 4
                call reada(lower(ifilter)%mask)
                !
              case('UPPER')
                !
                call readi(upper(ifilter)%i) ; upper(ifilter)%i = upper(ifilter)%i - 4
                call reada(upper(ifilter)%mask)
                !
              case default
                !
                call report ("Unrecognized unit name "//trim(w),.true.)
                !
              end select
              !
              call readu(w)
              !
            enddo
            !
            call read_line(eof) ; if (eof) exit
            call readu(w)
            !
          enddo
          !
          Nfilters = ifilter
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
           ! for J-dependent histograms a J-value is expected in the last column 
           !
           if (histogramJ)   call readi(intJvalue(1))
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
           ! for J-dependent histograms a J-value is expected in the last column 
           !
           if (histogramJ) then 
             call readi(intJvalue(i))
           endif 
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
            call report ("Unrecognized unit name in tranisions "//trim(w),.true.)
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
       case('GAUSSIAN','GAUSS','DOPPL','DOPPLER','RECT','BOX','BIN','STICKS','STICK','GAUS0','DOPP0',&
            'LOREN','LORENTZIAN','LORENTZ','MAX','VOIGT','PSEUDO','PSE-ROCCO','PSE-LIU','VOI-QUAD')
          !
          proftype = trim(w)
          !
          if (trim(w(1:5))=="LOREN") ioffset = 500
          if (trim(w(1:))=="VOI") ioffset = 500
          if (trim(w(1:3))=="PSE") ioffset = 500
          offset = 25.0_rk
          !
       case ("HWHM","HALFWIDTH")
          !
          call readf(halfwidth)
          !
       case ("IOFFSET")
          !
          call readi(ioffset)
          !
       case ("OFFSET")
          !
          call readf(offset)
          !
       case ("MASS")
          !
          call readf(meanmass)
          !
       case("SPECIES","BROADENER")
          !
          i = 0
          !
          call read_line(eof) ; if (eof) exit
          !
          call readu(w)
          !
          do while (trim(w)/="".and.trim(w)/="END")
            !
            i = i + 1
            !
            if (i>nspecies_max) call report ("Too many species, increase nspecies_max"//trim(w),.true.)
            !
            species(i)%name = trim(w)
            !
            do while (trim(w)/="".and.item<Nitems)
              !
              call readu(w)
              !
              select case(w)
                !
              case("GAMMA","GAMMA0")
                !
                call readf(species(i)%gamma)
                !
              case ("N")
                !
                call readf(species(i)%N)
                !
              case ("RATIO")
                !
                call readf(species(i)%ratio)
                !
              case ("P0")
                !
                call readf(species(i)%P0)
                !
              case ("T0")
                !
                call readf(species(i)%T0)
                !
              case ("FILENAME","FILE")
                !
                call reada(w)
                !
                species(i)%filename = trim(w)
                call upcase(w)
                !
              case ("MODEL")
                !
                call readu(species(i)%model)
                !
              case default 
                !
                call report ("Unrecognized unit name in sub-SPECIES "//trim(w),.true.)
                !
              end select 
                !
            enddo
            !
            call read_line(eof) ; if (eof) exit
            !
            call readu(w)
            !
          enddo
          !
          nspecies = i
          !
      case default
        !
        call report ("Unrecognized unit name "//trim(w),.true.)
        !
      end select 
      !
    enddo
    !
    !   half width for Doppler profiling
    if (proftype(1:3)=='VOI'.or.proftype(1:3)=='PSE'.or.proftype(1:3)=='LOR'.and.Nspecies>0) then
      !
      halfwidth = 0 
      !
      do i=1,Nspecies
        !
        halfwidth =  halfwidth + species(i)%ratio*species(i)%gamma*(species(i)%T0/Temp)**species(i)%N*pressure/species(i)%P0
        !
      enddo
      !
    endif 
    !
    !write(out,"('...done!'/)")
    !
  end subroutine ReadInput
  !
  subroutine intensity
   !
   use  input
   !
   integer(ik) :: info,ipoint,nlevels,i,itemp,enunit,tunit,sunit,bunit,j,ilog,ib,ie,j0,ilevelf,ileveli,indexi,indexf,iline,maxitems,kitem,l,nlines,iquad,ifilter
   integer(ik) :: indexf_,indexi_
   real(rk)    :: beta,ln2,ln22,dtemp,dfreq,temp0,beta0,intband,dpwcoef,x0,tranfreq,tranfreq_i,abscoef,dfreq_,xp,xm,de,lor,b,lor2,dfreq_2,halfwidth0,dnu_half,dxp
   real(rk)    :: cmcoef,emcoef,energy,energyf,energyi,elow,jf,ji,acoef,j0rk,bnorm,f,eta1,eta2,intens1,intens2,Va0,gammaV,a,wg,d
   real(rk)    :: sigma,alpha,gamma,y,x1,x2,voigt_,dx2,xi,L1,L2,acoef_
   integer(ik) :: Jmax,Jp,Jpp,Kp,Kpp,J_,ispecies
   real(rk)    :: gamma_,n_
   character(len=cl) :: ioname
   character(wl) :: string_tmp
   real(rk) :: Lorentz(nquadmax)
   !
   real(rk),allocatable :: freq(:),intens(:),jrot(:),pf(:,:),energies(:),Asum(:),weight(:),abciss(:),bnormq(:)
   integer(ik),allocatable :: gtot(:),indices(:)
   character(len=20),allocatable :: quantum_numbers(:,:)
   !
   integer(ik),allocatable :: nchars_quanta(:)  ! max number of characters used for each quantum number
   !
   character(len=9) :: npoints_fmt  !text variable containing formats for reads/writes
   character(len=9) :: b_fmt
   integer(hik):: Nintens(-60:60)
   !
   logical :: eof
   character(len=cl) :: w
   !
   integer(ik) :: imol
   real(rk)    :: gf,gi
   character(55) ch_q
   !
   ln2=log(2.0_rk)
   ln22 = ln2*2.0_rk
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
   ! open and count number of lines (levels) in the Energy files
   !
   if (trim(enrfilename)/="none".and..not.hitran) then
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
         if (nitems<4) then 
           call report ("In .states less than  than 4 columns, mixed up with .trans?"//trim(w),.true.)
         endif 
         !
         call readi(itemp)
         call readf(energy)
         !
         if (energy>enermax) cycle
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
      nlines = iline
      !
      if (verbose>=3) write(out,"('Max number of energy records = ',i9)") maxitems
      !
      open(unit=enunit,file=trim(enrfilename),action='read',status='old')
      !
      rewind(enunit)
      !
      nlevels = i
      !
      maxitems = maxitems-4
      !
      if (verbose>=3) print*,"nlevels: (total)",nlevels," (selected)",nlines
      !
      allocate(energies(nlines),Jrot(nlines),gtot(nlines),indices(nlevels),stat=info)
      call ArrayStart('energies',info,size(energies),kind(energies))
      call ArrayStart('Jrot',info,size(Jrot),kind(Jrot))
      call ArrayStart('gtot',info,size(gtot),kind(gtot))
      call ArrayStart('indices',info,size(indices),kind(indices))
      !
      if (trim(specttype)=='LIFETIME') THEN 
        !
        ! allocate the matrix for the sum of the A-coeffs 
        !
        allocate(Asum(nlines),stat=info)
        call ArrayStart('Asum',info,size(Asum),kind(Asum))
        !
        Asum = -1.0_rk
        !
        write(my_fmt,'(a)')  '(1x,i11,1x,f12.4,1x,f5.1,1x,es16.8,5x)'
        !
      endif
      !
      allocate(quantum_numbers(0:maxitems,iline),nchars_quanta(maxitems),stat=info)
      call ArrayStart('quantum_numbers',info,size(quantum_numbers),kind(quantum_numbers))
      call ArrayStart('nchars_quanta',info,size(nchars_quanta),kind(nchars_quanta))
      !
      ! default, min value of number of characters for quantum numbers outputs
      !
      nchars_quanta = 3
      !
      if (verbose>=3) call MemoryReport
      !
      quantum_numbers = ""
      !
      indices = 0
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
        if (trim(enrfilename)=="none") then 
           !
           write(out,"('PARTITION FUNCTION ERROR: no energy file is specified')") 
           stop "no energy file is specified"
           !
        endif
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
         call readf(energy)
         !
         if (energy>enermax) then
           i = i - 1
           cycle
         endif
         !
         energies(i) = energy
         !
         call readi(gtot(i))
         call readf(jrot(i))
         !
         if (mod(nint(jrot(i)*2),2)/=0) then 
           write(b_fmt,"(f5.1)") jrot(i)             ; quantum_numbers(0,i) = adjustl(b_fmt)
         else
           write(b_fmt,"(tl1,i5,tl1)") nint(jrot(i)) ; quantum_numbers(0,i) = adjustl(b_fmt)
         endif 
         !
         do kitem = 5,nitems
            !
            call reada(quantum_numbers(kitem-4,i))
            nchars_quanta(kitem-4) = max(len(trim(quantum_numbers(kitem-4,i))),nchars_quanta(kitem-4))
            !
         enddo
         !
         if (filter.and.i==1) then
           !
           do ifilter = 1,Nfilters
             !
             if (lower(ifilter)%i<0.or.lower(ifilter)%i>nitems-4.or.upper(ifilter)%i<0.or.upper(ifilter)%i>nitems-4) then
               !
               write(out,"('wrong filter indices, upper or lower: ',2i9)") upper(ifilter)%i+4,lower(ifilter)%i+4
               print*,quantum_numbers(:,i)
               stop 'wrong filter indices, upper or lower'
               !
             endif
             !
           enddo
           !
         endif 
         !
         indices(itemp) = i
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
   endif  ! end of processing the states files 
   !   
   if (verbose>=2.or.(partfunc_do.and.verbose>0)) print('(1x,a,1x,es16.8/)'),'! partition function value is',partfunc
   !
   ! prepare the quadratures (Gauss-Hermite)
   !
   if (trim(proftype(1:5))=='VOI-Q') then 
     write(out,"(/'Number of quadrature poitns (Gauss-Hermite) = ',i7/)") nquad
     !
     allocate(weight(nquad),abciss(nquad),bnormq(nquad),stat=info)
     call ArrayStart('weight-abciss',info,size(weight),kind(weight))
     call ArrayStart('weight-abciss',info,size(abciss),kind(abciss))
     call ArrayStart('weight-abciss',info,size(bnormq),kind(bnormq))
     !
     call gauher(abciss,weight,nquad)
     !
   endif
   !
   ! open and read broenning files
   !
   do i =1,Nspecies
     !
     if ( trim(species(i)%filename)/="" ) then
       !
       write(ioname, '(a)') 'brodenning file'
       !
       call IOstart(trim(ioname),bunit)
       open(unit=bunit,file=trim(species(i)%filename),action='read',status='old')
       !
       Jmax = 0 
       !
       do
         !
         ! Scan and find Jmax
         read(bunit,*,end=14) gamma_,n_,J_
         !
         Jmax = max(Jmax,J_)
         !
         cycle 
         14  exit
         !
       enddo
       !
       if (Jmax>0) then
         !
         allocate(species(i)%gammaQN(0:Jmax,-1:1),stat=info)
         call ArrayStart('gammaQN',info,size(species(i)%gammaQN),kind(species(i)%gammaQN))
         species(i)%gammaQN(:,:) = species(i)%gamma
         !
         allocate(species(i)%nQN(0:Jmax,-1:1),stat=info)
         call ArrayStart('nQN',info,size(species(i)%nQN),kind(species(i)%nQN))
         species(i)%nQN(:,:) = species(i)%N
         !
        else
         !
         write(out,"('Jmax = 0 in ',a)") trim(species(i)%filename)
         stop 'Jmax = 0 in a broadening file'
         !
       endif
       !
       rewind(bunit)
       !
       ! read in the QN-dependent broadening parameters
       !
       do
         !
         ! model specific read
         !
         select case ( trim(species(i)%model) )
           !
         case ('J','A0')
           !
           read(bunit,*,end=15) gamma_,n_,J
           species(i)%gammaQN(J,:) = gamma_
           species(i)%nQN(J,:) = n_
           !
         case ('JJ','A1')
           !
           read(bunit,*,end=15) gamma_,n_,Jpp,Jp
           !
           if (abs(Jp-Jpp)>1) stop 'Jp-Jpp in broadening file'
           !
           species(i)%gammaQN(Jpp,Jp-Jpp) = gamma_
           species(i)%nQN(Jpp,Jp-Jpp) = n_
           !
         case ('JJKK-X')
           !
           read(bunit,*,end=15) gamma_,n_,Jpp,Jp,Kp,Kpp
           !
           if (abs(Jp-Jpp)>1) stop 'Jp-Jpp in broadening file'
           if (Kp>Jp.or.Kpp>Jpp) stop 'K > J in broadening file'
           !
           !species(i)%gammaQN(Jpp,Jp-Jpp,Kp,Kpp) = gamma_
           !species(i)%nQN(Jpp,Jp-Jpp,Kp,Kpp) = n_
           !
         case default 
           !
           write(out,"('The broadening model',a,' is not implemented yet')") trim(species(i)%model)
           stop 'Illegal broadening model'
           !
         end select 
         !
         cycle 
         15  exit
         !
       enddo
       !
     endif
     !
     close(bunit)
     ! 
   enddo
   !
   ! Intensities
   !
   intens = 0.0
   intband = 0.0
   sunit = 0
   acoef_ = 0
   indexf_ = 0
   indexi_ = 0
   !
   Nintens = 0
   !
   write(ioname, '(a)') 'Transition file'
   call IOstart(trim(ioname),tunit)
   !
   select case (trim(proftype(1:5)))
       !
   case ('GAUSS','DOPPL','LOREN','GAUS0','DOPP0','VOIGT','PSEUD','PSE-R','PSE-L','VOI-Q')
       !
       if (verbose>=2) then
          !
          write(out,"(10x,/'Cross-sections using ',a,' profile with HWHM = ',f17.8)") trim(proftype),halfwidth
          write(out,"(10x,'Number of grid points = ',i8)") Npoints
          write(out,"(10x,'Range = ',f18.7,'-',f18.7)") freql,freqr
          write(out,"(10x,'Temperature = ',f18.7)") temp
          write(out,"(10x,'Partition function = ',f17.4)") partfunc
          write(out,"(10x,'Spectrum type = ',a/)") trim(specttype)
          if (offset<0) write(out,"(10x,'Line truncated at ',f9.2/)") offset
          if (Nspecies>0) then 
             write(out,"(10x,'Pressure = ',e18.7)") pressure
             write(out,"(10x,'Voigt parameters:  gamma       n         T0            P0')")
             do i =1,Nspecies
               write(out,"(21x,a,4f12.4)") trim(species(i)%name),species(i)%gamma,species(i)%n,species(i)%t0,species(i)%p0
             enddo 
          endif
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
   ! estimate memory 
   ! memory_limit-memory_now
   !
   if (hitran.and.partfunc<0.0) then
     write(out,"('For HITRAN partition function must be defined using PF or QSTAT keywords')")
     stop 'Undefined PF'
   endif
   !
   do i = 1,nfiles
     !
     open(unit=tunit,file=trim(intfilename(i)),action='read',status='old')
     loop_tran: do
        !
        if (histogram) then 
           !
           ! using pre-computed integrated intensities for a given Temperature and bin
           !
           read(tunit,*,end=20) tranfreq,abscoef
           !
           ! for pre-computed J-dependent integrated intensities 
           if (histogramJ) then 
             Ji = IntJvalue(i)
             Jf = Ji
           endif
           !
        elseif (hitran) then 
           !
           ! using pre-computed integrated intensities for a given Temperature and bin
           !
           read(tunit,"(i3,f12.6,e10.3e3,e10.3,10x,f10.4,12x,a55,23x,2f7.1)",end=20) imol,tranfreq,abscoef,acoef,energyi,ch_q,gf,gi
           !
           energyf = tranfreq + energyi
           !
           if (imol/=iso) cycle
           !
           if (tranfreq<small_) cycle
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
             abscoef=cmcoef*acoef*gf*exp(-beta*energyi)*(1.0_rk-exp(-beta*tranfreq))/(tranfreq**2*partfunc)
             !
           end select
           !
        else
           !
           !   read new line from intensities file
           !
           read(tunit,*,end=20) indexf,indexi,acoef
           !
           if (indexf>nlevels.or.indexi>nlevels) cycle
           !
           ilevelf = indices(indexf)
           ileveli = indices(indexi)
           !
           if (ilevelf==0.or.ileveli==0) cycle
           !
           energyf = energies(ilevelf)
           energyi = energies(ileveli)
           !
           if (energyf-energyi<-1e1) then 
             write(out,"(4i12,2x,3es16.8)"),ilevelf,ileveli,indexf,indexi,acoef,energyf,energyi
             stop 'wrong order of indices'
             cycle
           elseif (energyf-energyi<-small_) then
             cycle
           endif 
           !
           if (filter) then
             !
             do ifilter = 1,Nfilters
               !
               if (upper(ifilter)%mask/=trim(quantum_numbers(upper(ifilter)%i,ilevelf)).and.trim(upper(ifilter)%mask)/="") cycle loop_tran
               if (lower(ifilter)%mask/=trim(quantum_numbers(lower(ifilter)%i,ileveli)).and.trim(lower(ifilter)%mask)/="") cycle loop_tran
               !
             enddo
             !
           endif
           !
           jf = jrot(ilevelf)
           ji = jrot(ileveli)
           !
           tranfreq = energyf-energyi
           !
           if (tranfreq<small_) cycle
           !
           ! check for duplicates 
           !
           if (indexf==indexf_.and.indexi_==indexi.and.abs(acoef-acoef_)<small_) then
             !
             write(out,"('Duplicate: J = ',2f9.1,' E = ',2f18.6,' Ind = ',2i9,' A = ',e18.5)") jf,ji,energyf,energyi,indexf_,indexi_,acoef
             if (abs(acoef-acoef_)>small_) then
               write(out,"('but A-coefs are different:',2e17.5)") acoef_,acoef
             endif 
             stop 'Duplicates'
             !
           endif
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
            case ('LIFETIME')
             !
             if (Asum(ilevelf)<0) Asum(ilevelf) = 0 
             !
             Asum(ilevelf) = Asum(ilevelf) + acoef
             !
             !print*,ilevelf,indexf,indexi,acoef
             !
             cycle 
             !
           end select
           !
        endif
        !
        !   half width for Doppler profiling
        if (proftype(1:5)=='DOPPL') then
          !
          if (tranfreq<small_) cycle
          !
          halfwidth=dpwcoef*tranfreq
          !
        else
          !
          halfwidth = 0
          !
          do ispecies =1,Nspecies
            !
            if ( trim(species(ispecies)%filename)/="" ) then
              !
              Jpp = nint(Ji)
              Jp  = nint(Jf)
              !
              gamma_ = species(ispecies)%gammaQN(Jpp,Jp-Jpp)
              n_     = species(ispecies)%nQN(Jpp,Jp-Jpp)
              !
            else
              !
              gamma_ = species(ispecies)%gamma
              n_     = species(ispecies)%n
              !
            endif
            !
            halfwidth =  halfwidth + species(ispecies)%ratio*gamma_*(species(ispecies)%T0/Temp)**N_*pressure/species(ispecies)%P0
            !
          enddo
          !
        endif
        !
        if (abscoef<thresh) cycle
        !
        if (offset<0) offset = ioffset*halfwidth
        !
        x0 = sqrt(ln2)/halfwidth*dfreq*0.5_rk
        !
        !   half width for Doppler profiling
        if (proftype(1:3)=='VOI'.or.proftype(1:3)=='PSE') then
          if (tranfreq<small_) cycle
          halfwidth0=dpwcoef*tranfreq
          x0 = sqrt(ln2)/halfwidth0*dfreq*0.5_rk
        endif 
        !
        !   if transition frequency is out of selected range
        if (tranfreq>freqr+offset.or.tranfreq<freql-offset) cycle
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
           if (hitran) then 
             !
             write(out,"('HITRAN option has not been implemeneted to worj with STICK yet, try other options')")
             stop 'HITRAN is not working with STICK'
             !
           endif
           !
           !if (abscoef>thresh)  then 
              !
              ! write to .stick-file
              write(sunit,my_fmt) tranfreq,abscoef
              !
              write(out,my_fmt,advance="no"), &
              tranfreq,abscoef,jrot(ilevelf),energyf, jrot(ileveli),energyi 
              !
              !write(out,"(a4)",advance="no"), " <- "
              !
              do kitem = 1,maxitems 
                !
                l = len(trim(quantum_numbers(kitem,ilevelf)))
                !
                b_fmt = "(1x,a3)" ; if (l>3) b_fmt = "(1x,a8)"
                !
                write(b_fmt,"('(1x,a',i1,')')") nchars_quanta(kitem)
                !
                write(out,b_fmt,advance="no"), trim(quantum_numbers(kitem,ilevelf))
                !
              enddo
              !
              write(out,"(a3)",advance="no") " <-"
              !
              do kitem = 1,maxitems 
                !
                !l = len(trim(quantum_numbers(kitem,ileveli)))
                !
                !b_fmt = "(1x,a3)" ; if (l>3) b_fmt = "(1x,a8)"
                !
                write(b_fmt,"('(1x,a',i1,')')") nchars_quanta(kitem)
                !
                write(out,b_fmt,advance="no"), trim(quantum_numbers(kitem,ileveli))
                !
              enddo
              !
              write(out,"(a1)",advance="yes") " "
              !
           !endif 
           !
           cycle
        end if
        !
        if (abscoef<abscoef_thresh) cycle
        !
        ib =  max(nint( ( tranfreq-offset-freql)/dfreq )+1,1)
        ie =  min(nint( ( tranfreq+offset-freql)/dfreq )+1,npoints)
        !
        select case (trim(proftype(1:5)))
            !
        case ('GAUS0','DOPP0')
            !
            alpha = -ln2/halfwidth**2
            !
            abscoef=abscoef*sqrt(ln2/pi)/halfwidth
            !
            !$omp parallel do private(ipoint,dfreq_,de) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               dfreq_=freq(ipoint)-tranfreq
               !
               de = exp(alpha*dfreq_**2)
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
               dfreq_=freq(ipoint)-tranfreq
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
            ib =  max(nint( ( tranfreq-offset-freql)/dfreq )+1,1)
            ie =  min(nint( ( tranfreq+offset-freql)/dfreq )+1,npoints)
            !
            !abscoef=abscoef*sqrt(ln2/pi)/halfwidth
            !
            !lor = halfwidth/dfreq
            !b = 0.25_rk*pi*lor*( 4.0_rk+lor**2**2 )/( 2.0_rk+lor**2 )
            !
            !b = 1.0_rk
            !
            !lor = 0.5_rk/pi*halfwidth*abscoef*b
            !lor2 = 0.25_rk*halfwidth**2
            !
            dfreq_=freq(ib)-tranfreq
            !
            xm = atan( (freq(ib)-tranfreq)/halfwidth )
            xp = atan( (freq(ie)-tranfreq)/halfwidth )
            !
            b = 1.0_rk/(xp-xm)
            !
            dnu_half = dfreq*0.5_rk
            !
            !$omp parallel do private(ipoint,dfreq_,xp,xm,de) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               dfreq_=freq(ipoint)-tranfreq
               !
               xp = atan((dfreq_+dnu_half)/halfwidth)
               xm = atan((dfreq_-dnu_half)/halfwidth)
               !
               de = (xp)-(xm)
               !
               intens(ipoint)=intens(ipoint)+abscoef*de/dfreq*b
               !
               !dfreq_2=(freq(ipoint)-tranfreq)**2+lor2
               !
               !intens(ipoint)=intens(ipoint)+lor/dfreq_2
               !
            enddo
            !$omp end parallel do 
            !
        case ('PSEUD')
            !
            ib =  max(nint( ( tranfreq-offset-freql)/dfreq )+1,1)
            ie =  min(nint( ( tranfreq+offset-freql)/dfreq )+1,npoints)
            !
            !abscoef=abscoef*sqrt(ln2/pi)/halfwidth
            !
            !
            dfreq_=freq(ib)-tranfreq
            !
            xm = atan( (freq(ib)-tranfreq)/halfwidth )
            xp = atan( (freq(ie)-tranfreq)/halfwidth )
            !
            bnorm = 1.0_rk/(xp-xm)
            !
            f = ( halfwidth0**5 + 2.69269_rk*halfwidth0**4*halfwidth+2.42843_rk*halfwidth0**3*halfwidth**2+4.47163_rk*halfwidth0**2*halfwidth**3+0.07842_rk*halfwidth0*halfwidth**4+halfwidth**5)**0.2_rk
            !
            eta1 = 1.36603_rk*(halfwidth/f)-0.47719_rk*(halfwidth/f)**2+0.11116_rk*(halfwidth/f)**3
            !
            eta2 = 1.0_rk - eta1
            !
            dnu_half = dfreq*0.5_rk
            !
            !$omp parallel do private(ipoint,dfreq_,xp,xm,de,intens1,intens2) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               dfreq_=freq(ipoint)-tranfreq
               !
               xp = atan((dfreq_+dnu_half)/halfwidth)
               xm = atan((dfreq_-dnu_half)/halfwidth)
               !
               de = (xp)-(xm)
               !
               intens1 = de/dfreq*bnorm
               !
               xp = sqrt(ln2)/halfwidth0*(dfreq_)+x0
               xm = sqrt(ln2)/halfwidth0*(dfreq_)-x0
               !
               de = erf(xp)-erf(xm)
               !
               intens2 = 0.5_rk/dfreq*de
               !
               intens(ipoint)=intens(ipoint)+abscoef*(intens1*eta1+intens2*eta2)
               !
            enddo
            !$omp end parallel do
            !
        case ('PSE-R') ! pseudo_Voigt_Rocco_Cruzado_ActaPhysPol_2012.pdf
            !
            ib =  max(nint( ( tranfreq-offset-freql)/dfreq )+1,1)
            ie =  min(nint( ( tranfreq+offset-freql)/dfreq )+1,npoints)
            !
            dfreq_=freq(ib)-tranfreq
            !
            xm = atan( (freq(ib)-tranfreq)/halfwidth )
            xp = atan( (freq(ie)-tranfreq)/halfwidth )
            !
            bnorm = 1.0_rk/(xp-xm)
            !
            ! halfwidth0 is Doppler and halfwidth is Lorentzian
            !
            a = halfwidth/halfwidth0
            !
            wG = halfwidth0/ln2
            !
            Va0 = exp(a**2)*(1.0_rk-erf(a))/sqrt(pi)
            !
            gammaV = ( a+ln2*exp(-0.6055_rk*a+0.0718_rk*a**2-0.0049_rk*a**3+ 0.000136*a**4) )
            !
            eta1 = ( sqrt(pi)*Va0*gammaV-ln2 )/( sqrt(pi)*Va0*gammaV*( 1.0_rk - sqrt(pi*ln2) ) )
            !
            eta2 = 1.0_rk - eta1
            !
            dnu_half = dfreq*0.5_rk
            !
            !$omp parallel do private(ipoint,dfreq_,xp,xm,de,intens1,intens2) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               dfreq_=freq(ipoint)-tranfreq
               !
               xp = atan((dfreq_+dnu_half)/halfwidth)
               xm = atan((dfreq_-dnu_half)/halfwidth)
               !
               de = (xp)-(xm)
               !
               intens1 = de/dfreq*bnorm
               !
               xp = sqrt(ln2)/halfwidth0*(dfreq_)+x0
               xm = sqrt(ln2)/halfwidth0*(dfreq_)-x0
               !
               de = erf(xp)-erf(xm)
               !
               intens2 = 0.5_rk/dfreq*de
               !
               intens(ipoint)=intens(ipoint)+abscoef*(intens1*eta1+intens2*eta2)
               !
            enddo
            !$omp end parallel do       
            !
        case ('PSE-L') ! pseudo_Voigt_Liu_Lin_JOptSocAmB_2001.pdf
            !
            ib =  max(nint( ( tranfreq-offset-freql)/dfreq )+1,1)
            ie =  min(nint( ( tranfreq+offset-freql)/dfreq )+1,npoints)
            !
            dfreq_=freq(ib)-tranfreq
            !
            xm = atan( (freq(ib)-tranfreq)/halfwidth )
            xp = atan( (freq(ie)-tranfreq)/halfwidth )
            !
            bnorm = 1.0_rk/(xp-xm)
            !
            ! halfwidth0 is Doppler and halfwidth is Lorentzian
            !
            d = (halfwidth-halfwidth0)/(halfwidth+halfwidth0)
            !
            eta1 = 0.68188_rk+0.61293_rk*d-0.18384_rk*d**2-0.11568_rk*d**3
            !
            eta2 = 1.0_rk-eta1
            !
            dnu_half = dfreq*0.5_rk
            !
            !eta2 = 0.32460_rk-0.61825_rk*d+0.17681_rk*d**2+0.12109_rk*d**3
            !
            !$omp parallel do private(ipoint,dfreq_,xp,xm,de,intens1,intens2) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               dfreq_=freq(ipoint)-tranfreq
               !
               xp = atan((dfreq_+dnu_half)/halfwidth)
               xm = atan((dfreq_-dnu_half)/halfwidth)
               !
               de = (xp)-(xm)
               !
               intens1 = de/dfreq*bnorm
               !
               xp = sqrt(ln2)/halfwidth0*(dfreq_)+x0
               xm = sqrt(ln2)/halfwidth0*(dfreq_)-x0
               !
               de = erf(xp)-erf(xm)
               !
               intens2 = 0.5_rk/dfreq*de
               !
               intens(ipoint)=intens(ipoint)+abscoef*(intens1*eta1+intens2*eta2)
               !
            enddo
            !$omp end parallel do             
            !
        case ('VOI-Q') ! VOIGT-QUADRATURES
            !
            alpha = halfwidth0
            gamma = halfwidth
            !
            !call voi_quad(npoints,ib,ie,freq,abscoef,intens,dfreq,tranfreq,alpha,gamma,nquad,abciss,weight)
            !
            sigma = alpha/ln22
            !
            y = gamma/sigma
            !
            dx2 = dfreq*0.5_rk/sigma
            !
            do iquad=1,nquad
               !
               xi = abciss(iquad)
               !
               x1 = (freq(ib)-tranfreq)/sigma-xi
               x2 = (freq(ie)-tranfreq)/sigma-xi
               !
               L1 = atan( x1/y )
               L2 = atan( x2/y )
               !
               bnormq(iquad) = 1.0_rk/(L2-L1)*pi
               !
            enddo
            !
            Lorentz = 0
            !
            !$omp parallel do private(ipoint,xp,iquad,xi,x1,x2,dxp,L1,L2,Lorentz,voigt_) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               xp = (freq(ipoint)-tranfreq)/sigma
               !
               do iquad=1,nquad
                  !
                  xi = abciss(iquad)
                  !
                  dxp  = xp-xi
                  !
                  x1 = dxp-dx2
                  x2 = dxp+dx2
                  !
                  L1 = atan( x1/y )
                  L2 = atan( x2/y )
                  !
                  Lorentz(iquad) = (L2-L1)/pi
                  !
               enddo
               !
               voigt_ = sum(Lorentz(1:nquad)*weight(1:nquad)*bnormq(1:nquad))/dfreq
               !
               intens(ipoint)=intens(ipoint)+abscoef*voigt_
               !
            enddo
            !$omp end parallel do
            !
        case ('VOIGT')
            !
            ib =  max(nint( ( tranfreq-offset-freql)/dfreq )+1,1)
            ie =  min(nint( ( tranfreq+offset-freql)/dfreq )+1,npoints)
            !
            ! bnorm is due to to truncation of the wings, see Sharp and Burrows 2007 
            !
            !bnorm = 2.0_rk-2.0_rk/pi * atan(2.0_rk*offset/halfwidth)
            !
            !if (halfwidth<0.5_rk*dfreq**2) then
            !  !
            !  dfreq_2 = halfwidth/dfreq
            !  !
            !  bnorm = pi*0.25_rk*dfreq_2*(4.0_rk+dfreq_2**2)/(2.0_rk+dfreq_2**2)
            !  !
            !  tranfreq = min(max(nint( ( tranfreq)/dfreq )+1,1),npoints)
            !  !
            !endif
            !
            !abscoef = abscoef*bnorm
            !
            !$omp parallel do private(ipoint,tranfreq_i) shared(intens) schedule(dynamic)
            do ipoint=ib,ie
               !
               tranfreq_i = freq(ipoint)
               !
               !intens(ipoint)=intens(ipoint)+voigt_faddeeva(tranfreq_i,tranfreq,halfwidth0,halfwidth)*abscoef
               intens(ipoint)=intens(ipoint)+voigt_humlicek(tranfreq_i,tranfreq,halfwidth0,halfwidth)*abscoef
               !
            enddo
            !$omp end parallel do      
            !
        case ('MAX');
          !
          ipoint =  max(nint( ( tranfreq-freql)/dfreq )+1,1)
          !
          intens(ipoint)=max(intens(ipoint),abscoef)
          !
        case ('RECT','BOX');
          !
          !$omp parallel do private(ipoint) shared(intens) schedule(dynamic)
          do ipoint=ib,ie
             !
             intens(ipoint)=intens(ipoint)+abscoef/dfreq
             !
          enddo
          !$omp end parallel do
          !
        case ('BIN');
          !
          ipoint =  max(nint( ( tranfreq-freql)/dfreq )+1,1)
          !
          intens(ipoint)=intens(ipoint)+abscoef
          !
        end select
        !
        ! will be used to check duplicates 
        !
        indexf_ = indexf ; indexi_ = indexi ; acoef_ = acoef
        !
        cycle
     20  exit
     enddo loop_tran
     !
     close(tunit)
     !
   enddo 
   !
   call IOstop(trim(ioname))
   !
   if (trim(specttype)=='LIFETIME') THEN 
     !
     write(ioname, '(a)') 'Life times'
     call IOstart(trim(ioname),tunit)
     !
     open(unit=tunit,file=trim(output)//".life",action='write',status='replace')
     !
     do ilevelf = 1,iline
       !
       ! write to .life-file
       !
       write(tunit,my_fmt,advance="no") indices(ilevelf),energies(ilevelf),jrot(ilevelf),1.0_rk/Asum(ilevelf) 
       !
       do kitem = 1,maxitems 
         !
         !l = len(trim(quantum_numbers(kitem,ilevelf)))
         !
         !b_fmt = "(1x,a3)" ; if (l>3) b_fmt = "(1x,a8)"
         !
         write(b_fmt,"('(1x,a',i1,')')") nchars_quanta(kitem)
         !
         write(tunit,b_fmt,advance="no"), trim(quantum_numbers(kitem,ilevelf))
         !
       enddo
       !
       write(tunit,"(a1)",advance="yes") " "
       !
     enddo
     !
     close(tunit,status='keep')   
     !
   elseif (any( trim(proftype(1:3))==(/'DOP','GAU','REC','BIN','BOX','LOR','VOI','MAX','PSE'/)) ) then 
     !
     write(ioname, '(a)') 'Cross sections or intensities'
     call IOstart(trim(ioname),tunit)
     !
     open(unit=tunit,file=trim(output)//".xsec",action='write',status='replace')
     !
     write(tunit,'(2(1x,es16.8))'),(freq(ipoint),intens(ipoint),ipoint=1,npoints)
     !
     if (verbose>=2) print('(/"Total intensity  (sum):",es16.8," (int):",es16.8)'), intband,sum(intens)*dfreq
     !
     close(tunit,status='keep')
     !
   else
     !
     if (verbose>=2) print('("Total intensity = ",es16.8,f18.4)'), intband,temp
     !
   endif
   !
   !call IOStop(trim(ioname))
   !
   if (verbose>=2.and.any(Nintens(:)/=0) ) then 
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
   if (trim(enrfilename)/="none") call ArrayStop('energies')
   if (trim(enrfilename)/="none") call ArrayStop('Jrot')
   if (trim(enrfilename)/="none") call ArrayStop('gtot')
   if (specttype(1:4)=='PART') call ArrayStop('pf')
   !
  end subroutine intensity
  !



  subroutine voi_quad(npoints,ib,ie,freq,abscoef,intens,dfreq,tranfreq,alpha,gamma,nquad,abciss,weight)
     implicit none 
     !
     integer(ik),intent(in) :: npoints,nquad,ib,ie
     real(rk),intent(in)  :: abscoef,dfreq,alpha,gamma,tranfreq
     real(rk), intent(in) :: freq(npoints),abciss(nquad),weight(nquad)
     real(rk),intent(inout) :: intens(npoints)
     !
     real(rk) :: bnormq(nquad),Lorentz(nquad)
     integer(ik) :: iquad,ipoint
     real(rk)    :: sigma,y,dx2,xi,x1,x2,L1,L2,voigt_,ln22,xp
       !
       ln22 = log(2.0_rk)*2.0_rk
       !
       sigma = alpha/ln22
       !
       y = gamma/sigma
       !
       dx2 = dfreq*0.5_rk/sigma
       !
       do iquad=1,nquad
          !
          xi = abciss(iquad)
          !
          x1 = (freq(ib)-tranfreq)/sigma-xi
          x2 = (freq(ie)-tranfreq)/sigma-xi
          !
          L1 = atan( x1/y )
          L2 = atan( x2/y )
          !
          bnormq(iquad) = 1.0_rk/(L2-L1)*pi
          !
       enddo
       !
       Lorentz = 0
       !
       do ipoint=ib,ie
          !
          xp = (freq(ipoint)-tranfreq)/sigma
          !
          do iquad=1,nquad
             !
             xi = abciss(iquad)
             !
             x1 = xp-dx2-xi
             x2 = xp+dx2-xi
             !
             L1 = atan( x1/y )
             L2 = atan( x2/y )
             !
             Lorentz(iquad) = (L2-L1)/pi
             !
          enddo
          !
          voigt_ = sum(Lorentz(1:nquad)*weight(1:nquad)*bnormq(1:nquad))/dfreq
          !
          intens(ipoint)=intens(ipoint)+abscoef*voigt_
          !
       enddo
       !
  end subroutine voi_quad 




  !
  function voigt_faddeeva(nu,nu0,alpha,gamma) result(f)
   !
   real(rk),intent(in) :: nu,nu0,alpha,gamma
   real(rk) :: f,sigma,x,y
   integer(ik),parameter :: n=10
   real(rk) :: r(n),weight(n),ln2

   real(rk),parameter ::  RT2LN2 = 1.1774100225154747_rk     ! sqrt(2ln2)
   real(rk),parameter ::  RTPI   = 1.7724538509055159_rk     ! sqrt(pi)
   real(rk),parameter ::  RT2PI  = 2.5066282746310002_rk     ! sqrt(2pi)
   real(rk),parameter ::  RT2    = 1.4142135623730951_rk     ! sqrt(2)
   complex(rk) :: w,z
   !
   ln2=log(2.0_rk)
   !complex(rk) :: w_of_z
   !
   sigma = alpha / RT2LN2
   x = (nu0 - nu)*sqrt(2.0_rk*ln2)/alpha
   y = gamma / alpha*sqrt(2.0_rk*ln2)
   !
   call gauher(r,weight,n)
   !
   !r = r*alpha/RT2LN2
   !
   f = sum(weight(:)/( ( x-r(:) )**2+y**2   ))
   !
   f = f*gamma*2.0_rk*ln2/ alpha**2/pi/RTPI
   !
   f = sum(weight(:))/RTPI
   !
   !w = w_of_z(z)
   !
   !w = exp(-z**2)*(1.0_rk-erfc(z))
   !
  end function voigt_faddeeva

  function voigt_humlicek(nu,nu0,alpha,gamma) result(f)
   !
   real(rk),intent(in) :: nu,nu0,alpha,gamma
   real(rk) :: f,sigma,x,y

   real(rk),parameter ::  RT2LN2 = 1.1774100225154747_rk     ! sqrt(2ln2)
   real(rk),parameter ::  RTPI   = 1.7724538509055159_rk     ! sqrt(pi)
   real(rk),parameter ::  RT2PI  = 2.5066282746310002_rk     ! sqrt(2pi)
   real(rk),parameter ::  RT2    = 1.4142135623730951_rk     ! sqrt(pi)
   !
   sigma = alpha / RT2LN2
   x = (nu0 - nu) / sigma / RT2
   y = gamma / sigma / RT2
   !
   f = humlic(x, y) / sigma / RT2PI
   !
  end function voigt_humlicek

  !
  function humlic(x,y) result (f)
    real(rk),intent(in) :: x,y
    real(rk) :: f 
    complex(rk) :: t,humlic1,u
    real(rk) :: s
    !
    T = cmplx(y,-x)
    S = abs(x) + y
    !
    if (S >= 15) then
        ! Region I
        humlic1 = T*0.5641896_rk/(0.5_rk+T*T)
        !fprintf(stdout, "I") 
        !
        f = real(humlic1,rk)
        return
    endif
    ! 
    if (S >= 5.5) then 
        ! Region II
        U = T * T;
        humlic1 = T * (1.410474_rk + U*0.5641896_rk)/(0.75_rk + U*(3.0_rk+U))
        !
        !fprintf(stdout, "II");
        !
        f = real(humlic1,rk)
        return
        !
    endif
    !
    if (y >= 0.195_rk * abs(x) - 0.176_rk) then
        ! Region III
        humlic1 = (16.4955_rk+T*(20.20933_rk+T*(11.96482_rk &
                +T*(3.778987_rk+T*.5642236_rk)))) / (16.4955_rk+T*(38.82363_rk &
                +T*(39.27121_rk+T*(21.69274_rk+T*(6.699398_rk+T)))))
        !fprintf(stdout, "III");
        f = real(humlic1,rk)
        return
        !
    endif
    !
    ! Region IV
    U = T * T
    humlic1 = exp(U)-T*(36183.31_rk-U*(3321.9905_rk-U*(1540.787_rk-U*(219.0313_rk-U* &
       (35.76683_rk-U*(1.320522_rk-U*.56419_rk))))))/(32066.6_rk-U*(24322.84_rk-U* &
       (9022.228_rk-U*(2186.181_rk-U*(364.2191_rk-U*(61.57037_rk-U*(1.841439_rk-U)))))))
    !fprintf(stdout, "IV");
    !
    f = real(humlic1,rk)
    !
    end function humlic
    !
    subroutine gauher(x,w,n)
       !
       ! modified for w = exp(-x^2/2)/sqrt(2 pi)
       !
       integer(ik) :: n,MAXIT
       real(rk) :: w(n),x(n)
       real(rk) ::  EPS,PIM4
       parameter (EPS=3.D-16,PIM4=.75112554446494248285870300477622_rk,MAXIT=40)
       integer(ik) :: i,its,j,m
       real(rk) ::p1,p2,p3,pp,z,z1

       !
       !if (n==20) then 
       !  !
       !     x(1:20)  = &
       !       (/-7.619048541679757_rk,-6.510590157013656_rk,-5.578738805893203_rk,&
       !       -4.734581334046057_rk,-3.943967350657318_rk,-3.18901481655339_rk,&
       !       -2.458663611172367_rk,-1.745247320814127_rk,-1.042945348802751_rk,&
       !       -0.346964157081356_rk, 0.346964157081356_rk, 1.042945348802751_rk,&
       !        1.745247320814127_rk, 2.458663611172367_rk, 3.18901481655339_rk,&
       !        3.943967350657316_rk, 4.734581334046057_rk, 5.578738805893202_rk,&
       !        6.510590157013653_rk, 7.619048541679757_rk/)
       !   w(1:20) =  (/0.000000000000126_rk, 0.000000000248206_rk, 0.000000061274903_rk,&
       !        0.00000440212109_rk, 0.000128826279962_rk, 0.00183010313108_rk,&
       !        0.013997837447101_rk, 0.061506372063977_rk, 0.161739333984_rk,&
       !        0.260793063449555_rk, 0.260793063449555_rk, 0.161739333984_rk,&
       !        0.061506372063977_rk, 0.013997837447101_rk, 0.00183010313108_rk,&
       !        0.000128826279962_rk, 0.00000440212109_rk, 0.000000061274903_rk,&
       !        0.000000000248206_rk, 0.000000000000126_rk/)
       !   return
       !  !
       !endif 

       m=(n+1)/2
       do i=1,m
         !
         if(i.eq.1)then
           z=sqrt(float(2*n+1))-1.85575*(2*n+1)**(-.16667)
         else if(i.eq.2)then
           z=z-1.14*n**.426/z
         else if (i.eq.3)then
           z=1.86*z-.86*x(1)
         else if (i.eq.4)then
           z=1.91*z-.91*x(2)
         else
           z=2.*z-x(i-2)
         endif
         do its=1,MAXIT
           p1=PIM4
           p2=0.0_rk
           do j=1,n
             p3=p2
             p2=p1
             p1=z*sqrt(2.0_rk/j)*p2-sqrt(real(j-1,ark)/real(j,ark))*p3
           enddo
           pp=sqrt(2.0_rk*n)*p2
           z1=z
           z=z1-p1/pp
           if(abs(z-z1).le.EPS) exit
           !
         enddo
         if (its==MAXIT) stop 'too many iterations in gauher'
         x(i)=z
         x(n+1-i)=-z
         w(i)=2.0_rk/(pp*pp)
         w(n+1-i)=w(i)
       enddo
       !
       w = w / sqrt(pi)
       x = x*sqrt(2.0_rk)
       !
       return
    end subroutine gauher


    subroutine gauher_half(x,w,n)
       !
       ! modified for w = exp(-x^2/2)/sqrt(2 pi)
       !
       integer(ik) :: n,MAXIT
       real(rk) :: w(n),x(n)
       real(rk) ::  EPS,PIM4
       parameter (EPS=3.D-16,PIM4=.75112554446494248285870300477622_rk,MAXIT=40)
       integer(ik) :: i,its,j,m
       real(rk) ::p1,p2,p3,pp,z,z1

       !
       !if (n==20) then 
       !  !
       !     x(1:20)  = &
       !       (/-7.619048541679757_rk,-6.510590157013656_rk,-5.578738805893203_rk,&
       !       -4.734581334046057_rk,-3.943967350657318_rk,-3.18901481655339_rk,&
       !       -2.458663611172367_rk,-1.745247320814127_rk,-1.042945348802751_rk,&
       !       -0.346964157081356_rk, 0.346964157081356_rk, 1.042945348802751_rk,&
       !        1.745247320814127_rk, 2.458663611172367_rk, 3.18901481655339_rk,&
       !        3.943967350657316_rk, 4.734581334046057_rk, 5.578738805893202_rk,&
       !        6.510590157013653_rk, 7.619048541679757_rk/)
       !   w(1:20) =  (/0.000000000000126_rk, 0.000000000248206_rk, 0.000000061274903_rk,&
       !        0.00000440212109_rk, 0.000128826279962_rk, 0.00183010313108_rk,&
       !        0.013997837447101_rk, 0.061506372063977_rk, 0.161739333984_rk,&
       !        0.260793063449555_rk, 0.260793063449555_rk, 0.161739333984_rk,&
       !        0.061506372063977_rk, 0.013997837447101_rk, 0.00183010313108_rk,&
       !        0.000128826279962_rk, 0.00000440212109_rk, 0.000000061274903_rk,&
       !        0.000000000248206_rk, 0.000000000000126_rk/)
       !   return
       !  !
       !endif 
       !
       ! do only half 
       ! m=(n+1)/2
       !
       m=n
       !
       do i=1,m
         !
         if(i.eq.1)then
           z=sqrt(float(2*n+1))-1.85575*(2*n+1)**(-.16667)
         else if(i.eq.2)then
           z=z-1.14*n**.426/z
         else if (i.eq.3)then
           z=1.86*z-.86*x(1)
         else if (i.eq.4)then
           z=1.91*z-.91*x(2)
         else
           z=2.*z-x(i-2)
         endif
         do its=1,MAXIT
           p1=PIM4
           p2=0.0_rk
           do j=1,n
             p3=p2
             p2=p1
             p1=z*sqrt(2.0_rk/j)*p2-sqrt(real(j-1,ark)/real(j,ark))*p3
           enddo
           pp=sqrt(2.0_rk*n)*p2
           z1=z
           z=z1-p1/pp
           if(abs(z-z1).le.EPS) exit
           !
         enddo
         if (its==MAXIT) stop 'too many iterations in gauher'
         x(i)=z
         !
         ! exclude the symmetric point
         !x(n+1-i)=-z
         !
         w(i)=2.0_rk/(pp*pp)
         !
         ! exclude the symmetric point
         !w(n+1-i)=w(i)
         !
       enddo
       !
       w = w / sqrt(pi) * 2.0_rk
       x = x*sqrt(2.0_rk)
       !
       return
    end subroutine gauher_half



    !C  (C) Copr. 1986-92 Numerical Recipes Software 

  subroutine gauher_ark(xr,wr,n)
    !
    real(rk) :: xr(n),wr(n)
    !
    integer(ik) :: n,MAXIT
    real(ark) :: w(n),x(n)
    real(ark) ::  EPS,PIM4
    parameter (EPS=3.D-21,PIM4=.75112554446494248285870300477622_ark,MAXIT=40)
    integer(ik) :: i,its,j,m
    real(ark) ::p1,p2,p3,pp,z,z1
      !
      m=(n+1)/2
      do i=1,m
        !
        if(i.eq.1)then
          z=sqrt(float(2*n+1))-1.85575*(2*n+1)**(-.16667)
        else if(i.eq.2)then
          z=z-1.14*n**.426/z
        else if (i.eq.3)then
          z=1.86*z-.86*x(1)
        else if (i.eq.4)then
          z=1.91*z-.91*x(2)
        else
          z=2.*z-x(i-2)
        endif
        do its=1,MAXIT
          p1=PIM4
          p2=0.0_ark
          do j=1,n
            p3=p2
            p2=p1
            p1=z*sqrt(2.0_ark/j)*p2-sqrt(real(j-1,ark)/real(j,ark))*p3
          enddo
          pp=sqrt(2.0_ark*n)*p2
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS) exit
          !
        enddo
        if (its==MAXIT) stop 'too many iterations in gauher'
        x(i)=z
        x(n+1-i)=-z
        w(i)=2.0_ark/(pp*pp)
        w(n+1-i)=w(i)
      enddo
      !
      xr = real(x,rk)*sqrt(2.0_rk)
      !
      wr = real(w,rk) / sqrt(pi)
      !
      return
  end subroutine gauher_ark
  !C  (C) Copr. 1986-92 Numerical Recipes Software 


    !
  end module spectrum






