module spectrum
  !
  use accuracy
  use timer
  use VoigtKampff
  use Phoenix
  !
  implicit none
  !
  private
  public intensity,readinput,verbose
  !
  integer(ik),parameter   :: nfiles_max =1000, max_items = 1000, nspecies_max = 10, nquadmax = 101, filtermax = 100
  integer(ik),parameter   :: Ngrids_max = 100
  integer(ik),parameter :: HITRAN_max_ierr = 10            ! maximal number of QNs for error-specification in HITRAN
  integer(hik),parameter :: max_transitions_to_ram = 1000000000
  integer(ik) :: N_omp_procs=1
  !
  integer(ik)   :: GNS=1,npoints=1001,nchar=1,nfiles=1,ipartf=0,verbose=2,ioffset = 10,iso=1,imolecule =-1
  real(rk)      :: temp=298.0,partfunc=-1.0,partfunc_ref=-1.0,freql=-small_,freqr= 200000.0,thresh=1.0d-70
  real(rk)      :: voigt_gamma = 0.05, voigt_n = 0.44, offset = 25.0, pressure = 1.0_rk,halfwidth=1e-2,meanmass=1.0,maxtemp=10000.0
  real(rk)      :: lineshift=0
  real(rk)      :: enermax = 1e7, abscoef_thresh = 1.0d-50, abundance = 1.0d0, gf_factor = 1.0d0
  real(rk)      :: S_crit = 1e-29      ! cm/molecule, HITRAN cut-off paramater
  real(rk)      :: nu_crit = 2000.0d0  ! cm-1, HITRAN cut-off paramater
  real(rk)      :: resolving_power  = 1e6,resolving_f ! using resolving_power to set up grid
  character(len=cl) :: cutoff_model = "NONE"
  integer(ik)   :: nquad = 20      ! Number of quadrature points
  integer(hik)   :: N_to_RAM = -1000 ! Lines to keep in RAM
  !
  character(len=cl) :: specttype="ABSORPTION",proftype="DOPPL"
  character(len=wl) :: enrfilename="NONE",intfilename(nfiles_max),output="output"
  character(len=cl) :: pffilename="NONE"
  integer(ik)   :: intJvalue(nfiles_max)
  character(4) a_fmt
  character(9) b_fmt
  real(rk) :: temp_vib = -1.0, temp_ref = 296.0
  !
  !VoigtKampff parameters
  integer :: voigt_index=0
  !
  type selectT
    integer(ik)  :: i = 0
    character(len=cl) :: mask = ""
  end type selectT
  !
  type speciesT ! broadener
    !
    real(rk)       :: N = 0.5_rk      ! Voigt parameter N
    real(rk)       :: gamma = 0.0_rk  ! Voigt parameter gamma
    real(rk)       :: ratio = 1.0_rk  ! Ratio
    real(rk)       :: T0    = 296_rk  ! Reference T, K
    real(rk)       :: P0    = 1.0_rk  ! Reference P, bar
    character(len=cl) :: name = "NA"         ! Broadener name
    character(len=cl) :: filename=""  ! File name with QN-dependent broadening parameters
    character(len=cl) :: model="const"  ! Broadening model
    real(rk),pointer  :: gammaQN(:,:) ! Voigt parameter gamma as a function of QN
    real(rk),pointer  :: nQN(:,:)     ! Voigt parameter n as a function of QN
    real(rk)          :: delta = 0    ! Pressure shift induced by air, referred to p=1 atm
    logical           :: if_defined = .false. ! used to check if this broadenerer is taken fron SPECIES
    !
  end type speciesT

  type HitranErrorT ! broadener
    !
    integer(ik) :: iqn                ! the QN-index used for error
    integer(ik) :: error_vmax(0:6)=-1  ! errors vs range of qn
    integer(ik) :: N=0                ! Number of QN-entries
    integer(ik) :: ierr = 0           ! Number of QN-entries
    !
  end type HitranErrorT


  type QNT ! broadener
    !
    integer(ik) :: Kcol=1                ! Column with K 
    integer(ik) :: statecol=7            ! State column
    integer(ik) :: vibcol(2)=8           ! Range of columns with vib quanta 
    integer(ik) :: Nmodes= 1             ! Number of vib modes
    integer(ik) :: Nsym=1                ! Number of symmetries
    integer(ik) :: Rotcol(2)=(/5,6/)     ! Columns with Rot-QNs
    integer(ik) :: dens_col              ! column with number density 
    !
  end type QNT
  !
  type gridT
    integer(ik)  :: npoints = 10000 ! total number 
    integer(ik)  :: i1  = 1 ! first grid point as a global counter
    real(rk)     :: freql=-small_
    real(rk)     :: freqr= 200000.0
    real(rk)     :: dfreq
    real(rk)     :: offset = -100000.0_rk
  end type gridT
  !
  type(selectT),save :: upper(filtermax),lower(filtermax)
  type(QNT),save :: QN
  type(gridT),save :: grid(Ngrids_max)
  !
  integer(ik)    :: Nspecies = 1, Nfilters = 0, Ngrids = 0
  type(speciesT),save :: species(nspecies_max)
  type(HitranErrorT),target,save :: HITRAN_E(HITRAN_max_ierr),HITRAN_S(HITRAN_max_ierr)
  type(HitranErrorT),target,save :: HITRAN_Air,HITRAN_Self,HITRAN_N,HITRAN_Delta
  !
  integer(ik),allocatable,save  :: gamma_idx(:,:) !Used for indexing the gamma for the fast_voigt
  real(rk),allocatable,save  :: gamma_comb(:,:) !Used for indexing the gamma for the fast_voigt
  logical :: partfunc_do = .true., filter = .false., histogram = .false., hitran_do = .false.,  histogramJ = .false., &
             stick_hitran = .false.,stick_oxford = .false.,vibtemperature_do = .false., spectra_do = .false., &
             vibpopulation_do = .false.
  logical :: lineprofile_do = .false., use_width_offset = .false.
  logical :: microns = .false.
  logical :: use_resolving_power = .false.  ! using resolving for creating the grid
  logical :: ready_for_work = .false.
  logical :: accepted_work = .false.
  logical :: completed_work = .true.
  logical :: all_done = .false.
  
   type(VoigtKampffCollection),save :: fast_voigt
  
  
  !
  contains
  !
  subroutine ReadInput
    !
    use  input
    !
    implicit none
    !
    logical :: eof,if_halfwidth_defined = .false., if_species_defined = .false., if_QN_defined = .false.
    character(len=cl) :: w
    character(len=wl) :: vl
    integer(ik)   :: i,ifilter,iE,iS,npoints0,ierror,igrid
    type(HitranErrorT),pointer :: HITRAN
    real(rk) :: f_t
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
    write(my_fmt,'(a,a4,a,a4,a)')  '(1x,f16.6,1x,es16.8,1x,f5.1,1x,f12.4," <- ",f5.1,1x,f12.4,5x)'
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
          if (Nitems<2) then
              call report ("Temperatue value is not given"//trim(w),.true.)
          endif 
          !
          call readf(temp)
          !
          if (temp<small_) then
              call report ("Illegal Temperatue (must be positive>0)"//trim(w),.true.)
          endif 
          !
          !temp_ref = temp
          !
          if (Nitems>2) then 
             !
             call readu(w)
             !
             temp_vib = temp
             !
             if (w(1:3)=='VIB'.and.Nitems>3) then
                 !
                 call readf(temp_vib)
                 !
                 vibtemperature_do  = .true.
                 if (temp_vib<small_) call report ("Illegal Tvib"//trim(w),.true.)
                 !
             elseif(w(1:3)=='REF') then
                 !
                 call readf(temp_ref)
                 if (temp_ref<small_) call report ("Illegal Tref"//trim(w),.true.)
                 !
             else
                 call report ("Illegal key, expected vib or ref"//trim(w),.true.)
             endif 
             !
          endif
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
          !partfunc_ref = partfunc
          !
          if (Nitems>2) then 
             !
             call readu(w)
             !
             if (w(1:3)=='REF') then
                 !
                 call readf(partfunc_ref)
                 !
             else
                 call report ("Illegal key, expected ref"//trim(w),.true.)
             endif 
             !
          endif
         
          !
        case ("RESOLVING","R")
          !
          call readf(resolving_power)
          !
          resolving_f = log((resolving_power+1.0)/resolving_power)
          !
          use_resolving_power = .true.
          !
        case ("WINDOW","FREQUENCY","FREQUENCIES","WAVENUMBERS","RANGE")
          !
          call readf(freql)
          call readf(freqr)
          !
          if (Nitems>3) then
            !
            call readu(w)
            !
            select case(w)
              !
            case("UM","MICRON")
              !
              microns = .true.
              !
              if (freql<small_) then
                write(out,"('Error: Using zero micron is Illegal')")
                call report ("Using zero micron is Illegal"//trim(w),.true.)
              endif
              !
              ioffset = 0
              offset = 0
              !
            case ("CM-1")
              !
              microns = .false.
              !
            case default
              !
              call report ("Unrecognized unit name in units of RANGE "//trim(w),.true.)
              !
            end select
            !
          endif
          !
        case ("NPOINTS","NUMBER-OF-POINTS")
          !
          call readi(npoints)
          !
          ! must be an odd number
          !
          if (mod(npoints,2)==0) npoints = npoints + 1
          !
        case ("ABSORPTION","EMISSION","GF")
          !
          specttype = trim(w)
          !
       case ("MEM","MEMORY")
          !
         call readf(memory_limit)
         !
       case ("OMP_PROCS","NPROCS","OMP_NUM_PROCS")
          !
         call readi(N_omp_procs)
          !
        case ("ISO","ISOTOPE","IELION","IMOLECULE","IMOL")
          !
          call readi(imolecule)
          !
          if (nitems>2) then 
            call readi(iso)
          else
            iso = mod(imolecule,10)
            if (iso<1.or.iso>9) then 
              write(out,"('Input error: illegal iso/imol value = ',i7)") imolecule
              stop 'Input error: illegal iso/imol value'
            endif 
            imolecule = (imolecule-iso)/10
          endif
          !
        case ("ABUNDANCE")
          !
          call readf(abundance)
          !
        case ("GF_FACTOR")
          !
          call readf(gf_factor)
          !
        case ("VERBOSE")
          !
          call readi(verbose)
          !
        case ("PARTFUNC","PARTITION-FUNCTION","COOLING")
          !
          !specttype = trim(w)
          proftype = trim(w)
          !
          call read_line(eof) ; if (eof) exit
          call readu(w)
          !
          do while (trim(w)/="".and.trim(w)/="END")
            !
            select case(w)
            !
            case('TEMPMAX','MAXTEMP','MAX-TEMPERATURE',"TMAX")
              !
              call readf(maxtemp)
              !
              if (maxtemp<small_) call report ("Illegal Max Temperature maxtem=0 "//trim(w),.true.)
              !
            case ("NTEMPS","NPOINTS","NSTEPS","NUMBER-OF-TEMPERATURES","N")
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
             write (out,"('input: wrong last line in PARTFUNC =',a)") trim(w)
             stop 'input - illigal last line in PARTFUNC'
             !
          endif
          !
        case ("QN","QUNTUM-NUMBERS")
          !
          call read_line(eof) ; if (eof) exit
          call readu(w)
          !
          if_QN_defined = .true.
          !
          do while (trim(w)/="".and.trim(w)/="END")
            !
            select case(w)
            !
            case('NMODES')
              !
              call readi(QN%Nmodes)
              !
            case ("NSYM")
              !
              call readi(QN%Nsym)
              !
            case ("K")
              !
              call readi(QN%Kcol)
              !
            case ("ROT")
              !
              call readi(QN%rotcol(1))
              !
              if (Nitems>2) then 
                !
                call readi(QN%Rotcol(2))
                !
              else
                !
                QN%Rotcol(2) = QN%Rotcol(1)
                !
              endif
              !
            case ("DENSITY","DENS")
              !
              call readi(QN%dens_col)
              !
              vibpopulation_do = .true.
              !
            case ("VIB","VIBRATIONAL")
              !
              call readi(QN%Vibcol(1))
              !
              if (Nitems>2) then 
                !
                call readi(QN%Vibcol(2))
                !
              else
                !
                QN%Vibcol(2) = QN%Vibcol(1)
                !
              endif
              !
            case ("STATE","STATES")
              !
              call readi(QN%statecol)
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
             write (out,"('input: wrong last line in QN =',a)") trim(w)
             stop 'input - illigal last line in QN'
             !
          endif
          !
        case ("HISTOGRAM","SUPER-LINE","SUPER-LINES")
          !
          histogram = .true.
          !
        case ("HISTOGRAM-J")
          !
          histogram = .true.
          histogramJ = .true.
          !
        case ("NRAM","NLINES-TO-RAM","LINES-TO-CACHE","NCACHE")
          !
          call readf(f_t)
          N_to_RAM = int(f_t,hik)
          !
        case ("SPECTRA")
          !
          spectra_do = .true.
          !
        case ("HITRAN","OXFORD")
          !
          hitran_do = .true.
          !
          if (nitems>1) then
            !
            stick_hitran = .true.
            if (trim(w)=="OXFORD") then 
               stick_hitran = .false.
               stick_oxford = .true.
            endif
            !
            call readu(w)
            !
            if (trim(w)/="WRITE") then 
              call report ("Unrecognized unit name in HITRAN"//trim(w),.true.)
            endif
            !
            hitran_do = .false.
            proftype = 'STICK'
            specttype = 'ABSORPTION'
            !
            iE = 0
            iS = 0
            !
            call read_line(eof) ; if (eof) exit
            !
            call readu(w)
            !
            do while (trim(w)/="".and.trim(w)/="END")
               !
               select case(w)
                 !
               case("ERROR-E","ERROR-S")
                 !
                 if (iE>HITRAN_max_ierr.or.iS>HITRAN_max_ierr) call report &
                    ("Too many ERROR-E/S, increase HITRAN_max_irr"//trim(w),.true.)
                 !
                 if (trim(w)=="ERROR-E") then
                   iE = iE + 1
                   HITRAN => HITRAN_E(iE)
                   HITRAN_E(1)%N = iE
                 endif
                 if (trim(w)=="ERROR-S") then
                   iS = iS + 1
                   HITRAN => HITRAN_S(iS)
                   HITRAN_S(1)%N = iS
                 endif
                 !
                 do while (trim(w)/="".and.item<Nitems)
                   !
                   call readu(w)
                   !
                   select case(w)
                     !
                   case("QN")
                     !
                     call readi(HITRAN%iqn)
                     !
                     ! to make the QN strt from the 5th column in .states
                     HITRAN%iqn = HITRAN%iqn - 4
                     !
                   case ("IERR")
                     !
                     call readi(ierror)
                     !
                     if (ierror<0.or.ierror>6) &
                        call report("Illegal error code in HITRAN options, must be within [0..6] "//trim(w),.true.)
                     !
                     call readu(w)
                     !
                     if (trim(w)/="VMAX") call report("Illegal record, vmax is expected after ierr value in HITRAN"//trim(w),.true.)
                     !
                     call readi(hitran%error_vmax(ierror))
                     !
                   end select
                 enddo
                 !
               case("ERROR-AIR","ERROR-SELF","ERROR-N","ERROR-DELTA")
                 !
                 if (trim(w)=="ERROR-AIR")   HITRAN => HITRAN_Air
                 if (trim(w)=="ERROR-SELF")  HITRAN => HITRAN_Self
                 if (trim(w)=="ERROR-N")     HITRAN => HITRAN_N
                 if (trim(w)=="ERROR-DELTA") HITRAN => HITRAN_Delta
                 !
                 do while (trim(w)/="".and.item<Nitems)
                   !
                   call readu(w)
                   !
                   select case(w)
                     !
                   case ("IERR")
                     !
                     call readi(HITRAN%ierr)
                     !
                   case default
                     !
                     call report ("Unrecognized unit name in HITRAN-air .. options"//trim(w),.true.)
                     !
                   end select
                 enddo
                 !
               case default
                 !
                 !write(out,"('Please make sure that HITRAN is followed by a section body or an empty line')")
                 !call report ("Unrecognized unit name in HITRAN options"//trim(w),.true.)
                 !
               end select
               !
               call read_line(eof) ; if (eof) exit
               !
               call readu(w)
               !
            enddo
            !
          endif
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
             write (out,"('input: wrong last line in FILTER =',a)") trim(w)
             stop 'input - illigal last line in FILTER'
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
       case ("PFFILE")
          !
          call reada(pffilename)
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
         call reada(vl)
         vl = trim(vl)
         w = trim(vl)
         call upcase(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           i = i + 1
           !
           if (i>nfiles_max) call report ("Too many files, increase nfiles_max"//trim(w),.true.)
           !
           intfilename(i) = trim(vl)
           !
           ! for J-dependent histograms a J-value is expected in the last column
           !
           if (histogramJ) then
             call readi(intJvalue(i))
           endif
           !
           call read_line(eof) ; if (eof) exit
           call reada(vl)
           vl = trim(vl)
           w = trim(vl)
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
       case ("THRESHOLD","CUTOFF")
          !
          if (Nitems>2) then
            !
            call readu(w)
            !
            cutoff_model = trim(w)
            !
            select case(w)
              !
            case("HITRAN")
              !
              call readf(thresh)
              !
              S_crit = thresh
              !
              write(out,"('HITRAN intensity cut-off model is used, see HITRAN 2012 paper')")
              !
              call readu(w)
              !
              if (trim(w)/="NU_CRIT".and.trim(w)/="NUCRIT") &
                 call report ("Expected: keyword = NU_CRIT with the switching frequency"//trim(w),.true.)
              !
              call readf(nu_crit)
              !
            case ("EXP","EXP-WEAK","EXP-STRONG")
              !
              cutoff_model = "EXP"
              !
              call readf(thresh)
              !
              S_crit = thresh
              !
              write(out,"('EXP intensity cut-off model I0*exp(-nu/nucrit)')")
              !
              call readu(w)
              !
              if (trim(w)/="ALPHA".and.trim(w)/="NUCRIT") &
                 call report ("Expected: keyword = ALPHA with the switching frequency"//trim(w),.true.)
              !
              call readf(nu_crit)
              !
            case default
              !
              call report ("Unrecognized 2nd keyword in THRESHOLD "//trim(w),.true.)
              !
            end select
            !
          else
            !
            call readf(thresh)
            !
          endif
          !
       case ("ENERMAX")
          !
          call readf(enermax)
          !
       case('GAUSSIAN','GAUSS','DOPPL','DOPPLER','RECT','BOX','BIN','STICKS','STICK','GAUS0','DOPP0',&
            'LOREN','LORENTZIAN','LORENTZ','MAX','VOIGT','PSEUDO','PSE-ROCCO','PSE-LIU','VOI-QUAD','PHOENIX',&
            'LIFETIME','LIFETIMES','VOI-FAST','VOI-FNORM','VOI-916','T-LIFETIME','TRANS')
          !
          if (pressure<small_.and.(w(1:3)=='VOI'.or.w(1:3)=='PSE')) then
             !
             ! for pressure= 0 Voigt is repalced by Doppler
             !
             w = 'DOPPLER'
             !
             write(out,"('This is P=0 atm Voigt, which will be treated as Doppler')")
             !
          endif
          !
          if (any( w(1:3)==(/'VOI','PSE','LOR','PHO'/))) lineprofile_do = .true.
          !
          proftype = trim(w)
          !
          if (trim(w(1:5))=="LOREN") ioffset = 500
          if (trim(w(1:))=="VOI") ioffset = 500
          if (trim(w(1:3))=="PSE") ioffset = 500
          !offset = 25.0_rk
          !
          if (trim(w(1:5))=="STICK".and.Nitems>1) then
            call readu(w)
            if (trim(w)=="HITRAN") call report ("Illegal keyord for stick, HITRAN is no expected "//trim(w),.true.)
            !
            stick_hitran = .true.
            !
          endif
          !
          if (trim(w)=="LIFETIME".or.trim(w)=="LIFETIMES") then 
            specttype = "LIFETIME"
            proftype = "LIFETIME"
          endif
          !
          if (trim(w)=='T-LIFETIME') then 
            specttype = 'T-LIFETIME'
            proftype  = 'T-LIFETIME'
          endif
          !
          if (microns) then
            ioffset = 0
            offset = 0
          endif
          !
          if (nitems>1) then
            !
            call readu(w)
              !
              select case(w)
                !
              case("SAMPLING")
                !
                if (trim(proftype)=="DOPPLER")  proftype = 'DOPP0'
                if (trim(proftype)=="GAUSSIAN") proftype = 'GAUS0'
                if (trim(proftype)=="VOIGT")    proftype = 'VOIGT'
                !
              case ("NORM","NORMALIZED","AVERAGED")
                !
                if (trim(proftype)=="VOIGT") proftype = 'VOI-QUAD'
                !
              case default
                !
                call report ("Unrecognized 2nd keyword in Profile-type "//trim(w),.true.)
                !
              end select
          endif
          !
          if (trim(proftype)=='TRANS') then 
              write(my_fmt,'(a)')  '(i12,1x,i12,1x,es10.4,1x,f12.6)'
          endif
          !
       case ("HWHM","HALFWIDTH")
          !
          call readf(halfwidth)
          species(1)%gamma = halfwidth
          !
          if_halfwidth_defined = .true.
          !
       case ("IOFFSET")
          !
          call readi(ioffset)
          !
       case ("OFFSET")
          !
          call readf(offset) 
          if (nitems>2) then
            call readu(w)
            if (trim(w)/='HWHM') then
              write(out,"('illegal units of offset ',a,' must be HWHM')") trim(w)
              stop "Illegal units of offset"
            endif
            use_width_offset = .true.
          endif
          !
          if (offset<=0) stop 'offset cannot be negative or zero'
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
          if_species_defined = .true.
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
                species(i)%if_defined = .true.
                !
              case ("N")
                !
                call readf(species(i)%N)
                !
              case("DELTA")
                !
                call readf(species(i)%delta)
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
                species(i)%if_defined = .true.
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
        case ("GRID","GRIDS")
          !
          call read_line(eof) ; if (eof) exit
          !
          igrid = 0
          !
          call readu(w)
          !
          do while (trim(w)/="".and.trim(w)/="END")
            !
            igrid = igrid + 1
            !
            do while (trim(w)/="".and.item<Nitems)
              !
              select case(w)
              !
              case('RANGE')
                !
                call readf(grid(igrid)%freql)
                call readf(grid(igrid)%freqr)
                !
                if (igrid>1) then 
                  if ( abs(grid(igrid)%freql-grid(igrid-1)%freqr)>small_ ) then 
                     write (out,"('input: beginnning of the ',i5,'th grid is inconsistent with the grid on the left ',2f12.5)")&
                           igrid,grid(igrid-1)%freqr,grid(igrid)%freql
                     stop 'input - illigal last line in GRIDS'
                  endif
                endif
                !
              case('NPOINTS')
                !
                call readi(grid(igrid)%npoints)
                !
              case('OFFSET')
                !
                call readf(grid(igrid)%offset)
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
          Ngrids = igrid
          !
          if (trim(w)/="".and.trim(w)/="END") then
             !
             write (out,"('input: wrong last line in GRIDS =',a)") trim(w)
             stop 'input - illigal last line in GRIDS'
             !
          endif
          !
       case default
         !
         call report ("Unrecognized unit name "//trim(w),.true.)
         !
       end select
      !
    enddo
    !
    if (proftype(1:3)/='BIN'.and.microns) then
      write(out,"('Microns or um is currently only working with BIN, not ',a)") proftype(1:3)
      stop "Illegal profile for Microns"
    endif
    !
    if (proftype(1:3)=='BIN'.or.proftype(1:3)=='MAX') then
      offset = 0
    endif
    !
    if (proftype(1:3)/='BIN'.and.microns) then
      write(out,"('Microns or um is currently only working with BIN, not ',a)") proftype(1:3)
      stop "Illegal profile for Microns"
    endif
    !
    if (stick_hitran) then
      if (trim(proftype)/='STICK') then 
        write(out,"('For HITRAN WRITE use STICK or delete ',a)") trim(proftype)
        stop "Illegal profile for HITRAN WRITE, CHANGE TO STICK"
      endif
      !    
      if (trim(specttype)/='ABSORPTION') then 
        write(out,"('For HITRAN WRITE use ABSORPTION or delete TYPE spectral type ',a)") trim(specttype)
        stop "Illegal profile for HITRAN WRITE, CHANGE TO ABSORPTION"
      endif
      !    
      if (Nspecies/=2) then 
        write(out,"('For HITRAN WRITE the broadening parameters must be provided for two species, not ',i4)") nspecies
        stop "Illegal number of species for HITRAN WRITE, must be 2"
      endif
      !
    endif
    !
    if (stick_oxford) then
      if (trim(proftype)/='STICK') then 
        write(out,"('For OXFORD WRITE use STICK or delete ',a)") trim(proftype)
        stop "Illegal profile for OXFORD WRITE, CHANGE TO STICK"
      endif
      !    
      if (trim(specttype)/='ABSORPTION') then 
        write(out,"('For OXFORD WRITE use ABSORPTION or delete TYPE spectral type ',a)") trim(specttype)
        stop "Illegal profile for OXFORD WRITE, CHANGE TO ABSORPTION"
      endif
      !    
      if (Nspecies/=3) then 
        write(out,"('For OXFORD WRITE the broadening parameters must be provided for 3 species, not ',i4)") nspecies
        stop "Illegal number of species for OXFORD WRITE, must be 3"
      endif
      !
    endif
    !
    if (proftype(1:3)=='GAU'.and.halfwidth<small_) then 
      write(out,"('Error-Gaussian: HWHM is too small = ',g12.4)") halfwidth
      stop 'Error-Gaussian: Illegal value of HWHM'
    endif
    !
    if (hitran_do.and.(iso<0.or.imolecule<0)) then 
      write(out,"('Error-HITRAN/Oxford: iso or imolecule are illegal (negative) on undefined, use e.g. ISO 1 1  for 16 water')") 
      stop 'Error-HITRAN/Oxford: illegal iso'
    endif
    !
    if (proftype(1:3)=='BIN'.and.microns) then
      proftype = 'BIN-MICRON'
    endif
    !
    if (proftype(1:3)=='BIN'.and.use_resolving_power) then
      proftype = 'BIN-R'
    endif
    !
    if (use_resolving_power) then
      !
      if (trim(proftype) /= 'BIN-R' .and. trim(proftype) /= 'VOIGT') then 
        write(out,"(/'Warning: The RESOLVING option has not been tested for ',a,', only for BIN and VOIGT')") trim(proftype)
        !stop "The RESOLVING option can be used with BIN only"
      endif
      !
      if (Ngrids>0) then
        write(out,"('Error: resolving power (option R) is currently not working with multiple grids')")
        stop 'Error: resolving power (option R) cannot be used with multiple grids'
      endif
      !
      npoints0 = nint(real((log(freqr)-log(freql))/resolving_f,rk))+1
      !
      !if (npoints0>npoints) then
        write(out,"('For the resolving power of ',f15.1,' and range of',2f12.2)") resolving_power,freql,freqr
        write(out,"('Npoints is',i18)") npoints0
        !write(out,"('Consider increasing npoints > ',i15)") npoints
        write(out,"('npoints(max) = ln(nu2/nu1)/ln(1+1/R)+1')")
        !stop "Too small number of points for resolving_power and range given!"
      !endif
      !
      npoints = npoints0
      !
      if (freql<small_) then
        write(out,"('use_resolving_power cannot be used for range starting at zero',2f11.3)") freql,freqr
        stop "Illegal use_resolving_power with zero nu1"
      endif
      !
    endif
    !
    if (Ngrids>0.and.all( trim(proftype(1:5))/=(/ character(len=5) :: 'VOIGT','GAUS0','DOPP0','BIN'/))) then 
       write(out,"('Error: grids cannot be used with in combination with ',a)") trim(proftype(1:5))
       write(out,"('currently only with Voigt or Gaussian or Doppler or BIN SAMPLING')")
       stop 'Error: illegal use of grid and a profile type'
    endif
    !
    if (if_species_defined.and.if_halfwidth_defined) then
      write(out,"('Input Error: HWHM cannot be used together with SPECIES')")
      stop 'Input Error: HWHM cannot be used together with Species'
    endif
    !
    if (if_species_defined.and.trim(proftype)=="VOI-FAST") then
      write(out,"('Warning: VOI-FAST is now used together with SPECIES')")
      !stop 'Input Error: VOI-FAST cannot be used together with Species'
    endif
    !
    if (vibtemperature_do.and..not.if_QN_defined) then 
      write(out,"('Input Error: For non-LTE (Tvib/=Trot) QN must be defined')")
      stop 'Input Error: For non-LTE (Tvib/=Trot) QN must be defined'
    endif
    !
    !   half width for Doppler profiling
    if (any( proftype(1:3)==(/'VOI','PSE','LOR','PHO'/)).and.Nspecies>0) then
      !
      halfwidth = 0
      !
      lineprofile_do = .true.
      !
      if (if_species_defined.and..not.if_halfwidth_defined.and..not.hitran_do) then
        do i=1,Nspecies
          if (.not.species(i)%if_defined) then 
            write(out,"('Input Error: For ',i4,'-species gamma or .broad is undefined, add gamma 0.05 or name of .broad')") i
            stop 'Input Error: default gamma or .broad-file must be defined in SPECIES'
          endif
        enddo
      endif
      !
      if (.not.if_species_defined.and..not.if_halfwidth_defined.and..not.hitran_do.and..not.spectra_do) then
        write(out,"('Input Error: For ',a,' either HWHM or gamma@species must be defined, e.g. HWHM 0.05')") trim(proftype)
        stop 'Input Error: either  HWHM or gamma@SPECIES must be defined for this profile type'
      endif
      !
      do i=1,Nspecies
        !
        halfwidth =  halfwidth + species(i)%ratio*species(i)%gamma*(species(i)%T0/Temp)**species(i)%N*pressure/species(i)%P0
        lineshift = lineshift + species(i)%ratio*species(i)%delta*pressure
        !
      enddo
      !
    endif
    !
  end subroutine ReadInput
  !
  subroutine intensity
   !
   use  input
   !
   !
   real(rk)    :: hitran_Tref = 296_rk
   integer(ik) :: info,ipoint,ipoint_,nlevels,i,itemp,enunit,tunit,sunit,wunit,bunit,pfunit,j,j0,ilevelf,ileveli,indexi,indexf
   integer(ik) :: indexf_,indexi_,kitem,nlines,ifilter,k,igrid,maxitems,iline
   real(rk)    :: beta,ln2,ln22,dtemp,dfreq,temp0,beta0,intband,dpwcoef,tranfreq,abscoef,halfwidth0,tranfreq0,delta_air,beta_ref
   real(rk)    :: cmcoef,emcoef,energy,energyf,energyi,jf,ji,acoef,j0rk,gfcoef
   real(rk)    :: acoef_,cutoff,ndensity
   integer(ik) :: Jmax,Jp,Jpp,Noffset,Nspecies_,Nvib_states,ivib1,ivib2,ivib,JmaxAll,imin
   real(rk)    :: gamma_,n_,gamma_s,ener_vib,ener_rot,J_,pf_1,pf_2,t_1,t_2
   character(len=cl) :: ioname
   !
   real(rk),allocatable :: freq(:),intens(:),jrot(:),pf(:,:),energies(:),Asum(:),weight(:),abciss(:),bnormq(:)
   integer(ik),allocatable :: gtot(:),indices(:)
   character(len=20),allocatable :: quantum_numbers(:,:)
   !
   real(rk),allocatable :: acoef_RAM(:),abscoef_ram(:),nu_ram(:),intens_omp(:,:),gamma_ram(:)
   integer(ik),allocatable :: indexi_RAM(:),indexf_RAM(:)
   integer(ik),allocatable :: ileveli_RAM(:),ilevelf_RAM(:),gamma_idx_RAM(:)
   integer(ik) :: nswap_,nswap,iswap,iswap_,iomp,ichunk
   real(rk),allocatable :: energies_vib(:)
   integer(ik),allocatable :: ivib_state(:)
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
   integer(ik) :: imol,iostat_,isotope
   real(rk)    :: gf,gi,mem_t,temp_gamma_n,offset_
   character(55) ch_q,ch_broad
   !
   ln2=log(2.0_rk)
   ln22 = ln2*2.0_rk
   beta=c2/temp
   cmcoef=1.0_rk/(8.0_rk*pi*vellgt)
   emcoef=1.0_rk*planck*vellgt/(4.0_rk*pi)
   dfreq=(freqr-freql)/real(npoints-1,rk)
   !
   beta_ref = c2/temp_ref
   !
   gfcoef=1.34738d+21 ! vellgt**3*me/(8.0_rk*pi**2*e**2) e:=4.80320427e-10*Fr; me := .9109383560e-27*g ; Fr:=cm^(3/2)*sqrt(g)/s;
   !
   ! Let's include abundance intpo the constants
   !
   cmcoef = cmcoef*abundance
   emcoef = emcoef*abundance
   gfcoef = gfcoef*gf_factor
   !
   !   half width for Doppler profiling
   !
   dpwcoef=sqrt(2.0*ln2*boltz*avogno)/vellgt
   dpwcoef = dpwcoef*sqrt(temp/meanmass)
   !
   if (Ngrids>0) npoints = sum(grid(1:Ngrids)%npoints)+1
   !
   write(npoints_fmt,'(i9)') npoints
   !
   allocate(freq(npoints),stat=info)
   call ArrayStart('frequency',info,size(freq),kind(freq))
   !
   allocate(intens(npoints),stat=info)
   call ArrayStart('intens',info,size(intens),kind(intens))
   intens = 0
   !
   forall(ipoint=1:npoints) freq(ipoint)=freql+real(ipoint-1,rk)*dfreq
   !
   if (Ngrids>0) then 
     !
     if (use_resolving_power) then
       write(out,"('Error: resolving power (option R) is currently not working with multiple grids')")
       stop 'Error: resolving power (option R) is currently not working with multiple grids'
     endif
     !
     ipoint_ = 0
     freq(1) = grid(1)%freql
     do igrid = 1,Ngrids
       dfreq=(grid(igrid)%freqr-grid(igrid)%freql)/real(grid(igrid)%npoints,rk)
       grid(igrid)%i1 = ipoint_+1
       do ipoint = 1,grid(igrid)%npoints
          ipoint_ = ipoint_ + 1
          freq(ipoint_) = grid(igrid)%freql+real(ipoint-1,rk)*dfreq
       enddo
       grid(igrid)%dfreq = dfreq
       if (grid(igrid)%offset<0) grid(igrid)%offset = offset
     enddo
     ! last point 
     freq(npoints) = grid(Ngrids)%freqr
     !
     freql = grid(1)%freql
     freqr = grid(Ngrids)%freqr
     !
   endif
   !   
   if (use_resolving_power) then
     !
     do ipoint =  1,npoints
       !
       tranfreq0 = real(ipoint-1,rk)*resolving_f+log(freql)
       !
       tranfreq0 = exp(tranfreq0)
       freq(ipoint) = tranfreq0
     enddo
     !
   endif   
   !
   ! open and count number of lines (levels) in the Energy files
   !
   if (trim(enrfilename)/="NONE".and..not.hitran_do) then
      !
      if (verbose>=2) print('(/a,a,a/)'),'Reading energies from ',trim(enrfilename),'...'
      !
      write(ioname, '(a)') 'Energy file'
      call IOstart(trim(ioname),enunit)
      open(unit=enunit,file=trim(enrfilename),action='read',status='old')
      imin = huge(1)/2
      i = -imin
      iline = 0
      maxitems = 5
      Nvib_states = 0
      info = 0
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
         iline = iline + 1
         !
         i = max(i,itemp) ; imin = min(itemp,imin)
         !
         maxitems = max(nitems,maxitems)
         !
         if (vibtemperature_do) then 
            call readf(dtemp)
            call readf(ji)
            !
            if (info>0.and.int(ji)==0) then
              write(out,"('It is illegal to use not J-sorted states with vib-temperatures, E=',f14.6)") energy
              stop 'States file must be ordered by J to use with Tvib'
            endif
            !
            if (info>0) cycle
            !
            if (int(ji)==0) then 
              Nvib_states = Nvib_states + 1
            else 
              info = 1
            endif
            !
         endif
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
      if (verbose>=3) print*,"ID max:",nlevels,", min:",imin,"States selected:",nlines
      !
      allocate(energies(nlines),Jrot(nlines),gtot(nlines),indices(imin:nlevels),stat=info)
      call ArrayStart('energies',info,size(energies),kind(energies))
      call ArrayStart('Jrot',info,size(Jrot),kind(Jrot))
      call ArrayStart('gtot',info,size(gtot),kind(gtot))
      call ArrayStart('indices',info,size(indices),kind(indices))
      !
      if (vibtemperature_do) then
        write(out,"('Number of vibrational states = ',i6)") Nvib_states
        allocate(energies_vib(Nvib_states),ivib_state(nlines),stat=info)
        call ArrayStart('energies_vib',info,size(energies_vib),kind(energies_vib))
        call ArrayStart('ivib_state',info,size(ivib_state),kind(ivib_state))
        ivib_state = 1
        energies_vib = 0
        ivib1 = QN%vibcol(1)-4
        ivib2 = QN%vibcol(2)-4
      endif
      !
      if (trim(proftype)=='LIFETIME'.or.trim(proftype)=='T-LIFETIME') THEN
        !
        ! allocate the matrix for the sum of the A-coeffs
        !
        allocate(Asum(nlines),stat=info)
        call ArrayStart('Asum',info,size(Asum),kind(Asum))
        !
        Asum = -1.0_rk
        !
      endif
      !
      allocate(quantum_numbers(0:maxitems,nlines),nchars_quanta(maxitems),stat=info)
      call ArrayStart('quantum_numbers',info,size(quantum_numbers),kind(quantum_numbers))
      call ArrayStart('nchars_quanta',info,size(nchars_quanta),kind(nchars_quanta))
      !
      ! default, min value of number of characters for quantum numbers outputs
      !
      nchars_quanta = 3
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
      if (proftype(1:4)=='PART'.or.proftype(1:4)=='COOL') then
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
        if (proftype(1:4)=='PART') then
          !
          if (ipartf<0.or.ipartf>3) stop "illegal partfunc component, can be only 0,1,2,3"
          !
          if (verbose>=1) print("('0,1,2,3 stand for PF, 1st and 2d moments and Cp:')")
          !
          if (verbose>=1) print('("!",5x,a4,( 1x,'//npoints_fmt//'( 5x,"T=",f8.2,5x) ) )'),'  J ',(i*dtemp,i=1,npoints)
          !
        endif
        !
        allocate(pf(0:3,npoints),stat=info)
        call ArrayStart('pf',info,size(pf),kind(pf))
        !
        if (proftype(1:4)=='COOL') then
          !
          if (trim(specttype(1:5))/="EMISS") then
            write(out,"('Error: COOLING can work only with EMISSION, not',a)") trim(specttype)
            stop 'Error: COOLING can work only with EMISSION'
          endif
          !
          if (histogram) then
            write(out,"('Error: COOLING does not work for HISTOGRAM')")
            stop 'Error: COOLING does not work with HISTOGRAM'
          endif
          !
          !allocate(cooling(npoints),stat=info)
          !call ArrayStart('cooling',info,size(cooling),kind(cooling))
          !cooling = 0
          !
        endif
        !
        pf = 0
        !
      endif
      j0 = 0
      i = 0
      if (partfunc_do) partfunc = 0
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
         if (vibtemperature_do) then
           if (int(jrot(i))==0) then 
              energies_vib(i) = energy
              ivib_state(i) = i
              ener_vib = energy
              ener_rot = energy-ener_vib
           else
              loop_vib : do ivib  = 1,Nvib_states
                if (all(quantum_numbers(ivib1:ivib2,i)==quantum_numbers(ivib1:ivib2,ivib))) then
                  !
                  ivib_state(i) = ivib
                  ener_vib = energies_vib(ivib)
                  ener_rot = energy-ener_vib
                  exit loop_vib
                  !
                endif
              enddo loop_vib
           endif
         endif
         !
         if (partfunc_do) then
           if (j/=j0) then
              !
              if (proftype(1:4)=='PART') then
                if (verbose>=1.and.npoints<1000) print('("!",f8.1,3(1x,'//npoints_fmt//'es20.8))'),jrot(i),pf(ipartf,1:npoints)
              else
                if (verbose>=4) print('("|",f8.1,1x,es16.8)'),jrot(i),partfunc
              endif
              !
           endif
           !
           if (.not.vibtemperature_do) then
              partfunc=partfunc+gtot(i)*exp(-beta*energy)
           else 
              ! split PF into a product of vib and rot parts 
              partfunc=partfunc+gtot(i)*exp(-c2/temp*ener_rot)*exp(-c2/temp_vib*ener_vib)
           endif 
           !
           j0 = j
           !
         endif
         !
         if (proftype(1:4)=='PART'.or.proftype(1:4)=='COOL') then
           !
           do itemp = 1,npoints
             !
             if (energy>enermax) cycle
             !
             temp0 = real(itemp,rk)*dtemp
             !
             beta0 = c2/temp0
             !
             !beta0 = planck*vellgt/(boltz*temp0)
             !
             pf(0,itemp) = pf(0,itemp) + gtot(i)*exp(-beta0*energy)
             pf(1,itemp) = pf(1,itemp) + gtot(i)*exp(-beta0*energy)*(beta0*energy)
             pf(2,itemp) = pf(2,itemp) + gtot(i)*exp(-beta0*energy)*(beta0*energy)**2
             !
           enddo
           !
         endif
         !
      end do
      !
      if (proftype(1:4)=='PART'.or.proftype(1:4)=='COOL') then
         !
         pf(3,:) = ( pf(2,:)/pf(0,:)-(pf(1,:)/pf(0,:))**2 )*R_ + 2.5_rk*R_
         !
      endif
      !
      close(enunit)
      !
      ! print out the computed part. function objects and finish
      !
      if (proftype(1:4)=='PART'.or.proftype(1:4)=='COOL') then
        !
        write(ioname, '(a)') 'Partition function'
        call IOstart(trim(ioname),tunit)
        !
        open(unit=tunit,file=trim(output)//".pf",action='write',status='replace')
        !
        do itemp = 1,npoints
          !
          temp0 = real(itemp,rk)*dtemp
          !
          write(tunit,"(f8.1,1x,f15.4)") temp0,pf(0,itemp)
          !
        enddo
        !
        close(tunit,status='keep')
        !
      endif
      !
      if (proftype(1:4)=='PART') then
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
        write(ioname, '(a)') 'Partition function'
        call IOstart(trim(ioname),tunit)
        !
        open(unit=tunit,file=trim(output)//".pf",action='write',status='replace')
        !
        do itemp = 1,npoints
          !
          temp0 = real(itemp,rk)*dtemp
          !
          write(tunit,"(f8.1,1x,f15.4)") temp0,pf(0,itemp)
          !
        enddo
        !
        close(tunit,status='keep')
        !
        stop
        !
      else
        !
        if (verbose==3) print('("|",f8.1,1x,es16.8)'),real(j0-1)*0.5,partfunc
        !
      endif
      !
      if (verbose>=2) print('(1x,a/)'),'... reading energies done!'
      !
   endif  ! end of processing the states files
   !
   if (trim(pffilename)/='NONE'.and.proftype(1:4)/='PART') then 
     if (verbose>=2) print('(/a,a,a/)'),'Reading pf from ',trim(pffilename),'...'
     !
     write(ioname, '(a)') 'PF file'
     call IOstart(trim(ioname),pfunit)
     open(unit=pfunit,file=trim(pffilename),action='read',status='old')
     !
     pf_2 = 1.0_rk
     if (trim(enrfilename)/="NONE") pf_2 = gtot(1)
     T_1 = 0
     T_2 = 0
     do while (T_2<Temp)
      T_1 = T_2 ; pf_1 = pf_2
      read(pfunit,*,end=31) t_2,pf_2
      !
      cycle
      !
      31 continue
      !
      write(out,"('The tempreture requested ',f12.4,' is > Tmax=',f9.2,'K of ',a)") Temp,t_2,trim(pffilename)
      write(out,"('The PF-file is not applicable, consider providing PF value (keyword PF) explicitly ... ')")
      write(out,"('... or let exocross do re-compute pf using eneegies (no PF keyword needed).')")
      stop 'The PF-file is not applicable'
     enddo
     !
     ! From Numerical Recipes in Fortran 77
     if (abs(T_1-T_2)<small_) stop 'Duplicate temperature in partition file'
     partfunc = pf_1 + (temp - T_1) * (pf_2 - pf_1)/(T_2 - T_1)
     !
   endif
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
   !Prepare VoigtKampff
   if(trim(proftype(1:5))=='VOI-F') then
     !
     offset_ = offset
     if ( use_width_offset ) offset_ = max(offset*halfwidth,25._rk)
     !
     if(trim(proftype(1:6))=='VOI-FN') then
      call fast_voigt%construct(dpwcoef,offset_,dfreq,.true.)
      !call initalize_voigt_kampff(dpwcoef,dfreq,offset,voigt_index,1)
     else
      call fast_voigt%construct(dpwcoef,offset_,dfreq,.false.)
      !call initalize_voigt_kampff(dpwcoef,dfreq,offset,voigt_index,0)
     endif
     !call add_lorentzian(halfwidth,voigt_index,i)
   endif
   !
   ! open and read broadening files
   !
   if (verbose>=4.and.Nspecies>0) print('(1x,a/)'),'Reading broadening parameters.'
   !
   !Used for indexing the gamma for the fast_voigt
   JmaxAll = 0
   if (allocated(jrot)) JmaxAll = nint(maxval(jrot)+0.5_rk)
   !
   do i =1,Nspecies
     !
     if ( trim(species(i)%filename)/="" ) then
       !
       write(ioname, '(a)') 'broadening file'
       !
       call IOstart(trim(ioname),bunit)
       open(unit=bunit,file=trim(species(i)%filename),action='read',status='old')
       !
       Jmax = JmaxAll
       !
       do
         !
         ! Scan and find Jmax
         read(bunit,*,end=14) ch_broad(1:3),gamma_,n_,J_
         !
         if (all(trim(ch_broad(1:2))/=(/'a0','a1','A0','A1'/))) cycle 
         !
         Jmax = max(Jmax,nint(J_+0.5_rk))
         JmaxAll = max(nint(J_+0.5_rk),JmaxAll)
         !
         cycle
         14  exit
         !
       enddo
       !
       if (Jmax>0) then
         !
         allocate(species(i)%gammaQN(0:Jmax,-2:2),stat=info)
         call ArrayStart('gammaQN',info,size(species(i)%gammaQN),kind(species(i)%gammaQN))
         species(i)%gammaQN(:,:) = species(i)%gamma
         !
         allocate(species(i)%nQN(0:Jmax,-2:2),stat=info)
         !         
         call ArrayStart('nQN',info,size(species(i)%nQN),kind(species(i)%nQN))
         species(i)%nQN(:,:) = species(i)%N
         !
        else
         !
         write(out,"('Empty .broad file or wrong format (Jmax = 0) in ',a)") trim(species(i)%filename)
         stop 'Jmax = 0 in a broadening file, empty .broad or wrong format'
         !
       endif
       !
       rewind(bunit)
       !
       ! read in the QN-dependent broadening parameters
       !
       Jmax = 0
       !
       do
         !
         ! read the line and check the model
         !
         call read_line(eof,bunit) ; if (eof) exit
         !
         call readu(ch_broad)
         !
         ! Ignore all models if not implemented 
         if (all(trim(ch_broad(1:2))/=(/'A0','A1','J ','JJ'/))) cycle 
         !  write(out,"('Error, illegal model in .broad: ',a3,' inly a0 and a1 are implemented')") ch_broad(1:3)
         !  stop 'Error, illegal model in .broad: '
         !endif
         !
         ! model specific read
         !
         select case ( trim(ch_broad(1:3)) )
           !
         case ('a0','A0')
           !
           if (nitems<4) call report ("not enough data in .broad for a0"//trim(w),.true.)
           !
           call readf(gamma_)
           call readf(n_)
           call readf(Ji)
           !
           ! use integer value of J
           !
           j = int(Ji)
           Jmax = max(Jmax,J)
           !
           species(i)%gammaQN(J,:) = gamma_
           species(i)%nQN(J,:) = n_
           !
         case ('a1','A1')
           !
           if (nitems<5) call report ("not enough data in .broad for a1"//trim(w),.true.)
           !
           call readf(gamma_)
           call readf(n_)
           call readf(Ji)
           call readf(Jf)
           !
           ! use integer value of J
           !
           Jp  = int(Jf)
           Jpp = int(Ji)
           !
           Jmax = max(Jmax,min(Jp,Jpp))
           !
           if (abs(Jp-Jpp)>1) then 
             write(out,'("Jp-Jpp>1 in broadening file:",2i8)') Jp,Jpp
             stop 'Jp-Jpp>1 in broadening file'
           endif
           !
           species(i)%gammaQN(Jpp,Jp-Jpp) = gamma_
           species(i)%nQN(Jpp,Jp-Jpp) = n_
           !
         case default
           !
           write(out,"('The broadening model',a,' is not implemented yet or illegal')") trim(ch_broad(1:3))
           stop 'Illegal broadening model'
           !
         end select
         !
       enddo
       !
       ! using the gamma(Jmax) value for all gamms for J>Jmax
       !
       forall(k=-2:2) species(i)%gammaQN(Jmax+1:,k) = species(i)%gammaQN(Jmax,k)
       forall(k=-2:2) species(i)%nQN(Jmax+1:,k) = species(i)%nQN(Jmax,k)
       !
       close(bunit)
       !
    elseif (stick_hitran.or.stick_oxford) then
       !
       Jmax = JmaxAll
       !
       allocate(species(i)%gammaQN(0:Jmax,-2:2),stat=info)
       call ArrayStart('gammaQN',info,size(species(i)%gammaQN),kind(species(i)%gammaQN))
       species(i)%gammaQN(:,:) = species(i)%gamma
       !
       allocate(species(i)%nQN(0:Jmax,-2:2),stat=info)
       !         
       call ArrayStart('nQN',info,size(species(i)%nQN),kind(species(i)%nQN))
       species(i)%nQN(:,:) = species(i)%N
       !
     endif
     !
   enddo
   !
   allocate(gamma_idx(0:JmaxAll,-2:2),stat=info)
   call ArrayStart('gamma_idx',info,size(gamma_idx),kind(gamma_idx))
   allocate(gamma_comb(0:JmaxAll,-2:2),stat=info)
   call ArrayStart('gamma_comb',info,size(gamma_comb),kind(gamma_comb))
   !
   gamma_idx = 1
   gamma_comb = -1
   do i=0,JmaxAll
     do j= max(0,i-2),min(JmaxAll,i+2)
      gamma_comb(j,i-j) = get_Voigt_gamma_val(Nspecies,real(i,rk),real(j,rk))
     enddo
   enddo
   !
   if(proftype(1:5) == 'VOI-F') then
     !    
     gamma_idx = 1
     !
     !Combine all the gammas
     do i=0,JmaxAll
      do j= max(0,i-2),min(JmaxAll,i+2) 
       call fast_voigt%generate_indices(gamma_comb(j,i-j),gamma_idx(j,i-j),JmaxAll)
      enddo
     enddo  
     !   
     !gamma_idx = -1
     !
     write(out,"('Generated ',i8,' fast voigt grids')") fast_voigt%get_size()
     !stop "Not yet implemented for VOI-F"
     !Lets perform a precomputation of all species involved
     !
   endif
   !
   ! Intensities
   !
   intens = 0
   intband = 0
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
   case ('GAUSS','DOPPL','LOREN','GAUS0','DOPP0','VOIGT','PSEUD','PSE-R','PSE-L','VOI-Q','VOI-F','VOI-9')
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
          !
       endif
       !
   case ('RECT','BOX','MAX')
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
       !
       if (stick_hitran.or.stick_oxford) then 
         open(unit=sunit,file=trim(output)//".par",action='write',status='replace')
       else
         open(unit=sunit,file=trim(output)//".stick",action='write',status='replace')
       endif
       if (verbose>=2) then
          !
          if (cutoff_model=="HITRAN") then
             write(out,"(10x,/'Stick pectra of type ',a,' with the HITRAN cut-off model, Scrit = ',e18.5,' nu_crit = ',f16.4)") &
             trim(proftype),S_crit,nu_crit
          else
             write(out,"(10x,/'Stick spectra of type ',a,' stronger than ',e18.5)") trim(proftype),thresh
          endif
          write(out,"(10x,'Range = ',f18.7,'-',f18.7)") freql,freqr
          write(out,"(10x,'Temperature = ',f18.7)") temp
          write(out,"(10x,'Partition function = ',f17.4)") partfunc
          write(out,"(10x,'Spectrum type = ',a/)") trim(specttype)
          !
       endif
       !
   case ('TRANS')
       !
       if (verbose>=2) then
          !
          write(out,"(10x,a)") 'Re-print .trans file'
          !
          if (cutoff_model=="EXP") then
             write(out,"(10x,/'Strong/weak lines ',a,' with the EXP cut-off model, Scrit = ',e18.5,' nu_crit = ',f16.4)") &
             trim(proftype),S_crit,nu_crit
          endif
          !
       endif
       !
   end select
   !
   if (lineprofile_do) then
      write(out,"(10x,'Pressure = ',e18.7)") pressure
      write(out,"(10x,'Voigt parameters:   gamma       n         T0            P0     Ratio')")
      do i =1,Nspecies
        write(out,"(21x,a,4f12.4,f13.6)") trim(species(i)%name),species(i)%gamma,species(i)%n,species(i)%t0,&
                                          species(i)%p0,species(i)%ratio
      enddo
   endif
   !
   if ( trim(proftype)=='PHOENIX' ) then
     !
     if (verbose>=4) write(out,"('    iel          wv   ener     gf gamma1     n1 gamma2     n2')") 
     !
     if (trim(specttype)/="GF") then 
       write(out,"('Change the spectral type to GF (oscillator strength) for Phoenix'/)")
       specttype = "GF"
     endif
     !
   endif
   !
   if (hitran_do.and.partfunc<0.0) then
     write(out,"('For HITRAN/OXFORD partition function must be defined using PF or QSTAT keywords')")
     stop 'Undefined PF'
   endif
   !
   if (spectra_do.and.(partfunc<0.0.or.partfunc_ref<0.0)) then
     write(out,"('For SPECTRA pf and ref-pf must be defined using PF XXXX REF YYYY')")
     stop 'Undefined PF'
   endif
   !
   if (spectra_do.and.(temp_ref<0.0)) then
     write(out,"('For SPECTRA reference temperatures must be defined using temperature XXXX REF YYYY')")
     stop 'Undefined reference Temperature'
   endif
   !
   if (any( trim(proftype(1:3))==(/'DOP','GAU','REC','BIN','BOX','LOR','VOI','MAX','PSE','COO'/)) ) then
     allocate(intens_omp(npoints,N_omp_procs),stat=info)
     intens_omp = 0
     call ArrayStart('swap:intens_omp',info,size(intens_omp),kind(intens_omp))
   endif
   !
   ! estimate how many transitions we can out into RAM
   !
   if (N_to_RAM<0) then
     !
     mem_t = real(6*rk+5*ik,rk)/1024_rk**3
     !
     N_to_RAM = max(100,min( max_transitions_to_ram,int((memory_limit-memory_now)/mem_t,hik),hik))
     !
     if (verbose>=4) print('("We can cache ",i12," transitions into RAM")'),N_to_RAM
     !
   endif
   !
   ! allocating all different arrays to keep the data in RAM
   !
   allocate(acoef_RAM(N_to_RAM),indexf_RAM(N_to_RAM),indexi_RAM(N_to_RAM),stat=info)
   call ArrayStart('swap:acoef_RAM',info,size(acoef_RAM),kind(acoef_RAM))
   call ArrayStart('swap:indexf_RAM',info,size(indexf_RAM),kind(indexf_RAM))
   call ArrayStart('swap:indexi_RAM',info,size(indexi_RAM),kind(indexi_RAM))
   !
   allocate(abscoef_ram(N_to_RAM),nu_ram(N_to_RAM),stat=info)
   call ArrayStart('swap:abscoef_ram',info,size(abscoef_ram),kind(abscoef_ram))
   call ArrayStart('swap:nu_ram',info,size(nu_ram),kind(nu_ram))
   !
   allocate(ileveli_RAM(N_to_RAM),ilevelf_RAM(N_to_RAM),stat=info)
   call ArrayStart('swap:ileveli_RAM',info,size(ileveli_RAM),kind(ileveli_RAM))
   call ArrayStart('swap:ilevelf_RAM',info,size(ilevelf_RAM),kind(ilevelf_RAM))
   !
   if (lineprofile_do) then
     allocate(gamma_ram(N_to_RAM),stat=info)
     allocate(gamma_idx_RAM(N_to_RAM),stat=info)
     call ArrayStart('swap:gamma_ram',info,size(gamma_ram),kind(gamma_ram))
   endif
   !
   if (verbose>=4) call MemoryReport
   !
   nswap = 0
   Noffset = max(nint(offset/dfreq,ik),1)
   !
   if (verbose>=3) print('(/"Computing ",a," with ",a," ... "/)'), trim(specttype),trim(proftype)
   !
   if (verbose>=3) print('(/"Transition files ... "/)')
   !
   loop_file : do i = 1,nfiles
     !
     open(unit=tunit,file=trim(intfilename(i)),action='read',status='old',iostat=iostat_)
     !
     eof = .false. ! this will be used to control if the end of file is reached
     !
     if (verbose>=3) print('("  Processing ",a,"...")'), trim(intfilename(i))
     !
     if (iostat_/=0) then 
       print*,"Error opening .trans file (istat is not 0, ",iostat_,")  ",trim(intfilename(i))
       stop 'Error opening .trans file'
     endif
     !
     select case (proftype(1:5)) 
     case('TRANS')
       !
       write(ioname, '(a)') 'Trans file strong'
       call IOstart(trim(ioname),sunit)
       open(unit=sunit,file=trim(intfilename(i))//trim(output)//".s",action='write',status='replace')
       write(ioname, '(a)') 'Trans file weak'
       call IOstart(trim(ioname),wunit)
       open(unit=wunit,file=trim(intfilename(i))//trim(output)//".w",action='write',status='replace')
       !
     end select
     !
     ichunk = 0
     !
     loop_tran: do
        !
        ! finis this if the eof has been reached 
        !
        if (eof) exit loop_tran
        !
        nswap_ = N_to_RAM
        ichunk = ichunk + 1 ! counting the courent swap
        !
        if (histogram) then
           !
           ! using pre-computed integrated intensities for a given Temperature and bin
           !
           iswap_ = 0
           !
           do while(iswap_<N_to_RAM)
             !
             !   read new line from intensities file
             read(tunit,*,end=118) tranfreq,abscoef
             !
             if (tranfreq<small_) cycle
             !
             if (tranfreq<freql.or.tranfreq>freqr) cycle
             !
             iswap_ = iswap_ + 1
             !
             nu_ram(iswap_) = tranfreq
             abscoef_RAM(iswap_) = abscoef
             !
             if (lineprofile_do) then
               gamma_ram(iswap_)=halfwidth
               gamma_idx_RAM(iswap_) = 0
             endif
             !
             cycle
             !
           118 continue
             nswap = iswap_-1
             eof = .true.
             exit
           enddo
           !
           nswap = iswap_
           !
           intband = intband + sum(abscoef_RAM(1:nswap))
           !
           if (nswap<1) then
              if (verbose>=5.and.proftype(1:5)/="STICK") write(out,"('... Finish swap')")
              exit loop_tran
           endif
           !
           ! for pre-computed J-dependent integrated intensities
           if (histogramJ) then
             stop 'histogramJ is currently not working'
             Ji = IntJvalue(i)
             Jf = Ji
           endif
           !
        elseif (hitran_do) then
           !
           ! HITRAN format
           !
           if(any(trim(proftype(1:5))==(/'PHOEN'/))) then
             print('(a,2x,a)'),'Illegal proftype HITRAN',trim(proftype)
             stop 'Illegal proftype for HITRAN'
           endif
           !
           Nspecies_ = max(Nspecies,2)
           !
           if (stick_oxford) Nspecies_ = max(Nspecies,3)
           !
           iswap_ = 0
           !
           do while(iswap_<N_to_RAM)
             !
             !read(tunit,*,end=119) indexf_RAM(iswap),indexi_RAM(iswap),acoef_RAM(iswap)
             read(tunit,"(i2,i1,f12.6,e10.3,e10.3,f5.4,f5.3,f10.4,f4.2,f8.6,a55,23x,2f7.1)",end=119) &
                  imol,isotope,tranfreq,abscoef,acoef,gamma_,gamma_s,energyi,n_,delta_air,ch_q,gf,gi
             !
             tranfreq = tranfreq + delta_air * pressure
             energyf = tranfreq + energyi
             !
             if (imol/=imolecule.or.iso/=isotope) cycle
             !
             if (tranfreq<small_) cycle
             !
             if (tranfreq<freql.or.tranfreq>freqr) cycle
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
               ! the reference partition function is undefined, we use the standard expression 
               if (nint(partfunc_ref)==-1) then
                  !
                  abscoef = cmcoef*acoef*gf*exp(-beta*energyi)*(1.0_rk-exp(-beta*tranfreq))/(tranfreq**2*partfunc)
                  !
                else
                  !
                  abscoef=abscoef*exp(-beta    *energyi)*(1.0_rk-exp(-beta    *tranfreq))*partfunc_ref/&
                                 (exp(-beta_ref*energyi)*(1.0_rk-exp(-beta_ref*tranfreq))*partfunc   ) 
               endif
               !
             case ('EMISSION')
               !
               ! emission coefficient [Ergs/mol/Sr]
               !
               abscoef=emcoef*acoef*gf*exp(-beta*energyf)*tranfreq/(partfunc)
               !
             end select
             !
             !
             if (abscoef<thresh) cycle
             !
             iswap_ = iswap_ + 1
             !
             abscoef_ram(iswap_)=abscoef
             !
             nu_ram(iswap_) = tranfreq
             !
             ! estimate the line shape parameter
             !
             if (lineprofile_do) then
             !if (lineprofile_do.and..not.species(1)%if_defined) then
               !
               ! The current version of ExoCross does no know how to locate J-values in HITRAN-input.
               ! Therefore we can only use the static values of Voigt-parameters together with Nspecies
               !
               Jf = 0.0_rk
               Ji = 0.0_rk
               species(1)%gamma = gamma_
               species(1)%n = n_
               species(2)%gamma = gamma_s
               !species(2)%n = n_
               !
               gamma_ram(iswap_) = gamma_ * pressure * (hitran_Tref/temp)**n_*species(1)%ratio
               if (Nspecies>1) gamma_ram(iswap_) = gamma_ram(iswap_) + gamma_s * pressure * species(2)%ratio
               !
               !Below: probably not correct, but gamma_idx_RAM definitely needs
               !to be initialized somehow
               gamma_idx_RAM(iswap_) = fast_voigt%get_size()
               !
             endif
             !
             cycle
             !
           119 continue
             nswap = iswap_-1
             eof = .true.
             exit
           enddo
           !
           nswap = iswap_
           !
           intband = intband + sum(abscoef_RAM(1:nswap))
           !
           if (nswap<1) then
              if (verbose>=5) write(out,"('... Finish swap')")
              exit loop_tran
           endif
           !
        elseif (spectra_do) then
           !
           ! SPECTRA format (ioa tomsk database)
           !
           iswap_ = 0
           !
           do while(iswap_<N_to_RAM)
             !
             !read(tunit,*,end=119) indexf_RAM(iswap),indexi_RAM(iswap),acoef_RAM(iswap)
             read(tunit,*,end=120) &
                  imol,isotope,tranfreq,abscoef,energyi,gamma_,n_
             !
             energyf = tranfreq + energyi
             !
             if (imol/=imolecule.or.iso/=isotope) cycle
             !
             if (tranfreq<small_) cycle
             !
             if (tranfreq<freql.or.tranfreq>freqr) cycle
             !
             iswap_ = iswap_ + 1
             !
             abscoef_ram(iswap_)=abscoef*exp(-beta    *energyi)*(1.0_rk-exp(-beta    *tranfreq))*partfunc_ref/&
                                        (exp(-beta_ref*energyi)*(1.0_rk-exp(-beta_ref*tranfreq))*partfunc   ) 
             !
             nu_ram(iswap_) = tranfreq
             !
             ! estimate the line shape parameter
             !
             if (lineprofile_do) then
               !
               ! The current version of ExoCross does no know how to locate J-values in HITRAN-input.
               ! Therefore we can only use the static values of Voigt-parameters together with Nspecies
               !
               Jf = 0.0_rk
               Ji = 0.0_rk
               species(1)%gamma = gamma_
               species(1)%n = n_
               !
               gamma_ram(iswap_) = gamma_ * pressure * (temp_ref/temp)**n_

               !Below: probably not correct, but gamma_idx_RAM definitely needs
               !to be initialized somehow
               gamma_idx_RAM(iswap_) = fast_voigt%get_size()
               !
             endif
             !
             cycle
             !
           120 continue
             nswap = iswap_-1
             eof = .true.
             exit
           enddo
           !
           nswap = iswap_
           !
           intband = intband + sum(abscoef_RAM(1:nswap))
           !
           if (nswap<1) then
              if (verbose>=5) write(out,"('... Finish swap')")
              exit loop_tran
           endif
           !
        else ! ExoMol format of the trans-file
           !
           do iswap = 1,N_to_RAM
             !
             read(tunit,*,end=121) indexf_RAM(iswap),indexi_RAM(iswap),acoef_RAM(iswap)
             !
             cycle
             !
           121 continue
             !
             nswap_ = iswap-1
             eof = .true.
             exit
             !
           enddo
           !
           if (nswap_<1) then
              !
              if (verbose>=5) write(out,"('... Finish swap')")
              !
              exit loop_tran
              !
              stop 'Error nswap_=0'
           endif
           !
           if (verbose>=5.and.nswap_<N_to_RAM.and.proftype(1:5)/="STICK") print('(a,i12,a)'), &
                                                  "Now computing ... ",nswap_," transitions..."
           !
           !read(tunit,*,end=20) indexf,indexi,acoef
           iswap_ = 1
           !
           ileveli_ram = -1
           abscoef_ram = 0
           !
           !$omp  parallel do private(iswap,indexf,indexi,acoef,ilevelf,ileveli,energyf,energyi,ifilter,&
           !$omp& ivib,ener_vib,ener_rot,jf,ji,tranfreq,tranfreq0,offset,abscoef,cutoff)&
           !$omp& schedule(static) shared(ilevelf_ram,ileveli_ram,abscoef_ram,acoef_ram,nu_ram,Asum,gamma_ram)
           loop_swap : do iswap = 1,nswap_
             !
             indexf = indexf_RAM(iswap)
             indexi = indexi_RAM(iswap)
             acoef = acoef_RAM(iswap)
             !
             if (indexf>nlevels.or.indexi>nlevels) cycle loop_swap
             !
             ilevelf = indices(indexf)
             ileveli = indices(indexi)
             !
             if (ilevelf==0.or.ileveli==0) cycle loop_swap
             !
             energyf = energies(ilevelf)
             energyi = energies(ileveli)
             !
             if (filter) then
               !
               do ifilter = 1,Nfilters
                 !
                 if (upper(ifilter)%mask/=trim(quantum_numbers(upper(ifilter)%i,ilevelf)).and.&
                     trim(upper(ifilter)%mask)/="") cycle loop_swap
                 if (lower(ifilter)%mask/=trim(quantum_numbers(lower(ifilter)%i,ileveli)).and.&
                     trim(lower(ifilter)%mask)/="") cycle loop_swap
                 !
               enddo
               !
             endif
             !
             if (vibtemperature_do) then
                !
                if (trim(specttype)=='ABSORPTION') then 
                  ivib = ivib_state(ileveli)
                  ener_vib = energies_vib(ivib)
                  ener_rot = energyi-ener_vib
                else
                  ivib = ivib_state(ilevelf)
                  ener_vib = energies_vib(ivib)
                  ener_rot = energyf-ener_vib
                endif
             endif
             !
             jf = jrot(ilevelf)
             ji = jrot(ileveli)
             !
             tranfreq = energyf-energyi
             !
             tranfreq = tranfreq + lineshift
             !
             tranfreq0 = tranfreq
             !
             if (microns) then
               offset = offset/tranfreq**2
               tranfreq0 = 10000.0_rk/(tranfreq+small_)
             endif
             !
             if (tranfreq0<freql.or.tranfreq0>freqr) cycle
             !
             if (energyf-energyi<-1e1) then
               write(out,"('Error Ei>Ef+10: i,f,indi,indf,Aif,Ef,Ei = ',4i12,2x,3es16.8)") & 
                    ilevelf,ileveli,indexf,indexi,acoef,energyf,energyi
               stop 'wrong order of indices'
               cycle loop_swap
             elseif (energyf-energyi<small_) then
               cycle loop_swap
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
               !abscoef=cmcoef*acoef*gtot(ilevelf)*exp(-beta*energyi)*(1.0_rk-exp(-beta*tranfreq))/(tranfreq**2*partfunc)
               !
               if (.not.vibtemperature_do) then
                  abscoef=cmcoef*acoef*gtot(ilevelf)*exp(-beta*energyi)*(1.0_rk-exp(-beta*tranfreq))/(tranfreq**2*partfunc)
               elseif(vibpopulation_do) then
                  ! apply vibrational population
                  read(quantum_numbers(QN%dens_col-4,ileveli),*) ndensity
                  abscoef=cmcoef*acoef*gtot(ilevelf)*exp(-c2/temp*ener_rot)*ndensity*&
                          (1.0_rk-exp(-c2/ener_rot*tranfreq))/(tranfreq**2*partfunc)
                  ! split into a product of vib and rot parts 
               else
                  abscoef=cmcoef*acoef*gtot(ilevelf)*exp(-c2/temp*ener_rot)*exp(-c2/temp_vib*ener_vib)*&
                          (1.0_rk-exp(-c2/temp_vib*tranfreq))/(tranfreq**2*partfunc)
               endif 
               !
             case ('EMISSION')
               !
               ! emission coefficient [Ergs/mol/Sr]
               !
               if (.not.vibtemperature_do) then
                  abscoef=emcoef*acoef*gtot(ilevelf)*exp(-beta*energyf)*tranfreq/(partfunc)

               elseif(vibpopulation_do) then
                  ! apply vibrational population
                  read(quantum_numbers(QN%dens_col-4,ileveli),*) ndensity
                  !
                  abscoef=emcoef*acoef*gtot(ilevelf)*exp(-c2/temp*ener_rot)*ndensity*&
                          tranfreq/(partfunc)
                  !
               else 
                  ! split into a product of vib and rot parts 
                  abscoef=emcoef*acoef*gtot(ilevelf)*exp(-c2/temp*ener_rot)*&
                          exp(-c2/temp_vib*ener_vib)*tranfreq/(partfunc)

                  !
               endif 
               !
             case ('GF')
               !
               ! oscillator strength gf; abscoef is to check agains the intensity thresholds 
               !
               abscoef=cmcoef*acoef*gtot(ilevelf)*exp(-beta*energyi)*(1.0_rk-exp(-beta*tranfreq))/(tranfreq**2*partfunc)
               acoef=gfcoef*acoef*gtot(ilevelf)/(vellgt*tranfreq)**2
               !
             case ('LIFETIME','T-LIFETIME')
               !
               ilevelf_ram(iswap) = ilevelf
               ileveli_ram(iswap) = ileveli
               abscoef_ram(iswap) = 0
               acoef_ram(iswap) = acoef
               !
               nu_ram(iswap) = 0
               !
               cycle loop_swap
               !
             end select
             !
             !if (abscoef<abscoef_thresh) cycle loop_swap
             !
             cutoff = thresh
             !
             select case (trim(cutoff_model))
                !
             case ("EXP")
                !
                continue
                !
             case ("HITRAN") 
                !
                if (tranfreq<=nu_crit) then
                  cutoff = S_crit*tranfreq/nu_crit*tanh(0.5_rk*c2/temp*tranfreq)
                else
                  cutoff = S_crit
                endif
                !
                if (abscoef<cutoff) then
                  cycle loop_swap
                endif
                !
             case default
               !
               if (abscoef<cutoff) then
                 cycle loop_swap
               endif
               !
             end select
             !
             ilevelf_ram(iswap) = ilevelf
             ileveli_ram(iswap) = ileveli
             abscoef_ram(iswap) = abscoef
             acoef_ram(iswap) = acoef
             !
             nu_ram(iswap) = tranfreq
             !
             if (lineprofile_do) then
                gamma_ram(iswap)=get_Voigt_gamma_n(Nspecies,Jf,Ji)
                temp_gamma_n = get_Voigt_gamma_n(Nspecies,Jf,Ji,gamma_idx_RAM(iswap))
             endif
             !
             !do do_log_intensity
             !
           enddo loop_swap
           !$omp end parallel do
           !
           iswap_ = 1
           !
           do iswap = 1,nswap_
             !
             if (ileveli_ram(iswap)>0) then 
                ilevelf_ram(iswap_) = ilevelf_ram(iswap)
                ileveli_ram(iswap_) = ileveli_ram(iswap)
                abscoef_ram(iswap_) = abscoef_ram(iswap)
                acoef_ram(iswap_)   = acoef_ram(iswap)
                !
                nu_ram(iswap_) = nu_ram(iswap)
                !
                if (lineprofile_do) then
                  gamma_ram(iswap_)=gamma_ram(iswap)
                  gamma_idx_RAM(iswap_)=gamma_idx_RAM(iswap)
                endif
                !
                iswap_ = iswap_ + 1
                !
             endif
             !
           enddo
           nswap = iswap_-1
           !
           intband = intband + sum(abscoef_ram(1:nswap))
           !
        endif
        !
        !if (offset<0) offset = ioffset*halfwidth
        !
        !   if transition frequency is out of selected range
        !
        call TimerStart('Calc')
        !
        select case (trim(proftype(1:5)))
            !
        case ('STICK')
            !
            if (hitran_do) then 
               !
               do iomp = 1,N_omp_procs
                 !
                 do iswap = iomp,nswap,N_omp_procs
                   !
                   abscoef = abscoef_ram(iswap)
                   tranfreq = nu_ram(iswap)
                   !
                   ! write to .stick-file
                   write(sunit,my_fmt) tranfreq,abscoef
                   !
                 enddo
                 !
               enddo
               !
            else
               !
               call do_stick(sunit,nswap,nlines,maxitems,energies,Jrot,gtot,ilevelf_ram,ileveli_ram,abscoef_ram,acoef_ram,&
                            nchars_quanta,quantum_numbers)
               !
            endif
            !
            call TimerStop('Calc')
            cycle loop_tran
            !
        case ('VOIGT')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,halfwidth) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                halfwidth = gamma_ram(iswap)
                !
                call do_Voigt(tranfreq,abscoef,use_resolving_power,freq,halfwidth,dpwcoef,offset,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do            
            !
        case ('GAUS0')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                !
                call do_gauss0(tranfreq,abscoef,dfreq,freq,halfwidth,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('DOPPL')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,halfwidth0) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                halfwidth0=dpwcoef*tranfreq
                if (halfwidth0<100.0_rk*small_) cycle
                !
                call do_gauss(tranfreq,abscoef,dfreq,freq,halfwidth0,offset,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('DOPP0')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,halfwidth0) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                halfwidth0=dpwcoef*tranfreq
                if (halfwidth0<100.0_rk*small_) cycle
                !
                call do_gauss0(tranfreq,abscoef,dfreq,freq,halfwidth0,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('GAUSS')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                !
                call do_gauss(tranfreq,abscoef,dfreq,freq,halfwidth,offset,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('LOREN')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,halfwidth) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                halfwidth = gamma_ram(iswap)
                !
                call do_lorentz(tranfreq,abscoef,dfreq,freq,halfwidth,offset,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('LIFET')
            !
            do iswap = 1,nswap
              !
              ilevelf = ilevelf_ram(iswap)
              energyf = energies(ilevelf)
              !
              if (energyf>enermax) cycle
              !
              acoef = acoef_ram(iswap)
              !
              if (Asum(ilevelf)<0) Asum(ilevelf) = 0
              !
              Asum(ilevelf) = Asum(ilevelf) + acoef
              !
            enddo
            !
        case ('T-LIF')
            !
            do iswap = 1,nswap
              !
              ilevelf = ilevelf_ram(iswap)
              energyf = energies(ilevelf)
              !
              ileveli = ileveli_ram(iswap)
              energyi = energies(ileveli)
              !
              if (energyf>enermax.or.energyi>energy) cycle
              !
              acoef = acoef_ram(iswap)
              !
              if (Asum(ilevelf)<0) Asum(ilevelf) = 0
              !
              Asum(ilevelf) = Asum(ilevelf) + acoef+acoef/(exp(c2/temp*tranfreq)-1.0_rk)
              !
              if (Asum(ileveli)<0) Asum(ileveli) = 0
              !
              Asum(ileveli) = Asum(ileveli) + acoef/(exp(c2/temp*tranfreq)-1.0_rk)*gtot(ileveli)/gtot(ilevelf)
              !
            enddo
            !
        case ('TRANS')
            !
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                acoef = acoef_ram(iswap)
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                !
                ilevelf = ilevelf_ram(iswap)
                ileveli = ileveli_ram(iswap)
                !
                cutoff = S_crit*exp(-tranfreq/nu_crit)
                !
                if (abscoef>=cutoff) then
                    write(sunit,my_fmt) ilevelf,ileveli,acoef,tranfreq
                else
                    write(wunit,my_fmt) ilevelf,ileveli,acoef,tranfreq
                endif
                !
              enddo
              !
            enddo
            !
            call TimerStop('Calc')
            cycle loop_tran
            !
        case ('PHOEN')
            !
            if (Nspecies/=2) then 
              write(out,"('Phoenix-Error: illegal Nspecies must be 2 not ',i4)") Nspecies
              stop 'Phoenix-Error: illegal Nspecies '
            endif
            !
            call do_gf_oscillator_strength_Phoenix(i,ichunk,iso,nswap,nlines,energies,Jrot,ilevelf_ram,ileveli_ram,&
                                                   acoef_ram,abscoef_ram,jmax,species(1)%gammaQN,&
                                                   species(1)%nQN,species(2)%gammaQN,species(2)%nQN,output)
            !
            call TimerStop('Calc')
            cycle loop_tran
            !
        case ('COOLI')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                ilevelf = ilevelf_ram(iswap)
                energyf = energies(ilevelf)
                abscoef = abscoef_ram(iswap)
                !
                do itemp = 1,npoints
                  !
                  temp0 = real(itemp,rk)*dtemp
                  !
                  beta0 = c2/temp0
                  !
                  intens_omp(itemp,iomp) = intens_omp(itemp,iomp) + abscoef*exp(-(beta0-beta)*energyf)*partfunc/(pf(0,itemp))
                  !
                enddo
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('PSEUD')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,halfwidth) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                halfwidth = gamma_ram(iswap)
                !
                call do_pseud(tranfreq,abscoef,dfreq,freq,halfwidth,offset,freql,dpwcoef,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
            !
        !case ('PSE-R') ! pseudo_Voigt_Rocco_Cruzado_ActaPhysPol_2012.pdf
            !
            !
        !case ('PSE-L') ! pseudo_Voigt_Liu_Lin_JOptSocAmB_2001.pdf
            !
            !
        case ('VOI-Q') ! VOIGT-QUADRATURES
            !
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,halfwidth) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                halfwidth = gamma_ram(iswap)
                !
                call do_Voi_Q(tranfreq,abscoef,dfreq,freq,abciss,weight,halfwidth,offset,freql,dpwcoef,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('VOI-9') ! VOIGT-Alg-916
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,halfwidth) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                halfwidth = gamma_ram(iswap)
                !
                call do_Voigt_916(tranfreq,abscoef,dfreq,freq,halfwidth,dpwcoef,offset,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('VOI-F')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,halfwidth) shared(intens_omp,fast_voigt) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                halfwidth = gamma_ram(iswap)
                !
                call do_Voigt_Fast(tranfreq,abscoef,use_resolving_power, &
                     freq,gamma_idx_RAM(iswap),dpwcoef,offset,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('MAX')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,cutoff,ipoint) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                !
                cutoff =  apply_HITRAN_cutoff(tranfreq)
                !
                if (abscoef<cutoff) cycle
                !
                call get_grid_ipoint(tranfreq,freq,ipoint)
                !
                !ipoint =  max(nint( ( tranfreq-freql)/dfreq )+1,1)
                !
                intens_omp(ipoint,iomp) = max(intens_omp(ipoint,iomp),abscoef)
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('RECT','BOX')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,ipoint) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                !
                call do_rect(tranfreq,abscoef,dfreq,freq,offset,freql,intens_omp(:,iomp))
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('BIN')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,ipoint) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                !
                !cutoff =  apply_HITRAN_cutoff(tranfreq)
                !
                !if (abscoef<cutoff) cycle
                !
                call get_grid_ipoint(tranfreq,freq,ipoint)
                !
                !ipoint =  max(nint( ( tranfreq-freql)/dfreq )+1,1)
                !
                intens_omp(ipoint,iomp) = intens_omp(ipoint,iomp)+abscoef
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        case ('BIN-MICRON');
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,ipoint) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                !
                tranfreq = 10000.0_rk/tranfreq
                !
                ipoint =  max(nint( ( tranfreq-freql)/dfreq )+1,1)
                !
                intens_omp(ipoint,iomp) = intens_omp(ipoint,iomp)+abscoef
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
            !
        case ('BIN-R')
            !
            !$omp parallel do private(iomp,iswap,abscoef,tranfreq,ipoint) shared(intens_omp) schedule(dynamic)
            do iomp = 1,N_omp_procs
              !
              do iswap = iomp,nswap,N_omp_procs
                !
                abscoef = abscoef_ram(iswap)
                tranfreq = nu_ram(iswap)
                !
                ipoint =  nint(real((log(tranfreq)-log(freql))/resolving_f,rk))+1
                if (ipoint<1.or.ipoint>npoints) cycle
                !
                intens_omp(ipoint,iomp) = intens_omp(ipoint,iomp)+abscoef
                !
              enddo
              !
            enddo
            !$omp end parallel do
            !
        end select
        !       
        call TimerStop('Calc')
        !
        ! will be used to check duplicates
        !
        !!!indexf_ = indexf ; indexi_ = indexi ; acoef_ = acoef
        !
        !cycle loop_tran
        !!
        !333 print*,iostat_
        !!
        !stop
        !
     enddo loop_tran
     !
     select case (proftype(1:5)) 
     case('TRANS')
       close(tunit)
       close(wunit)
     end select
     !
     if (nswap_<nswap) cycle loop_file
     !
   enddo loop_file
   !
   !Do all the summation at the end
   !
   if (any( trim(proftype(1:3))==(/'DOP','GAU','REC','BIN','BOX','LOR','VOI','PSE','COO'/)) ) then
     !
     do i=1,N_omp_procs
      intens(:) = intens(:) + intens_omp(:,i)
     enddo
     !
   else if (trim(proftype(1:3))== 'MAX') then
     do iomp = 1,N_omp_procs
       do ipoint = 1,npoints
         intens(ipoint) = max(intens_omp(ipoint,iomp),intens(ipoint))
       enddo
     enddo  
   endif  
   !
   call IOstop(trim(ioname))
   !
   if (trim(proftype)=='LIFETIME'.or.trim(proftype)=='T-LIFETIME') THEN
     !
     write(ioname, '(a)') 'Life times'
     call IOstart(trim(ioname),tunit)
     !
     open(unit=tunit,file=trim(output)//".life",action='write',status='replace')
     !
     if ( mod(nint(2.0_rk*jrot(1)),2)==0 ) then 
       write(my_fmt,'(a)')  '(1x,i11,1x,f12.6,1x,i6,1x,i7,1x,es12.4,3x)'
     else
       write(my_fmt,'(a)')  '(1x,i11,1x,f12.6,1x,i6,1x,f7.1,1x,es12.4,3x)'
     endif
     !
     do ilevelf = 1,nlines
       !
       ! write to .life-file
       !
       if ( mod(nint(2.0_rk*jrot(1)),2)==0 ) then 
          write(tunit,my_fmt,advance="no") indices(ilevelf),energies(ilevelf),gtot(ilevelf),nint(jrot(ilevelf)),1.0_rk/Asum(ilevelf)
       else
          write(tunit,my_fmt,advance="no") indices(ilevelf),energies(ilevelf),gtot(ilevelf),jrot(ilevelf),1.0_rk/Asum(ilevelf)
       endif
       !
       do kitem = 1,maxitems
         !
         !l = len(trim(quantum_numbers(kitem,ilevelf)))
         !
         !b_fmt = "(1x,a3)" ; if (l>3) b_fmt = "(1x,a8)"
         !
         if (nchars_quanta(kitem)<10) then
            write(b_fmt,"('(1x,a',i1,')')") nchars_quanta(kitem)
         else
            write(b_fmt,"('(1x,a',i2,')')") nchars_quanta(kitem)
         endif
         !
         write(tunit,b_fmt,advance="no") trim(quantum_numbers(kitem,ilevelf))
         !
       enddo
       !
       write(tunit,"(a1)",advance="yes") " "
       !
     enddo
     !
     close(tunit,status='keep')
     !
   elseif (trim(proftype)=='COOLING') then
     !
     write(ioname, '(a)') 'Cooling'
     call IOstart(trim(ioname),tunit)
     !
     open(unit=tunit,file=trim(output)//".cooling",action='write',status='replace')
     !
     do itemp = 1,npoints
       !
       if (energy>enermax) cycle
       !
       temp0 = real(itemp,rk)*dtemp
       !
       write(tunit,"(1x,f12.3,1x,es20.8)") temp0,intens(itemp)
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
     if (microns) then
       !
       ! change from microns to wavenumbers
       !
       do ipoint =  npoints,1,-1
          !
         write(tunit,'(2(1x,es16.8))') 10000.0_rk/freq(ipoint),intens(ipoint)
         !
       enddo
       !
     elseif (use_resolving_power.and.trim(proftype(1:3))=='BIN') then
       !
       do ipoint =  1,npoints
         !
         tranfreq0 = real(ipoint-1,rk)*resolving_f+log(freql)
         !
         tranfreq0 = exp(tranfreq0)
         !
         write(tunit,'(2(1x,es16.8))') tranfreq0,intens(ipoint)
         !
       enddo
       !
     else
       !
       write(tunit,'(2(1x,es16.8E3))') (freq(ipoint),intens(ipoint),ipoint=1,npoints)
       !
     endif
     !
     ! no need to scale with dfreq for bin or max
     if (any( trim(proftype(1:3))==(/'BIN','MAX'/)) ) dfreq = 1.0d0
     !
     if (verbose>=2) print('(/"Total intensity  (sum):",es16.8," (int):",es16.8)'), intband,sum(intens)*dfreq
     !
     close(tunit,status='keep')
     !
   elseif (trim(proftype)=='PHOENIX') then
     !
     !write(gfunit,*) int(-1,kind=2),int(-1,kind=4),int(-1,kind=2),int(-1,kind=2),int(-1,kind=2),int(-1,kind=2),int(-1,kind=2)
     !close(gfunit,status='keep')
     !
   else
     !
     if (verbose>=2) print('(/"Total intensity = ",es16.8,", T = ",f18.4)'), intband,temp
     !
   endif
   !
   !call IOStop(trim(ioname))
   !
   if (verbose>=-2.and.any(Nintens(:)/=0) ) then
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
        write(tunit,'(i4,2x,i12)') i,Nintens(i)
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
   if (allocated(energies)) call ArrayStop('energies')
   if (allocated(Jrot)) call ArrayStop('Jrot')
   if (allocated(gtot)) call ArrayStop('gtot')
   if (specttype(1:4)=='PART') call ArrayStop('pf')
   !
  end subroutine intensity
  !
  !
  function apply_HITRAN_cutoff(tranfreq) result(cutoff)
      !
      implicit none
      !
      real(rk),intent(in) :: tranfreq
      real(rk)  :: cutoff
      !
      ! HITRAN cutoff model as an intensity threshold
      if (trim(cutoff_model)=="HITRAN") then
         if (tranfreq<=nu_crit) then
           cutoff = S_crit*tranfreq/nu_crit*tanh(0.5_rk*c2/temp*tranfreq)
         else
           cutoff = S_crit
         endif
      else
         cutoff = thresh
      endif
      !
   end function apply_HITRAN_cutoff

   subroutine do_log_intensity(abscoef,Nintens)
     !
     implicit none
     !
     !
     real(rk),intent(in) :: abscoef
     integer(ik),intent(out) :: Nintens(:)
     integer(ik) :: ilog
       !
       ilog=log10(abscoef)
       !
       ilog = max(min(ubound(Nintens,dim=1),ilog),lbound(Nintens,dim=1))
       !
       Nintens(ilog) = Nintens(ilog)+1
       !
   end subroutine do_log_intensity


  subroutine do_gauss0(tranfreq,abscoef0,dfreq,freq,halfwidth,freql,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef0,dfreq,halfwidth,freql
     real(rk),intent(in) :: freq(npoints)
     real(rk),intent(inout) :: intens(npoints)
     real(rk) :: alpha,abscoef,dfreq_,de,ln2,offset_
     integer(ik) :: ib,ie,ipoint
     !
     if (halfwidth<small_) return
     !
     offset_ = offset
     if ( use_width_offset ) offset_ = offset*halfwidth
     !
     call get_ipoint_ranges(tranfreq,freq,offset_,ib,ie)
     !
     ln2=log(2.0_rk)
     alpha = -ln2/halfwidth**2
     !
     abscoef=abscoef0*sqrt(ln2/pi)/halfwidth
     !
     do ipoint=ib,ie
        !
        dfreq_=freq(ipoint)-tranfreq
        !
        de = exp(alpha*dfreq_**2)
        !
        intens(ipoint)=intens(ipoint)+abscoef*de
        !
     enddo
     !
  end subroutine do_gauss0
  !
  !
  ! Obtain left and right grid points around the current frequency
  subroutine get_ipoint_ranges(tranfreq,freq,offset,ib,ie)
     !
     real(rk),intent(in) :: tranfreq,offset
     real(rk),intent(in) :: freq(npoints)
     integer(ik),intent(out) :: ib,ie
     !
     real(rk) :: dfreq_t,freql_t,offset_
     integer(ik) :: igrid
     !
     ib = 1
     ie = Npoints
     !
     if (use_resolving_power) then
         !
         dfreq_t = log(freq(2)) - log(freq(1))
         ib = nint((log(tranfreq - offset) - log(freql))/dfreq_t) + 1
         ib = max(ib, 1)
         !
         ie = nint((log(tranfreq + offset) - log(freql))/dfreq_t) + 1
         ie = min(ie,npoints)
         !
         return
         !
     endif
     !
     if (Ngrids>0) then 
       igrid = 1
       do while(tranfreq>grid(igrid+1)%freql.and.igrid<Ngrids)
         igrid = igrid + 1
       enddo
       !
       dfreq_t = grid(igrid)%dfreq
       freql_t = grid(igrid)%freql
       offset_ = grid(igrid)%offset
       !
       if (tranfreq-offset_>=grid(igrid)%freql) then
         ib =  max(nint( ( tranfreq-offset_-freql_t)/dfreq_t )+1,1)
         ib = grid(igrid)%i1+ib-1
       elseif(igrid>1) then
         ib =  max(nint( ( tranfreq-offset_-grid(igrid-1)%freql)/grid(igrid-1)%dfreq )+1,1)
         ib = grid(igrid-1)%i1+ib-1
       endif
       !
       if (tranfreq+offset_<=grid(igrid)%freqr) then
         ie =  min(nint( ( tranfreq+offset_-freql_t)/dfreq_t )+1,npoints)
         ie = grid(igrid)%i1+ie-1
       elseif(igrid<Ngrids) then
         ie =  min(nint( ( tranfreq+offset_-grid(igrid+1)%freql)/grid(igrid+1)%dfreq )+1,npoints)
         ie = grid(igrid+1)%i1+ie-1
       endif
       !
     else
       !
       freql_t = freq(1)
       dfreq_t = freq(2) - freq(1)
       !
       ib =  max(nint( ( tranfreq-offset-freql_t)/dfreq_t )+1,1)
       ie =  min(nint( ( tranfreq+offset-freql_t)/dfreq_t )+1,npoints)
       !
     endif
     !
  end subroutine get_ipoint_ranges
  !
  ! Obtain the ipoint for the current frequency
  subroutine get_grid_ipoint(tranfreq,freq,ipoint)
     !
     real(rk),intent(in) :: tranfreq
     real(rk),intent(in) :: freq(npoints)
     integer(ik),intent(out) :: ipoint
     !
     real(rk) :: dfreq_t,freql_t
     integer(ik) :: igrid
     !
     if (use_resolving_power) then
         !
         dfreq_t = log(freq(2)) - log(freq(1))
         freql_t = freq(1)
         !
         ipoint =  nint(real((log(tranfreq)-log(freql_t))/dfreq_t,rk))+1
         !
         return
         !
     endif
     !
     if (Ngrids>0) then 
       igrid = 1
       do while(tranfreq>grid(igrid+1)%freql.and.igrid<Ngrids)
         igrid = igrid + 1
       enddo
       !
       dfreq_t = grid(igrid)%dfreq
       freql_t = grid(igrid)%freql
       !
       ipoint =  min(grid(igrid)%i1+max(nint( ( tranfreq-freql_t)/dfreq_t )+1,1),npoints)
       !
     else
       !
       freql_t = freq(1)
       dfreq_t = freq(2) - freq(1)
       !
       ipoint =  min(max(nint( ( tranfreq-freql_t)/dfreq_t )+1,1),npoints)
       !
     endif
     !
  end subroutine get_grid_ipoint
  !
  !
  subroutine do_gauss(tranfreq,abscoef,dfreq,freq,halfwidth,offset,freql,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,dfreq,halfwidth,offset,freql
     real(rk),intent(in) :: freq(npoints)
     real(rk),intent(inout) :: intens(npoints)
     real(rk) :: dfreq_,de,xp,xm,ln2,x0,offset_
     integer(ik) :: ib,ie,ipoint
     !
     if (halfwidth<small_) return 
     !
     ln2=log(2.0_rk)
     !
     x0 = sqrt(ln2)/halfwidth*dfreq*0.5_rk
     !
     offset_ = offset
     if ( use_width_offset ) offset_ = offset*halfwidth
     !
     ib =  max(nint( ( tranfreq-offset_-freql)/dfreq )+1,1)
     ie =  min(nint( ( tranfreq+offset_-freql)/dfreq )+1,npoints)
     !
     !omp parallel do private(ipoint,dfreq_,xp,xm,de) shared(intens) schedule(dynamic)
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
     !omp end parallel do

  end subroutine do_gauss
  !
  !
  subroutine do_lorentz(tranfreq,abscoef,dfreq,freq,halfwidth,offset,freql,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,dfreq,halfwidth,offset,freql
     real(rk),intent(in) :: freq(npoints)
     real(rk),intent(inout) :: intens(npoints)
     real(rk) :: dfreq_,de,xp,xm,b,dnu_half,offset_
     integer(ik) :: ib,ie,ipoint
     !
     offset_ = offset
     if ( use_width_offset ) offset_ = offset*halfwidth
     !
     ib =  max(nint( ( tranfreq-offset_-freql)/dfreq )+1,1)
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
  end subroutine do_lorentz

  subroutine do_pseud(tranfreq,abscoef,dfreq,freq,halfwidth,offset,freql,dpwcoef,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,dfreq,halfwidth,offset,freql,dpwcoef
     real(rk),intent(in) :: freq(npoints)
     real(rk),intent(inout) :: intens(npoints)
     real(rk) :: dfreq_,de,xp,xm,ln2,dnu_half,bnorm,f,eta1,eta2,intens2,halfwidth0,x0,intens1,offset_
     integer(ik) :: ib,ie,ipoint

      !
      ln2=log(2.0_rk)
      halfwidth0=dpwcoef*tranfreq
      x0 = sqrt(ln2)/halfwidth0*dfreq*0.5_rk
      !
      offset_ = offset
      if ( use_width_offset ) offset_ = offset*max(halfwidth0,halfwidth)
      !
      ib =  max(nint( ( tranfreq-offset_-freql)/dfreq )+1,1)
      ie =  min(nint( ( tranfreq+offset_-freql)/dfreq )+1,npoints)
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
      f = ( halfwidth0**5 + 2.69269_rk*halfwidth0**4*halfwidth+2.42843_rk*halfwidth0**3*halfwidth**2+&
          4.47163_rk*halfwidth0**2*halfwidth**3+0.07842_rk*halfwidth0*halfwidth**4+halfwidth**5)**0.2_rk
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


  end subroutine do_pseud

  subroutine do_pse_R(tranfreq,abscoef,dfreq,freq,halfwidth,offset,freql,dpwcoef,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,dfreq,halfwidth,offset,freql,dpwcoef
     real(rk),intent(in) :: freq(npoints)
     real(rk),intent(inout) :: intens(npoints)
     real(rk) :: dfreq_,de,xp,xm,ln2,dnu_half,bnorm,eta1,eta2,intens2,halfwidth0,x0,a,wg,va0
     real(rk) :: gammaV,intens1,offset_
     integer(ik) :: ib,ie,ipoint
     !
     ln2=log(2.0_rk)
     halfwidth0=dpwcoef*tranfreq
     x0 = sqrt(ln2)/halfwidth0*dfreq*0.5_rk
     !
     offset_ = offset
     if ( use_width_offset ) offset_ = offset*max(halfwidth0,halfwidth)
     !
     ib =  max(nint( ( tranfreq-offset_-freql)/dfreq )+1,1)
     ie =  min(nint( ( tranfreq+offset_-freql)/dfreq )+1,npoints)
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
  end subroutine do_pse_R
      !
  subroutine  do_pse_L(tranfreq,abscoef,dfreq,freq,halfwidth,offset,freql,dpwcoef,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,dfreq,halfwidth,offset,freql,dpwcoef
     real(rk),intent(in) :: freq(npoints)
     real(rk),intent(inout) :: intens(npoints)
     real(rk) :: dfreq_,de,xp,xm,ln2,bnorm,eta1,eta2,intens2,halfwidth0
     real(rk) :: d,x0,intens1,dnu_half,offset_
     integer(ik) :: ib,ie,ipoint
     !
     ln2=log(2.0_rk)
     halfwidth0=dpwcoef*tranfreq
     x0 = sqrt(ln2)/halfwidth0*dfreq*0.5_rk
     !
     offset_ = offset
     if ( use_width_offset ) offset_ = offset*max(halfwidth0,halfwidth)
     !
     ib =  max(nint( ( tranfreq-offset_-freql)/dfreq )+1,1)
     ie =  min(nint( ( tranfreq+offset_-freql)/dfreq )+1,npoints)
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
  end subroutine  do_pse_L


  subroutine do_Voi_Q(tranfreq,abscoef,dfreq,freq,abciss,weight,halfwidth,offset,freql,dpwcoef,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,dfreq,halfwidth,offset,freql,dpwcoef
     real(rk),intent(in) :: freq(npoints),abciss(nquad),weight(nquad)
     real(rk),intent(inout) :: intens(npoints)
     real(rk) :: alpha,ln2,halfwidth0,gamma,xp
     real(rk) :: ln22,y,dx2,x0,sigma,xi,x1,x2,l1,l2,bnormq(1:nquad),Lorentz(nquad),dxp,voigt_,offset_
     integer(ik) :: ib,ie,ipoint,iquad
     !
     halfwidth0=dpwcoef*tranfreq
     if (halfwidth0<100.0_rk*small_) return
     !
     offset_ = offset
     if ( use_width_offset ) offset_ = offset*max(halfwidth0,halfwidth)
     !
     ln2=log(2.0_rk)
     ln22 = ln2*2.0_rk
     x0 = sqrt(ln2)/halfwidth0*dfreq*0.5_rk
     !
     ib =  max(nint( ( tranfreq-offset_-freql)/dfreq )+1,1)
     ie =  min(nint( ( tranfreq+offset_-freql)/dfreq )+1,npoints)
     !
     if (ie<=ib) return
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
  end subroutine  do_Voi_Q
  !
  subroutine  do_Voigt(tranfreq,abscoef,use_resolving_power,freq,halfwidth_Lorentz,dpwcoef,offset,freql,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,halfwidth_Lorentz,offset,freql,dpwcoef
     real(rk),intent(in) :: freq(:)
     real(rk),intent(out) :: intens(:)
     real(rk) :: tranfreq_i,halfwidth_doppler,offset_
     integer(ik) :: ib,ie,ipoint
     logical, intent(in) :: use_resolving_power
      !
      halfwidth_doppler=dpwcoef*tranfreq
      !
      offset_ = offset
      if ( use_width_offset ) offset_ = offset*halfwidth_Lorentz
      !
      call get_ipoint_ranges(tranfreq,freq,offset_,ib,ie)
      !
      if (halfwidth_doppler<100.0_rk*small_) return
      !
      !if (use_resolving_power) then
      !   d_ln_freq = log(freq(2)) - log(freq(1))
      !   ib = nint((log(tranfreq - offset) - log(freql))/d_ln_freq) + 1
      !   ib = max(ib, 1)
      !   
      !   ie = nint((log(tranfreq + offset) - log(freql))/d_ln_freq) + 1
      !   ie = min(ie, npoints)
      !else
      !   d_freq = freq(2) - freq(1)
      !   ib =  max(nint( ( tranfreq-offset-freql)/d_freq )+1,1)
      !   ie =  min(nint( ( tranfreq+offset-freql)/d_freq )+1,npoints)
      !   !
      !endif
      !
      if (ie<=ib) return
      !
      !$omp parallel do private(ipoint,tranfreq_i) shared(intens) schedule(dynamic)
      do ipoint = ib, ie
         !
         tranfreq_i = freq(ipoint)
         !
         !intens(ipoint)=intens(ipoint)+voigt_faddeeva(tranfreq_i,tranfreq,halfwidth_doppler,halfwidth_Lorentz)*abscoef
         intens(ipoint)=intens(ipoint)+voigt_humlicek(tranfreq_i,tranfreq,halfwidth_doppler,halfwidth_Lorentz)*abscoef
         !
      enddo
      !$omp end parallel do
      !
  end subroutine  do_Voigt
  !

  subroutine  do_Voigt_Fast(tranfreq,abscoef,use_resolving_power,freq,index_Lorentz,dpwcoef,offset,freql,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,offset,freql,dpwcoef
     real(rk),intent(in) :: freq(:)
     integer,intent(in)  :: index_Lorentz
     real(rk),intent(out) :: intens(:)
     real(rk) :: halfwidth_doppler,d_freq,d_ln_freq
     integer(ik) :: ib,ie
     logical :: use_resolving_power
      !
      halfwidth_doppler=dpwcoef*tranfreq
      if (halfwidth_doppler<100.0_rk*small_) return
      !
      if (use_resolving_power) then
         d_ln_freq = log(freq(2)) - log(freq(1))
         ib = nint((log(tranfreq - offset) - log(freql))/d_ln_freq) + 1
         ib = max(ib, 1)
         
         ie = nint((log(tranfreq + offset) - log(freql))/d_ln_freq) + 1
         ie = min(ie, npoints)
      else
         d_freq = freq(2) - freq(1)
         ib =  max(nint( ( tranfreq-offset-freql)/d_freq )+1,1)
         ie =  min(nint( ( tranfreq+offset-freql)/d_freq )+1,npoints)
      endif
      !
      call fast_voigt%compute(freq,intens, abscoef,ib,ie,freql,tranfreq,index_Lorentz)
      !
  end subroutine  do_Voigt_Fast
  !
  !
  subroutine  do_Voigt_916(tranfreq,abscoef,dfreq,freq,halfwidth_Lorentz,dpwcoef,offset,freql,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,dfreq,halfwidth_Lorentz,offset,freql,dpwcoef
     real(rk),intent(in) :: freq(:)
     real(rk),intent(out) :: intens(:)
     real(rk) :: tranfreq_i,halfwidth_doppler
     integer(ik) :: ib,ie,ipoint
      !
      stop 'do_Voigt_916 is not working yet'
      !
      halfwidth_doppler=dpwcoef*tranfreq
      if (halfwidth_doppler<100.0_rk*small_) return
      ib =  max(nint( ( tranfreq-offset-freql)/dfreq )+1,1)
      ie =  min(nint( ( tranfreq+offset-freql)/dfreq )+1,npoints)
      !
      !$omp parallel do private(ipoint,tranfreq_i) shared(intens) schedule(dynamic)
      do ipoint=ib,ie
         !
         tranfreq_i = freq(ipoint)
         !
         intens(ipoint)=intens(ipoint)+voigt_humlicek(tranfreq_i,tranfreq,halfwidth_doppler,halfwidth_Lorentz)*abscoef
         !
      enddo
      !$omp end parallel do
      !
  end subroutine  do_Voigt_916
  !

  subroutine  do_rect(tranfreq,abscoef,dfreq,freq,offset,freql,intens)
     !
     implicit none
     !
     real(rk),intent(in) :: tranfreq,abscoef,dfreq,offset,freql
     real(rk),intent(in) :: freq(:)
     real(rk),intent(out) :: intens(:)
     real(rk) :: tranfreq_i
     integer(ik) :: ib,ie,ipoint
      !
      ib =  max(nint( ( tranfreq-offset-freql)/dfreq )+1,1)
      ie =  min(nint( ( tranfreq+offset-freql)/dfreq )+1,npoints)
      !
      !$omp parallel do private(ipoint,tranfreq_i) shared(intens) schedule(dynamic)
      do ipoint=ib,ie
         !
         tranfreq_i = freq(ipoint)
         !
         intens(ipoint)=intens(ipoint)+abscoef/dfreq
         !
      enddo
      !$omp end parallel do

  end subroutine  do_rect


  !
  subroutine do_stick(sunit,nswap,nlines,maxitems,energies,Jrot,gtot,ilevelf_ram,ileveli_ram,&
                      abscoef_ram,acoef_ram,nchars_quanta,quantum_numbers)

     !
     implicit none
     !
     integer(ik),intent(in) :: sunit,nswap,maxitems,nlines,ilevelf_ram(nswap),ileveli_ram(nswap)
     real(rk),intent(in) :: abscoef_ram(nswap),acoef_ram(nswap),jrot(nlines),energies(nlines)
     real(rk)  :: Jf,Ji,abscoef,acoef,tranfreq,energyf,energyi
     integer(ik),intent(in) :: nchars_quanta(maxitems),gtot(nlines)
     character(len=20),intent(in) :: quantum_numbers(0:maxitems,nlines)
     integer(ik) :: nchars_,kitem,ierror_nu,iE,ierror_i,ierror,ierror_S,l,qni,qnf,ierror_f,&
                    ilevelf,ileveli,iswap,nchars_tot
     character(9) :: b_fmt
     character(len=cl) :: h_fmt
     integer(ik) :: Jpp,Jp,i
     real(rk)    :: gamma1,gamma2,gamma3,n1,n2,n3
     !
     !
     do iswap = 1,nswap
       !
       ilevelf = ilevelf_ram(iswap)
       ileveli = ileveli_ram(iswap)
       abscoef = abscoef_ram(iswap)
       acoef   = acoef_ram(iswap)
       energyf = energies(ilevelf)
       energyi = energies(ileveli)
       tranfreq = energyf - energyi
       !
       jf = jrot(ilevelf)
       ji = jrot(ileveli)
       !
       !
       if (abscoef>1.0D-99) then 
          write(h_fmt,"(a)") "(i2,i1,f12.6,es10.3,e10.3,f5.4,f5.3,f10.4,f4.2,f8.6)"
       else
          write(h_fmt,"(a)") "(i2,i1,f12.6,es10.2E3,e10.3,f5.4,f5.3,f10.4,f4.2,f8.6)"
       endif
       !
       if (stick_hitran) then
         !
         Jp  = int(Jf)
         Jpp = int(Ji)
         !
         gamma1 = species(1)%gammaQN(Jpp,Jp-Jpp)
         gamma2 = species(2)%gammaQN(Jpp,Jp-Jpp)
         n1 = species(1)%nQN(Jpp,Jp-Jpp)
         !
         write(sunit,h_fmt,advance="no") &
                     imolecule,iso,tranfreq,abscoef,acoef,&
                     gamma1,gamma2,&
                     energyi,n1,species(1)%delta
         !
         nchars_ = 0
         nchars_tot = 0
         !
         do kitem = 3,maxitems
           !
           l = len(trim(quantum_numbers(kitem,ilevelf)))
           !
           write(b_fmt,"('(a',i1,')')") nchars_quanta(kitem)
           !
           nchars_  = nchars_ + nchars_quanta(kitem)
           !
           if (nchars_>21) cycle
           nchars_tot = nchars_tot + nchars_quanta(kitem)
           !
           write(sunit,b_fmt,advance="no") trim(quantum_numbers(kitem,ilevelf))
           !
         enddo
         !
         nchars_ = 0
         !
         do kitem = 3,maxitems
           !
           !l = len(trim(quantum_numbers(kitem,ileveli)))
           !
           !b_fmt = "(1x,a3)" ; if (l>3) b_fmt = "(1x,a8)"
           !
           write(b_fmt,"('(a',i1,')')") nchars_quanta(kitem)
           !
           nchars_  = nchars_ + nchars_quanta(kitem)
           !
           if (nchars_>21) cycle
           nchars_tot = nchars_tot + nchars_quanta(kitem)
           !
           write(sunit,b_fmt,advance="no") trim(quantum_numbers(kitem,ileveli))
           !
         enddo
         !
         ! topup to 44
         do kitem = nchars_tot+1,42
           !
           write(sunit,'(1x)',advance="no")
           !
         enddo
         !
         if ( mod(nint(2*Jf),2)==0 ) then
           !
           write(sunit,'(i4,a2,a2)',advance="no") nint(Jf),trim(quantum_numbers(1,ilevelf)),trim(quantum_numbers(2,ilevelf))
           write(sunit,'(1x,i4,a2,a2)',advance="no") nint(Ji),trim(quantum_numbers(1,ileveli)),&
                                                           trim(quantum_numbers(2,ileveli))
           !
         else
           !
           write(sunit,'(f5.1,a1,a2)',advance="no") Jf,trim(quantum_numbers(1,ilevelf)),trim(quantum_numbers(2,ilevelf))
           write(sunit,'(f5.1,a1,a2)',advance="no") Ji,trim(quantum_numbers(1,ileveli)),trim(quantum_numbers(2,ileveli))
           !
         endif
         !
         ierror_nu = 0
         !
         do iE = 1,HITRAN_E(1)%N
            !
            read(quantum_numbers(HITRAN_E(iE)%iqn,ileveli),*) qni
            read(quantum_numbers(HITRAN_E(iE)%iqn,ilevelf),*) qnf
            !
            ierror_i = 0
            do ierror = 0,6
              if (qni<=HITRAN_E(iE)%error_vmax(ierror)) ierror_i = ierror
              if (qnf<=HITRAN_E(iE)%error_vmax(ierror)) ierror_f = ierror
            enddo
            !
            ierror_nu = max(ierror_i,ierror_f)
            !
         enddo
         !
         ierror_S = 0
         !
         do iE = 1,HITRAN_S(1)%N
            !
            read(quantum_numbers(HITRAN_S(iE)%iqn,ileveli),*) qni
            read(quantum_numbers(HITRAN_S(iE)%iqn,ilevelf),*) qnf
            !
            ierror_i = 0
            do ierror = 0,6
              if (qni<=HITRAN_S(iE)%error_vmax(ierror)) ierror_i = ierror
              if (qnf<=HITRAN_S(iE)%error_vmax(ierror)) ierror_f = ierror
            enddo
            !
            ierror_S = max(ierror_i,ierror_f)
            !
         enddo
         !
         write(sunit,'(1x,i1)',advance="no") ierror_nu
         write(sunit,'(i1)',advance="no") ierror_S
         write(sunit,'(i1)',advance="no") HITRAN_Air%ierr
         write(sunit,'(i1)',advance="no") HITRAN_Self%ierr
         write(sunit,'(i1)',advance="no") HITRAN_n%ierr
         write(sunit,'(i1)',advance="no") HITRAN_Delta%ierr
         !
         write(sunit,"(' 0 0 0 0 0 0 ',f7.1,f7.1)",advance="yes") real(gtot(ilevelf)),real(gtot(ileveli))
         !
       elseif (stick_oxford) then
         !
         Jp  = int(Jf)
         Jpp = int(Ji)
         !
         gamma1 = species(1)%gammaQN(Jpp,Jp-Jpp)
         gamma2 = species(2)%gammaQN(Jpp,Jp-Jpp)
         gamma3 = species(3)%gammaQN(Jpp,Jp-Jpp)
         n1 = species(1)%nQN(Jpp,Jp-Jpp)
         n2 = species(2)%nQN(Jpp,Jp-Jpp)
         n3 = species(3)%nQN(Jpp,Jp-Jpp)
         !
         write(sunit,"(i3,f12.6,e10.3,3(f5.4,f4.2),f10.4)",advance="no") &
                     iso,tranfreq,abscoef,&
                     gamma1,n1,gamma2,n2,gamma3,n3,&
                     energyi
         !
         nchars_tot = 1
         nchars_ = 1
         !
         do i = QN%Vibcol(1),QN%Vibcol(2)
           !
           kitem=i-4
           !
           !l = len(trim(quantum_numbers(kitem,ilevelf)))
           !
           write(b_fmt,"('(a',i1,')')") nchars_quanta(kitem)
           !
           nchars_  = nchars_ + nchars_quanta(kitem)
           !
           if (nchars_>22) cycle
           nchars_tot = nchars_tot + nchars_quanta(kitem)
           !
           write(sunit,b_fmt,advance="no") trim(quantum_numbers(kitem,ilevelf))
           !
         enddo
         !
         !write(sunit,'(1x)',advance="no")
         !
         nchars_tot = nchars_tot + 1
         !
         nchars_ = 1
         do i = QN%Vibcol(1),QN%Vibcol(2)
           !
           kitem=i-4
           !
           write(b_fmt,"('(a',i1,')')") nchars_quanta(kitem)
           !
           nchars_  = nchars_ + nchars_quanta(kitem)
           !
           if (nchars_>22) cycle
           nchars_tot = nchars_tot + nchars_quanta(kitem)
           !
           write(sunit,b_fmt,advance="no") trim(quantum_numbers(kitem,ileveli))
           !
         enddo
         !
         ! topup to 42
         do kitem = nchars_tot+1,44
           !
           write(sunit,'(1x)',advance="no")
           !
         enddo
         !
         write(sunit,'(1x)',advance="no")
         !
         if ( mod(nint(2*Jf),2)==0 ) then
           !
           write(sunit,'(i4,1x)',advance="no") nint(Jf)
           nchars_ = 5
           do i = 1,2
             kitem=QN%Rotcol(i)-4
             !
             write(b_fmt,"('(a',i1,')')") nchars_quanta(kitem)
             nchars_  = nchars_ + nchars_quanta(kitem)
             !
             if (nchars_>10) cycle
             nchars_tot = nchars_tot + nchars_quanta(kitem)
             !
             write(sunit,b_fmt,advance="no") trim(quantum_numbers(kitem,ilevelf))
             !
           enddo
           !
           write(sunit,'(i4,1x)',advance="no") nint(Ji)
           nchars_ = 5
           do i = 1,2
             kitem=QN%Rotcol(i)-4
             !
             write(b_fmt,"('(a',i1,')')") nchars_quanta(kitem)
             nchars_  = nchars_ + nchars_quanta(kitem)
             !
             if (nchars_>10) cycle
             nchars_tot = nchars_tot + nchars_quanta(kitem)
             !
             write(sunit,b_fmt,advance="no") trim(quantum_numbers(kitem,ileveli))
             !
           enddo
             !
         else
           !
           write(sunit,'(f5.1,1x,a1,a1)',advance="no") Jf,trim(quantum_numbers(1,ilevelf)),trim(quantum_numbers(2,ilevelf))
           write(sunit,'(f5.1,1x,a1,a1)',advance="no") Ji,trim(quantum_numbers(1,ileveli)),trim(quantum_numbers(2,ileveli))
           !
         endif
         !
         ierror_nu = 0
         !
         do iE = 1,HITRAN_E(1)%N
            !
            read(quantum_numbers(HITRAN_E(iE)%iqn,ileveli),*) qni
            read(quantum_numbers(HITRAN_E(iE)%iqn,ilevelf),*) qnf
            !
            ierror_i = 0
            do ierror = 0,6
              if (qni<=HITRAN_E(iE)%error_vmax(ierror)) ierror_i = ierror
              if (qnf<=HITRAN_E(iE)%error_vmax(ierror)) ierror_f = ierror
            enddo
            !
            ierror_nu = max(ierror_i,ierror_f)
            !
         enddo
         !
         ierror_S = 0
         !
         do iE = 1,HITRAN_S(1)%N
            !
            read(quantum_numbers(HITRAN_S(iE)%iqn,ileveli),*) qni
            read(quantum_numbers(HITRAN_S(iE)%iqn,ilevelf),*) qnf
            !
            ierror_i = 0
            do ierror = 0,6
              if (qni<=HITRAN_S(iE)%error_vmax(ierror)) ierror_i = ierror
              if (qnf<=HITRAN_S(iE)%error_vmax(ierror)) ierror_f = ierror
            enddo
            !
            ierror_S = max(ierror_i,ierror_f)
            !
         enddo
         !
         write(sunit,'(1x,i1)',advance="no") ierror_nu
         write(sunit,'(i1)',advance="no") ierror_S
         write(sunit,'(i1)',advance="no") HITRAN_Air%ierr
         write(sunit,'(i1)',advance="no") HITRAN_Self%ierr
         write(sunit,'(i1)',advance="no") HITRAN_n%ierr
         write(sunit,'(i1)',advance="yes") HITRAN_Delta%ierr
         !
       else
          !
          ! Not HITRAN: the main detailed stick-output is written to out and 
          ! also a simple two-column-version into .stick 
          !
          ! write to .stick-file
          write(sunit,my_fmt) tranfreq,abscoef
          !
          write(out,my_fmt,advance="no") &
          tranfreq,abscoef,jrot(ilevelf),energyf, jrot(ileveli),energyi
          !
          ! printing the line width
          !
          if (lineprofile_do) then
            !
            write(out,"(f12.5)",advance="no") halfwidth
            !
          endif
          !
          !write(out,"(a4)",advance="no"), " <- "
          !
          do kitem = 1,maxitems
            !
            l = len(trim(quantum_numbers(kitem,ilevelf)))
            !
            b_fmt = "(1x,a3)" ; if (l>3) b_fmt = "(1x,a8)"
            !
            if (nchars_quanta(kitem)<10) then
              write(b_fmt,"('(1x,a',i1,')')") nchars_quanta(kitem)
            else
              write(b_fmt,"('(1x,a',i2,')')") nchars_quanta(kitem)
            endif
            !
            write(out,b_fmt,advance="no") trim(quantum_numbers(kitem,ilevelf))
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
            if (nchars_quanta(kitem)<10) then
              write(b_fmt,"('(1x,a',i1,')')") nchars_quanta(kitem)
            else
              write(b_fmt,"('(1x,a',i2,')')") nchars_quanta(kitem)
            endif
            !
            !write(out,"(2i7,e17.4)",advance="no") ilevelf,ileveli,acoef
            !
            write(out,b_fmt,advance="no") trim(quantum_numbers(kitem,ileveli))
            !
          enddo
          !
          write(out,"(a1)",advance="yes") " "
          !
       endif
       !
    enddo

  end subroutine do_stick



  subroutine voi_quad(npoints,ib,ie,freq,abscoef,intens,dfreq,tranfreq,alpha,gamma,nquad,abciss,weight)
     !
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

  function get_Voigt_gamma_n(Nspecies,Jf,Ji,idx) result(f)
     !
     implicit none
     !
     integer(ik),intent(in) :: Nspecies
     real(rk),intent(in) :: Jf,Ji
     real(rk) :: halfwidth,gamma_,n_,f
     integer(ik) :: ispecies,Jpp,Jp
     integer(ik),optional,intent(inout) :: idx
     Jpp = nint(Ji)
     Jp  = nint(Jf)
     !
     if (species(1)%if_defined) then
       !
       if(present(idx)) then
          !
          idx = gamma_idx(Jpp,Jp-Jpp)
          f   = gamma_comb(Jpp,Jp-Jpp)
          return
          !
       endif
       !     
       f = gamma_comb(Jpp,Jp-Jpp)
       return
       !
     else
       !
       halfwidth = 0
       !
       do ispecies =1,Nspecies
           !
           gamma_ = species(ispecies)%gamma
           n_     = species(ispecies)%n
           !
           halfwidth =  halfwidth + species(ispecies)%ratio*gamma_*(species(ispecies)%T0/Temp)**N_*pressure/species(ispecies)%P0
           !
       enddo
       !
       f = halfwidth
       !
     endif
     !
  end function get_Voigt_gamma_n
 
  function get_Voigt_gamma_val(Nspecies,Jf,Ji) result(f)
     !
     implicit none
     !
     integer(ik),intent(in) :: Nspecies
     real(rk),intent(in) :: Jf,Ji
     real(rk) :: halfwidth,gamma_,n_,f
     integer(ik) :: ispecies,Jpp,Jp
     !
     halfwidth = 0
     Jpp = nint(Ji)
     Jp  = nint(Jf)
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
     f = halfwidth
     !
  end function get_Voigt_gamma_val
  !
  function voigt_faddeeva(nu,nu0,alpha,gamma) result(f)
   !
   implicit none
   !
   real(rk),intent(in) :: nu,nu0,alpha,gamma
   real(rk) :: f,sigma,x,y
   integer(ik),parameter :: n=10
   real(rk) :: r(n),weight(n),ln2

   real(rk),parameter ::  RT2LN2 = 1.1774100225154747_rk     ! sqrt(2ln2)
   real(rk),parameter ::  RTPI   = 1.7724538509055159_rk     ! sqrt(pi)
   real(rk),parameter ::  RT2PI  = 2.5066282746310002_rk     ! sqrt(2pi)
   real(rk),parameter ::  RT2    = 1.4142135623730951_rk     ! sqrt(2)
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
  !
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
    integer(ik) :: n,MAXIT
    real(rk) :: xr(n),wr(n)
    !
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
