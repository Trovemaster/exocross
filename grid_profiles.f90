! User-supplied, temperature- and band-symmetry-dependent line profiles.
module grid_profiles
  use accuracy, only: rk, ik, wl, out
  use input, only: read_line, reada, readu, reread, nitems
  use timer, only: ArrayStart, ArrayStop
  use symmetry, only: sym, irrep_index, gamma_lookup
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: grid_profiles_do, read_grid_profiles, configure_grid_profiles, require_grid_pair_profiles
  public :: load_grid_profiles, free_grid_profiles, grid_label_index, grid_profile_extent
  public :: do_grid_sampling, grid_profile_count

  integer, parameter :: label_length = 20 ! Matches spectrum%quantum_numbers.
  type :: profile_file_t
    real(rk) :: temperature
    character(len=label_length) :: label
    character(len=wl) :: filename
  end type
  type(profile_file_t), allocatable :: files(:)
  character(len=label_length), allocatable :: labels(:)
  integer, allocatable :: file_index(:,:)
  real(rk), allocatable :: nu_profile(:,:,:), f_profile(:,:,:)
  real(rk), allocatable :: profile_inverse_step(:,:)
  real(rk), allocatable :: profile_weights(:,:,:)
  integer(ik), allocatable :: missing_band(:,:)
  integer :: n_profile_grid = 0, n_profile_t = 0, n_profile_sym = 0
  logical :: grid_profiles_do = .false.
  real(rk), parameter :: normalization_tolerance = 1.0e-3_rk

contains

  subroutine grid_error(message)
    character(len=*), intent(in) :: message
    write(out,'(a)') 'GRID: '//trim(message)
    error stop 1
  end subroutine

  logical function same_temperature(a,b)
    real(rk), intent(in) :: a,b
    same_temperature = abs(a-b) <= 1.0e-10_rk*max(1.0_rk,abs(a),abs(b))
  end function

  subroutine read_grid_profiles()
    type(profile_file_t) :: entry
    character(len=wl) :: word
    logical :: eof
    integer :: ios

    if (allocated(files)) call grid_error('Only one PROFILES block is allowed')
    allocate(files(0))
    do
      call read_line(eof)
      if (eof) call grid_error('PROFILES requires END')
      if (nitems == 0) cycle
      call readu(word)
      if (trim(word) == 'END') then
        if (nitems /= 1) call grid_error('Unexpected fields after END in PROFILES')
        exit
      endif
      if (nitems /= 3) call grid_error('Expected: temperature label filename in PROFILES')
      read(word,*,iostat=ios) entry%temperature
      if (ios /= 0) call grid_error('Invalid profile temperature: '//trim(word))
      if (.not.ieee_is_finite(entry%temperature)) call grid_error('Profile temperature must be finite')
      if (entry%temperature <= 0) call grid_error('Profile temperature must be positive')
      call reada(word) ! Canonicalize band labels later; preserve filenames as supplied.
      if (len_trim(word) == 0.or.len_trim(word) > label_length) &
        call grid_error('Profile label must contain 1 to 20 characters')
      entry%label = word
      call reada(entry%filename)
      if (len_trim(entry%filename) == 0) call grid_error('Empty profile filename')
      files = [files,entry]
    enddo
    if (size(files) == 0) call grid_error('PROFILES is empty')
    grid_profiles_do = .true.
  end subroutine

  ! Validate the full T x label product after all input blocks have been read.
  ! No temperature interpolation or extrapolation is performed.
  subroutine configure_grid_profiles(temperatures)
    real(rk), intent(in) :: temperatures(:)
    integer :: i,j,it,is,info,irrep
    character(len=wl) :: message

    if (.not.allocated(files)) call grid_error('GRID requires a PROFILES block')
    if (sym%Nirreps == 0) call grid_error('Define SYMMETRY C2v, C3v or Cs for band profiles')
    do i = 1,size(files)
      irrep = irrep_index(files(i)%label)
      if (irrep == 0) call grid_error('Unknown band irrep '//trim(files(i)%label)//' for '//trim(sym%group))
      files(i)%label = sym%label(irrep) ! Canonicalize aliases before duplicate checks.
    enddo
    n_profile_t = size(temperatures)
    if (n_profile_t < 1) call grid_error('Empty temperature list')
    do it = 1,n_profile_t
      if (.not.ieee_is_finite(temperatures(it))) call grid_error('Temperature must be finite')
      if (temperatures(it) <= 0) call grid_error('Temperature must be positive')
      do j = 1,it-1
        if (same_temperature(temperatures(it),temperatures(j))) &
          call grid_error('Duplicate requested temperature')
      enddo
    enddo
    allocate(labels(size(files)))
    labels = ''
    n_profile_sym = 0
    do i = 1,size(files)
      if (grid_label_index(files(i)%label) /= 0) cycle
      n_profile_sym = n_profile_sym+1
      labels(n_profile_sym) = files(i)%label
    enddo
    allocate(file_index(n_profile_t,n_profile_sym),stat=info)
    call ArrayStart('GRID:file_index',info,size(file_index),kind(file_index))
    file_index = 0
    do i = 1,size(files)
      is = grid_label_index(files(i)%label)
      it = 0
      do j = 1,n_profile_t
        if (same_temperature(files(i)%temperature,temperatures(j))) it = j
      enddo
      if (it == 0) then
        write(message,'(a,f15.6,a,a)') 'Unrequested temperature ',files(i)%temperature,' for ',trim(files(i)%label)
        call grid_error(message)
      endif
      if (file_index(it,is) /= 0) then
        write(message,'(a,f15.6,1x,a)') 'Duplicate (temperature,label): ',temperatures(it),trim(labels(is))
        call grid_error(message)
      endif
      file_index(it,is) = i
    enddo
    do is = 1,n_profile_sym
      do it = 1,n_profile_t
        if (file_index(it,is) /= 0) cycle
        write(message,'(a,f15.6,1x,a)') 'Missing (temperature,label): ',temperatures(it),trim(labels(is))
        call grid_error(message)
      enddo
    enddo
    write(out,'(a,i0,a,i0)') 'GRID: temperatures = ',n_profile_t,', labels = ',n_profile_sym
    call configure_band_weights()
  end subroutine

  subroutine configure_band_weights()
    integer :: up,low,gamma,is,info
    integer(ik) :: multiplicity(sym%Nirreps)
    allocate(profile_weights(n_profile_sym,sym%Nirreps,sym%Nirreps),stat=info)
    call ArrayStart('GRID:profile_weights',info,size(profile_weights),kind(profile_weights))
    allocate(missing_band(sym%Nirreps,sym%Nirreps),stat=info)
    call ArrayStart('GRID:missing_band',info,size(missing_band),kind(missing_band))
    profile_weights = 0
    missing_band = 0
    do low = 1,sym%Nirreps
      do up = 1,sym%Nirreps
        multiplicity = gamma_lookup(up,low)
        do gamma = 1,sym%Nirreps
          if (multiplicity(gamma) == 0) cycle
          is = grid_label_index(sym%label(gamma))
          if (is == 0) then
            missing_band(up,low) = gamma
          else
            ! Provisional equal average over irrep components, NOT over dimensions.
            ! In C3v, E x E therefore has weights 1/3,1/3,1/3 for A1,A2,E.
            profile_weights(is,up,low) = real(multiplicity(gamma),rk)/real(sum(multiplicity),rk)
          endif
        enddo
      enddo
    enddo
    if (any(sum(sym%product,dim=1) > 1)) &
      write(out,'(a)') 'GRID: reducible products use a provisional equal-component profile average'
  end subroutine

  ! Check only products used by retained transitions. Never renormalize a partial
  ! set of components when a required band profile is missing.
  subroutine require_grid_pair_profiles(up,low)
    integer(ik), intent(in) :: up,low
    integer :: gamma
    if (min(up,low) < 1.or.max(up,low) > sym%Nirreps) call grid_error('Unknown state irrep')
    gamma = missing_band(up,low)
    if (gamma == 0) return
    call grid_error('Missing band profile '//trim(sym%label(gamma))//' required by '// &
      trim(sym%label(up))//' x '//trim(sym%label(low)))
  end subroutine

  integer function grid_label_index(label) result(is)
    character(len=*), intent(in) :: label
    integer :: i
    is = 0
    do i = 1,n_profile_sym
      if (label /= labels(i)) cycle
      is = i
      return
    enddo
  end function

  integer function grid_profile_count() result(n)
    n = n_profile_sym
  end function

  ! Uses the existing input tokeniser (including parenthesised comments).
  subroutine read_profile_point(unit,filename,x,y,eof)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename
    real(rk), intent(out) :: x,y
    logical, intent(out) :: eof
    character(len=wl) :: word
    integer :: ios
    do
      call read_line(eof,unit)
      if (eof) return
      if (nitems == 0) cycle
      if (nitems /= 2) call grid_error('Expected two columns in '//trim(filename))
      call reada(word)
      read(word,*,iostat=ios) x
      if (ios /= 0) call grid_error('Invalid wavenumber in '//trim(filename))
      call reada(word)
      read(word,*,iostat=ios) y
      if (ios /= 0) call grid_error('Invalid profile value in '//trim(filename))
      if (.not.ieee_is_finite(x).or..not.ieee_is_finite(y)) &
        call grid_error('Non-finite data in '//trim(filename))
      if (y < 0) call grid_error('Negative profile value in '//trim(filename))
      return
    enddo
  end subroutine

  subroutine load_grid_profiles()
    integer :: it,is,i,j,n,unit,ios,info
    real(rk) :: x,y,area
    logical :: eof
    character(len=wl) :: filename
    do is = 1,n_profile_sym
      do it = 1,n_profile_t
        j = file_index(it,is)
        filename = files(j)%filename
        open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
        if (ios /= 0) call grid_error('Cannot open '//trim(filename))
        n = 0
        do
          call read_profile_point(unit,filename,x,y,eof)
          if (eof) exit
          n = n+1
        enddo
        if (n < 2) call grid_error('At least two profile points required in '//trim(filename))
        if (.not.allocated(nu_profile)) then
          n_profile_grid = n
          allocate(nu_profile(n,n_profile_t,n_profile_sym),stat=info)
          call ArrayStart('GRID:nu_profile',info,size(nu_profile),kind(nu_profile))
          allocate(f_profile(n,n_profile_t,n_profile_sym),stat=info)
          call ArrayStart('GRID:f_profile',info,size(f_profile),kind(f_profile))
          allocate(profile_inverse_step(n_profile_t,n_profile_sym),stat=info)
          call ArrayStart('GRID:profile_inverse_step',info,size(profile_inverse_step),kind(profile_inverse_step))
        endif
        if (n /= n_profile_grid) call grid_error('Profile point counts differ: '//trim(filename))
        rewind(unit)
        do i = 1,n
          call read_profile_point(unit,filename,x,y,eof)
          if (eof) call grid_error('Unexpected EOF in '//trim(filename))
          if (i > 1) then
            if (x <= nu_profile(i-1,it,is)) call grid_error('Grid must strictly increase: '//trim(filename))
          endif
          nu_profile(i,it,is) = x
          f_profile(i,it,is) = y
        enddo
        close(unit)
        if (nu_profile(1,it,is) >= 0.or.nu_profile(n,it,is) <= 0) &
          call grid_error('Offset grid must straddle zero: '//trim(filename))
        profile_inverse_step(it,is) = uniform_inverse_step(nu_profile(:,it,is))
        area = sum(0.5_rk*(f_profile(1:n-1,it,is)+f_profile(2:n,it,is))* &
                         (nu_profile(2:n,it,is)-nu_profile(1:n-1,it,is)))
        write(out,'(a,f15.6,1x,a,1x,a,1x,es20.12)') 'GRID: ',files(j)%temperature,trim(labels(is)), &
          trim(filename)//' integral =',area
        if (.not.ieee_is_finite(area)) call grid_error('Non-finite integral in '//trim(filename))
        if (area <= 0.or.abs(area-1.0_rk) > normalization_tolerance) &
          call grid_error('Profile integral must be within 1e-3 of unity: '//trim(filename))
        ! Remove small finite-support/roundoff errors, after reporting the original area.
        f_profile(:,it,is) = f_profile(:,it,is)/area
      enddo
    enddo
  end subroutine

  real(rk) function grid_profile_extent() result(offset)
    offset = maxval(abs(nu_profile))
  end function

  ! Identify uniform input grids once at load time, allowing only roundoff.
  ! Zero signals a nonuniform grid, for which searches remain necessary.
  pure real(rk) function uniform_inverse_step(values) result(inverse_step)
    real(rk), intent(in) :: values(:)
    real(rk) :: step,tolerance
    integer :: i,n
    n = size(values)
    step = (values(n)-values(1))/real(n-1,rk)
    tolerance = 32.0_rk*epsilon(1.0_rk)*max(1.0_rk,abs(values(1)),abs(values(n)))
    inverse_step = 0
    do i = 2,n-1
      if (abs(values(i)-(values(1)+real(i-1,rk)*step)) > tolerance) return
    enddo
    inverse_step = 1.0_rk/step
  end function

  ! Same lower-bound semantics as the binary search, with direct indexing when
  ! spacing is known. Neighbour checks retain exact bracketing at rounded knots.
  pure integer function grid_lower_bound(values,x,inverse_step) result(lo)
    real(rk), intent(in) :: values(:),x,inverse_step
    integer :: n
    if (inverse_step <= 0) then
      lo = lower_bound(values,x)
      return
    endif
    n = size(values)
    if (x <= values(1)) then
      lo = 1
      return
    elseif (x > values(n)) then
      lo = n+1
      return
    endif
    ! x is within the grid, so INT acts as FLOOR. Clamp before conversion.
    lo = 1+int(min(real(n-1,rk),max(0.0_rk,(x-values(1))*inverse_step)))
    if (values(lo) < x) lo = min(n,lo+1)
    if (lo > 1) then
      if (values(lo-1) >= x) lo = lo-1
    endif
    ! A fallback also protects unusually ill-conditioned floating-point grids.
    if (values(lo) < x) then
      lo = lower_bound(values,x)
    elseif (lo > 1) then
      if (values(lo-1) >= x) lo = lower_bound(values,x)
    endif
  end function

  ! First index whose value is >= x; n+1 when all values are smaller.
  pure integer function lower_bound(values,x) result(lo)
    real(rk), intent(in) :: values(:),x
    integer :: hi,mid
    lo = 1
    hi = size(values)+1
    do while (lo < hi)
      mid = lo+(hi-lo)/2
      if (values(mid) < x) then
        lo = mid+1
      else
        hi = mid
      endif
    enddo
  end function

  ! Sampling of the piecewise-linear profile on the actual output grid.
  ! All saved profile data are read-only here; intens must be thread-local.
  subroutine do_grid_sampling(tranfreq,abscoef,freq,cutoff,it,irrep_up,irrep_low,intens)
    real(rk), intent(in) :: tranfreq,abscoef,freq(:),cutoff
    integer(ik), intent(in) :: it,irrep_up,irrep_low
    real(rk), intent(inout) :: intens(:)
    integer :: is
    real(rk) :: weight
    do is = 1,n_profile_sym
      weight = profile_weights(is,irrep_up,irrep_low)
      if (weight == 0) cycle
      ! Evaluate each component on its own grid before forming the average.
      call do_grid_component_sampling(tranfreq,abscoef*weight,freq,cutoff,it,is,intens)
    enddo
  end subroutine

  subroutine do_grid_component_sampling(tranfreq,abscoef,freq,cutoff,it,is,intens)
    real(rk), intent(in) :: tranfreq,abscoef,freq(:),cutoff
    integer, intent(in) :: it,is
    real(rk), intent(inout) :: intens(:)
    integer :: ipoint,ib,ie,j,n
    real(rk) :: left,right,x,fraction,value,output_inverse_step,inverse_step
    n = n_profile_grid
    left = max(nu_profile(1,it,is),-cutoff)
    right = min(nu_profile(n,it,is),cutoff)
    if (right < left) return
    ! GRID currently requires a uniform output grid (checked in ReadInput).
    output_inverse_step = real(size(freq)-1,rk)/(freq(size(freq))-freq(1))
    inverse_step = profile_inverse_step(it,is)
    ib = grid_lower_bound(freq,tranfreq+left,output_inverse_step)
    ie = min(size(freq),grid_lower_bound(freq,tranfreq+right,output_inverse_step))
    if (ib > ie) return
    j = max(1,min(n-1,grid_lower_bound(nu_profile(:,it,is),freq(ib)-tranfreq,inverse_step)-1))
    do ipoint = ib,ie
      x = freq(ipoint)-tranfreq
      if (x < left.or.x > right) cycle
      if (x > nu_profile(j+1,it,is)) then
        ! Reuse the neighbouring interval for finely sampled output. Jump over
        ! larger gaps instead of walking through every skipped profile point.
        j = min(j+1,n-1)
        if (x > nu_profile(j+1,it,is)) &
          j = max(1,min(n-1,grid_lower_bound(nu_profile(:,it,is),x,inverse_step)-1))
      endif
      fraction = (x-nu_profile(j,it,is))/(nu_profile(j+1,it,is)-nu_profile(j,it,is))
      value = (1.0_rk-fraction)*f_profile(j,it,is)+fraction*f_profile(j+1,it,is)
      intens(ipoint) = intens(ipoint)+abscoef*value
    enddo
  end subroutine

  subroutine free_grid_profiles()
    if (allocated(profile_weights)) then
      deallocate(profile_weights,missing_band)
      call ArrayStop('GRID:profile_weights')
      call ArrayStop('GRID:missing_band')
    endif
    if (allocated(nu_profile)) then
      deallocate(nu_profile,f_profile,profile_inverse_step)
      call ArrayStop('GRID:nu_profile')
      call ArrayStop('GRID:f_profile')
      call ArrayStop('GRID:profile_inverse_step')
    endif
    if (allocated(file_index)) then
      deallocate(file_index)
      call ArrayStop('GRID:file_index')
    endif
    if (allocated(labels)) deallocate(labels)
    if (allocated(files)) deallocate(files)
    n_profile_sym = 0
    n_profile_t = 0
    n_profile_grid = 0
    grid_profiles_do = .false.
  end subroutine
end module grid_profiles
