! Minimal irrep algebra, following TROVE's SymmetryT/SymmetryInitialize interface.
! Product multiplicities describe the algebra only; profile weights belong in GRID.
module symmetry
  use accuracy, only: ik, cl, out
  implicit none
  private
  public :: SymmetryT, sym, SymmetryInitialize, SymmetryClear, irrep_index, gamma_lookup

  type :: SymmetryT
    character(len=cl) :: group = ''
    integer(ik) :: Nirreps = 0
    character(len=20), allocatable :: label(:)
    integer(ik), allocatable :: degen(:)
    ! product(gamma, upper, lower) is the multiplicity of gamma in upper x lower.
    integer(ik), allocatable :: product(:,:,:)
  end type
  type(SymmetryT), save, protected :: sym

contains

  function uppercase(text) result(word)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: word
    integer :: i,k
    word = adjustl(text)
    do i = 1,len_trim(word)
      k = iachar(word(i:i))
      if (k >= iachar('a').and.k <= iachar('z')) word(i:i) = achar(k-32)
    enddo
  end function

  subroutine symmetry_error(message)
    character(len=*), intent(in) :: message
    write(out,'(a)') 'SYMMETRY: '//trim(message)
    error stop 1
  end subroutine

  subroutine SymmetryInitialize(group)
    character(len=*), intent(in) :: group
    integer(ik) :: i,j,n
    integer(ik), parameter :: c2v_table(4,4) = reshape([ &
      1,2,3,4, 2,1,4,3, 3,4,1,2, 4,3,2,1], [4,4])
    integer(ik), parameter :: cs_table(2,2) = reshape([1,2,2,1], [2,2])

    if (sym%Nirreps /= 0) call symmetry_error('The group has already been defined')
    select case(trim(uppercase(group)))
    case ('C2V','C2V(M)')
      sym%group = 'C2v'
      n = 4
    case ('C3V','C3V(M)')
      sym%group = 'C3v'
      n = 3
    case ('CS','CS(M)')
      sym%group = 'Cs'
      n = 2
    case default
      call symmetry_error('Unsupported group '//trim(group)//'; expected C2v, C3v or Cs')
    end select
    sym%Nirreps = n
    allocate(sym%label(n),sym%degen(n),sym%product(n,n,n))
    sym%degen = 1
    sym%product = 0

    select case(trim(sym%group))
    case ('C2v')
      sym%label = [character(len=20) :: 'A1','A2','B1','B2']
      do j = 1,n
        do i = 1,n
          sym%product(c2v_table(i,j),i,j) = 1
        enddo
      enddo
    case ('C3v')
      sym%label = [character(len=20) :: 'A1','A2','E']
      sym%degen = [1,1,2]
      do i = 1,n
        sym%product(i,1,i) = 1 ! A1 x gamma = gamma
        sym%product(i,i,1) = 1
      enddo
      sym%product(1,2,2) = 1 ! A2 x A2 = A1
      sym%product(3,2,3) = 1 ! A2 x E = E
      sym%product(3,3,2) = 1
      sym%product(:,3,3) = [1,1,1] ! E x E = A1 + A2 + E
    case ('Cs')
      sym%label = [character(len=20) :: "A'",'A"']
      do j = 1,n
        do i = 1,n
          sym%product(cs_table(i,j),i,j) = 1
        enddo
      enddo
    end select

    ! Catch incomplete or dimensionally inconsistent tables when adding groups.
    do j = 1,n
      do i = 1,n
        if (any(sym%product(:,i,j) < 0)) call symmetry_error('Negative product multiplicity')
        if (any(sym%product(:,i,j) /= sym%product(:,j,i))) call symmetry_error('Noncommutative irrep product')
        if (sum(sym%product(:,i,j)*sym%degen) /= sym%degen(i)*sym%degen(j)) &
          call symmetry_error('Inconsistent product dimensions')
      enddo
    enddo
    write(out,'(a,a,a,i0)') 'SYMMETRY: ',trim(sym%group),', Nirreps = ',sym%Nirreps
  end subroutine

  integer(ik) function irrep_index(label) result(irrep)
    character(len=*), intent(in) :: label
    character(len=len(label)) :: word
    integer(ik) :: i
    word = uppercase(label)
    ! Accept two apostrophes as an alternative spelling of the Cs double prime.
    if (trim(sym%group) == 'Cs'.and.trim(word) == "A''") word = 'A"'
    irrep = 0
    do i = 1,sym%Nirreps
      if (word /= sym%label(i)) cycle
      irrep = i
      return
    enddo
  end function

  function gamma_lookup(irrep_up,irrep_low) result(multiplicity)
    integer(ik), intent(in) :: irrep_up,irrep_low
    integer(ik) :: multiplicity(sym%Nirreps)
    if (sym%Nirreps == 0) call symmetry_error('Define a symmetry group before gamma_lookup')
    if (min(irrep_up,irrep_low) < 1.or.max(irrep_up,irrep_low) > sym%Nirreps) &
      call symmetry_error('Irrep index outside the group in gamma_lookup')
    multiplicity = sym%product(:,irrep_up,irrep_low)
  end function

  subroutine SymmetryClear()
    if (allocated(sym%label)) deallocate(sym%label,sym%degen,sym%product)
    sym%group = ''
    sym%Nirreps = 0
  end subroutine
end module symmetry
