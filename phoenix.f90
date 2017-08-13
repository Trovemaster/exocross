module Phoenix
  !
  use accuracy
  use timer

  implicit none


  private
  public do_gf_oscillator_strength_Phoenix

  !
  ! Phoenix type
  type gauss_linesT
   sequence
   integer(kind=4) :: iwl
   integer(kind=2) :: ielion,ielo,igflog,igr,igs,igw
  end type gauss_linesT

  integer(ik),parameter :: verbose = 2

  contains


  subroutine do_gf_oscillator_strength_Phoenix(itrans,ichunk,iso,nswap,nlines,energies,Jrot,ilevelf_ram,ileveli_ram,gf_ram,&
             abscoeff_ram,jmax,gamma_1,n_1,gamma_2,n_2,output)
     !
     implicit none
     !
     integer(ik),intent(in) :: itrans,ichunk,iso,nswap,nlines,ilevelf_ram(nswap),ileveli_ram(nswap),jmax
     real(rk),intent(in) :: gf_ram(nswap),jrot(nlines),energies(nlines)
     real(rk),intent(in) :: abscoeff_ram(nswap)
     real(rk),intent(in) :: gamma_1(0:jmax,-1:1),gamma_2(0:jmax,-1:1),n_1(0:jmax,-1:1),n_2(0:jmax,-1:1)
     character(len=cl),intent(in) :: output
     !
     real(rk)  :: gf,tranfreq,energyf,energyi,gamma1,gamma2,n1,n2,ji,jf,abscoeff
     integer(ik) :: gfunit
     integer(ik) :: i,ileveli,ilevelf,jp,jpp,wl2ind,il,ib,current_blk,irec_length_ph
     real(rk) :: lambda,ratiolog
     !
     integer(kind=2) :: iloggf,ilog_n1,ilog_n2,ilog_gamma1,ilog_gamma2,ielion,iener
     !
     integer, parameter :: blksize=65536
     type(gauss_linesT) :: oneline(blksize)
     character(len=9) :: c_fmt,t_fmt
     integer(ik) :: iblksize_3,nlines_3,nblock,ibsu,istat
     character(len=cl) :: ioname
     !
     write(ioname, '(a)') 'Phoenix'
     call IOstart(trim(ioname),gfunit)
     !
     ielion = iso
     !
     ! for Phoemix conversion a new .bin will be used for each .trans
     !
     inquire(iolength=irec_length_ph) oneline
     !
     write(t_fmt,'(i5.5)') itrans
     write(c_fmt,'(i5.5)') ichunk
     !
     current_blk = 1
     iblksize_3 = 65536
     nlines_3 = nswap
     nblock =  nswap / iblksize_3
     if (nswap>nblock*iblksize_3) nblock = nblock + 1
     !
     ibsu = 262144 ! value for PIOFS ...
     !
     ratiolog = log(1.0_rk+1._rk/2000000.0_rk)
     !
     !open(13,file=trim(outfile),form='unformatted',access='direct',recl=irec_length_3,action='write')
     open(unit=gfunit,file=trim(output)//trim(t_fmt)//"_"//trim(c_fmt)//".bin",action='write',status='replace',&
                      form='unformatted',access='direct',recl=irec_length_ph)
     !
     write(gfunit,iostat=istat,rec=1) nlines_3,iblksize_3,nblock,ibsu
     !
     if (verbose>=3) write(out,"('lines = ',i8,' chunks = ',i8,' block = ',i5)") nswap,nswap/iblksize_3,iblksize_3
     !
     do ib = 1,nswap,blksize
       !
       oneline(:)%iwl = -1
       oneline(:)%ielion = -1
       oneline(:)%ielo = -1
       oneline(:)%igflog = -1
       oneline(:)%igr = -1
       oneline(:)%igs = -1
       oneline(:)%igw = -1
       !oneline(:)%ign = -1
       !
       if (verbose>=4) write(out,"(' block = ',i8)") current_blk
       !
       loop_block :  do il = 1,blksize
         !
         if (ib+il==2) cycle
         !
         ! the wavenumbers are sorted in increasing order, therefore we start from the last 
         ! to get the libe list sorted with wavelength increasing 
         !
         i = nswap - (ib-1+il)
         !
         if (i<1) exit
         !
         ilevelf = ilevelf_ram(i)
         ileveli = ileveli_ram(i)
         gf   = gf_ram(i)
         abscoeff   = abscoeff_ram(i)
         energyf = energies(ilevelf)
         energyi = energies(ileveli)
         tranfreq = energyf - energyi
         !
         lambda = 1e7/tranfreq ! wavelength in nm 
         !
         wl2ind = int(log(lambda)/ratiolog+0.5_rk)
         !
         !g_l*fij = c2 * Aji * g_u / (c*wn)**2
         !gf = 1.34738e+21*acoef*g_u / (c*tranfreq)**2
         !
         iloggf = int(log(gf)*(1.0_rk/(log(10.0_rk)*0.001_rk)),kind=2) + 2**14
         iener  = int(log(energyi)*(1.0_rk/(log(10.0_rk)*0.001_rk)),kind=2) + 2**14
         !
         jf = jrot(ilevelf)
         ji = jrot(ileveli)
         !
         !gamma1 = species(1)%gamma
         !n1     = species(1)%n
         !
         Jpp = nint(Ji)
         Jp  = nint(Jf)
         !
         !if ( trim(species(1)%filename)/="" ) then
         gamma1 = gamma_1(Jpp,Jp-Jpp)
         n1     = n_1(Jpp,Jp-Jpp)
         !endif
         !
         ilog_gamma1 = int(log(gamma1)*(1.0_rk/(log(10.0_rk)*0.001_rk)),kind=2) + 2**14
         ilog_n1 = int(log(n1)*(1.0_rk/(log(10.0_rk)*0.001_rk)),kind=2) + 2**14
         !
         !gamma2 = voigt_gamma
         !
         gamma2 = gamma_2(Jpp,Jp-Jpp)
         n2     = n_2(Jpp,Jp-Jpp)
         !
         ilog_gamma2 = int(log(gamma2)*(1.0_rk/(log(10.0_rk)*0.001_rk)),kind=2) + 2**14
         ilog_n2 = int(log(n2)*(1.0_rk/(log(10.0_rk)*0.001_rk)),kind=2) + 2**14
         !
         oneline(il)%iwl    = wl2ind
         oneline(il)%ielion = ielion
         oneline(il)%ielo   = iener
         oneline(il)%igflog = iloggf
         oneline(il)%igr    = ilog_gamma1
         oneline(il)%igs    = ilog_n1
         oneline(il)%igw    = ilog_gamma2
         !oneline(il)%ign = ilog_n2
         !
         if (verbose>=4) write(out,"(i7,i12,1x,6(1x,i7),f16.6,1x,e16.6,1x,f16.6,1x,f16.6)") ielion,wl2ind,iener,& 
                                    iloggf,ilog_gamma1,ilog_n1,&
                                    ilog_gamma2,ilog_n2,tranfreq,abscoeff,energyi,lambda
         !
         !write(out,"(f16.6,1x,e16.6,1x,f16.6)") tranfreq,abscoeff,energyi
         !
       enddo loop_block
       !
       current_blk = current_blk+1
       !
       write(gfunit,rec=current_blk) oneline
       !
     enddo
     !
     oneline(:)%iwl = -1
     oneline(:)%ielion = -1
     oneline(:)%ielo = -1
     oneline(:)%igflog = -1
     oneline(:)%igr = -1
     oneline(:)%igs = -1
     oneline(:)%igw = -1
     !
     current_blk = current_blk+1
     !
     write(gfunit,rec=current_blk) oneline
     !
     !write(gfunit,rec=ib) int(-1,kind=2),int(-1,kind=4),int(-1,kind=2),int(-1,kind=2),int(-1,kind=2),int(-1,kind=2),int(-1,kind=2)
     !
     close(gfunit,status='keep')
     !
     !ielion, wl, gf-value, xi, igr, igs, igw
     !where ielion: internal code of the molecule for Phoenix,
     !wl:wavelength,
     !xi:lower state energy,
     !gf-value: log gf coded in integer format,
     !igr: natural width constant,
     !igs: Stark broadening constant,
     !igw: van der Waals damping constants.
     !igr for gamma0 and igs for n (Voigt parameters)
     !
  end subroutine do_gf_oscillator_strength_Phoenix


end module Phoenix
