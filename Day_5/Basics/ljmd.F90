!***********************************************************************
!
!     CCP5 Summer School Lennard Jones Molecular Dynamics Program
!
!     copyright Daresbury Laboratory
!
!     originally written by W Smith
!     modified and updated by J Purton
!     modernised by A Elena July 2022
!***********************************************************************

program lj_nve
  use iso_fortran_env, only : wp=>real64, iostat_end
  implicit none

  integer, parameter :: matms = 4
  integer, parameter :: ngrid = 128
  real(kind=wp) ,parameter ::  xn(4) = [-0.25_wp,0.25_wp,-0.25_wp,0.25_wp] ,&
    yn(4) = [-0.25_wp,0.25_wp,0.25_wp,-0.25_wp],&
    zn(4) = [-0.25_wp,-0.25_wp,0.25_wp,0.25_wp]
  real(kind=wp), parameter ::  pi = 4.0_wp*atan(1.0_wp)
  real(kind=wp), parameter :: bol =1.5_wp
  integer, parameter :: X = 1, Y = 2, Z = 3

  logical :: ltraj, lpos
  integer :: i, j, k, l, n, istep, nstep, neql, nav, isamp, nrep
  integer :: natms
  integer :: ur, ut, us

  real(kind=wp), allocatable ::  r(:,:),v(:,:),f(:,:)
  real(kind=wp) :: gofr(ngrid),m(3)
  real(kind=wp) :: sid,  scl, tst, rct,  time,  nrm
  real(kind=wp) :: tot, totav, pot, vir, kin, tmp, tmpav, den, tem, dtm, vol
  real(kind=wp) :: prs, prsav, box, dr

  integer :: is, ip
  real(kind=wp) :: tf, ts

  call cpu_time(ts)
  natms = 4 * matms * matms * matms

  allocate(r(3,natms), v(3,natms), f(3,natms))
  write(*,'(a)')'CCP5 Summer School Lennard Jones MD Program'

  ! read control data
  open(newunit=ur, file="MDIN", action="read",status="old")
  read(ur,*) nstep
  read(ur,*) neql
  read(ur,*) isamp
  read(ur,*) den
  read(ur,*) tem
  read(ur,*) dtm
  read(ur,*) ltraj
  read(ur,*) lpos
  read(ur,*) is,ip
  close(ur)
  write(*,'(a32,1x,i0)')'Number of time steps', nstep
  write(*,'(a32,1x,i0)')'Number of equilibration steps', neql
  write(*,'(a32,1x,i0)')'Properties sampling interval', isamp
  write(*,'(a32,1x,i0)')'Number of atoms', natms
  write(*,'(a32,1x,es16.8)')'System density', den
  write(*,'(a32,1x,es16.8)')'System temperature', tem
  write(*,'(a32,1x,es16.8)')'Simulation time step', dtm
  write(*,'(a32,1x,l4)')'Trajectory file switch', ltraj
  write(*,'(a32,1x,l4)')'Trajectory file content', lpos
  if (lpos) then
    write(*,'(a32,1x,a)')'Trajectory file', 'positions'
  else
    write(*,'(a32,1x,a)')'Trajectory file', 'velocities'
  end if
  write(*,'(a32,1x,2(i0,1x))')'Random seed', is,ip

!set up simulation
  nav = 0
  tot = 0.0_wp
  totav = 0.0_wp
  prsav = 0.0_wp
  tmpav = 0.0_wp
  tst = dtm * 0.5_wp
  vol = real(natms,wp) / den
  box = vol**(1.0_wp / 3.0_wp)
  sid = box/real(matms,wp)
  rct = 2.5_wp
  dr = rct / real(ngrid,wp)
  nrep = max(1,nstep / 100)

  write(*,'(a32,1x,es16.8)')'Volume', vol
  write(*,'(a32,1x,es16.8)')'Box side', box
  write(*,'(a32,1x,es16.8)')'cutoff', rct
  write(*,'(a32,1x,es16.8)')'bin size', dr
  write(*,'(a32,1x,i0)')'print frequency', nrep
  write(*,'(a)')"check these parameters are what you want !!!"

  !set initial atom positions
  k = 1
  do n = 1,4
    do i = 1,matms
      do j = 1,matms
        do l = 1,matms
          r(X,k) = (real(l,wp) - 0.5_wp + xn(n)) * sid - 0.5_wp * box
          r(Y,k) = (real(j,wp) - 0.5_wp + yn(n)) * sid - 0.5_wp * box
          r(Z,k) = (real(i,wp) - 0.5_wp + zn(n)) * sid - 0.5_wp * box
          k = k + 1
        enddo
      enddo
    enddo
  enddo

  !initialise the rdf array
  gofr = 0.0_wp

! set random seed
  call init_random_seed(is,ip)

  ! set initial velocities
  call init_velocities(natms,tem,v)

  !open statistics file
  open(newunit = us,file = "STA",action="write",status='unknown')

  !open trajectory file
  if(ltraj)open(newunit=ut,file = "TRJ",action="write",status="unknown")

  !calculate initial forces
  call forces(natms, pot, vir, rct, box, dr, r, f, gofr,.false.)
  kin = compute_kinetic_energy(natms,v)

  tmp = temperature(kin,natms)
  tot = kin + pot
  prs = pressure(kin,vir,vol)

  !velocity verlet algorithm
  call write_xyz("initial.xyz",r,natms,box)
  write( * ,'(/,10x,a4,9x,a7,12x,a4,11x,a5)')'time','tot_engy','temp','press'
  write( * ,'(1x,1p,4(e16.8,1x))')0.0, tot, tmp, prs

  do istep = 1, nstep

    ! first half step motion
    do k = 1,natms
      v(:,k) = v(:,k) + tst * f(:,k)

      r(:,k) = r(:,k) + dtm * v(:,k)
    enddo

    ! calculate new forces
    call forces(natms, pot, vir, rct, box, dr, r, f, gofr,istep > neql)

    ! second half step motion
    do k = 1, natms
      v(:,k) = v(:,k) + tst * f(:,k)
    enddo

    kin = compute_kinetic_energy(natms,v)

    !restore periodic boundary
    do k = 1,natms
      r(:,k) = r(:,k) - box * nint(r(:,k) / box)
    enddo

    ! system properties
    tmp = temperature(kin,natms)
    tot = kin + pot
    prs = pressure(kin,vir,vol)

    if(istep > neql)then
      nav = nav + 1
      tmpav = running_average(tmpav, tmp, nav)
      totav = running_average(totav, tot, nav)
      prsav = running_average(prsav, prs, nav)
    end if

    ! checkpoint report
    time = dtm * istep
    if(mod(istep,nrep) == 0) then
      write( * ,'(1x,1p,4(e16.8,1x))')time, tot, tmp, prs
    end if
    if(istep == neql) write( * ,'(a)')'  End of equlibration'

    ! write statistics file
    write(us,'(1p,4(e16.8,1x))')time, tot, tmp, prs

    if (mod(istep,isamp) == 0) then
      if(ltraj)then
        call write_traj(ut,r,v,natms,time,box,lpos)
      endif

      ! temperature scaling
      if(istep <= neql) then
        scl = sqrt(tem / tmp)
        do k = 1, natms
          v(:,k) = scl * v(:,k)
        enddo
      endif
    endif
  enddo

  !final average energy temperature and pressure
  write( * ,'(a)')'final average energy temperature and pressure:'
  write( * ,'(1x,i10,1p,3(e16.8,1x))')nav, totav, tmpav, prsav

  call write_xyz("final.xyz",r,natms,box)

!    calculate final rdf
  nrm = 1.0_wp / (2.0_wp * pi * dr * den * real(nstep - neql,wp) * real(natms - 1,wp))
  call write_rdf("RDF",ngrid,gofr,nrm,dr)
  !close files
  if(ltraj)close(ut)
  close(us)

  deallocate(r,v,f)
  call cpu_time(tf)
  write( * ,'(a,f16.8,a)')'Job done in ',tf-ts,'s'
contains

  subroutine init_random_seed(is,ip)
    integer, intent(in) :: is,ip

    integer, allocatable :: seed(:)
    integer :: i,p

    call random_seed(size=p)
    allocate(seed(p))
    seed = ip*[(i-is,i=1,p)]
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_random_seed

  subroutine write_rdf(filename,ngrid,gofr,nrm,dr)
    integer, intent(in) :: ngrid
    real(kind=wp), intent(in) :: nrm,dr
    real(kind=wp), intent(inout) :: gofr(:)
    character(len=*), intent(in) :: filename

    integer :: ux,i
    real(kind=wp) :: d

    open(newunit=ux,file = trim(filename), action="write", status="unknown")
    do i = 1, ngrid
      d = (real(i,wp) - 0.5_wp) * dr
      gofr(i) = nrm * gofr(i) / (d**2 + dr**2 / 12.0_wp)
      write(ux,'(2(f16.8,1x))')d, gofr(i)
    enddo
    close(ux)
  end subroutine write_rdf

  subroutine write_traj(ut,r,v,natoms,time,box,lpos)
    integer, intent(in) :: ut,natoms
    real(kind=wp), intent(in) :: time,box,r(:,:),v(:,:)
    logical,intent(in) :: lpos

    integer :: i

    write(ut,'(i0,2(f16.8,1x))') natms, time, box
    if(lpos) then
      do i = 1,natoms
        write(ut,'(3(f16.8,1x))') r(:,i)
      enddo
    else
      do i = 1,natoms
        write(ut,'(3(f16.8,1x))') v(:,i)
      enddo
    endif

  end subroutine write_traj

  pure real(kind=wp) function pressure(kin,vir,vol)
    real(kind=wp),intent(in) :: kin,vir,vol

    pressure = (2.0_wp * kin - vir) / (3.0_wp * vol)
  end function pressure

  pure real(kind=wp) function temperature(kin,natoms)
    real(kind=wp),intent(in) :: kin
    integer, intent(in) :: natoms

    temperature = kin/(bol*real(natoms,wp))

  end function temperature

  pure real(kind=wp) function running_average(cav, c ,n)
    real(kind=wp), intent(in) :: cav, c
    integer, intent(in) :: n

    running_average = ((cav * (n - 1)) / n) + (c / n)

  end function running_average

  subroutine write_xyz(filename,r,natoms,b)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: r(:,:),b
    integer, intent(in) :: natoms

    integer :: ux
    real(wp) :: z

    z = 0.0_wp

    open(newunit=ux,file = trim(filename),action="write",status="unknown")
    write(ux,'(i0)')natms
    write(ux,'(a,9(g0,1x),a)')'Lattice="', &
                                   b,z,z, &
                                   z,b,z, &
                                   z,z,b, &
                                   '" Properties=species:S:1:pos:R:3 pbc="T T T"'

    do i=1,natoms
      write(ux,'(*(a4,1x,3(f16.8,1x)))')'Ar  ',r(:,i)
    end do
    close (ux)

  end subroutine write_xyz

  pure subroutine forces(natms, pot, vir, rct, box, dr, r, f, gofr, lrdf)
! * **********************************************************************
!
!     program ljmd  -  calculate lennard jones forces
!
!     copyright Daresbury Laboratory
!     author W. Smith
!
! * **********************************************************************
    real(kind=wp), intent(in) :: r(:,:), box, rct, dr
    real(kind=wp), intent(out) :: f(:,:), pot, vir
    real(kind=wp), intent(inout) :: gofr(:)
    logical,intent(in) :: lrdf
    integer, intent(in) :: natms

    integer :: i, j, ix
    real(kind=wp) :: rc2, alpha, beta, dd(3)
    real(kind=wp) :: sg6, sg12, gamma, rrs, rr0, rr3, d, rsq

    pot = 0.0_wp
    vir = 0.0_wp
    rc2 = rct * rct
    alpha = 24.0_wp * (2.0_wp * rct**(-12) - rct**(-6)) / rct
    beta = 4.0_wp * (rct**(-12) - rct**(-6)) + alpha * rct

    f = 0.0_wp

    do i = 1,natms - 1
      do j = i + 1,natms

        dd(:) = r(:,j) - r(:,i)

        !minimum image
        dd = dd - box * nint(dd/box)

        !interatomic distance
        rsq = dd(X) * dd(X) + dd(Y) * dd(Y) + dd(Z) * dd(Z)

        !spherical cutoff
        if (rsq < rc2)then
          d = sqrt(rsq)
          !rdf calculation
          if(lrdf)then
            ix  =  int(d / dr) + 1
            gofr(ix) = gofr(ix) + 1.0_wp
          endif

          !potential calculation
          rrs = 1.0_wp / rsq
          rr0 = d * rrs
          rr3 = rr0 * rrs
          sg6 = rr3 * rr3
          sg12 = sg6 * sg6
          pot = pot + 4.0_wp * (sg12 - sg6) + alpha * d - beta

          !force calculation
          gamma = 24.0_wp * rrs * (2.0_wp * sg12 - sg6) - alpha * rr0
          ! if(gamma > 1.0e6_wp)gamma = 1.0e6_wp
          vir = vir - gamma * rsq
          f(:,i) = f(:,i) - gamma * dd(:)
          f(:,j) = f(:,j) + gamma * dd(:)

        endif
      enddo
    enddo

  end subroutine forces

  pure real(kind=wp) function compute_kinetic_energy(natms,v) result(kin)

    integer, intent(in) :: natms
    real(kind=wp), intent(in) :: v(:,:)

    integer :: k

    kin = 0.0_wp
    do k = 1, natms
      kin = kin + v(X,k) * v(X,k) + v(Y,k) * v(Y,k) + v(Z,k) * v(Z,k)
    enddo

    kin = 0.5_wp * kin
  end function compute_kinetic_energy

  subroutine init_velocities(natoms,temp,v)
    integer, intent(in) :: natoms
    real(kind=wp),intent(in) :: temp
    real(kind=wp),intent(inout) :: v(:,:)


    real(kind=wp) :: m(3),tmp,scl
    integer :: k

    call random_number(v)
    m = 0.0_wp
    do k = 1,natoms
      v(:,k) = v(:,k) - 0.5_wp

      m(:) = m(:) + v(:,k)
    enddo

    m = m / real(natoms,wp)

    ! remove net system momentum
    do k = 1,natoms
      v(:,k) = v(:,k) - m(:)
    enddo
    kin = compute_kinetic_energy(natoms,v)

    ! temperature scaling
    tmp = temperature(kin,natoms)
    scl = sqrt(temp / tmp)

    do k = 1,natoms
      v(:,k) = scl * v(:,k)
    enddo
  end subroutine init_velocities
end program lj_nve
