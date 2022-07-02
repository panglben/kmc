! A fortran95 program for G95
! By yl
module outer_data
    implicit none
    ! This module is used to read and save parameters form 'input_tot'. All parameters is read-only.
    
    integer(kind=1), protected :: restart  
    integer(kind=8), protected :: nLoop 
    
    real(kind=8), protected :: s0_CO_edge, s0_CO_facet, s0_O2_edge, s0_O2_facet
    real(kind=4), protected :: EdiffCO, EdiffO2, EdiffO
    real(kind=8), protected :: Eads_CO_a, Eads_CO_b, Eads_O_a, Eads_O_b, Eads_O2_a, Eads_O2_b
    real(kind=8), protected :: BEP1_a, BEP1_b
    real(kind=8), protected :: BEP2_a, BEP2_b
    real(kind=8), protected :: Ea3
    real(kind=4), protected :: dimx, dimy, dimz
    
    character(len=2) :: elem
    real(kind=4), protected :: latt_para, Temp, p_tot, ppCO, ppO2
    real(kind=4), protected :: Eco_co, Eo2_o2, Eo_o
    real(kind=4), protected :: Eco_o2, Eco_o, Eo2_o

    contains 
    
    subroutine read_data
        open(1, file='input_tot', action='read')
            read(1, *) restart ! 0: start this_one, 1: start last_one
            read(1, *) nLoop ! steps
            read(1, *) dimx, dimy, dimz ! boxsizex, boxsizey, boxsizez
            read(1, *) s0_CO_edge, s0_CO_facet ! sticking coffecient of CO
            read(1, *) s0_O2_edge, s0_O2_facet ! sticking coffecient of O
            read(1, *) EdiffCO, EdiffO, EdiffO2 ! diffusion barrier
            read(1, *) Eads_CO_a, Eads_CO_b ! Eads = a*gcn + b
            read(1, *) Eads_O2_a, Eads_O2_b
            read(1, *) Eads_O_a, Eads_O_b
            read(1, *) BEP1_a, BEP1_b ! BEP coeffient for CO* + O2* -- OCOO*
            read(1, *) BEP2_a, BEP2_b ! BEP coeffient for OCOO* -- CO2 + O*
            read(1, *) Ea3 ! Ea for CO + O* -- CO2

            read(1, *) elem
            read(1, *) latt_para
            read(1, *) Temp, p_tot
            read(1, *) ppCO, ppO2 ! partial pressure of CO and O2
            read(1, *) Eco_co, Eo2_o2, Eo_o ! Interaction between CO-CO; O2-O2; O-O
            read(1, *) Eco_o2, Eco_o, Eo2_o ! Interaction between CO-O2; CO-O; O2-O
        close(1)
    end subroutine read_data    
end module outer_data
    
module struct_data
    use outer_data, only : elem, latt_para, dimx, dimy, dimz
    implicit none
    ! This module is used to read and record structure data of NPs
    
    integer(kind=4), private :: i, j, k
    
    integer(kind=4), protected :: nbulk ! number of atoms in bulk
    integer(kind=2), private :: nlx, nly, nlz
    real(kind=4), protected :: clx, cly, clz
    real(kind=4), protected :: bx, by, bz
    real(kind=4), protected :: xxx(100000), yyy(100000), zzz(100000)  ! coordinate of the grid
    
    integer(kind=4), protected :: natoms ! number of atoms in NPs
    real(kind=4), protected :: x(100000), y(100000), z(100000)  ! coordinate of the NP
    
    real(kind=8), protected :: gcn(100000)
    integer(kind=4), protected :: cn(100000), cn_2(100000), effcn(100000)
    integer(kind=4), protected :: nnsite(12, 100000), nnnsite(54, 100000)
    integer(kind=4), public :: cov_type(100000) ! 0-none; 1-CO; 2-O2; 3-O; 4-OCOO_CO; -4-OCOO_OO; 5-bulk
    integer(kind=4), public :: site_type(100000) ! 1-np; 2-surface site; 3-outer grid
    integer(kind=4), public :: co_ads_site(100000) ! record the # the coads site
    
    contains
    
    subroutine calc_para
        clx = latt_para
        cly = latt_para
        clz = latt_para
        nlx = int(dimx/clx)
        nly = int(dimy/cly)
        nlz = int(dimz/clz)
        bx = (nlx-0.5)*clx
        by = (nly-0.5)*cly
        bz = (nlz-0.5)*clz
    end subroutine calc_para
    
    subroutine read_str_new
        real(kind=4) :: x3_sum, y3_sum, z3_sum
        real(kind=4) :: x3cent, y3cent, z3cent, dx3a, dy3a, dz3a
        real(kind=4) :: xcent, ycent, zcent, dxa, dya, dza
        real(kind=4) :: x_sum, y_sum, z_sum
        real(kind=4) :: dx, dy, dz
        
        call calc_para
        ! generate grid
        nbulk = 0
        do i = 1, nlx
            do j = 1, nly
                do k = 1, nlz
                    nbulk = nbulk + 1
                    xxx(nbulk) = i*clx
                    yyy(nbulk) = j*cly
                    zzz(nbulk) = k*clz
                    nbulk = nbulk + 1
                    xxx(nbulk) = (0.5+i)*clx
                    yyy(nbulk) = (0.5+j)*cly
                    zzz(nbulk) = k*clz
                    nbulk = nbulk + 1
                    xxx(nbulk) = i*clx
                    yyy(nbulk) = (0.5+j)*cly
                    zzz(nbulk) = (0.5+k)*clz
                    nbulk = nbulk + 1
                    xxx(nbulk) = (i+0.5)*clx
                    yyy(nbulk) = j*cly
                    zzz(nbulk) = (0.5+k)*clz
                end do
            end do
        end do
        x3_sum = 0
        y3_sum = 0
        z3_sum = 0
        do i = 1, nbulk
            x3_sum = x3_sum + xxx(i)
            y3_sum = y3_sum + yyy(i)
            z3_sum = z3_sum + zzz(i)
        end do
        x3cent = real(NINT(x3_sum/nbulk))
        y3cent = real(NINT(y3_sum/nbulk))
        z3cent = real(NINT(z3_sum/nbulk))
        ! center bulk sites
        do i = 1, nbulk
            dx3a = abs(xxx(i) - x3cent)
            dy3a = abs(yyy(i) - y3cent)
            dz3a = abs(zzz(i) - z3cent)
            if(dx3a.lt.3.0 .and. dy3a.lt.3.0 .and. dz3a.lt.3.0) then
                x3cent = xxx(i)
                y3cent = yyy(i)
                z3cent = zzz(i)
                exit
            end if
        end do
        ! output grid to 'bulk.xyz'
        open(1, file='bulk.xyz')
            write(1, *) nbulk
            write(1, *)
            do i=1,nbulk
                xxx(i)=xxx(i)-x3cent
                yyy(i)=yyy(i)-y3cent
                zzz(i)=zzz(i)-z3cent
                write(1, *) elem, xxx(i), yyy(i), zzz(i)
                site_type(i) = 3  ! fake site
                cov_type(i) = 5   ! do not adsorb  
                co_ads_site(i) = 0 
            end do
        close(1)
        
        ! read ini.xyz file and get mass center
        x_sum = 0
        y_sum = 0
        z_sum = 0
        open(1, file='ini.xyz')
            read(1, *) natoms
            read(1, *)
            do i = 1, natoms
            read(1, *) elem, x(i), y(i), z(i)
            x_sum = x_sum + x(i)
            y_sum = y_sum + y(i)
            z_sum = z_sum + z(i)
            end do
        close(1)
        xcent = real(NINT(x_sum/natoms))
        ycent = real(NINT(y_sum/natoms))
        zcent = real(NINT(z_sum/natoms))
        ! center atom site
        do i = 1, natoms
            dxa = abs(x(i) - xcent)
            dya = abs(y(i) - ycent)
            dza = abs(z(i) - zcent)
            if(dxa.lt.3.0 .and. dya.lt.3.0 .and. dza.lt.3.0) then
                xcent = x(i)
                ycent = y(i)
                zcent = z(i)
                exit
            end if
        end do
        ! output new initial structure
        open(1, file='ini_new.xyz')
            write(1, *) natoms
            write(1, *)
            do i = 1, natoms
            x(i) = x(i) - xcent
            y(i) = y(i) - ycent
            z(i) = z(i) - zcent
            write(1, *) elem, x(i), y(i), z(i)
            end do
        close(1)
        ! record NPs on grid (site_type = 1)
        do i = 1,natoms
            do j = 1, nbulk
                dx = abs(x(i)-xxx(j))
                dy = abs(y(i)-yyy(j))
                dz = abs(z(i)-zzz(j))
                if(dx.le.0.001 .and. dy.le.0.001 .and. dz.le.0.001) then
                    site_type(j) = 1
                    exit
                end if
            end do
        end do
    end subroutine read_str_new
    
    ! generate grids and record them on "bulk.xyz "
    subroutine read_str_last
        real(kind=4) :: dx, dy, dz
        integer(kind=4) :: cov_type_ini(100000)
        
        call calc_para
        ! read bulk
        open(11, file = 'bulk.xyz')
        read(11, *) nbulk
        read(11, *)
        do i = 1, nbulk
            read(11, *) elem, xxx(i), yyy(i), zzz(i)
            site_type(i) = 3
            cov_type(i) = 5
        end do
        close(11)

        ! read ini file and cov_type
        open(11, file = 'last_one.xyz')
        read(11, *) natoms
        read(11, *)
        do i = 1, natoms
            read(11, *) elem, x(i), y(i), z(i), cov_type_ini(i), co_ads_site(i)
        end do
        close(11)

        ! bulk sites and covs classification
        do i = 1,natoms
            do j = 1, nbulk
            dx = abs(x(i)-xxx(j))
            dy = abs(y(i)-yyy(j))
            dz = abs(z(i)-zzz(j))
            if(dx.le.0.001 .and. dy.le.0.001 .and. dz.le.0.001) then
                site_type(j) = 1
                cov_type(j) = cov_type_ini(i)
                exit
            end if
            end do
        end do
    end subroutine read_str_last
    
    subroutine count_cn
        ! record cn (nnsite) and scn (nnnsite) of each atom
        real(kind=4) :: rx, ry, rz, dr
        integer(kind=4) :: nn_j
        cn = 0
        cn_2 = 0
        effcn = 0
        nnsite = 0
        nnnsite = 0
        do i = 1, nbulk
            do j = i + 1, nbulk
                rx = abs(xxx(i) - xxx(j))
                ry = abs(yyy(i) - yyy(j))
                rz = abs(zzz(i) - zzz(j))
                ! extended boundary
                if(rx * 2 .gt. bx) rx = (bx - rx) + 0.5 * clx
                if(ry * 2 .gt. by) ry = (by - ry) + 0.5 * cly
                if(rz * 2 .gt. bz) rz = (bz - rz) + 0.5 * clz
                dr = sqrt(rx**2 + ry**2 + rz**2)
                if(dr.lt.3.0) then
                    cn(i) = cn(i) + 1
                    cn(j) = cn(j) + 1
                    nnsite(cn(i), i) = j
                    nnsite(cn(j), j) = i
                    if(site_type(i).lt.3 .and. site_type(j).lt.3) then
                        effcn(i) = effcn(i) + 1
                        effcn(j) = effcn(j) + 1
                    end if
                end if
                if(dr.lt.5.8) then
                    cn_2(i) = cn_2(i) + 1
                    cn_2(j) = cn_2(j) + 1
                    nnnsite(cn_2(i), i) = j
                    nnnsite(cn_2(j), j) = i
                end if
            end do
        end do
        ! calculate and record gcn of each atom
        gcn = 0
        do i = 1, nbulk
            if(site_type(i).lt.3) then
                do j = 1, 12
                    nn_j = nnsite(j, i)
                    gcn(i) = gcn(i) + effcn(nn_j)
                end do
                gcn(i) = gcn(i)/12.0
            end if
        end do
    end subroutine count_cn
end module struct_data
    
subroutine init_random_seed()
    integer(kind=4) :: i, n ,clock
    integer, dimension(:), allocatable :: seed  ! ArrayList

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count = clock)

    seed = clock + 37*(/(i-1, i=1, n)/)
    call random_seed(put=seed)

    deallocate(seed)
end subroutine

program main
    use outer_data
    use struct_data
    implicit none

    integer(kind=4), parameter :: n_process = 13 ! number of processes
    ! 2,3: CO ads & dis;   4,5: O2 ads & dis;   6,7,8: CO,O2,O diff;   9: OCOO diff
    ! 10: CO*+O2*--OCOO*;   11: OCOO*--CO2+O*;   12: CO*+O*--CO2;   13: O*+O*--O2*
    integer(kind=4), parameter :: n_singleP = 6 ! number of single-site-process
    integer(kind=8), parameter :: record_int = 10000  ! interval of record 'ijk_rec.dat'
    integer, DIMENSION(n_singleP) :: singleP = (/2, 3, 4, 5, 11, 12/) ! # of single-site-process
    integer, DIMENSION(n_singleP) :: kArr

    character(len=120) :: filename1, filename2, filename3, char_seq
    character(len=120) :: filename4, filename5

    integer(kind=4) :: i, j, k, co_ads

    real(kind=8) :: r_site, rpoint, rtot, rsite(100000), revent(n_process, 100000), rneis(12, n_process, 100000)
    real(kind=8) :: Ea_r(2, 12, 100000)

    integer(kind=4) :: j_cn ! # of nnsite 
    integer(kind=4) :: ipick, jpick, kpick, nn_i, nn_j, nf_1
    integer(kind=8) :: n_step_tot, n_step_each(n_process)
    real(kind=8) :: ctime, dtime, randum, sum_rsite, sum_revent, sum_rneis
    
    ! get parameter from 'input_tot'
    call read_data 

    ! generate grid and read initial structure
    if(restart .eq. 0) then
        call read_str_new
    else if(restart .eq. 1) then
        call read_str_last
    end if 
    ! calculate cn, cnn, gcn
    call count_cn
    ! define surface atom
    do i = 1, nbulk
        if(effcn(i).gt.0 .and. effcn(i).lt.12) then
            if(restart .eq. 0) then
                site_type(i) = 2
                cov_type(i) = 0
            else if(restart .eq. 1) then
                if(effcn(i).gt.0 .and. effcn(i).lt.12) site_type(i) = 2
            end if
        end if
    end do

    ! calculate rtot, rsite, rneis, revent
    rtot = 0.0
    rsite = 0.0
    revent = 0.0
    rneis = 0.0
    Ea_r = 0.0
    do i = 1, nbulk
        if(site_type(i).eq.2) then
            do k = 2, n_process
                kArr = k
                if(any(kArr .eq. singleP)) then
                    r_site = 0
                    call rijk(r_site, i, k, 0, Ea_r)
                    revent(k, i) = r_site
                else
                    do j = 1, 12
                        r_site = 0
                        call rijk(r_site, i, k, j, Ea_r)
                        rneis(j, k, i) = r_site
                        revent(k, i) = revent(k, i) + r_site
                    end do
                end if
                rsite(i) = rsite(i) + revent(k, i)
            end do
            rtot = rtot + rsite(i)
        end if
    end do

    call init_random_seed()
    n_step_tot = 0
    n_step_each = 0
    ctime = 0.0
    ! new files
    nf_1 = 0
    write(char_seq,'(i2.2)') nf_1
    write(filename1, *) 'atom_str_',trim(char_seq),'.xyz'  !real atoms coordination
    write(filename2, *) 'ijk_rec_',trim(char_seq),'.dat'   !ctime, ipick, jpick, kick of every step
    !(filename3, *) 'Ea_atom_',trim(char_seq),'.dat'       !energy barrier of atom jumping, cn, n_tot_step
    write(filename4, *) 'Ea_COO_',trim(char_seq),'.dat'    !energy barrier of CO oxdaition, cn, n_tot_step
    write(filename5, *) 'step_rec_',trim(char_seq),'.dat'  !steps of every event

    open(10, file=filename1, status='replace', action='write')
    open(11, file=filename2, status='replace', action='write')
    !open(12, file=filename3, status='replace', action='write')
    open(13, file=filename4, status='replace', action='write')
    open(14, file=filename5, status='replace', action='write')
    ! record first structural coordination
    write(10, *) natoms
    write(10, *) 
    do i = 1, nbulk
        if(site_type(i).lt.3) write(10, *) elem, xxx(i), yyy(i), zzz(i), cov_type(i)
    end do
    write(14, "(es22.12, 14I12)") ctime, n_step_tot, n_step_each

    ! KCM cycle
    do while (n_step_tot .le. nLoop)
        !do while (ctime.le. tend)
        ! Circular recording
        if(n_step_tot.gt.0 .and. mod(n_step_tot, 20000000).eq.0) then
            nf_1 = nf_1 + 1
            write(char_seq,'(i2.2)') nf_1
            !write(filename1, *) 'atom_str_',trim(char_seq),'.xyz'  ! real atoms coordination
            write(filename2, *) 'ijk_rec_',trim(char_seq),'.dat'   ! ctime, ipick, jpick, kick of every step
            !write(filename3, *) 'Ea_atom_',trim(char_seq),'.dat'   ! energy barrier of atom jumping, cn, n_tot_step
            write(filename4, *) 'Ea_COO_',trim(char_seq),'.dat'    ! energy barrier of CO oxdaition, cn, n_tot_step
            write(filename5, *) 'step_rec_',trim(char_seq),'.dat'  ! steps of every event

            !open(10, file=filename1, status='new', action='write')
            open(11, file=filename2, status='new', action='write')
            !open(12, file=filename3, status='new', action='write')
            open(13, file=filename4, status='new', action='write')
            open(14, file=filename5, status='new', action='write')
        end if

        ! cal time
        call random_number(randum)
        dtime = -dlog(randum)/rtot
        ctime = ctime + dtime
            
        ! ipick
        !do while (1 .gt. 0)
        call random_number(randum)
        rpoint = rtot * randum
        sum_rsite = 0.0
        do i = 1, nbulk
            if(site_type(i).eq.2) then
                sum_rsite = sum_rsite + rsite(i)
                if(sum_rsite.ge.rpoint) exit
            end if
        end do
        !end do
        ipick = i

        ! event pick
        !do while (1 .gt. 0)
        call random_number(randum)
        rpoint = rsite(ipick) * randum
        sum_revent = 0
        do k = 2, n_process
            sum_revent = sum_revent + revent(k, ipick)
            if(sum_revent.ge.rpoint) exit
        end do
        !end do
        kpick = k
        !write(*, *) kpick
        
        ! j_cn & jpick
        kArr = kpick
        if(all(kArr .ne. singleP)) then
            !do while(1.gt.0)
            call random_number(randum)
            rpoint = revent(kpick, ipick) * randum
            sum_rneis = 0
            do j = 1, 12
                sum_rneis = sum_rneis + rneis(j, kpick, ipick)
                if(sum_rneis.ge.rpoint) exit
            end do
            !end do
            j_cn = j
            jpick = nnsite(j_cn, ipick)
        end if

        n_step_tot = n_step_tot + 1

        if(kpick.eq.10) then
            write(13, "(2f15.6, 2I5, 2f11.4, 4x, I0)") &
            Ea_r(1, j_cn, ipick), Ea_r(2, j_cn, ipick), effcn(ipick), effcn(jpick), gcn(ipick), gcn(jpick), n_step_tot
        end if

        SELECT CASE(kpick)
            CASE(2) ! CO adsorption
                cov_type(ipick) = 1
                n_step_each(2) = n_step_each(2) + 1

            CASE(3) ! CO desorption
                cov_type(ipick) = 0
                n_step_each(3) = n_step_each(3) + 1

            CASE(4) ! O2 adsorption
                cov_type(ipick) = 2
                n_step_each(4) = n_step_each(4) + 1

            CASE(5) ! O2 desorption
                cov_type(ipick) = 0
                n_step_each(5) = n_step_each(5) + 1

            CASE(6) ! CO diffusion
                cov_type(ipick) = 0
                cov_type(jpick) = 1
                n_step_each(6) = n_step_each(6) + 1

            CASE(7) ! O2 diffusion
                cov_type(ipick) = 0
                cov_type(jpick) = 2
                n_step_each(7) = n_step_each(7) + 1

            CASE(8) ! O diffusion
                cov_type(ipick) = 0
                cov_type(jpick) = 3
                n_step_each(8) = n_step_each(8) + 1

            CASE(9) ! OCOO diffusion
                co_ads = co_ads_site(ipick)
                cov_type(jpick) = cov_type(ipick)
                co_ads_site(co_ads) = jpick
                co_ads_site(jpick) = co_ads
                n_step_each(9) = n_step_each(9) + 1

            CASE(10) ! CO+O2 -- OCOO -- CO2 + O2
                cov_type(ipick) = 4
                cov_type(jpick) = -4
                co_ads_site(ipick) = jpick
                co_ads_site(jpick) = ipick
                n_step_each(10) = n_step_each(10) + 1
                
                cov_type(ipick) = 0
                cov_type(jpick) = 3
                co_ads_site(ipick) = 0
                co_ads_site(jpick) = 0
                n_step_each(11) = n_step_each(11) + 1
                ! n_step_tot = n_step_tot + 1
                write(11, *) ctime, n_step_tot, kpick, kpick+1

            CASE(11) ! OCOO--CO2 + O2
                co_ads = co_ads_site(ipick)
                if(cov_type(ipick) .eq. 4) then
                    cov_type(ipick) = 0
                    cov_type(co_ads) = 3
                else
                    cov_type(ipick) = 3
                    cov_type(co_ads) = 0
                end if
                co_ads_site(ipick) = 0
                co_ads_site(co_ads) = 0
                n_step_each(11) = n_step_each(11) + 1
                write(11, *) ctime, n_step_tot, kpick

            CASE(12) ! CO+O reaction
                cov_type(ipick) = 0
                n_step_each(12) = n_step_each(12) + 1
                write(11, *) ctime, n_step_tot, kpick
        
            CASE(13) ! O+O reaction
                cov_type(ipick) = 2
                cov_type(jpick) = 0
                write(11, *) ctime, n_step_tot, kpick
                n_step_each(13) = n_step_each(13) + 1
        END SELECT
            
        if(mod(n_step_tot, record_int).eq.0) then
            write(14, "(es22.12, 14I12)") ctime, n_step_tot, n_step_each
        end if

        ! initial and cal new relative r
        ! update rsite, revent, and rneis of all atoms in nnnsite of ipick
        do i = 1, 54
            nn_i = nnnsite(i, ipick)
            if(nn_i .gt. 0) then
                rtot = rtot - rsite(nn_i)
                rsite(nn_i) = 0
                do k = 2, n_process
                    revent(k, nn_i) = 0
                    kArr = k
                    if(all(kArr .ne. singleP)) then
                        do j = 1, 12
                            rneis(j, k, nn_i) = 0
                        end do
                    end if
                end do
                if(site_type(nn_i).eq.2) then
                    do k = 2, n_process
                        kArr = k
                        if(any(kArr .eq. singleP)) then
                            r_site = 0
                            call rijk(r_site, nn_i, k, 0, Ea_r)
                            
                            revent(k, nn_i) = revent(k, nn_i) + r_site
                        else
                            do j = 1, 12
                                r_site = 0
                                call rijk(r_site, nn_i, k, j, Ea_r)
                                
                                rneis(j, k, nn_i) = r_site
                                revent(k, nn_i) = revent(k, nn_i) + r_site
                            end do
                        end if
                        rsite(nn_i) = rsite(nn_i) + revent(k, nn_i)                        
                    end do
                    rtot = rtot + rsite(nn_i)
                end if
            end if
        end do

        kArr = kpick
        if(any(kArr .eq. singleP)) then
        ! update rsite, revent, and rneis of ipick
            rtot = rtot - rsite(ipick)
            rsite(ipick) = 0
            do k = 2, n_process
                ! initial
                revent(k, ipick) = 0
                kArr = k
                if(all(kArr .ne. singleP)) then
                    do j = 1, 12
                        rneis(j, k, ipick) = 0
                    end do
                end if
                ! cal new r
                if(site_type(ipick).eq.2) then
                    kArr = k
                    if(any(kArr .eq. singleP))then
                        r_site = 0
                        call rijk(r_site, ipick, k, 0, Ea_r)

                        revent(k, ipick) = revent(k, ipick) + r_site
                    else
                        do j = 1, 12
                            r_site = 0
                            call rijk(r_site, ipick, k, j, Ea_r)
                            
                            rneis(j, k, ipick) = r_site
                            revent(k, ipick) = revent(k, ipick) + r_site
                        end do
                    end if
                    rsite(ipick) = rsite(ipick) + revent(k, ipick)                    
                end if
            end do
            rtot = rtot + rsite(ipick)
        else
        ! update rsite, revent, and rneis of all atoms in nnnsite of jpick
            do i = 1, 54
                nn_i = nnnsite(i, jpick)
                if(nn_i.gt.0) then
                    ! initial
                    rtot = rtot - rsite(nn_i)
                    rsite(nn_i) = 0
                    do k = 2, n_process
                        revent(k, nn_i) = 0
                        kArr = k
                        if(any(kArr .eq. singleP)) then
                            do j = 1, 12
                                rneis(j, k, nn_i) = 0
                            end do
                        end if
                    end do
                    ! cal new r
                    if(site_type(nn_i).eq.2) then
                        do k = 2, n_process
                            kArr = k
                            if(any(kArr .eq. singleP))then
                                r_site = 0
                                call rijk(r_site, nn_i, k, 0, Ea_r)

                                revent(k, nn_i) = revent(k, nn_i) + r_site
                            else
                                do j = 1, 12
                                    r_site = 0
                                    call rijk(r_site, nn_i, k, j, Ea_r)
                            
                                    rneis(j, k, nn_i) = r_site
                                    revent(k, nn_i) = revent(k, nn_i) + r_site
                                end do
                            end if
                            rsite(nn_i) = rsite(nn_i) + revent(k, nn_i)                           
                        end do
                        rtot = rtot + rsite(nn_i)
                    end if
                end if
            end do
        end if
    end do
    
    ! record 14
    write(14, "(es22.12, 14I12)") ctime, n_step_tot, n_step_each

    close(10)
    close(11)
    !close(12)
    close(13)
    close(14)

    write(*, *) "====== END ======="
    ! write last_one (xyz and cov_type)
    open(11, file = 'last_one.xyz')
        write(11, *) natoms
        write(11, *) 
        do i = 1,nbulk
            if(site_type(i).lt.3) write(11, *) elem, xxx(i), yyy(i), zzz(i), cov_type(i), co_ads_site(i)
        end do
    close(11)
    
    !pause
end

subroutine rijk(r_site, ipick, kpick, j_cn, Ea_r)
    use outer_data
    use struct_data
    
    implicit none

    real(kind=8), intent(OUT) :: r_site
    integer(kind=4), intent(IN):: ipick, kpick, j_cn
    real(kind=8), intent(INOUT) :: Ea_r(2, 12, 100000)

    integer(kind=4) :: ni_ads, nj_ads
    integer(kind=4) :: i, j, k, co_ads
    integer(kind=4) :: jpick
    real(kind=8) :: rco_ads, ro2_ads, ro_ads
    real(kind=8) :: gcn_ipick, gcn_jpick

    real(kind=8), parameter :: pi = 3.141592654
    real(kind=8), parameter :: kb = 8.6173324D-05
    real(kind=8), parameter :: h = 4.1356676D-15
    real(kind=8), parameter :: eV2J = 1.60217662D-19  ! 1 eV = 1.60217662D-19 J
    real(kind=8), parameter :: Na = 6.0221409D23
    real(kind=8), parameter :: Asite = (10D-10)**2
    real(kind=8), parameter :: p0 = 100000  ! p0 = 100kPa
    
    real(kind=8) :: s0CO, s0O2, s0O, mCO, mO, Sco_0, So2_0
    real(kind=8) :: Eads_i, Eads_j, dE, Ea, Ea2
    real(kind=8) :: dS_CO, dS_O2, r_K_eq
    real(kind=4) :: pCO, pO2

    mCO = 28.01D-3/Na   ! kg/atom
    mO  = 15.999D-3/Na

    Sco_0 = 85.142*(Temp**0.14709)/(Na*eV2J)
    So2_0 = 89.655*(Temp**0.14489)/(Na*eV2J)
    
    pCO = ppCO * p_tot
    pO2 = ppO2 * p_tot

    ni_ads = cov_type(ipick)
    gcn_ipick = gcn(ipick)
    if(j_cn .ne. 0) then
        jpick = nnsite(j_cn, ipick)
        nj_ads = cov_type(jpick)
        gcn_jpick = gcn(jpick)
    end if
    
    SELECT CASE (kpick)
        CASE(2) ! CO adsorption
            if(ni_ads.eq.0) then
                s0CO = s0_CO_facet
                if(gcn_ipick.ge.0 .and. gcn_ipick.lt.5.33) s0CO = s0_CO_edge
                if(effcn(ipick).ge.10) s0CO = 0.0
                r_site = (s0CO * pCO * Asite)/sqrt(2. * pi * mCO * kb * eV2J * Temp)
            end if
    
        CASE(3) ! CO desorption
            if(ni_ads.eq.1) then
                s0CO = s0_CO_facet        
                ! CO ads
                if(gcn_ipick.ge.0 .and. gcn_ipick.lt.5.33) s0CO = s0_CO_edge
                if(effcn(ipick).ge.10) s0CO = 0.0
                rco_ads = (s0CO * pCO * Asite)/sqrt(2. * pi * mCO * kb * eV2J * Temp)
                ! CO equilibrium
                dS_CO = 0 - (Sco_0 - kb*log(pCO/p0))
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                r_K_eq = exp(-(Eads_i - Temp*dS_CO)/(kb*Temp))
                ! r of CO des
                r_site = rco_ads/(pCO*r_K_eq)
            end if
            
        CASE(4) ! O2 adsorption
            if (ni_ads.eq.0) then
                s0O2 = s0_O2_facet
                if(gcn_ipick.ge.0 .and. gcn_ipick.lt.5.33) s0O2 = s0_O2_edge
                if(effcn(ipick).ge.10) s0O2 = 0.0
                r_site = (s0O2 * pO2 * Asite)/sqrt(2. * pi * mO * kb * eV2J * Temp)
            end if

        CASE(5) ! O2 desorption
            if(ni_ads.eq.2) then
                ! O2 ads
                s0O2 = s0_O2_facet
                if(gcn_ipick.ge.0 .and. gcn_ipick.lt.5.33) s0O2 = s0_O2_edge
                if(effcn(ipick).ge.10) s0O2 = 0.0
                ro2_ads = (s0O2 * pO2 * Asite)/sqrt(2. * pi * mO * kb * eV2J * Temp)
                ! O2 equilibrium
                dS_O2 = 0 - (So2_0 - kb*log(pO2/p0))
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                r_K_eq = exp(-(Eads_i - Temp*dS_O2)/(kb*Temp))
                ! O2 des
                r_site = ro2_ads/(pO2*r_K_eq)
            end if
            
        CASE(6) ! CO diffusion
            if(ni_ads.eq.1 .and. nj_ads.eq.0 .and. effcn(jpick).lt.10) then
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                cov_type(ipick) = 0
                cov_type(jpick) = 1
                Eads_j = 0.0
                call Eads_site(Eads_j, jpick, 1)
                dE = Eads_j - Eads_i
                if(dE.le.0) dE = 0.0
                Ea = dE + EdiffCO
                r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
                ! recovery
                cov_type(ipick) = 1
                cov_type(jpick) = 0
            end if

        CASE(7) ! O2 diffusion
            if(ni_ads.eq.2 .and. nj_ads.eq.0) then
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                cov_type(ipick) = 0
                cov_type(jpick) = 2
                Eads_j = 0.0
                call Eads_site(Eads_j, jpick, 1)
                dE = Eads_j - Eads_i
                if(dE.le.0) dE = 0.0
                Ea = dE + EdiffO
                r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
                ! recovery
                cov_type(ipick) = 2
                cov_type(jpick) = 0
            end if

        CASE(8) ! O diffusion
            if(ni_ads.eq.3 .and. nj_ads.eq.0) then
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                cov_type(ipick) = 0
                cov_type(jpick) = 3
                Eads_j = 0.0
                call Eads_site(Eads_j, jpick, 1)
                dE = Eads_j - Eads_i
                if(dE.le.0) dE = 0.0
                Ea = dE + EdiffO2
                r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
                ! recovery
                cov_type(ipick) = 3
                cov_type(jpick) = 0
            end if
        
        CASE(9) ! OCOO diffusion
            if(abs(ni_ads).eq.4 .and. nj_ads.eq.0) then
                r_site = 0
            end if

        CASE(10) ! CO + O2 -- OCOO -- CO2 + O
            if(ni_ads.eq.1 .and. nj_ads.eq.2) then
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                Eads_j = 0.0
                call Eads_site(Eads_j, jpick, 1)
                Ea = BEP1_a*(Eads_i + Eads_j) + BEP1_b
                Ea2 = BEP2_a*(Eads_i + Eads_j) + BEP2_b
                if(Ea .le. 0) Ea = 0.0
                if(Ea2 .le. 0) Ea2 = 0.0
                Ea_r(1, j_cn, ipick) = Ea
                Ea_r(2, j_cn, ipick) = Ea2
                Ea = Ea + Ea2
                r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
            end if
        
        CASE(11) ! OCOO - CO2 + O reaction
            if(ni_ads .eq. 4) then 
                r_site = 0
            end if
            
        CASE(12) ! CO + O - CO2 reaction
            if(ni_ads .eq. 3) then
                Ea = Ea3
                r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
            end if
        
        CASE(13) ! O + O - O2 desorption
            if(ni_ads.eq.3 .and. nj_ads.eq.3) then
                r_site = 0
            end if            
    END SELECT
end subroutine

subroutine Eads_site(Eads_s, site, num)
    use outer_data
    use struct_data
    implicit none

    integer(kind=4), intent(IN) :: site, num
    real(kind=8), intent(OUT) :: Eads_s

    integer(kind=4) :: i, j
    integer(kind=4) :: nn_s, nn_j
    integer(kind=4) :: n_CO, n_O, n_O2, effcn_ni, effcn_nj

    n_CO = 0
    n_O = 0
    n_O2 = 0

    if(num.eq.1) then
        effcn_ni = effcn(site)
        if(effcn_ni.gt.6) then
            do i = 1, 12
                nn_s = nnsite(i, site)
                if(effcn(nn_s).gt.6) then
                    if(cov_type(nn_s).eq.1) n_CO = n_CO + 1
                    if(cov_type(nn_s).eq.2) n_O2 = n_O2 + 1
                    if(cov_type(nn_s).eq.3) n_O = n_O + 1
                    !if(cov_type(nn_s).eq.4) n_CO = n_CO + 1
                    !if(cov_type(nn_s).eq.-4) n_O2 = n_O2 + 1
                end if
            end do
        end if
    else if(num.eq.2) then
        effcn_ni = 0
        do i = 1, 12
            nn_s = nnsite(i, site)
            if(site_type(nn_s).lt.3) effcn_ni = effcn_ni + 1
        end do
        if(effcn_ni.gt.6) then
            do i = 1, 12
                nn_s = nnsite(i, site)
                if(cov_type(nn_s).eq.1 .or. cov_type(nn_s).eq.2 .or. cov_type(nn_s).eq.3) then
                    effcn_nj = 0
                    do j = 1, 12
                        nn_j = nnsite(j, nn_s)
                        if(site_type(nn_j).lt.3) effcn_nj = effcn_nj + 1
                    end do
                    if(effcn_nj.gt.6) then
                        if(cov_type(nn_s).eq.1) n_CO = n_CO + 1
                        if(cov_type(nn_s).eq.2) n_O2 = n_O2 + 1
                        if(cov_type(nn_s).eq.3) n_O = n_O + 1
                        !if(cov_type(nn_s).eq.4) n_CO = n_CO + 1
                        !if(cov_type(nn_s).eq.-4) n_O2 = n_O2 + 1                        
                    end if
                end if
            end do
        end if
    end if
    if(cov_type(site).eq.1) Eads_s = (Eads_CO_a*gcn(site) + Eads_CO_b) - (n_CO*Eco_co + n_O2*Eco_o2 + n_O*Eco_o)
    if(cov_type(site).eq.2) Eads_s = (Eads_O2_a*gcn(site) + Eads_O2_b) - (n_CO*Eco_o2 + n_O2*Eo2_o2 + n_O*Eo2_o)
    if(cov_type(site).eq.3) Eads_s = (Eads_O_a*gcn(site) + Eads_O_b) - (n_CO*Eco_o + n_O2*Eo2_o + n_O*Eo_o)
end subroutine
