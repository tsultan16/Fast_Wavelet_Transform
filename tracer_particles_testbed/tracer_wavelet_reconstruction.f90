PROGRAM tracer_reconstruction


!USE MPI
USE OMP_LIB

IMPLICIT NONE



TYPE :: patch_scalar

    INTEGER :: nparticles = 0
    REAL(4), ALLOCATABLE :: patch_array(:,:,:,:)
    REAL(4), ALLOCATABLE :: tracer_wavelets(:,:,:,:)

END TYPE patch_scalar

INTEGER, PARAMETER :: dump_num = 5
INTEGER, PARAMETER :: nx = 32, ny = 32, nz = 32
INTEGER, PARAMETER :: nranks_x = 3, nranks_y = 3, nranks_z = 3
INTEGER, PARAMETER :: nwavelets_per_bin = 5 ! # of wavelets per tracer per bin
INTEGER, PARAMETER :: bin_neighborhood = 0
INTEGER, PARAMETER :: nscalars = 3   ! # of scalars

INTEGER, PARAMETER :: daub_num = 20 ! set this to either 4, 12 or 20
REAL*4,  PARAMETER :: TWOPI  = 8.d0*ATAN(1.d0)
LOGICAL, PARAMETER :: output_basis = .FALSE.

CHARACTER(LEN=200), PARAMETER :: output_filepath = '/data/uchu/tanzid/Wombat_Sprint_5/WAVELETS/wombat/build/Snapshots/TRACER/'


INTEGER :: i, j, k, nxtot, nytot, nztot
TYPE(patch_scalar), ALLOCATABLE :: dat_rec(:,:,:)
REAL(8) :: t1, t2, p(3)
INTEGER ::  nwavelets_per_scalar, nwavelets_max
    
! max number of wavelets per tracer
nwavelets_per_scalar =  nwavelets_per_bin * (1+2*bin_neighborhood)**3
nwavelets_max = nscalars * nwavelets_per_scalar

! world grid size
nxtot = nranks_x*nx
nytot = nranks_y*ny
nztot = nranks_z*nz


! allocate memory for patch arrays
CALL create_patch_arrays()


! compute the wavelet transform
t1 = OMP_GET_WTIME() ! MPI_WTIME()


CALL reconstruct()

! save to file
CALL file_output()



t2 = OMP_GET_WTIME() ! MPI_WTIME()




CALL destroy_patch_arrays()

PRINT*,''
PRINT*,'WT time (sec) = ',t2-t1
PRINT*,''

PRINT*,'Done.'

CONTAINS


SUBROUTINE create_patch_arrays()

    INTEGER :: i, j, k, nparticles

    PRINT*,''
    PRINT*,'Constructing patch arrays and reading tracer data from files...'


    ALLOCATE(dat_rec(0:nranks_x-1,0:nranks_y-1,0:nranks_z-1))
    
    DO k = 0, nranks_z-1
        DO j = 0, nranks_y-1
            DO i = 0, nranks_x-1
    
                CALL read_metadata(i, j, k, nparticles)             
                dat_rec(i,j,k)%nparticles = nparticles

                ALLOCATE(dat_rec(i,j,k)%patch_array(1:nx,1:ny,1:nz,1:nscalars))
                dat_rec(i,j,k)%patch_array = 0.D0
    
                ! allocate memory for tracer data
                ALLOCATE(dat_rec(i,j,k)%tracer_wavelets(1:nparticles,nwavelets_max,nscalars,2))
                dat_rec(i,j,k)%tracer_wavelets = 0.0
    
                ! now read tracer data                 
                CALL read_tracers(i, j, k)
    
            END DO
        END DO
    END DO
   

END SUBROUTINE create_patch_arrays


SUBROUTINE destroy_patch_arrays()

    INTEGER :: i, j, k
    
    DO k = 0, nranks_z-1
        DO j = 0, nranks_y-1
            DO i = 0, nranks_x-1
    
                DEALLOCATE(dat_rec(i,j,k)%patch_array)
                DEALLOCATE(dat_rec(i,j,k)%tracer_wavelets)
    
            END DO
        END DO
    END DO
    
    DEALLOCATE(dat_rec)
   


END SUBROUTINE destroy_patch_arrays


SUBROUTINE collapse_index(ii, i, j, k)

    INTEGER, INTENT(IN) :: ii
    INTEGER, INTENT(OUT) :: i, j, k

    
    i = MOD(ii-1, nranks_x)
    j = MOD(ii-1, nranks_x*nranks_y) / nranks_x 
    k = (ii-1) / (nranks_x*nranks_y) 
              
END SUBROUTINE collapse_index



! this subroutine assumes non-asynchornous WOMBAT I/O, i.e. each MPI ranks output's it's own file
SUBROUTINE reconstruct()

    CHARACTER(LEN=100) :: filename, buffer, label, value, tmp
    CHARACTER(LEN=6) :: uniti, rx, ry, rz
    INTEGER :: ii, i, j, k, ios, line, pos, nranks_tot
    REAL :: time, dt, cmplt
    
    cmplt = 0.D0    
    
    nranks_tot = nranks_x * nranks_y * nranks_z
    
    ! loop over MPI ranks
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k)
    DO ii = 1, nranks_tot
    
        CALL collapse_index(ii,i,j,k)

        
        PRINT*,''
        PRINT*,'Reconstructing for rank: ',i,j,k
  
        CALL tracer_wavelet_reconstruction(i, j, k)     
   
        !$OMP ATOMIC
        cmplt = cmplt + 1
   
        PRINT*,'Done reconstructing for rank: ',i,j,k
        PRINT*,'Percent complete = ',100.0*REAL(cmplt)/REAL(nranks_tot)
                
    END DO
   !$OMP END PARALLEL DO


END SUBROUTINE reconstruct


SUBROUTINE read_metadata(ri, rj, rk, nparticles)

    INTEGER, INTENT(IN) :: ri, rj, rk
    INTEGER, INTENT(INOUT) :: nparticles
    CHARACTER(LEN=100) :: filename, buffer, label, value, tmp
    CHARACTER(LEN=6) :: uniti, rx, ry, rz
    INTEGER :: i, j, k, ios, line, pos
    REAL :: time, dt, origin_zone(3)
    
    
    IF(dump_num<10) THEN
        WRITE(uniti,'(I1.1)') dump_num
    ELSE IF(dump_num>=10 .and. dump_num<100) THEN
        WRITE(uniti,'(I2.2)') dump_num
    ELSE IF(dump_num>=100 .and. dump_num<1000) THEN
        WRITE (uniti,'(I3.3)') dump_num
    ELSE IF(dump_num>=1000 .and. dump_num<10000) THEN
        WRITE (uniti,'(I4.3)') dump_num
    ELSE IF(dump_num>=10000 .and. dump_num<100000) THEN
        WRITE (uniti,'(I5.3)') dump_num  
    END IF

    WRITE(rx,'(I1.1)') ri
    WRITE(ry,'(I1.1)') rj
    WRITE(rz,'(I1.1)') rk
    
    filename = TRIM(output_filepath)//TRIM('tracer_meta_')//TRIM(rx)//TRIM(ry)// & 
       TRIM(rz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.txt')        
    
    
    OPEN(UNIT=10, FILE = filename)
    
    ios = 0
    line = 0

    DO WHILE(ios .EQ. 0)
        READ(10, FMT='(A)', IOSTAT = ios) buffer
        !PRINT*,buffer
        
        IF(ios .EQ. 0) THEN
            line = line + 1
            
            ! if the first character is '#', skip this line
            IF(buffer(1:1) .EQ. '#')THEN
                !PRINT*,'Skipping comment at line ',line                        
                CYCLE                        
            END IF
            
            ! find the first instance of a whitespace
            pos = SCAN(buffer, ' ')
            
            label  = ADJUSTL(TRIM(buffer(1:pos)))
            value  = ADJUSTL(TRIM(buffer(pos+1:)))
            
            !PRINT*,'label =',label
            !PRINT*,'value =',value
            
            SELECT CASE(label)
            
                CASE('time')
            
                    READ(value,*,iostat=ios) time
             
                    !PRINT*,'Time value read from file: ',time
            
                CASE('dt')
            
                    READ(value,*,iostat=ios) dt
             
                    !PRINT*,'dt value read from file: ',dt
            
                CASE('origin_zone')

                    tmp = value
                    pos = SCAN(tmp, ' ')
                    value = ADJUSTL(TRIM(tmp(1:pos)))
                    tmp = ADJUSTL(TRIM(tmp(pos+1:)))                                
                    READ(value,*,iostat=ios) origin_zone(1)

                    pos = SCAN(tmp, ' ')
                    value = ADJUSTL(TRIM(tmp(1:pos)))
                    tmp = ADJUSTL(TRIM(tmp(pos+1:)))
                    READ(value,*,iostat=ios) origin_zone(2)

                    pos = SCAN(tmp, ' ')
                    value = ADJUSTL(TRIM(tmp(1:pos)))
                    tmp = ADJUSTL(TRIM(tmp(pos+1:)))
                    READ(value,*,iostat=ios) origin_zone(3)


                    !PRINT*,'origin_zone value read from file: ',origin_zone

                CASE('nparticles')

                    READ(value,*,iostat=ios) nparticles
             
                    !PRINT*,'nparticles value read from file: ',nparticles

                CASE DEFAULT
            
                    !PRINT*,'Skipping invalid label at line ',line
            
            END SELECT
            
            
        END IF
        
        
    END DO

    CLOSE(UNIT=10)
   
END SUBROUTINE read_metadata



SUBROUTINE read_tracers(ri, rj, rk)

    INTEGER, INTENT(IN) :: ri, rj, rk
    CHARACTER(LEN=100) :: filename, buffer, label, value, tmp
    CHARACTER(LEN=6) :: uniti, rx, ry, rz
    INTEGER :: i, j, k, isc, item_bytes, byte_offset,  nwavelets_per_scalar, nwavelets_max
    REAL :: time, dt, pid, init_rank, x(3)    
       
    
    ! max number of wavelets per tracer
    nwavelets_per_scalar =  nwavelets_per_bin * (1+2*bin_neighborhood)**3
    nwavelets_max = nscalars * nwavelets_per_scalar

    ! byte size per tracer data item
    item_bytes = 4
        
    
    IF(dump_num<10) THEN
        WRITE(uniti,'(I1.1)') dump_num
    ELSE IF(dump_num>=10 .and. dump_num<100) THEN
        WRITE(uniti,'(I2.2)') dump_num
    ELSE IF(dump_num>=100 .and. dump_num<1000) THEN
        WRITE (uniti,'(I3.3)') dump_num
    ELSE IF(dump_num>=1000 .and. dump_num<10000) THEN
        WRITE (uniti,'(I4.3)') dump_num
    ELSE IF(dump_num>=10000 .and. dump_num<100000) THEN
        WRITE (uniti,'(I5.3)') dump_num  
    END IF

    WRITE(rx,'(I1.1)') ri
    WRITE(ry,'(I1.1)') rj
    WRITE(rz,'(I1.1)') rk
    
    filename = TRIM(output_filepath)//TRIM('tracer_')//TRIM(rx)//TRIM(ry)// & 
       TRIM(rz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')        
    
   
    OPEN(UNIT=10, FILE = filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM')
   
    byte_offset = 1
    
    DO i = 1, dat_rec(ri,rj,rk)%nparticles  
    
        !PRINT*,'PARTICLE# ',i
    
        READ(10, POS = byte_offset) pid 
        READ(10, POS = byte_offset + item_bytes) init_rank 
        READ(10, POS = byte_offset + 2*item_bytes) x 
        
        byte_offset = byte_offset + 5*item_bytes

        !PRINT*,'Particle id, init_rank = ',pid, init_rank
        !PRINT*,'Particle position = ',x
    
        IF(INT(x(1)) .LT. 0 .OR. INT(x(1)) .GT. nx+1 .OR. INT(x(2)) .LT. 0 .OR. INT(x(2)) .GT. ny+1 .OR. INT(x(3)) .LT. 0 .OR. INT(x(3)) .GT. nz+1) THEN
            PRINT*,'ERROR!! Invalid particle coordinates: ',x
            STOP
        END IF
    
        DO isc = 1, nscalars
   
            !PRINT*,'Scalar# ',isc
      
            DO j = 1, nwavelets_per_scalar
            
                DO k = 1, 2
                    READ(10, POS = byte_offset) dat_rec(ri,rj,rk)%tracer_wavelets(i,j,isc,k)
                    byte_offset = byte_offset + item_bytes            
                END DO  


                
            END DO
            
        END DO
        
    END DO

    CLOSE(10)

   
   

END SUBROUTINE read_tracers


! partial reconstruction using wavelets picked up by the tracers
SUBROUTINE tracer_wavelet_reconstruction(ri, rj, rk)

    INTEGER, INTENT(IN) :: ri, rj, rk
    REAL(4) :: fbasis(nx,ny,nz), amp
    INTEGER :: i, j, k, ii, ix, iy, iz, isc, nn

    nn = 0

    ! loop over scalars
    DO isc = 1, nscalars
        
        ! loop over particles
        DO i = 1, dat_rec(ri,rj,rk)%nparticles
        
            !PRINT*,'Particle# ',i
            
            ! loop over wavelets
            DO j = 1, nwavelets_per_scalar
                
                IF(INT(dat_rec(ri,rj,rk)%tracer_wavelets(i,j,isc,1)) .LE. 0) EXIT
        
                ! flattened wavelet index
                ii = dat_rec(ri,rj,rk)%tracer_wavelets(i,j,isc,1)
                
                
                ! collapse flattened wavelet index
                ix = 1 + MOD(ii-1,nx)
                iy = 1 + MOD(ii-1,nx*ny)/nx 
                iz = 1 + (ii-1) / (nx*ny) 
                               
                                                           
                IF(ix .LT. 1 .OR. ix .GT. nx .OR. iy .LT. 1 .OR. iy .GT. ny .OR. iz .LT. 1 .OR. iz .GT. nz) THEN
                    PRINT*,'ERROR!! Invalid wavelet index: ',ix,iy,iz
                    STOP
                END IF       
                                               
                ! compute basis wavelet
                CALL compute_basis_wavelet(ix,iy,iz,fbasis)
       
                ! add wavelet contribution to the reconstructed function
                dat_rec(ri,rj,rk)%patch_array(:,:,:,isc) = dat_rec(ri,rj,rk)%patch_array(:,:,:,isc) + dat_rec(ri,rj,rk)%tracer_wavelets(i,j,isc,2)*fbasis(:,:,:) 
        
                nn = nn + 1
        
            END DO
        
        
        END DO
     
    END DO
 
    PRINT*,'Total number of wavelets used in the reconstruction  =',nn
 
 
END SUBROUTINE tracer_wavelet_reconstruction




! computes basis wavelet
SUBROUTINE compute_basis_wavelet(i,j,k,fx)

    INTEGER, INTENT(IN) :: i,j,k
    REAL, INTENT(INOUT) :: fx(nx,nx,nx)
    REAL :: fw(nx,nx,nx)
    INTEGER :: ix, iy, iz
    
    fw = 0.0
    fw(i,j,k) = 1.0

    ! inverse dwt
    CALL fwt_3d(fw,fx,-1)
    
    !fx = fw
    !CALL wtn(fx,3,-1)


END SUBROUTINE compute_basis_wavelet



! fully separable 2d wavelet transform
! Note: The output of a fully separable fwt has a diffferent layout than the usual Mallat approach. Be warned.. 
SUBROUTINE fwt_3d(fx, fw, sgn)

    REAL*4, INTENT(IN) :: fx(1:nx,1:nx,1:nx) 
    REAL*4, INTENT(INOUT) :: fw(1:nx,1:nx,1:nx) 
    INTEGER, INTENT(IN) :: sgn ! 1: forward dwt, -1: inverse dwt
    REAL*4 :: buffer(1:nx)
    INTEGER :: kx, ky, kz, i, ix, iy, iz
  
    
    !PRINT*,' Doing x pass..'    
        
    ! DWT in x-direction

    DO iz = 1, nx
        DO iy = 1, nx
                                
            ! copy strips-along x into 1d buffer    
            DO ix = 1, nx
                buffer(ix) = fx(ix,iy,iz)
            END DO
        
            ! perform 1D dwt 
            CALL fwt_1d(buffer,nx,sgn)
        
            ! copy back
            DO ix = 1, nx
                fw(ix,iy,iz) = buffer(ix)  
            END DO

        END DO       
    END DO       
    
    !PRINT*,' Doing y pass..'    


    ! DWT in y-direction
    DO iz = 1, nx 
        DO ix = 1, nx 
                
            ! copy strips-along y into 1d buffer    
            DO iy = 1, nx
                buffer(iy) = fw(ix,iy,iz)
            END DO
        
            ! perform 1D dwt 
            CALL fwt_1d(buffer,nx,sgn)
                    
            ! copy back
            DO iy = 1, nx
                fw(ix,iy,iz) = buffer(iy)  
            END DO            
                
        END DO
    END DO

    ! DWT in z-direction

    DO iy = 1, nx 
        DO ix = 1, nx 
                
            ! copy strips-along y into 1d buffer    
            DO iz = 1, nx
                buffer(iz) = fw(ix,iy,iz)
            END DO
        
            ! perform 1D dwt 
            CALL fwt_1d(buffer,nx,sgn)
                    
            ! copy back
            DO iz = 1, nx
                fw(ix,iy,iz) = buffer(iz)  
            END DO            
                
        END DO
    END DO

    
END SUBROUTINE fwt_3d



! This subroutine computes the (inverse) wavelet transform of the input data vector of length 'n' for sgn = (-1) 1 
! n has to be power of 2
SUBROUTINE fwt_1d(a,n,sgn)

    REAL*4, INTENT(INOUT) :: a(:)    ! input data vector
    INTEGER, INTENT(IN) :: sgn, n 
    INTEGER :: nn
    
    IF(n .LT. 4) RETURN
      
    ! compute the wavelet transform  
    IF(sgn .GE. 0) THEN
        
        nn = n ! start at largest hierarchy
        
        DO WHILE(nn .GE. 4) 
            CALL pwt(a,nn,daub_num,sgn)
            nn = nn/2
        END DO
        
    ELSE ! inverse transform

        nn = 4 ! start at smallest

        DO WHILE(nn .LE. n)
            CALL pwt(a,nn,daub_num,sgn)
            nn = nn*2        
        END DO
        
    END IF    


END SUBROUTINE fwt_1d


! initilaization for wavelet filter co-efficients
SUBROUTINE wtset(ncof, ioff, joff, cc, cr)

    INTEGER, INTENT(IN) :: ncof
    INTEGER, INTENT(INOUT) :: ioff, joff
    REAL*4, INTENT(INOUT) :: cc(:), cr(:)
    REAL*4 :: c4(4), c12(12), c20(20)
    INTEGER :: sig, k
    
    ! DAUB4 filter co-efficients
    c4 = (/ 0.4829629131445341, 0.8365163037378079, 0.2241438680420134,-0.1294095225512604 /)    
        
    ! DAUB12 filter co-efficients
    c12 = (/.111540743350, .494623890398, .751133908021, .315250351709,-.226264693965,-.129766867567, &
            .097501605587, .027522865530,-.031582039318, .000553842201, .004777257511,-.001077301085/)
    
    ! DAUB20 filter co-efficients
    c20 = (/ .026670057901, .188176800078, .527201188932, .688459039454, .281172343661,-.249846424327, &
            -.195946274377, .127369340336, .093057364604, -.071394147166,-.029457536822, .033212674059, &
             .003606553567,-.010733175483, .001395351747, .001992405295,-.000685856695,-.000116466855, &
             .000093588670,-.000013264203 /)

    sig = -1
    
    DO k = 1, ncof
    
        IF(ncof .EQ. 4)THEN
            cc(k) = c4(k)
        ELSE IF(ncof .EQ. 12) THEN
            cc(k) = c12(k)
        ELSE IF(ncof .EQ. 20) THEN
            cc(k) = c20(k)
        ELSE
            PRINT*,'Unimplemented value ncof in pwtset. Need to choose from 4, 12 and 20.'
            STOP
        END IF
        
        cr(ncof+1-k) = sig * cc(k)
        sig = -sig

    END DO
    
    joff = -ncof/2 ! center for wavelet function support
    ioff = -ncof/2


END SUBROUTINE wtset


! partial wavelet transform (i.e. multiplication by the wavelet coefficient matrix followed by a permutation that rearrages the resulting vector
! so that all smooth components occupy the foirst half followed by the detail coefficients)
SUBROUTINE pwt(a,n,filter,sgn)

    REAL*4, INTENT(INOUT) :: a(:)    ! input data vector
    INTEGER, INTENT(IN) :: sgn, n, filter
    INTEGER, PARAMETER :: nmax = 1024 ! maximum allowed value of n    
    INTEGER, PARAMETER :: ncmax = 50
    REAL*4 :: wksp(nmax), cc(ncmax), cr(ncmax) 
    INTEGER :: ncof, ioff, joff
    INTEGER :: i, ii, j, jf, jr, k, n1, ni, nj, nh, nmod
    REAL*4 :: ai, ai1

    IF(n .LT. 4) RETURN

    IF(n .GT. nmax) THEN 
        PRINT*,'nmax too small in daub4...'
        STOP
    END IF

    ! set filter co-efficients
    
    ncof = filter
    CALL wtset(ncof, ioff, joff, cc, cr)

    nmod = ncof*n        ! A positive constant equal to zero mod n.
    n1 = n-1             ! Mask of all bits, since n a power of 2.
    nh = n/2
    DO j=1,n
        wksp(j) = 0.
    END DO

    ! apply filter
    IF(sgn .GT. 0) THEN
        
        ii = 1
        
        DO i = 1, n, 2
        
            ni = i + nmod + ioff ! Pointer to be incremented and wrapped-around.
            nj = i + nmod + joff
           
            DO k = 1, ncof
                jf = IAND(n1,ni+k) ! We use bitwise and to wrap-around the pointers.
                jr = IAND(n1,nj+k)
                wksp(ii) = wksp(ii) + cc(k) * a(jf+1)        ! these are the smooth coefficients (stored in first half of array)
                wksp(ii+nh) = wksp(ii+nh) + cr(k) * a(jr+1)  ! these are the detail coefficients (stored in the second half of array)
            END DO

            ii = ii + 1
        
        END DO

    ELSE ! inverse transform

        ii = 1
        
        DO i = 1, n, 2
        
            ai = a(ii)
            ai1 = a(ii+nh)
            ni = i + nmod + ioff ! See comments above.
            nj = i + nmod + joff
            
            DO k = 1, ncof
                jf = IAND(n1,ni+k) + 1
                jr = IAND(n1,nj+k) + 1
                wksp(jf) = wksp(jf) + cc(k) * ai
                wksp(jr) = wksp(jr) + cr(k) * ai1
            END DO

            ii=ii+1

        END DO
                
    END IF    

    ! copy from buffer array into input array
    a(1:n) = wksp(1:n)
 

END SUBROUTINE pwt

! n-dimensional discrete wavelet transform (from Numerical Recipes)
SUBROUTINE wtn(a, ndim, sgn)

    INTEGER, INTENT(IN) :: sgn, ndim
    REAL, INTENT(INOUT) ::  a(*)
    INTEGER :: nn(ndim)
    INTEGER, PARAMETER :: NMAX  =1024
    INTEGER i1,i2,i3,idim,k,n,nnew,nprev,nt,ntot
    REAL wksp(NMAX)
    
    nn(:) = nx
    ntot=1
    
    DO idim=1,ndim
        ntot = ntot * nn(idim)
    END DO
    
    nprev=1
    
    DO idim  = 1,ndim ! Main loop over the dimensions.
    
        n = nn(idim)
        nnew = n*nprev
        
        IF(n .gt. 4) THEN
        
            DO i2 = 0, ntot-1, nnew
             
                DO i1 = 1, nprev
                
                    i3=i1+i2
                    
                    DO k = 1, n         ! Copy the relevant row or column or etc. into
                        wksp(k) = a(i3) !workspace.
                        i3 = i3+nprev
                    END DO
                    
                    IF (sgn.ge.0) THEN !Do one-dimensional wavelet transform.
                        
                        nt = n
                        DO WHILE (nt .ge. 4) 
                            CALL pwt(wksp,nt,daub_num,sgn)
                            nt = nt/2
                        END DO
                        
                    ELSE !Or inverse transform.
                        
                        nt = 4
                        DO WHILE(nt .le. n)
                            CALL pwt(wksp,nt,daub_num,sgn)
                            nt=nt*2
                        END DO
                        
                    endif
                    
                    i3 = i1 + i2
                        
                    DO k = 1,n ! Copy back from workspace.
                        a(i3)=wksp(k)
                        i3=i3+nprev
                    END DO
                
                END DO
          
            END DO
                
        END IF
       
        nprev = nnew
       
    END DO


END SUBROUTINE wtn


SUBROUTINE file_output()

    INTEGER :: i, j, k, isc, ix, iy, iz
    CHARACTER(LEN=100) :: filename
    CHARACTER(LEN=6) :: uniti
   
    IF(dump_num<10) THEN
        WRITE(uniti,'(I1.1)') dump_num
    ELSE IF(dump_num>=10 .and. dump_num<100) THEN
        WRITE(uniti,'(I2.2)') dump_num
    ELSE IF(dump_num>=100 .and. dump_num<1000) THEN
        WRITE (uniti,'(I3.3)') dump_num
    ELSE IF(dump_num>=1000 .and. dump_num<10000) THEN
        WRITE (uniti,'(I4.3)') dump_num
    ELSE IF(dump_num>=10000 .and. dump_num<100000) THEN
        WRITE (uniti,'(I5.3)') dump_num  
    END IF
    
    
    filename =TRIM('Output/reconstruction')//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')                           
    OPEN(UNIT = 10, FILE = filename, FORM = 'UNFORMATTED', ACCESS='STREAM')

    DO k = 0, nranks_z-1
        DO j = 0, nranks_y-1
            DO i = 0, nranks_x-1
    
     
                DO isc = 1, nscalars
                DO iz = 1, nz
                DO iy = 1, ny
                DO ix = 1, nx
    
                WRITE(10) dat_rec(i,j,k)%patch_array(ix,iy,iz,isc)
        
                END DO
                END DO
                END DO
                END DO
        
            END DO
        END DO
    END DO
   
    CLOSE(UNIT=10) 

   
END SUBROUTINE file_output







END PROGRAM tracer_reconstruction








