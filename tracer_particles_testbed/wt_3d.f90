! Implementation of 3D Fast Discrete Wavelet Transform (Mallat's algorithm)
! (via 1D wavelet transforms along each coordinate direction)
! Reference: Numerical Recipes Fortran 77, Ch. 13

PROGRAM wavelet_3d

!USE MPI
USE OMP_LIB

IMPLICIT NONE


INTEGER, PARAMETER :: max_per_bin = 200

TYPE leveltype
    
    INTEGER :: nxl
    REAL, ALLOCATABLE :: C(:,:)
    
END TYPE leveltype

TYPE wbin

    INTEGER :: wcount = 0
    INTEGER :: wavelet_index(max_per_bin) = 0
    REAL :: amplitude(max_per_bin) = 0.0

END TYPE wbin 

INTEGER, PARAMETER :: nx = 64
INTEGER, PARAMETER :: bin_width = 2
INTEGER, PARAMETER :: max_wavelets  = 100000
INTEGER, PARAMETER :: bin_neighborhood  = 2
INTEGER, PARAMETER :: levels = INT(LOG(SNGL(nx))/LOG(2.0))-1
INTEGER, PARAMETER :: daub_num = 20 ! set this to either 4 or 12 or 20
REAL*4,  PARAMETER :: TWOPI  = 8.d0*ATAN(1.d0)
LOGICAL, PARAMETER :: output_basis = .FALSE.
LOGICAL, PARAMETER :: output_loc = .TRUE.
LOGICAL, PARAMETER :: output_bins = .TRUE.
INTEGER :: i, j, k
REAL*4, ALLOCATABLE :: dat_in(:,:,:), dat_out(:,:,:), dat_rec(:,:,:), fwb_loc(:,:,:,:)
REAL*8 :: t1, t2, p(3)
REAL*4 :: wave_x(3), wave_x2(3)
INTEGER :: x_bins, x(nx,nx,nx,3), tracer_loc(3), tracer_bin(3), tracer_wavelet(max_wavelets,2), wcount
TYPE(wbin), ALLOCATABLE :: wavelet_bin(:,:,:) 

x_bins = nx/bin_width

ALLOCATE(dat_in(nx,nx,nx), dat_out(nx,nx,nx), dat_rec(nx,nx,nx), fwb_loc(nx,nx,nx,3), wavelet_bin(1:x_bins,1:x_bins,1:x_bins)) 

! set up a test function
PRINT*,'Setting up test function...'

dat_in = 0.d0
dat_out = 0.d0
dat_rec = 0.d0

!CALL compute_basis_wavelet(8,8,13,dat_in)
!dat_in = dat_in*1.e5
wave_x = (/ 0.25*nx, 0.25*nx, 0.75*nx/) 
wave_x2 = (/ 0.6*nx, 0.6*nx, 0.75*nx/) 

DO k = 1, nx
DO j = 1, nx
DO i = 1, nx
    dat_in(i,j,k) = 100.0 * SIN(8.0*(TWOPI/nx)*i) * EXP( -(SNGL(k-wave_x(3))/SNGL(nx/16))**2 - (SNGL(j-wave_x(2))/SNGL(nx/16))**2 - (SNGL(i-wave_x(1))/SNGL(nx/16))**2 ) +&
                    200.0 * SIN(8.0*(TWOPI/nx)*i) * EXP( -(SNGL(k-wave_x2(3))/SNGL(nx/16))**2 - (SNGL(j-wave_x2(2))/SNGL(nx/16))**2 - (SNGL(i-wave_x2(1))/SNGL(nx/16))**2 ) 
END DO
END DO
END DO

! assign test tracer particle positions

CALL RANDOM_NUMBER(p)
 
tracer_loc(1) =  wave_x2(1) + 2 * p(1)
tracer_loc(2) =  wave_x2(2) + 2 * p(2)
tracer_loc(3) =  wave_x2(3) + 2 * p(3) 
tracer_bin(1) = 1 + (tracer_loc(1)-1)/bin_width
tracer_bin(2) = 1 + (tracer_loc(2)-1)/bin_width
tracer_bin(3) = 1 + (tracer_loc(3)-1)/bin_width 
wcount = 0

PRINT*,'Tracer location = ',tracer_loc
PRINT*,'Tracer bin = ',tracer_bin

! compute the wavelet transform
!PRINT*,'Computing wavelet transform...'
t1 = OMP_GET_WTIME() ! MPI_WTIME()

CALL fwt_3d(dat_in, dat_out, 1)

PRINT*,'Computing wavelet basis functions...'
!CALL generate_wavelet_basis()
!CALL predict_loc(dat_out, x)

CALL get_wavelet_locations(x)


t2 = OMP_GET_WTIME() ! MPI_WTIME()



! save to file
!PRINT*,'Saving to file...'
CALL file_output()


DEALLOCATE(dat_in, dat_out, dat_rec, fwb_loc, wavelet_bin)

PRINT*,''
PRINT*,'WT time (sec) = ',t2-t1
PRINT*,''

PRINT*,'Done.'

CONTAINS

SUBROUTINE get_wavelet_locations(x)

    INTEGER, INTENT(INOUT) :: x(nx,nx,nx,3)
    REAL :: fx(nx,nx,nx), fw(nx,nx,nx)
    INTEGER :: i, j, k, ibin, jbin, kbin, ix, nt, nbasis, ntot, cmplt1, cmplt2, wc
    CHARACTER(LEN=100) :: filename1, filename2
    CHARACTER(LEN=6) :: uniti, unitj, unitk

    IF(output_bins) THEN
        filename1 = TRIM('Output/fwb_3d_bins.dat')
        OPEN(UNIT = 11, FILE = filename1, FORM = 'UNFORMATTED', ACCESS='STREAM')
    END IF
    
    IF(output_loc) THEN
        filename2 = TRIM('Output/fwb_3d_loc.dat')
        OPEN(UNIT = 12, FILE = filename2, FORM = 'UNFORMATTED', ACCESS='STREAM')
    END IF

    nbasis = 1
    ntot = nx**3
    cmplt1 = 100*(REAL(nbasis)/REAL(ntot))
    cmplt2 = cmplt1
    wc = 0
    
    ! loop over all basis wavelets and get their locations
    DO k = 1, nx
        DO j = 1, nx
            DO i = 1, nx
                

                cmplt1 = 100*(REAL(nbasis)/REAL(ntot))

                IF(cmplt1 .NE. cmplt2) PRINT*,'Generating wavelet basis. % complete = ',cmplt1
           
                nbasis = nbasis + 1
                cmplt2 = cmplt1
                
                
                fw(:,:,:) = 0.0
                fw(i,j,k) = 1.0

                ! perform inverse dwt to get the basis wavelet
                CALL fwt_3d(fw,fx,-1)   

                ! assign a "location" to this basis wavelet
                x(i,j,k,:) = 0.5*(MAXLOC(fx)+MINLOC(fx)) ! MAXLOC(fx) ! 

                ! put this wavelet into the appropriate bin
                ix = i + (j-1) * nx + (k-1) * nx * nx  ! flattened index
                
                ibin = 1 + (x(i,j,k,1) - 1) / bin_width 
                jbin = 1 + (x(i,j,k,2) - 1) / bin_width  
                kbin = 1 + (x(i,j,k,3) - 1) / bin_width
 
                wavelet_bin(ibin,jbin,kbin)%wcount = wavelet_bin(ibin,jbin,kbin)%wcount + 1
                
                IF(wavelet_bin(ibin,jbin,kbin)%wcount .LE. max_per_bin) THEN
                
                    ! store the flattened index
                    wavelet_bin(ibin,jbin,kbin)%wavelet_index(wavelet_bin(ibin,jbin,kbin)%wcount) = ix
                ELSE
                    PRINT*,'Error! Wavelet bin ran out of space. Need to make max_per_bin bigger...'
                    STOP
                END IF
 

                ! add the wavelet to tracer reconstruction
                IF((ABS(tracer_bin(1)-ibin) .LE. bin_neighborhood) .AND. (ABS(tracer_bin(2)-jbin) .LE. bin_neighborhood) &
                 .AND. (ABS(tracer_bin(3)-kbin) .LE. bin_neighborhood)) THEN
                
                    wc = wc+1
                    tracer_wavelet(wc,1) = ix
                    tracer_wavelet(wc,2) = dat_out(i,j,k)

                    dat_rec(:,:,:) = dat_rec(:,:,:) + dat_out(i,j,k)*fx(:,:,:)                    
                
                    PRINT*,'Picking up wavelet # ',i,j,k
                    PRINT*,'# of wavelets picked up by tracer = ',wc
                
                END IF
 
                !dat_rec(:,:,:) = dat_rec(:,:,:) + dat_out(i,j,k)*fx(:,:,:)                             
 
            END DO
        END DO
    END DO


    CALL get_wavelet_amplitudes(dat_out)

    PRINT*,'Saving to file..'

    IF(output_bins) THEN
        DO kbin = 1, x_bins 
        DO jbin = 1, x_bins 
        DO ibin = 1, x_bins
        DO ix = 1, max_per_bin
            WRITE(11)  REAL(wavelet_bin(ibin,jbin,kbin)%wavelet_index(ix)),wavelet_bin(ibin,jbin,kbin)%amplitude(ix) 
        END DO
        END DO
        END DO
        END DO
    END IF


    IF(output_loc) THEN
        DO k = 1, nx 
        DO j = 1, nx 
        DO i = 1, nx
            WRITE(12) REAL(x(i,j,k,1)),REAL(x(i,j,k,2)),REAL(x(i,j,k,3))
        END DO
        END DO
        END DO
    END IF
   
    IF(output_bins) CLOSE(UNIT = 11)
    IF(output_loc) CLOSE(UNIT = 12)

END SUBROUTINE get_wavelet_locations



SUBROUTINE get_wavelet_amplitudes(fw)
 
    REAL, INTENT(IN) :: fw(nx,nx,nx)
    INTEGER :: ibin, jbin, kbin, n, ix, i, j, k
    
    PRINT*,'Acquiring wavelet amplitudes..'
    
    ! loop over wavelet bins     
    DO kbin = 1, x_bins
        DO jbin = 1, x_bins 
            DO ibin = 1, x_bins
                DO n = 1, wavelet_bin(ibin,jbin,kbin)%wcount            
        
                    ix = wavelet_bin(ibin,jbin,kbin)%wavelet_index(n)
                    
                    ! collapse flattened wavelet index
                    i = 1 + MOD(ix-1,nx)
                    j = 1 + MOD(ix-1,nx*nx)/nx !1 + (1 + MOD(ix-1,nx*nx) -1)/nx   !1 + (ix-1) / nx
                    k = 1 + (ix-1) / (nx*nx) 
       
                    IF((i .GT. nx) .OR. (j .GT. nx) .OR. (k .GT. nx)) THEN
                        PRINT*,'Error! Invalid index: ix,i,j,k = ',ix,i,j,k
                        STOP
                    END IF        
                   
        
                    ! get the amplitude
                    wavelet_bin(ibin,jbin,kbin)%amplitude(n) = fw(i,j,k)
        
                END DO
            END DO
        END DO
    END DO


END SUBROUTINE get_wavelet_amplitudes


SUBROUTINE generate_wavelet_basis()

    INTEGER :: i, j, k, ix, iy, iz, nbasis, ntot, cmplt1, cmplt2
    REAL :: fx(nx,nx,nx), xmin(3)
    CHARACTER(LEN=100) :: filename1, filename2
    CHARACTER(LEN=6) :: uniti, unitj, unitk

    IF(output_basis) THEN
        filename1 = TRIM('Output/fwb_3d.dat')
        OPEN(UNIT = 11, FILE = filename1, FORM = 'UNFORMATTED', ACCESS='STREAM')   
    END IF

    IF(output_loc) THEN
        filename2 = TRIM('Output/fwb_3d_loc.dat')
        OPEN(UNIT = 12, FILE = filename2, FORM = 'UNFORMATTED', ACCESS='STREAM')
    END IF

    nbasis = 1
    ntot = nx**3
    cmplt1 = 100*(REAL(nbasis)/REAL(ntot))
    cmplt2 = cmplt1
    
    DO k = 1, nx 
    DO j = 1, nx 
    DO i = 1, nx 

        cmplt1 = 100*(REAL(nbasis)/REAL(ntot))

        IF(cmplt1 .NE. cmplt2) PRINT*,'Generating wavelet basis. % complete = ',cmplt1
   
        nbasis = nbasis + 1
        cmplt2 = cmplt1
       
       
        ! compute basis wavelet funciton
        CALL compute_basis_wavelet(i,j,k,fx)


        ! find wavelet's "location"
        !xmin = 0.5*(MINLOC(fx)+MAXLOC(fx)) 
        xmin = MAXLOC(fx) 
        x(i,j,k,:) = xmin(:)

        ! save to file
        IF(output_basis) THEN
            DO iz = 1, nx
                DO iy = 1, nx
                    DO ix = 1, nx
                        WRITE(11) fx(ix,iy,iz)
                    END DO
                END DO
            END DO   
        END IF
        
    END DO
    END DO
    END DO
   


    IF(output_loc) THEN
        DO k = 1, nx 
        DO j = 1, nx 
        DO i = 1, nx
            WRITE(12) x(i,j,k,1),x(i,j,k,2),x(i,j,k,3)
        END DO
        END DO
        END DO
    END IF
   
    
    IF(output_basis) CLOSE(UNIT = 11)
    IF(output_loc) CLOSE(UNIT = 12)

END SUBROUTINE generate_wavelet_basis


! computes "location" of the test wave using the wavelets coefficients
SUBROUTINE predict_loc(fw, x)

    REAL :: fw(nx,nx,nx), x(nx,nx,nx,3), locmax(3), locmin(3)
    INTEGER :: xmax(3), xmin(3)
    
    ! array indices for largest wavelet amplitudes
    xmin  = MINLOC(fw)
    xmax  = MAXLOC(fw)

    PRINT*,''
    PRINT*,'xmin = ',xmin
    PRINT*,'xmax = ',xmax
    PRINT*,''

    locmax(:) = x(xmax(1),xmax(2),xmax(3),:)/nx
    locmin(:) = x(xmin(1),xmin(2),xmin(3),:)/nx

    PRINT*,''
    PRINT*,'Actual Location = ',wave_x/nx
    PRINT*,''
    PRINT*,'Wavelet predicted location # 1 = ', locmin
    PRINT*,'Wavelet predicted location # 2 = ', locmax
    PRINT*,''

END SUBROUTINE predict_loc



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
    !$OMP PARALLEL PRIVATE(buffer,ix,iy,iz) SHARED(fx,fw)
    !$OMP DO
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
    !$OMP END PARALLEL    
    
    !PRINT*,' Doing y pass..'    


    ! DWT in y-direction
    !$OMP PARALLEL PRIVATE(buffer,ix,iy,iz) SHARED(fx,fw)
    !$OMP DO
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
    !$OMP END PARALLEL

    ! DWT in z-direction
    !$OMP PARALLEL PRIVATE(buffer,ix,iy,iz) SHARED(fx,fw)
    !$OMP DO
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
    !$OMP END PARALLEL
    
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

    INTEGER :: i, j, k, l
    CHARACTER(LEN=100) :: filename
    CHARACTER(LEN=6) :: uniti
    
    filename =TRIM('Output/wavelet_3d.dat')                        
    OPEN(UNIT = 10, FILE = filename, FORM = 'UNFORMATTED', ACCESS='STREAM')

    DO k = 1, nx
    DO j = 1, nx
    DO i = 1, nx
        WRITE(10) dat_in(i,j,k),dat_out(i,j,k),dat_rec(i,j,k)
    END DO
    END DO
    END DO

    CLOSE(UNIT=10) 

   
END SUBROUTINE file_output



END PROGRAM wavelet_3d