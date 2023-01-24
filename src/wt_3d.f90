! Implementation of 3D Fast Discrete Wavelet Transform (Mallat's algorithm)
! (via 1D wavelet transforms along each coordinate direction)
! Reference: Numerical Recipes Fortran 77, Ch. 13

PROGRAM wavelet_3d

USE MPI

IMPLICIT NONE


TYPE gridtype

    INTEGER :: nn
    REAL*4, ALLOCATABLE :: Dh(:,:,:), Dv(:,:,:), Dd(:,:,:) 

END type gridtype

TYPE leveltype
    
    INTEGER :: nxl
    REAL, ALLOCATABLE :: C(:,:)
    
END TYPE leveltype


INTEGER, PARAMETER :: nx = 32
INTEGER, PARAMETER :: levels = INT(LOG(SNGL(nx))/LOG(2.0))-1
INTEGER, PARAMETER :: daub_num = 20 ! set this to either 4 or 12 or 20
REAL*4,  PARAMETER :: TWOPI  = 8.d0*ATAN(1.d0)
INTEGER :: i, j, k
REAL*4, ALLOCATABLE :: dat_in(:,:,:), dat_out(:,:,:), fwb_loc(:,:,:,:)
REAL*8 :: t1, t2
TYPE(gridtype) :: grid(levels)


! allocate memory for detail coefficient arrays for each grid level
DO i = 1, levels

    j = nx/(2**i) 
    grid(i)%nn = j    
    ALLOCATE(grid(i)%Dh(j,j,j))
    ALLOCATE(grid(i)%Dv(j,j,j))
    ALLOCATE(grid(i)%Dd(j,j,j))

END DO

ALLOCATE(dat_in(nx,nx,nx), dat_out(nx,nx,nx), fwb_loc(nx,nx,nx,3)) 

! set up a test function
PRINT*,'Setting up test function...'

dat_in = 0.d0
dat_out = 0.d0

!CALL compute_basis_wavelet(8,8,13,dat_in)
!dat_in = dat_in*1.e5

!DO k = 1, nx
!DO j = 1, nx
!DO i = 1, nx
!    dat_in(i,j,k) = 100.0 * SIN(26.0*(TWOPI/nx)*i) * EXP( -(SNGL(k-0.5*nx)/SNGL(nx/20))**2 - (SNGL(j-0.65*nx)/SNGL(nx/20))**2 - (SNGL(i-0.65*nx)/SNGL(nx/20))**2 )
!END DO
!END DO
!END DO


! compute the wavelet transform
!PRINT*,'Computing wavelet transform...'
t1 =  MPI_WTIME()

!CALL fwt_3d(dat_in, dat_out, 1)

!dat_out = dat_in
!CALL wtn(dat_out,3,1)


! decompose detail coefficients at each level
!CALL decompose_detailcoeff(dat_out)

! save to file
!PRINT*,'Saving to file...'
!CALL file_output()

PRINT*,'Computing wavelet basis functions...'
CALL generate_wavelet_basis()


t2 =  MPI_WTIME()


DEALLOCATE(dat_in, dat_out, fwb_loc)

DO i = 1, levels
    DEALLOCATE(grid(i)%Dh,grid(i)%Dv,grid(i)%Dd)
END DO

PRINT*,''
PRINT*,'WT time (sec) = ',t2-t1
PRINT*,''

PRINT*,'Done.'

CONTAINS


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



SUBROUTINE generate_wavelet_basis()

    INTEGER :: i, j, k, ix, iy, iz
    REAL :: fx(nx,nx,nx), x(3)
    CHARACTER(LEN=100) :: filename1, filename2
    CHARACTER(LEN=6) :: uniti, unitj, unitk

    filename1 = TRIM('Output/fwb_3d.dat')
    filename2 = TRIM('Output/fwb_3d_loc.dat')
    OPEN(UNIT = 11, FILE = filename1, FORM = 'UNFORMATTED', ACCESS='STREAM')
    OPEN(UNIT = 12, FILE = filename2, FORM = 'UNFORMATTED', ACCESS='STREAM')

    DO k = 1, nx 
    DO j = 1, nx 
    DO i = 1, nx 

        PRINT*,'Wavelet # ',i,j,k

        ! compute basis wavelet funciton
        CALL compute_basis_wavelet(i,j,k,fx)


        ! find wavelet center by computing the first moment of it's magnitude
        !x(:) = 0.0
        !DO iz = 1, nx
        !    DO iy = 1, nx
        !        DO ix = 1, nx
        !            x(1) = x(1) + i * ABS(fx(ix,iy,iz))
        !            x(2) = x(2) + j * ABS(fx(ix,iy,iz))
        !            x(3) = x(3) + k * ABS(fx(ix,iy,iz))
        !        END DO
        !    END DO
        !END DO
        !x(:) = x(:) / SUM(SQRT(fx**2))
        
        x = MINLOC(fx) 
        

        WRITE(12) x(1),x(2),x(3)

        ! save to file
        DO iz = 1, nx
            DO iy = 1, nx
                DO ix = 1, nx
                    WRITE(11) fx(ix,iy,iz)
                END DO
            END DO
        END DO
        
    END DO
    END DO
    END DO

    CLOSE(UNIT = 11)
    CLOSE(UNIT = 12)

END SUBROUTINE generate_wavelet_basis


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
            DO kx = 1, nx
                fw(kx,iy,iz) = buffer(kx)  
            END DO

        END DO       
    END DO       
        
    !PRINT*,' Doing y pass..'    


    ! DWT in y-direction
    DO iz = 1, nx 
        DO kx = 1, nx 
                
            ! copy strips-along y into 1d buffer    
            DO iy = 1, nx
                buffer(iy) = fw(kx,iy,iz)
            END DO
        
            ! perform 1D dwt 
            CALL fwt_1d(buffer,nx,sgn)
                    
            ! copy back
            DO ky = 1, nx
                fw(kx,ky,iz) = buffer(ky)  
            END DO            
                
        END DO
    END DO
    

    ! DWT in z-direction
    DO ky = 1, nx 
        DO kx = 1, nx 
                
            ! copy strips-along y into 1d buffer    
            DO iz = 1, nx
                buffer(iz) = fw(kx,ky,iz)
            END DO
        
            ! perform 1D dwt 
            CALL fwt_1d(buffer,nx,sgn)
                    
            ! copy back
            DO kz = 1, nx
                fw(kx,ky,kz) = buffer(kz)  
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


! subroutine for decomposing the detail co-efficients at each scale level from the main wt matrix
SUBROUTINE decompose_detailcoeff(fw)

    REAL*4, INTENT(IN) :: fw(nx,nx,nx) 
    INTEGER :: i, j, k, l, ix, iy
    
    !PRINT*,'Detail levels = ',levels
    
    ! extract horizontal, vertical and diagonal detail coefficients for each level
    DO l = 1, levels
    
        !PRINT*,'level, nn =',l,grid(l)%nn
    
        ix = grid(l)%nn
        iy = grid(l)%nn
        
        ! horizontal detail coefficients
        DO k = 1, ix 
            DO j = 1, ix 
                DO i = ix+1, 2*ix
                    grid(l)%Dh(i-ix,j,k) = fw(i,j,k)     
                END DO
            END DO
        END DO

        ! vertical detail coefficients
        DO k = 1, ix 
            DO j = iy+1, 2*iy 
                DO i = 1, iy
                    grid(l)%Dv(i,j-iy,k) = fw(i,j,k)     
                END DO
            END DO
        END DO

        ! diagonal detail coefficients
        DO k = 1, nx
            DO j = ix+1, 2*iy 
                DO i = ix+1, 2*ix
                    grid(l)%Dd(i-ix,j-iy,k) = fw(i,j,k)     
                END DO
            END DO
        END DO


    END DO

    
END SUBROUTINE decompose_detailcoeff


SUBROUTINE file_output()

    INTEGER :: i, j, k, l
    CHARACTER(LEN=100) :: filename
    CHARACTER(LEN=6) :: uniti
    
    filename =TRIM('Output/wavelet_3d.dat')                        
    OPEN(UNIT = 10, FILE = filename, FORM = 'UNFORMATTED', ACCESS='STREAM')

    DO k = 1, nx
    DO j = 1, nx
    DO i = 1, nx
        WRITE(10) dat_in(i,j,k), dat_out(i,j,k)
    END DO
    END DO
    END DO

    CLOSE(UNIT=10) 

    GO TO 101
    
    DO l = 1, levels

        IF(l<10) THEN
            WRITE(uniti,'(I1.1)') l
        ELSE IF(l>=10 .and. l<100) THEN
            WRITE(uniti,'(I2.2)') l
        ELSE IF(l>=100 .and. l<1000) THEN
            WRITE (uniti,'(I3.3)') l
        ELSE IF(l>=1000 .and. l<10000) THEN
            WRITE (uniti,'(I4.3)') l
        ELSE IF(l>=10000 .and. l<100000) THEN
            WRITE (uniti,'(I5.3)') l  
        END IF
        
        filename = TRIM('Output/d_level_')//TRIM(uniti)//TRIM('.dat')
        OPEN(UNIT = 11, FILE = filename, FORM = 'UNFORMATTED', ACCESS='STREAM')

        DO k = 1, grid(l)%nn
            DO j = 1, grid(l)%nn
                DO i = 1, grid(l)%nn
                    WRITE(11) grid(l)%Dh(i,j,k),grid(l)%Dv(i,j,k),grid(l)%Dd(i,j,k)
                END DO
            END DO
        END DO

        CLOSE(UNIT = 11)

    END DO

    101 CONTINUE

END SUBROUTINE file_output



END PROGRAM wavelet_3d