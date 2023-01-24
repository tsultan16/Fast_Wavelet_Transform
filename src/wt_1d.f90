! Implementation of 1D Fast Discrete Wavelet Transform (Mallat's algorithm)
! Reference: Numerical Recipes Fortran 77, Ch. 13

PROGRAM wavelet

IMPLICIT NONE


INTEGER, PARAMETER :: nx = 256
INTEGER, PARAMETER :: daub_num = 20 ! set this to either 4 or 12 or 20
REAL*4,  PARAMETER :: TWOPI  = 8.d0*ATAN(1.d0)
INTEGER :: i, j
REAL*4 :: dat_in(1:nx), dat_out(1:nx), fftin_re(1:nx), fftin_im(1:nx), fftout_re(1:nx), fftout_im(1:nx)
CHARACTER(LEN=100) :: filename
CHARACTER(LEN=6) :: uniti



! set up a test function
dat_in = 0.d0

DO i = 1, nx

    dat_in(i) = 100.0 * SIN(25.0*(TWOPI/nx)*i) * EXP(-(SNGL(i-nx/2)/SNGL(nx/50))**2)


    GO TO 1000
    IF(i .LT. 0.45*nx) THEN
        dat_in(i) = SIN(2.0*(TWOPI/nx)*i)
    ELSE IF(i .GE. 0.45*nx .AND. i .LT. 0.58*nx) THEN
        dat_in(i) = 10.*SIN(64.0*(TWOPI/nx)*i)
    ELSE IF(i .GE. 0.58*nx .AND. i .LT. 0.85*nx) THEN
        dat_in(i) = 9.*SIN(6.0*(TWOPI/nx)*i)
    ELSE IF(i .GE. 0.85*nx .AND. i .LT. 0.95*nx) THEN
        dat_in(i) = 12.*SIN(80.0*(TWOPI/nx)*i)
    ELSE IF(i .GE. 0.95*nx) THEN
        dat_in(i) = 16.*SIN(256.0*(TWOPI/nx)*i)
    END IF
    1000 CONTINUE
END DO

dat_out = dat_in

! compute the wavelet transform
CALL fwt(dat_out, nx, 1)

! compute FFT
fftin_re = dat_in
fftin_im = 0.d0
CALL cfft_1d(fftin_re, fftin_im, fftout_re, fftout_im)

filename =TRIM('Output/wavelet_1d.dat')
                    
OPEN(UNIT = 10, FILE = filename, FORM = 'UNFORMATTED', ACCESS='STREAM')

DO i = 1, nx
    WRITE(10) dat_in(i), dat_out(i), fftout_re(i), fftout_im(i)
END DO

CLOSE(UNIT=10) 


! compute basis wavelets
CALL generate_basis_functions()


PRINT*,'Done.'

CONTAINS


SUBROUTINE generate_basis_functions()
    
    INTEGER :: i,j
    CHARACTER(LEN=100) :: filename

  
    filename =TRIM('Output/fwb_1d.dat')
                    
    OPEN(UNIT = 10, FILE = filename, FORM = 'UNFORMATTED', ACCESS='STREAM')

  
    DO j = 1, nx

        dat_in = 0.d0
        dat_in(j) = 1.d0
        dat_out = dat_in

        ! compute the basis wavelet function
        CALL fwt(dat_out, nx, -1)

        ! compute FFT of basis wavelet function
        fftin_re = dat_out
        fftin_im = 0.d0
        CALL cfft_1d(fftin_re, fftin_im, fftout_re, fftout_im)

        DO i = 1, nx
            WRITE(10) dat_out(i), fftout_re(i), fftout_im(i) 
        END DO


    END DO

    CLOSE(UNIT=10)


END SUBROUTINE generate_basis_functions



! This subroutine computes the (inverse) wavelet transform of the input data vector of length 'n' for sgn = (-1) 1 
! n has to be power of 2
SUBROUTINE fwt(a,n,sgn)

    REAL*4, INTENT(INOUT) :: a(:)    ! input data vector
    INTEGER, INTENT(IN) :: sgn, n 
    INTEGER :: nn
    
    IF(n .LT. 4) RETURN
      
    ! compute the wavelet transform  
    IF(sgn .GE. 0) THEN
        
        nn = n ! start at largest hierarchy
        
        DO WHILE(nn .GE. 4) 
            !CALL daub4(a,nn,sgn)  ! work towards smallest
            CALL pwt(a,nn,daub_num,sgn)
            nn = nn/2
        END DO
        
    ELSE ! inverse transform

        nn = 4 ! start at smallest

        DO WHILE(nn .LE. n)
            !CALL daub4(a,nn,sgn)  ! work towards largest
            CALL pwt(a,nn,daub_num,sgn)
            nn = nn*2        
        END DO
        
    END IF    


END SUBROUTINE fwt


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


! wavelet filter instance
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
                wksp(ii) = wksp(ii) + cc(k) * a(jf+1)
                wksp(ii+nh) = wksp(ii+nh) + cr(k) * a(jr+1)
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




! DAUB4 wavelet filter, i.e. multiplication of the input vector with the wavelet transform matrix (or its transpose for the inverse transform case), see eqn 13.10.2 in Numerical Recipes
SUBROUTINE daub4(a,n,sgn)

    REAL*4, INTENT(INOUT) :: a(:)    ! input data vector
    INTEGER, INTENT(IN) :: sgn, n
    INTEGER, PARAMETER :: nmax = 1024 ! maximum allowed value of n    
    REAL*4, PARAMETER :: c0 = 0.4829629131445341, c1 = 0.8365163037378079, &   ! Daub-4 filter co-efficients
                         c2 = 0.2241438680420134, c3 = -0.1294095225512604  
    
    REAL*4 :: wksp(nmax) 
    INTEGER :: nh, nh1, i, j

    IF(n .LT. 4) RETURN

    IF(n .GT. nmax) THEN 
        PRINT*,'nmax too small in daub4...'
        STOP
    END IF


    nh = n/2
    nh1 = nh + 1

    ! apply filter
    IF(sgn .GT. 0) THEN
        
        i = 1
        
        DO j = 1, n-3, 2 
            wksp(i) = c0*a(j) + c1*a(j+1) + c2*a(j+2) + c3*a(j+3)
            wksp(i+nh) = c3*a(j) - c2*a(j+1) + c1*a(j+2) - c0*a(j+3)
            i = i + 1 
        END DO
        
        wksp(i) = c0*a(n-1) + c1*a(n) + c2*a(1) + c3*a(2)
        wksp(i+nh) = c3*a(n-1) - c2*a(n) + c1*a(1) - c0*a(2)

    ELSE ! inverse transform

        wksp(1) = c2*a(nh) + c1*a(n) + c0*a(1) + c3*a(nh1)
        wksp(2) = c3*a(nh) - c0*a(n) + c1*a(1) - c2*a(nh1)
        j=3
        
        DO i = 1, nh-1 
            wksp(j) = c2*a(i) + c1*a(i+nh) + c0*a(i+1) + c3*a(i+nh1)
            wksp(j+1) = c3*a(i) - c0*a(i+nh) + c1*a(i+1) - c2*a(i+nh1)
            j = j+2
        END DO
        
    END IF    

    ! copy from buffer array into input array
    a(1:n) = wksp(1:n)
 

END SUBROUTINE daub4

! 1D Daniel-Lanczos FFT algorithm for complex  input
SUBROUTINE cfft_1d(in_re, in_im, out_re, out_im)

    REAL(4), INTENT(IN) :: in_re(1:nx), in_im(1:nx)  
    REAL(4), INTENT(INOUT) :: out_re(1:nx), out_im(1:nx)  
    INTEGER :: sgn = 1   ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    INTEGER :: ix, kx
 

    REAL(4) :: buffer(1:2*nx)   ! input array gets replaced by output
    REAL(4) :: tempr, tempi, theta, wi, wr, wpi, wpr, wtemp
    INTEGER :: i, j, n, m, mmax, istep
    INTEGER :: i1, i2, i3, i4, n2p3
    REAL*4 :: c1, c2, h1i, h1r, h2i, h2r, wis, wrs
    
 
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
    
    
    ! load input array into work buffer
    ix = 1
    DO i = 1, nx
        buffer(ix)   = in_re(i)   
        buffer(ix+1) = in_im(i) 
        ix = ix+2        
    END DO
  
  
    !***************************************
    ! Sort input array in bit-reversed order  
    !***************************************    
    n = 2*nx
    j = 1
    
    DO i = 1, n, 2
    
        IF(j .GT. i) THEN  ! swap the two complex numbers
            tempr = buffer(j)
            tempi = buffer(j+1)
            buffer(j) = buffer(i)
            buffer(j+1) = buffer(i+1)
            buffer(i) = tempr
            buffer(i+1) = tempi                    
        END IF
    
        m = n/2
        DO WHILE ((m .GT. 2) .AND. (j .GT. m))
            j = j - m 
            m = m / 2
        END DO
        j = j + m
    END DO
      
    !********************************************************************************
    ! Using Danielson-Laczos lemma, compute the DFT by summing up the 1-pt base DFT's
    !********************************************************************************     
    mmax = 2
    
    DO WHILE(n .GT. mmax) 
    
        ! initialize for trigonometric recurrence
        istep = 2 * mmax
        theta = TWOPI / (sgn*mmax)  
        wpr = -2.d0 * SIN(0.5d0 * theta)**2 
        wpi =  SIN(theta)
        wr = 1.d0
        wi = 0.d0
        
        DO m = 1, mmax, 2 
            DO i = m, n, istep
       
                j = i + mmax
                
                ! apply Danielson-Lanczos lemma
                tempr = wr*buffer(j) - wi*buffer(j+1)
                tempi = wr*buffer(j+1) + wi*buffer(j)
                buffer(j) = buffer(i) - tempr 
                buffer(j+1) = buffer(i+1) - tempi 
                buffer(i) = buffer(i) + tempr 
                buffer(i+1) = buffer(i+1) + tempi 
                
            END DO
            
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
            
        END DO
        mmax = istep
        
    END DO
      
    ix = 1
    DO kx = 1, nx
        out_re(kx) =  buffer(ix)
        out_im(kx) =  buffer(ix+1) 
        ix = ix + 2
    END DO


END SUBROUTINE cfft_1d


END PROGRAM wavelet