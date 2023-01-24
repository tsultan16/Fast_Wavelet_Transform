PROGRAM quicksort

IMPLICIT NONE

INTEGER, PARAMETER :: arrsize = 50
REAL :: p, arr(1:arrsize) 
INTEGER :: i

DO i = 1, arrsize

 CALL RANDOM_NUMBER(p)
 arr(i) = 500.0*p 

END DO

PRINT*,''
PRINT*,'UNSORTED ARRAY= ',arr
PRINT*,''


CALL sort(arr,1,arrsize)

PRINT*,''
PRINT*,'SORTED ARRAY= ',arr
PRINT*,''

CALL sort_check(arr,1,arrsize)

CONTAINS


RECURSIVE SUBROUTINE sort(arr, low, hi)

    REAL, INTENT(INOUT) :: arr(:)
    INTEGER, INTENT(IN) :: low, hi
    INTEGER :: ip

    IF(low .LT. hi) THEN
    
        ! find partition index
        CALL partition(arr, low, hi, ip)
    
        ! apply quicksort to each of the two partitions
        CALL sort(arr, low, ip-1) ! lower than pivot
        CALL sort(arr, ip+1, hi)  ! higer than pivot
    
    END IF

END SUBROUTINE sort




SUBROUTINE partition(arr, low, hi, ip)

    REAL, INTENT(INOUT) :: arr(:)
    INTEGER, INTENT(IN) :: low, hi
    INTEGER, INTENT(OUT) :: ip
    REAL :: pivot, tmp
    INTEGER :: i, j
    
    ! choose last element as pivot
    pivot = arr(hi)
    
    ! index of smaller element    
    i = low - 1

    ! traverse through all elements
    DO j = low, hi-1
    
        IF(arr(j) .LT. pivot) THEN
        
            i = i + 1
            
            ! swap arr(i) and arr(j)
            IF(i .NE. j) THEN
                tmp = arr(i)
                arr(i) = arr(j)
                arr(j) = tmp
            END IF
                
        END IF
    
    END DO
    
    ! swap arr(i+1) and arr(hi)
    tmp = arr(i+1)
    arr(i+1) = arr(hi)
    arr(hi) = tmp
    
    ! final position of pivot
    ip = i+1

END SUBROUTINE partition


SUBROUTINE sort_check(arr, low, hi)
   
    REAL, INTENT(IN) :: arr(:)
    INTEGER, INTENT(IN) :: low, hi
    INTEGER :: i
    
    DO i = low, hi-1
        IF(arr(i) .GT. arr(i+1)) THEN
            PRINT*,'ERROR*** Array is unsorted! Aborting...'
            STOP
        END IF
    END DO

    PRINT*,'Array is sorted!'
    PRINT*,''

END SUBROUTINE sort_check


END PROGRAM quicksort