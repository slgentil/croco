
       recursive subroutine QSort(a,na)
       implicit none
! DUMMY ARGUMENTS
       integer, intent(in) :: nA
       real(kind=8), dimension(nA), intent(in out) :: A
! LOCAL VARIABLES
       integer :: left, right
       real(kind=8) :: random
       real(kind=8) :: pivot
       real(kind=8) :: temp
       integer :: marker
  
       if (nA > 1) then
!   call random_number(random)
         random=0.5
         pivot = A(int(random*real(nA-1))+1) ! random pivor (not best performance, but avoids worst-case)
         left = 0
         right = nA + 1
         do while (left < right)
           right = right - 1
           do while (A(right) > pivot)
             right = right - 1
           end do
           left = left + 1
           do while (A(left) < pivot)
             left = left + 1
           end do
           if (left < right) then
             temp = A(left)
             A(left) = A(right)
             A(right) = temp
           end if
         end do
         if (left == right) then
           marker = left + 1
         else
           marker = left
         end if
         call QSort(A(:marker-1),marker-1)
         call QSort(A(marker:),nA-marker+1)
       end if
       end subroutine QSort
 

