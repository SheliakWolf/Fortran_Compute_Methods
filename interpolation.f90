
! Author: Jiarui Xiong(Sheliak Wolf)
! Last update: 2024/10/09 Wed
module Interpolation

        implicit none


        private
        public          :: newton, lagrange, linear
        public          :: my

        type            :: my_type
        contains
                procedure,nopass        :: newton
                procedure,nopass        :: lagrange
                procedure,nopass        :: linear
        end type my_type
        type(my_type)   :: my

        contains
                subroutine newton(n,xi,yi,m,input,output) !n:插值点数
                        implicit none

                        integer                 :: n
                        integer                 :: m
                        real                    :: xi(n), yi(n)
                        real                    :: input(m), output(m)
                        real                    :: interpolation(n+1,n)
                        integer                 :: i, j

!                        print *, "Newton Interpolation:"
                        interpolation = 0
                        interpolation(1,:) = yi(:)
                        do i=2,n
                                do j=1,n+1-i
                                        interpolation(i,j) = (interpolation(i-1,j+1) - interpolation(i-1,j))/(xi(j+i-1)-xi(j))
                                enddo
                                print *, interpolation(i,:)
                        enddo
     
!                        print "(A9,$)", "\phi_n = "
                        do i=2,n
!                                print "(A1,$)", "("
                        enddo
!                        print "(F10.6,A8,F5.2,A4,F10.6,A1,$)", interpolation(n,1), " * (x - ", xi(n-1), ") + ",&
!                                interpolation(n-1,1), ")"
                        do i=2,n-2
!                                print "(A6,F5.2,A4,F10.6,A1,$)", "*(x - ", xi(n-i), ") + ", interpolation(n-i,1), ")"
                        enddo
!                        print "(A6,F5.2,A4, F5.2,A1)", "*(x - ", xi(1), ") + ", yi(1), ")"

                        !Using Qingjiushao method to calculate
                        do i=1,m
                                output(i) = interpolation(n,1) * (input(i) - xi(n-1)) + interpolation(n-1,1)
                                do j=2,n-2

                                        output(i) = output(i)* (input(i) - xi(n-j)) + interpolation(n-j,1)
                        
                                enddo
                                output(i) = output(i) * (input(i) - xi(1)) + yi(1)
                        !

                        enddo

                        print *, input
                        print *, output


                end subroutine

        !Lagrange
                subroutine Lagrange(n,xi,yi,m,input,output)

                        implicit none
                        integer                 :: n, m
                        real                    :: xi(n), yi(n)
                        real                    :: input(m), output(m)
                        integer                 :: i, j, k
                        real                    :: interpolation(n,m)

                        interpolation = 1
                        output = 0
                        do i=1,m
                                do j=1,n
                                        do k=1,n
                                                if (k == j) cycle
                                                interpolation(j,i) = interpolation(j,i) * 1.0 * (input(i) - xi(k)) / &
                                                        (1.0*xi(j) - 1.0*xi(k))
                                        enddo
                                enddo
                        enddo
                        
!                        print "(A9,$)", "\phi_n = "
                        do i=1,n
!                                print "(F8.3,A5,I3,A3,$)", yi(i), " * l_", i, "(x)" 
                        enddo

                        do i=1,m
                                do j=1,n
                                        output(i) = output(i) + yi(j) * interpolation(j,i) * 1.0
                                enddo
                        enddo
!                        print *, ""
!                        print *, input
!                        print *, output

                end subroutine
                
                subroutine linear(n,x,y,m,input,output)

                        implicit none
                        integer                 :: n,m,i,j,temp
                        real                    :: x(n),y(n),input(m),output(m)

                        i = 0
                        j = 0
                        temp = 0
                        do i=1,n-1
                                do j=1,i
                                        if (x(j) > x(j+1)) then
                                                temp = x(j)
                                                x(j) = x(j+1)
                                                x(j+1) = x(j)
                                                temp = y(j)
                                                y(j) = y(j+1)
                                                y(j+1) = y(j)
                                        endif
                                 end do
                        enddo
                        do i=1,m
                                do j=1,n-1
                                        if (abs(input(i) - x(j)) < abs(x(j) - x(j+1))) then
                                                output(i) = y(j)
                                                exit
                                        endif
                                enddo
                                if (input(i) < x(1)) then
                                        output(i) = (input(i) - x(1)) * &
                                                (y(2) - y(1))/(x(2)-x(1)) + y(1)
                                else if (input(i) > x(n-1)) then
                                        output(i) = (input(i) - x(n)) *&
                                                (y(n) - y(n-1))/(x(n)-x(n-1)) + y(n)
                                else
                                        do j=1,n-1
                                                if((input(i) > x(j) .and. input(i) < x(j+1)) .eqv. .true.) then
                                                        output(i) = &
                 (input(i) - x(j)) * (y(j+1) - y(j))/(x(j+1) - x(j))*1.0 + y(j)*1.0
                                                        exit
                                                endif
                                        enddo
                                endif
                                !print *, input(i), output(i)
                        enddo
                end subroutine

end module

