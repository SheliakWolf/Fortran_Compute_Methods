! Environment: Ubuntu 22.04 + gfortran
! Create: 2024/10/09 Wed
! Author: Jiarui Xiong(Sheliak Wolf)
! Last update: 2024/10/14 Mon

! gfortran -c interpolation.f90
! 说明：本程序依次包含牛顿法多项式插值、拉格朗日插值、分段线性插值、厄米（Hermite）插值、分段三次厄米插值。
module Interpolation

    implicit none


    private
    public      :: newton, lagrange, linear, hermite, cubic_hermite, spline
    public      :: my

    type        :: my_type
    contains
        procedure,nopass    :: newton
        procedure,nopass    :: lagrange
        procedure,nopass    :: linear
        procedure,nopass    :: hermite
        procedure,nopass    :: cubic_hermite
        procedure,nopass    :: spline
    end type my_type
    type(my_type)   :: my

    contains
        subroutine newton(n,xi,yi,m,input,output) 
!用法：n：数据点数，xi：数据点的x坐标，yi：数据点的y坐标，m：你希望插值的点数量，input：插值点的x坐标，output：程序输出，你只需声明数组
            implicit none

            integer         :: n
            integer         :: m
            real            :: xi(n), yi(n)
            real            :: input(m), output(m)
            real            :: interpolation(n+1,n)
            integer         :: i, j

            interpolation = 0
            interpolation(1,:) = yi(:)
            do i=2,n
                do j=1,n+1-i
                    interpolation(i,j) = (interpolation(i-1,j+1) - interpolation(i-1,j))/(xi(j+i-1)-xi(j))
                enddo
                print *, interpolation(i,:)
            enddo


            !Using Qingjiushao method to calculate
            do i=1,m
                output(i) = interpolation(n,1) * (input(i) - xi(n-1)) + interpolation(n-1,1)
                do j=2,n-2

                    output(i) = output(i)* (input(i) - xi(n-j)) + interpolation(n-j,1)
            
                enddo
                output(i) = output(i) * (input(i) - xi(1)) + yi(1)
            enddo

        end subroutine

    !Lagrange
        subroutine Lagrange(n,xi,yi,m,input,output)
!同上
            implicit none
            integer         :: n, m
            real            :: xi(n), yi(n)
            real            :: input(m), output(m)
            integer         :: i, j, k
            real            :: interpolation(n,m)

            interpolation = 1
            output = 0
            do i=1,m
                do j=1,n
                    do k=1,n
                        if (k == j) cycle
                        interpolation(j,i) = interpolation(j,i) * 1.0 * (input(i) - xi(k)) / (1.0*xi(j) - 1.0*xi(k))
                    enddo
                enddo
            enddo

            do i=1,m
                do j=1,n
                    output(i) = output(i) + yi(j) * interpolation(j,i) * 1.0
                enddo
            enddo


        end subroutine
        
        subroutine linear(n,x,y,m,input,output)
!线性插值，用法同上
            implicit none
            integer         :: n,m,i,j,temp
            real            :: x(n),y(n),input(m),output(m)

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
                        y(j+1) = temp
                    endif
                 end do
            enddo

            do i=1,m
                do j=1,n-1
                    if (abs(input(i) - x(j)) < 0.001 * abs(x(j) - x(j+1))) then
                        output(i) = y(j)
                        exit
                    endif
                enddo
                if (input(i) < x(1)) then
                    output(i) = (input(i) - x(1)) * (y(2) - y(1))/(x(2)-x(1)) + y(1)
                else if (input(i) > x(n-1)) then
                    output(i) = (input(i) - x(n)) * (y(n) - y(n-1))/(x(n)-x(n-1)) + y(n)
                else
                    do j=1,n-1
                        if((input(i) > x(j) .and. input(i) < x(j+1)) .eqv. .true.) then
                            output(i) = (input(i) - x(j)) * (y(j+1) - y(j))/(x(j+1) - x(j))*1.0 + y(j)*1.0
                            exit
                        endif
                    enddo
                endif
            enddo

        end subroutine
        

        subroutine hermite(n,x,y,yp,m,input,output)
!厄米插值，需要注意，yp为函数在数据点处的导数值。
            implicit none
            integer         :: n,m,i,j,k,temp
            real            :: x(n), y(n), yp(n), input(m), output(m)
            real            :: lagrange(n,m), hermite1(2,n,m)

            
            temp = 0
            do i=1,n-1
                do j=1,i
                    if (x(j) > x(j+1)) then
                        temp = x(j)
                        x(j) = x(j+1)
                        x(j+1) = temp
                        temp = y(j)
                        y(j) = y(j+1)
                        y(j+1) = temp
                        temp = yp(j)
                        yp(j) = yp(j+1)
                        yp(j+1) = temp
                    endif
                 end do
            enddo

            lagrange = 1.0
            hermite1 = 0.0
            output = 0.0
            do i=1,m
                do j=1,n
                    do k=1,n
                        if (k == j) cycle
                        lagrange(j,i) = lagrange(j,i) * (1.0*input(i) - 1.0*x(k)) / (1.0*x(j) - 1.0*x(k))
                    enddo
                enddo
            enddo
            
            do i=1,n
                do j=1,n
                    do k=1,m
                        if (i .ne. j) hermite1(1,i,k) = hermite1(1,i,k) + 1.0 / (x(j) - x(i))
                    enddo
                enddo
            enddo
            do i=1,n
                do j=1,m
                    hermite1(2,i,j) = (lagrange(i,j)**2.0) * (input(j) - x(i))
                    hermite1(1,i,j) = (1.0 + 2.0 * (input(j) - x(i)) * hermite1(1,i,j)) * lagrange(i,j)**2.0
                enddo
            enddo
            do i=1,n
                do j=1,m
                    output(j) = output(j) +  y(i) * hermite1(1,i,j) + yp(i) * hermite1(2,i,j)
                enddo
            enddo
            do i=1,m
                do j=1,n-1
                    if (abs(input(i) - x(j)) < (0.001 * (x(j+1) - x(j)))) then
                        output(i) = y(j)
                        exit
                    endif
                enddo
            enddo

        end subroutine

        subroutine cubic_hermite(n,x,y,yp,m,input,output)
!同上
            implicit none

            integer         :: n,m,i,j,k,temp
            real            :: x(n), y(n), yp(n), input(m), output(m)
            real            :: lagrange(n,m), hermite1(2,n,m)


            i = 0
            j = 0
            temp = 0
            do i=1,n-1
                do j=1,i
                    if (x(j) > x(j+1)) then
                        temp = x(j)
                        x(j) = x(j+1)
                        x(j+1) = temp
                        temp = y(j)
                        y(j) = y(j+1)
                        y(j+1) = temp
                        temp = yp(j)
                        yp(j) = yp(j+1)
                        yp(j+1) = temp
                    endif
                 end do
            enddo

            do i=1,m
                if (input(i) < x(2) ) then
                    call hermite(2,x(1:2),y(1:2),yp(1:2),1,input(i),output(i))
                    cycle
                endif
                if (input(i) > x(n-1)) then
                    call hermite(2,x(n-1:n),y(n-1:n),yp(n-1:n),1,input(i),output(i))
                    cycle
                endif
                do j=2,n-2
                    if (input(i) > x(j) .and. input(i) < x(j+1)) then
                        call hermite(2,x(j:j+1),y(j:j+1),yp(j:j+1),1,input(i),output(i))
                        exit
                    endif
                enddo
            enddo
            do i=1,m
                do j=1,n-1
                    if (abs(input(i) - x(j)) < (0.001 * (x(j+1) - x(j)))) then
                        output(i) = y(j)
                        exit
                    endif
                enddo
            enddo

        end subroutine



        subroutine spline(n,x,y,m,input,output,BC_type,bc)

            implicit none
            integer             :: BC_type, n, m, i, j, k, temp
            real                :: x(0:n-1), y(0:n-1), input(m), output(m),bc(2)
            real                :: mi(0:n-1), alpha(0:n-1), beta(0:n-1), A(0:n-1), B(0:n-1),h(0:n-1)

            do i=0,n-2
                do j=0,i
                    if (x(j) > x(j+1)) then
                        temp = x(j)
                        x(j) = x(j+1)
                        x(j+1) = temp
                        temp = y(j)
                        y(j) = y(j+1)
                        y(j+1) = temp
                    endif
                 end do
            enddo


            if ((BC_type == 1) .or. (BC_type == 2) .eqv. .false.) then
                print *, char(27) // '[31m' // "错误！非法的输入（BC_type）边界条件必须取1或2" 
                stop
            endif
            
            do i=0,n-2
                h(i) = x(i+1) - x(i)
            enddo
            if(BC_type == 1) then
                do i=1,n-2
                    alpha(i) = h(i-1)/(h(i-1) + h(i))
                    beta(i) = 3.0 * ((1.0 - alpha(i)) * (y(i) - y(i-1)) / h(i-1) + alpha(i) * (y(i+1) - y(i)) / h(i))
                enddo
                A(0) = 0.0
                B(0) = bc(1)
                do i=1,n-1
                    A(i) = - alpha(i) / ((1.0 - alpha(i)) * A(i-1) + 2)
                    B(i) = (beta(i) - (1.0 - alpha(i)) * B(i-1)) / ((1.0 - alpha(i)) * A(i-1) + 2)
                enddo
                mi(0) = bc(1)
                mi(n-1) = bc(2)
                do i=n-2,1,-1
                    mi(i) = A(i) * mi(i+1) + B(i)
                enddo
                call cubic_hermite(n,x,y,mi,m,input,output)
            endif

            if (BC_type == 2) then
                alpha(0) = 1.0
                beta(0) = (3.0 / h(0)) * (y(1) - y(0)) - bc(1) * h(0) / 2.0
                alpha(n-1) = 0.0
                beta(n-1) = (3.0 / h(n-2)) * (y(n-1) - y(n-2)) + bc(2) * h(n-2) / 2.0
                do i=1,n-2
                    alpha(i) = h(i-1) / (h(i) + h(i-1))
                    beta(i) = 3.0 * ((1.0 - alpha(i)) * (y(i) - y(i-1)) / h(i-1) + (y(i+1) - y(i)) * alpha(i) / h(i))
                enddo
                A(0) = -0.5 * alpha(0)
                B(0) = 0.5 * beta(0)
                do i=1,n-1
                    A(i) = -alpha(i) / ((1.0 - alpha(i)) * A(i-1) + 2.0)
                    B(i) = (beta(i) - (1.0 - alpha(i)) * B(i-1)) / ((1 - alpha(i) * A(i-1)) + 2.0)
                enddo
                mi(n-1) = B(n-1)
                do i=n-2,0,-1
                    mi(i) = A(i) * mi(i+1) + B(i)
                enddo
                call cubic_hermite(n,x,y,mi,m,input,output)
            endif
            do i=1,m
                do j=1,n-1
                    if (abs(input(i) - x(j)) < (0.001 * (x(j+1) - x(j)))) then
                        output(i) = y(j)
                        exit
                    endif
                enddo
            enddo
            
        end subroutine

end module

