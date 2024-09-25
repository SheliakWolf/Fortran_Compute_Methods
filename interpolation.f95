! Environment: Linux(Ubuntu 24.04), Windows11
! Complier: gfortran, Intel Visual Fortran
! Author: Jiarui Xiong(Sheliak Wolf)
! Last update: 2024/09/25 Wed

! 注意！注释由GPT-4自动生成！！
program Interpolation
    ! 定义主要插值程序，包含不同插值方法的子程序

    contains
        ! 牛顿插值法子程序
        ! 参数说明:
        ! n: 插值点数
        ! xi: 自变量的插值点数组
        ! yi: 对应的因变量值数组
        ! m: 待求的输入点数
        ! input: 输入需要插值的点
        ! output: 输出插值后的值
        subroutine newtonian(n,xi,yi,m,input,output)
            implicit none

            ! 声明变量类型
            integer                 :: n  ! 插值点数
            integer                 :: m  ! 待求点数
            real                    :: xi(n), yi(n)  ! 自变量和因变量
            real                    :: input(m), output(m)  ! 输入点和输出结果
            real                    :: interpolation(n+1,n)  ! 存储插值多项式系数
            integer                 :: i, j  ! 循环变量

            ! 初始化插值矩阵，存储因变量
            interpolation = 0
            interpolation(1,:) = yi(:)

            ! 计算差分表，逐步计算牛顿插值系数
            do i=2,n
                do j=1,n+1-i
                    interpolation(i,j) = (interpolation(i-1,j+1) - interpolation(i-1,j)) / (xi(j+i-1)-xi(j))
                enddo
            enddo

            ! 使用秦九韶算法（霍纳法则）计算插值值
            do i=1,m
                ! 计算牛顿插值多项式
                output(i) = interpolation(n,1) * (input(i) - xi(n-1)) + interpolation(n-1,1)
                do j=2,n-2
                    output(i) = output(i) * (input(i) - xi(n-j)) + interpolation(n-j,1)
                enddo
                output(i) = output(i) * (input(i) - xi(1)) + yi(1)
            enddo

            ! 输出输入点和插值后的结果
            print *, input
            print *, output

        end subroutine

        ! 拉格朗日插值法子程序
        ! 参数说明:
        ! n: 插值点数
        ! xi: 自变量的插值点数组
        ! yi: 对应的因变量值数组
        ! m: 待求的输入点数
        ! input: 输入需要插值的点
        ! output: 输出插值后的值
        subroutine Lagrange(n,xi,yi,m,input,output)
            implicit none

            ! 声明变量类型
            integer                 :: n, m  ! 插值点数和输入点数
            real                    :: xi(n), yi(n)  ! 自变量和因变量
            real                    :: input(m), output(m)  ! 输入点和输出结果
            integer                 :: i, j, k  ! 循环变量
            real                    :: interpolation(n,m)  ! 存储插值基函数值

            ! 初始化插值基函数
            interpolation = 1
            output = 0

            ! 计算拉格朗日基函数
            do i=1,m
                do j=1,n
                    do k=1,n
                        if (k == j) cycle  ! 跳过自身节点
                        interpolation(j,i) = interpolation(j,i) * 1.0 * (input(i) - xi(k)) / &
                                             (1.0*xi(j) - 1.0*xi(k))
                    enddo
                enddo
            enddo

            ! 计算插值结果
            do i=1,m
                do j=1,n
                    output(i) = output(i) + yi(j) * interpolation(j,i) * 1.0
                enddo
            enddo

        end subroutine

end program
