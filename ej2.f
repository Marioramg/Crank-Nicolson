      program CrankNicolson
      implicit none
C
      complex*16, dimension(:), allocatable :: psi, omega, e, f
C----- Funcion de onda y auxiliares
      real*8, dimension (:), allocatable :: pot, x
C----- Potencial y pos discretizadas
      complex*16 i  !Numero complejo
C
      real*8 x0, dx, dt, sig, k0   !Parametros iniciales
      real*8 l, pi, t, prob, norma   !lambda, pi, tiempo y probabilidad
      real*8 v !potencial auxiliar para graficas
C
      integer h, j, k    !Contadores para los bucles
      integer M, N   !Caja 0<x<1, M div espaciales y N temporales
      integer resto
C
      open(11,file='prob_gif.dat')
      open(12,file='pot.dat')
      open(13,file='prob_png.dat')
C
      M=10000
      N=150
C
      allocate(omega(1:M-1), f(1:M-1), e(1:M-1))
      allocate(psi(0:M))
      allocate(pot(0:M), x(0:M))

      i=(0.d0,1.d0)
      pi=4.d0*atan(1.d0)
C
      dx=1.d-4
      dt=5.d-7
      sig=1.d-3
      x0=4.d-1
C
      l=2.d0*dx*dx/dt
C
      do j=0,M
         x(j)=j*dx
      enddo
C----------------- Definir el Potencial --------------------------------
      do j=0,M
C----- Particula libre  (a y b)
         pot(j)=0.
         v=0.
C----- Pared  (c)
c         if(x(j).ge.0.6d0) pot(j)=1.d6
c         if(x(j).ge.0.6d0) v=0.0025
C----- Barranco  (d)
c         if(x(j).ge.0.6d0) pot(j)=-1.d6
c         if(x(j).ge.0.6d0) v=-0.0025
C----- Barrera  (e)
c         if((x(j).ge.0.6d0).and.(x(j).le.0.7d0)) pot(j)=1.25d5
c         if((x(j).ge.0.6d0).and.(x(j).le.0.7d0)) v=0.002
C----- Pozo  (f)
         if((x(j).ge.0.6d0).and.(x(j).le.0.7d0)) pot(j)=-1.25d5
         if((x(j).ge.0.6d0).and.(x(j).le.0.7d0)) v=-0.002
         write(12,*) x(j), v
      enddo
C-----------------------------------------------------------------------
C-------------------Inicializar la funcion de onda----------------------
c      k0=0.d0  ! apartado (a)
      k0=4.75d2   !para el resto
      norma=0.d0
      do j=1,M-1
         psi(j)=exp(-(x(j)-x0)**2/sig)*exp(i*k0*x(j))
         prob=abs(psi(j))**2
         norma=norma+prob
      enddo
      norma=dsqrt(norma)
      psi(M)=(0.d0,0.d0)
      psi(0)=(0.d0,0.d0)
      do j=1,M-1
         psi(j)=psi(j)/norma    !Normalizar la funci¢n de onda
         prob=abs(psi(j))**2
         write(11,*) x(j), prob
         write(13,*) x(j), prob
      enddo
      write(11,*)
      write(11,*)
      write(13,*)
      write(13,*)
C-------------------Evolucion temporal----------------------------------
      t=0.d0
      do k=1,N     !Se escribe la probabilidad N veces para hacer el gif
        do h=1,100  !se escribe cada 100 pasos temporales (5ú10E-5 s)
C------------------Calculo de las funciones auxiliares------------------
          omega(1)=-psi(1)+2.d0*(i*l+dx*dx*pot(1)+1.d0)*psi(1)-psi(0)
          e(1)=2.d0*(1.d0+dx*dx*pot(1)-i*l)
          f(1)=omega(1)
C
          do j=2,M-1
             omega(j)=-psi(j+1)-psi(j-1)+2.d0*(i*l+dx*dx*pot(j)
     &        +1.d0)*psi(j)
             e(j)=2.d0*(1.d0+dx*dx*pot(j)-i*l)-1.d0/e(j-1)
             f(j)=omega(j)+f(j-1)/e(j-1)
          enddo
C
          psi(M-1)=-f(M-1)/e(M-1)
          do j=M-2,1,-1
             psi(j)=(psi(j+1)-f(j))/e(j)
          enddo
          norma=0.d0
          do j=1,M-1
             prob=abs(psi(j))**2
             norma=norma+prob
          enddo
          norma=dsqrt(norma)
          do j=0,M
             psi(j)=psi(j)/norma
          enddo
          t=t+dt
        enddo
        do j=0,M
           prob=abs(psi(j))**2
           write(11,*) x(j), prob
        enddo
        write(11,*)
        write(11,*)
        resto=mod(k,20)
        if(resto.eq.0) then
           write(13,*) '#', t, 'seg'
           do j=0,M
              prob=abs(psi(j))**2
              write(13,*) x(j), prob
           enddo
           write(13,*)
           write(13,*)
        endif
      enddo
C
      deallocate(omega, f, e, psi)
      deallocate(pot, x)
C
      close(11)
      close(12)
      close(13)
C
      stop
      end

