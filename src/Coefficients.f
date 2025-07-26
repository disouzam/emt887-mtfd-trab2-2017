*******************************************************
* FILE: Coefficients.f
* Author: Dickson Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      module Coefficients
        use Properties

        implicit none

      contains

      subroutine coeff2D(Aw,Ae,An,As,Ap,b,H,Tf,Delta_t,stage)



* Calculates matrix coefficients for 2D problems
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: stage

*       Stores heat transfer coefficient for all faces
*       H(1): left face
*       H(2): right face
*       H(3): bottom face
*       H(4): top face
        double precision :: H(1:4,1:10)

*       Stores temperature of fluid for all faces in order to calculate the coefficients
*       T(1): left face
*       T(2): right face
*       T(3): bottom face
*       T(4): top face
        double precision :: Tf(1:4,1:10)

        double precision :: Delta_t

        double precision :: Np(1:999,1:999,1:2), NI(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        double precision :: dens
        common /density/ dens

        double precision :: DelMinus, DelPlus, Ti, Tip1, Ap0, Cp, hcalc

        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J, Nx,Ny

        double precision :: Del(1:999,1:999,1:2), Delta(1:999,1:999,1:2)
        common / mesh/ Del, Delta

        double precision :: Aw(1:999,1:999), Ae(1:999,1:999)
        double precision :: An(1:999,1:999), As(1:999,1:999)
        double precision :: Ap(1:999,1:999), b(1:999,1:999)

        character(LEN=50) :: section

        double precision :: dummy, maxDumm

        Nx = N(1)
        Ny = N(2)

        dens = 7934.115D0

        maxDumm = 1.0D4

*        if (debugmode .EQV. .TRUE.) then
*          print*,"Density = ", dens
*        end if

        section = "Left Face"
*       Equations for points at left face (excluding corners at top and bottom)
        I = 1
        do J = 2, Ny - 1, 1

          Aw(I,J) = 0.0D0
          dummy = Aw(I,J)

          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)
          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)
          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)
          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)
          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)
          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)
          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          Cp = Cp_T(Tnew(I,J))
          hcalc = h_Total(H(1,stage),Tf(1,stage),Tnew(I,J))
          Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,2)
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + hcalc * Delta(I,J,2) * Tf(1,stage)
        end do

        section = "Right Face"
*       Equations for points at right face (excluding corners at top and bottom)
        I = Nx
        do J = 2, Ny - 1, 1

          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)
          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)
          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          Ae(I,J) = 0.0D0
          dummy = Ae(I,J)

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)
          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)
          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)
          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)
          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          Cp = Cp_T(Tnew(I,J))
          hcalc = h_Total(H(2,stage),Tf(2,stage),Tnew(I,J))
          Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,2)
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + hcalc * Delta(I,J,2) * Tf(2,stage)
        end do


        section = "Bottom Face"
*       Equations for points at bottom face (excluding left and right corners)
        J = 1
        do I = 2, Nx - 1, 1

          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)
          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)
          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)
          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)
          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          As(I,J) = 0.0D0
          dummy = As(I,J)

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)
          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)
          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          Cp = Cp_T(Tnew(I,J))
          hcalc = h_Total(H(3,stage),Tf(3,stage),Tnew(I,J))
          Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,1)
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + hcalc * Delta(I,J,1) * Tf(3,stage)
        end do


        section = "Top Face"
*       Equations for points at top face (excluding left and right corners)
        J = Ny
        do I = 2, Nx - 1, 1

          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)
          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)
          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)
          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)
          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            read*
          end if

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)
          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)
          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          An(I,J) = 0.0D0
          dummy = An(I,J)

          Cp = Cp_T(Tnew(I,J))
          hcalc = h_Total(H(4,stage),Tf(4,stage),Tnew(I,J))
          Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,1)
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J)  +  hcalc * Delta(I,J,1) * Tf(4,stage)
        end do

        section = "Top Left Corner"
        J = Ny
*       Equation for top left corner
        do I = 1, 1, 1
          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)
          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)
          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)
          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)
          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          Aw(I,J) = 0.0D0
          dummy = Aw(I,J)
          An(I,J) = 0.0D0
          dummy = An(I,J)

          Cp = Cp_T(Tnew(I,J))
          Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
          Ap(I,J) = Ap(I,J) + h_Total(H(1,stage),Tf(1,stage),Tnew(I,J)) * Delta(I,J,2)
          Ap(I,J) = Ap(I,J) + h_Total(H(4,stage),Tf(4,stage),Tnew(I,J)) * Delta(I,J,1)
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + h_Total(H(1,stage),Tf(1,stage),Tnew(I,J)) * Delta(I,J,2) * Tf(1,stage)
          b(I,J) = b(I,J) + h_Total(H(4,stage),Tf(4,stage),Tnew(I,J)) * Delta(I,J,1) * Tf(4,stage)
        end do

        section = "Top Right Corner"
        J = Ny
*       Equation for top right corner
        do I = Nx, Nx, 1
          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)
          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)
          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)
          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)
          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          Ae(I,J) = 0.0D0
          dummy = Ae(I,J)
          An(I,J) = 0.0D0
          dummy = An(I,J)

          Cp = Cp_T(Tnew(I,J))
          Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
          Ap(I,J) = Ap(I,J) + h_Total(H(2,stage),Tf(2,stage),Tnew(I,J)) * Delta(I,J,2)
          Ap(I,J) = Ap(I,J) + h_Total(H(4,stage),Tf(4,stage),Tnew(I,J)) * Delta(I,J,1)
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + h_Total(H(2,stage),Tf(2,stage),Tnew(I,J)) * Delta(I,J,2) * Tf(2,stage)
          b(I,J) = b(I,J) + h_Total(H(4,stage),Tf(4,stage),Tnew(I,J)) * Delta(I,J,1) * Tf(4,stage)
        end do

        section = "Bottom Left Corner"
        J = 1
*       Equation for bottom left corner
        do I = 1, 1, 1
          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)
          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)
          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)
          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)
          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          As(I,J) = 0.0D0
          dummy = As(I,J)
          Aw(I,J) = 0.0D0
          dummy = Aw(I,J)

          Cp = Cp_T(Tnew(I,J))
          Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
          Ap(I,J) = Ap(I,J) + h_Total(H(1,stage),Tf(1,stage),Tnew(I,J)) * Delta(I,J,2)
          Ap(I,J) = Ap(I,J) + h_Total(H(3,stage),Tf(3,stage),Tnew(I,J)) * Delta(I,J,1)
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + h_Total(H(1,stage),Tf(1,stage),Tnew(I,J)) * Delta(I,J,2) * Tf(1,stage)
          b(I,J) = b(I,J) + h_Total(H(3,stage),Tf(3,stage),Tnew(I,J)) * Delta(I,J,1) * Tf(3,stage)
        end do

        section = "Bottom Right Corner"
        J = 1
*       Equation for bottom right corner
        do I = Nx, Nx, 1
          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)
          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)
          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*            *              read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)
          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)
          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          As(I,J) = 0.0D0
          dummy = As(I,J)
          Ae(I,J) = 0.0D0
          dummy = Ae(I,J)

          Cp = Cp_T(Tnew(I,J))
          Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
          Ap(I,J) = Ap(I,J) + h_Total(H(2,stage),Tf(2,stage),Tnew(I,J)) * Delta(I,J,2)
          Ap(I,J) = Ap(I,J) + h_Total(H(3,stage),Tf(3,stage),Tnew(I,J)) * Delta(I,J,1)
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
*           read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + h_Total(H(2,stage),Tf(2,stage),Tnew(I,J)) * Delta(I,J,2) * Tf(2,stage)
          b(I,J) = b(I,J) + h_Total(H(3,stage),Tf(3,stage),Tnew(I,J)) * Delta(I,J,1) * Tf(3,stage)
        end do

        section = "Inner Points"
*       Equations for inner points
        do I = 2, Nx - 1, 1
          do J = 2, Ny - 1, 1

            DelMinus = Ni(I,J,1) - Np(I-1,J,1)
            DelPlus = Np(I,J,1) - Ni(I,J,1)
            Ti = Tnew(I-1,J)
            Tip1 = Tnew(I,J)
            Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)
            dummy = Aw(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
*             read*
            end if

            DelMinus = Ni(I+1,J,1) - Np(I,J,1)
            DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
            Ti = Tnew(I,J)
            Tip1 = Tnew(I+1,J)
            Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)
            dummy = Ae(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
*             read*
            end if

            DelMinus = Ni(I,J,2) - Np(I,J-1,2)
            DelPlus = Np(I,J,2) - Ni(I,J,2)
            Ti = Tnew(I,J-1)
            Tip1 = Tnew(I,J)
            As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)
            dummy = As(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
*             read*
            end if

            DelMinus = Ni(I,J+1,2) - Np(I,J,2)
            DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
            Ti = Tnew(I,J)
            Tip1 = Tnew(I,J+1)
            An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)
            dummy = An(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
*             read*
            end if

            Cp = Cp_T(Tnew(I,J))
            Ap0 = (dens * Cp * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

            Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap0
            dummy = Ap(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
*            read*
            end if

            b(I,J) = Ap0 * Told(I,J)

          end do
        end do


      end subroutine coeff2D

      end module Coefficients
