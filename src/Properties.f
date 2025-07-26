*******************************************************
* FILE: Properties.f
* Author: Dickson Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      module Properties
            implicit none

      contains

      function k_Tn(Tn)
*       Implementation for AISI 1018 steel according to
*         KIM, Y.;FAROUK, B.; KEVERIAN, J. A mathematical model for Thermal Analysis of Thin Strip Casting of Low Carbon Steel
*         Journal of Engineering for Industry, February 1991, v. 113
        implicit none

        double precision :: k_Tn
        double precision :: Tn

*       Correct formula
        k_Tn = 65.214D0
        k_Tn = k_Tn + 2.715D-2 * Tn
        k_Tn = k_Tn - 1.628D-4 * Tn * Tn
        k_Tn = k_Tn + 1.39D-7 * Tn * Tn * Tn
        k_Tn = k_Tn - 3.041D-11 * Tn * Tn * Tn * Tn

*        k_Tn = 65.0D0

*       Testing formula
*        k_Tn = 425.0D0
*        k_Tn = k_Tn - 0.3D0 * Tn
*        k_Tn = k_Tn - 1.628D-4 * Tn * Tn
*        k_Tn = k_Tn + 1.39D-7 * Tn * Tn * Tn
*        k_Tn = k_Tn - 3.041D-11 * Tn * Tn * Tn * Tn

      end function





      function k_Inter(Ti,Tip1,DelMinus,DelPlus)
        implicit none

        double precision :: k_Inter
        double precision :: Ti, Tip1,DelMinus,DelPlus
        double precision :: dummy,DelTotal

        DelTotal = DelMinus + DelPlus

        dummy = (DelMinus / DelTotal) * (1.0D0 / k_Tn(Ti))
        dummy = dummy + (DelPlus / DelTotal) * (1.0D0 / k_Tn(Tip1))

        dummy = 1.0D0 / dummy

        k_Inter = dummy

      end function k_Inter





      function Cp_T(Tn)

*       Implementation for AISI 1018 steel according to
*       KIM, Y.;FAROUK, B.; KEVERIAN, J. A mathematical model for Thermal Analysis of Thin Strip Casting of Low Carbon Steel
*       Journal of Engineering for Industry, February 1991, v. 113
        implicit none

        double precision :: Cp_T
        double precision :: Tn

        Cp_T = 1000.0D0

        if (Tn .LT. 1033.0D0) then

          Cp_T = 2.368D0
          Cp_T = Cp_T - 1.492D-2 * Tn
          Cp_T = Cp_T + 4.107D-5 * Tn * Tn
          Cp_T = Cp_T - 4.696D-8 * Tn * Tn * Tn
          Cp_T = Cp_T + 1.953D-11 * Tn * Tn * Tn * Tn
          Cp_T = 1000.0D0 * Cp_T ! Conversion from KJ/Kg.K to J/Kg.K

        elseif (Tn .LT. 1200.0D0) then

          Cp_T = 7.802D0
          Cp_T = Cp_T - 5.278D-3 * Tn
          Cp_T = Cp_T - 3.676D-6 * Tn * Tn
          Cp_T = Cp_T + 1.388D-9 * Tn * Tn * Tn
          Cp_T = Cp_T + 1.031D-12 * Tn * Tn * Tn * Tn
          Cp_T = 1000.0D0 * Cp_T ! Conversion from KJ/Kg.K to J/Kg.K

        elseif (Tn .LT. 1776.0D0) then

          Cp_T = 0.703D0
          Cp_T = 1000.0D0 * Cp_T ! Conversion from KJ/Kg.K to J/Kg.K

        end if

      end function Cp_T




      function h_Total(h,Tf,Tp)
        implicit none

        double precision :: h_Total
        double precision :: h,Tf, Tp
        double precision :: sigma, eps, p1, p2

        sigma = 5.67D-8
        eps = emissiv(Tp)
        p1 = Tf + Tp
        p2 = Tf * Tf + Tp * Tp

        h_Total = h + sigma * eps * p1 * p2

      end function h_Total




      function emissiv(Tp)
        implicit none

        double precision :: emissiv
        double precision :: Tp
        double precision :: y1, y2, x1, x2, m

        if (Tp .LT. 590.0D0) then
          emissiv = 0.69D0
        elseif (Tp .LT. 755.0D0) then
          y1 = 0.69D0
          x1 = 590.0D0
          y2 = 0.72D0
          x2 = 755.0D0
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 920.0D0) then
          y1 = 0.72D0
          x1 = 755.0D0
          y2 = 0.76D0
          x2 = 920.0D0
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 1090.0D0) then
          y1 = 0.76D0
          x1 = 920.0D0
          y2 = 0.79D0
          x2 = 1090.0D0
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 1255.0D0) then
          y1 = 0.79D0
          x1 = 1090.0D0
          y2 = 0.82D0
          x2 = 1255.0D0
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        else
          emissiv = 0.82D0
        end if


      end function emissiv




      end module Properties
