*******************************************************
* FILE: main.f
* Author: Dickson Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      program Trab2

*Modules
        use Coefficients
        use Geometry
        use Results
        use Solver
        use Properties

*Settings
        implicit none

*Variable declaration
        logical :: debugmode
        common /dbgMode/ debugmode

*       Lth: Lenght along each direction
*          1,2 indices refer to x and y respectively
        double precision :: Lth(1:2)
        common /domainsize/ Lth

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        integer :: maxX, maxY, Xi(9), Yi(9)
        common /storeLimits/ maxX, maxY, Xi, Yi

*       Del: array containing distance between neighboring nodes in all axis(to be calculated)
*       Index ranges from 1 to N(1) - 1 in X direction, from 1 to N(2) - 1 in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Delta: array containing distance between parallel interfaces of a volume element in all axis (to be calculated)
*       Index ranges from 1 to N(1)  in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        double precision :: Del(1:999,1:999,1:2),Delta(1:999,1:999,1:2)
        common /mesh/ Del, Delta

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        double precision :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

*       Tf stores the temperatures of the neighbor fluid. If the temperature at some boundary should be kept constant,
*         just set a very high heat transfer coefficient (in variable H)
*       H stores the heat transfer coefficient at 4 sides of domain
*       Tf(1): temperature of the fluid at left boundary (West)
*       H(1): heat transfer coefficient at left face(just convection, radiation is calculated according to node temperature)

*       Tf(2): temperature of the fluid at right boundary (East)
*       H(2): heat transfer coefficient at right face(just convection, radiation is calculated according to node temperature)

*       Tf(3): temperature of the fluid at bottom boundary (South)
*       H(3): heat transfer coefficient at bottom face(just convection, radiation is calculated according to node temperature)

*       Tf(4): temperature of the fluid at top boundary (North)
*       H(4): heat transfer coefficient at top face(just convection, radiation is calculated according to node temperature)
        double precision :: Tf(1:4,1:10), H(1:4,1:10)
        common /boundCond/ Tf, H

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg

*       First index: 9 nodes x 9 nodes = 81 selected positions to be saved
*       Second index: up to 1000 time steps
*       Third index:
*          1: curTime
*          2: X position
*          3: Y position
*          4: Temperature
*          5: Specific heat
*          6: Thermal conductivity
        double precision :: TmRs(1:81,0:5000,1:6)
        common /transient/ TmRs

        double precision :: simTime(1:10), curTime, num_steps(1:10), timeStep(1:10)
        integer :: num_stage, timeChoice
        common /time/ simTime, curTime, num_steps, num_stage, timeStep, timeChoice

        integer :: meth,relax_m, ADI
        double precision :: tol, alpha
        common /solMethod/ meth, relax_m, ADI, tol, alpha

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        integer :: stage

        double precision :: dens
*       Tolerance value for iterative methods - it defines the stopping criteria
*       It is an absolute value, i.e. it is intented to be the maximum value
*         for the difference between temperatures at successive iterations
*       alpha: A constant to be used as a factor for relaxation in Gauss-Seidel and TDMA methods

        common /density/ dens

*       In debug mode (=.TRUE.), some info is displayed in screen
        debugmode = .TRUE.

*       Header of program
        call printHeader()

* Instructions
        call readGeometry()

        print*
        print*
        read*

        call readInitTemp()

        print*
        print*
        read*

        call readTimeInfo()

        print*
        print*
        read*

        call readSolverSel()
        read*

        Tnew = Told
        call SaveResults(0,0)
        call defineLimits()
        call storeTimeSteps(0,0.0D0)

        curTimeStep = 0
        stage = 1
        totalIter = 0
        do while (stage .LE. num_stage)

          print*
          print*

          call readBoundCond(stage)

          print*
          print*
          print*,"==================================================="
          print*,"Solving Stage ", stage
          print*,"==================================================="
          print*

          call solveStage(stage)

          stage = stage + 1

        end do

        call saveTimeSteps(curTimeStep)
        print*,"Simulation has finished!"
        print*,"Total iterations: ", totalIter
        read*

      end program Trab2



      subroutine readGeometry()
*       Read geometry parameters and generate mesh through subroutine Grid2DUni
        use Geometry

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

*       Lth: Lenght along each direction
*          1,2 indices refer to x and y respectively
        double precision :: Lth(1:2)
        common /domainsize/ Lth

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        double precision :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        print*,"What are the dimensions of rectangular plate (in meters)?"
        print*,"Lenght in X direction (in meters):"
        if (debugmode .EQV. .TRUE.) then
          Lth(1) = 1.0D0
          print*,"Lth(1) = ", Lth(1), " m"
        else
          read*, Lth(1)
        end if

        print*
        print*,"Lenght in Y direction (in meters):"
        if (debugmode .EQV. .TRUE.) then
          Lth(2) = 0.20D0
          print*,"Lth(2) = ", Lth(2), " m"
        else
          read*, Lth(2)
        end if

        print*
        print*
        print*,"What are the number of nodes in each direction?"
        print*,"Nodes in X direction:"
        if (debugmode .EQV. .TRUE.) then
          N(1) = 241
          print*,"N(1) = ", N(1), " nodes"
        else
          read*, N(1)
        end if

        print*
        print*,"Nodes in Y direction:"
        if (debugmode .EQV. .TRUE.) then
          N(2) = 241
          print*,"N(2) = ", N(2), " nodes"
        else
          read*, N(2)
        end if

*       Generate a 2D grid for solution
        call Grid2DUni(Lth,Np,Ni)

      end subroutine readGeometry




      subroutine readBoundCond(stage)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: stage

*       Tf stores the temperatures of the neighbor fluid. If the temperature at some boundary should be kept constant,
*         just set a very high heat transfer coefficient (in variable H)
*       H stores the heat transfer coefficient at 4 sides of domain
*       Tf(1): temperature of the fluid at left boundary (West)
*       H(1): heat transfer coefficient at left face(just convection, radiation is calculated according to node temperature)

*       Tf(2): temperature of the fluid at right boundary (East)
*       H(2): heat transfer coefficient at right face(just convection, radiation is calculated according to node temperature)

*       Tf(3): temperature of the fluid at bottom boundary (South)
*       H(3): heat transfer coefficient at bottom face(just convection, radiation is calculated according to node temperature)

*       Tf(4): temperature of the fluid at top boundary (North)
*       H(4): heat transfer coefficient at top face(just convection, radiation is calculated according to node temperature)
        double precision :: Tf(1:4,1:10), H(1:4,1:10)
        common /boundCond/ Tf, H

        print*,"=================================================="
        print*,"         Boundary conditions for Stage ", stage
        print*,"=================================================="


        print*,"Specify temperatures at boundaries (or of the fluid where convection takes place) (in degrees Kelvin):"
        print*,"Temperature T1 (left border) (Kelvin):"
        if (debugmode .EQV. .TRUE.) then
          if (stage .EQ. 1) then
            Tf(1,stage) = 950.0D0 + 273.15D0
            print*,"Tf(1) = ", Tf(1,stage), " K"
          end if
          if (stage .EQ. 2) then
            Tf(1,stage) = 1150.0D0 + 273.15D0
            print*,"Tf(1) = ", Tf(1,stage), " K"
          end if
          if (stage .EQ. 3) then
            Tf(1,stage) = 1240.0D0 + 273.15D0
            print*,"Tf(1) = ", Tf(1,stage), " K"
          end if
          if (stage .EQ. 4) then
            Tf(1,stage) = 1230.0D0 + 273.15D0
            print*,"Tf(1) = ", Tf(1,stage), " K"
          end if
          if (stage .EQ. 5) then
            Tf(1,stage) = 1120.0D0 + 273.15D0
            print*,"Tf(1) = ", Tf(1,stage), " K"
          end if
        else
          read*, Tf(1,stage)
        end if

        print*
        print*,"Temperature T2 (right border) (Kelvin):"
        if (debugmode .EQV. .TRUE.) then
          if (stage .EQ. 1) then
            Tf(2,stage) = 950.0D0 + 273.15D0
            print*,"Tf(2) = ", Tf(2,stage), " K"
          end if
          if (stage .EQ. 2) then
            Tf(2,stage) = 1150.0D0 + 273.15D0
            print*,"Tf(2) = ", Tf(2,stage), " K"
          end if
          if (stage .EQ. 3) then
            Tf(2,stage) = 1240.0D0 + 273.15D0
            print*,"Tf(2) = ", Tf(2,stage), " K"
          end if
          if (stage .EQ. 4) then
            Tf(2,stage) = 1230.0D0 + 273.15D0
            print*,"Tf(2) = ", Tf(2,stage), " K"
          end if
          if (stage .EQ. 5) then
            Tf(2,stage) = 1120.0D0 + 273.15D0
            print*,"Tf(2) = ", Tf(2,stage), " K"
          end if
        else
          read*, Tf(2,stage)
        end if

        print*
        print*,"Temperature T3 (bottom border) (Kelvin):"
        if (debugmode .EQV. .TRUE.) then
          if (stage .EQ. 1) then
            Tf(3,stage) = 950.0D0 + 273.15D0
            print*,"Tf(3) = ", Tf(3,stage), " K"
          end if
          if (stage .EQ. 2) then
            Tf(3,stage) = 1150.0D0 + 273.15D0
            print*,"Tf(3) = ", Tf(3,stage), " K"
          end if
          if (stage .EQ. 3) then
            Tf(3,stage) = 1240.0D0 + 273.15D0
            print*,"Tf(3) = ", Tf(3,stage), " K"
          end if
          if (stage .EQ. 4) then
            Tf(3,stage) = 1230.0D0 + 273.15D0
            print*,"Tf(3) = ", Tf(3,stage), " K"
          end if
          if (stage .EQ. 5) then
            Tf(3,stage) = 1120.0D0 + 273.15D0
            print*,"Tf(3) = ", Tf(3,stage), " K"
          end if
        else
          read*, Tf(3,stage)
        end if

        print*
        print*,"Temperature T4 (top border) (Kelvin):"
        if (debugmode .EQV. .TRUE.) then
          if (stage .EQ. 1) then
            Tf(4,stage) = 950.0D0 + 273.15D0
            print*,"Tf(4) = ", Tf(4,stage), " K"
          end if
          if (stage .EQ. 2) then
            Tf(4,stage) = 1150.0D0 + 273.15D0
            print*,"Tf(4) = ", Tf(4,stage), " K"
          end if
          if (stage .EQ. 3) then
            Tf(4,stage) = 1240.0D0 + 273.15D0
            print*,"Tf(4) = ", Tf(4,stage), " K"
          end if
          if (stage .EQ. 4) then
            Tf(4,stage) = 1250.0D0 + 273.15D0
            print*,"Tf(4) = ", Tf(4,stage), " K"
          end if
          if (stage .EQ. 5) then
            Tf(4,stage) = 1160.0D0 + 273.15D0
            print*,"Tf(4) = ", Tf(4,stage), " K"
          end if
        else
          read*, Tf(4,stage)
        end if

        print*
        print*
        print*,"Specify heat transfer coefficient at boundaries (specify a very high value where constant temperature is desired):"
        print*,"Heat transfer coefficient H1 (left border) (W/sq-m*K):"
        if (debugmode .EQV. .TRUE.) then
          H(1,stage) = 50.0D0
          print*,"H(1) = ", H(1,stage), " W/sq-m*K"
        else
          read*, H(1,stage)
        end if

        print*
        print*,"Heat transfer coefficient H2 (right border) (W/sq-m*K):"
        if (debugmode .EQV. .TRUE.) then
          H(2,stage) = 50.0D0
          print*,"H(2) = ", H(2,stage), " W/sq-m*K"
        else
          read*, H(2,stage)
        end if

        print*
        print*,"Heat transfer coefficient H3 (bottom border) (W/sq-m*K):"
        if (debugmode .EQV. .TRUE.) then
          H(3,stage) = 8.0D0
          print*,"H(3) = ", H(3,stage), " W/sq-m*K"
        else
          read*, H(3,stage)
        end if

        print*
        print*,"Heat transfer coefficient H4 (top border) (W/sq-m*K):"
        if (debugmode .EQV. .TRUE.) then
          H(4,stage) = 150.0D0
          print*,"H(4) = ", H(4,stage), " W/sq-m*K"
        else
          read*, H(4,stage)
        end if

      end subroutine readBoundCond




      subroutine printHeader()
        implicit none

        print*,"==============================================================="
        print*,"        Program for calculation of a 2D thermal profile        "
        print*,"     in a rectangular plate in unsteady state conditions.      "
        print*,"                                                               "
        print*," Boundary conditions can be constant temperature in borders    "
        print*,"      T1: temperature at left border                           "
        print*,"      T2: temperature at right border                          "
        print*,"      T3: temperature at bottom border                         "
        print*,"      T4: temperature at top border                            "
        print*,"                                                               "
        print*," or a known heat transfer coefficient at borders               "
        print*,"      h1: heat transfer coefficient at left border             "
        print*,"      h2: heat transfer coefficient at right border            "
        print*,"      h3: heat transfer coefficient at bottom border           "
        print*,"      h4: heat transfer coefficient at top border              "
        print*,"                                                               "
        print*," or mixing boundary conditions, just setting a very high       "
        print*," heat transfer coefficient in faces where temperature should be"
        print*," constant at interface                                         "
        print*,"                                                               "
        print*,"==============================================================="
        print*," Author: Dickson Alves de Souza                                "
        print*,"                                                               "
        print*," Based on lectures by professor Roberto Parreiras Tavares      "
        print*,"                                                               "
        print*,"      and book Numerical Heat Transfer and Fluid Flow          "
        print*,"      by Suhas V. Patankar (1980)                              "
        print*,"                                                               "
        print*," Federal University of Minas Gerais                            "
        print*," October 19th, 2017                                            "
        print*,"==============================================================="
        print*,"                                                               "
        print*,"                                                               "
        print*,"                                                               "

      end subroutine printHeader




      subroutine readInitTemp()
        use Solver
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        print*,"Specify the initial temperature of slab (in Kelvin):"
        if (debugmode .EQV. .TRUE.) then
          Told = 300.0D0
          print*,"Initial slab temperature = ", Told(1,1), " K"
        else
          read*, Told
        end if

        Told(0,:) = 0.0D0
        Told(:,0) = 0.0D0
        Told(N(1)+1,:) = 0.0D0
        Told(:,N(2)+1) = 0.0D0

        call calcOldStats()

      end subroutine readInitTemp




      subroutine readTimeInfo()
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: I

        double precision :: simTime(1:10), curTime, num_steps(1:10), timeStep(1:10)
        integer :: num_stage, timeChoice
        common /time/ simTime, curTime, num_steps, num_stage, timeStep, timeChoice

        double precision :: dummy

        num_stage = 100

        do while (num_stage .GT. 10)
          print*,"What is the number of stages you want to simulate?"
          print*,"Each stage can have different boundary conditions (heat transfer coefficient and fluid temperatures)."
          if (debugmode .EQV. .TRUE.) then
            num_stage = 5
            print*,"Number of stages = ", num_stage, " stage(s)"
            print*,"Stage 1 - Non-firing zone"
            print*,"Stage 2 - Charging zone"
            print*,"Stage 3 - Pre-Heating zone"
            print*,"Stage 4 - Heating zone"
            print*,"Stage 5 - Soaking zone"
          else
            read*, num_stage

            if (num_stage .GT. 10) then
              print*," ERROR: The maximum number of stages allowed is 10."
              print*," Select a number between 1 and 10."
              print*
              print*
              print*
            end if
          end if
        end do

        do I = 1, num_stage, 1

          print*,"=================================================="
          print*,"             STAGE ", I,"/", num_stage
          print*,"=================================================="

          print*,"What is the duration for stage ", I, " (in seconds)?"
          if (debugmode .EQV. .TRUE.) then
            simTime(1) = 1140.0D0
            simTime(2) = 1920.0D0
            simTime(3) = 1920.0D0
            simTime(4) = 3420.0D0
            simTime(5) = 10000.0D0
            print*,"Simulation time = ", simTime(I), " seconds"
          else
            read*, simTime(I)
          end if

          print*
          print*

          if (I .EQ. 1) then

            print*,"Which one would you like to define: number of steps in a given stage or directly the time step?"
            print*,"1 for define number of steps."
            print*,"2 for define directly the time step."

            if (debugmode .EQV. .FALSE.) then
              read*, timeChoice
            else
              timeChoice = 2
            end if

          end if

          if (debugmode .EQV. .TRUE.) then
*            num_steps(1) = 20
*            timeStep(1) = simTime(1)/ num_steps(1)
*
*            num_steps(2) = 34
*            timeStep(2) = simTime(2)/ num_steps(2)
*
*            num_steps(3) = 34
*            timeStep(3) = simTime(3)/ num_steps(3)
*
*            num_steps(4) = 60
*            timeStep(4) = simTime(4)/ num_steps(4)
*
*            num_steps(5) = 48
*            timeStep(5) = simTime(5)/ num_steps(5)
*
*            print*,"Number of steps = ", num_steps(I)
*            print*,"Time step =", timeStep(I)
            dummy = 1.0D0
            timeStep(1)= dummy
            num_steps(1) = ceiling(simTime(1) / timeStep(1))

            timeStep(2)= dummy
            num_steps(2) = ceiling(simTime(2) / timeStep(2))

            timeStep(3)= dummy
            num_steps(3) = ceiling(simTime(3) / timeStep(3))

            timeStep(4)= dummy
            num_steps(4) = ceiling(simTime(4) / timeStep(4))

            timeStep(5)= dummy
            num_steps(5) = ceiling(simTime(5) / timeStep(5))

            print*,"Number of steps = ", num_steps(I)
            print*,"Time step =", timeStep(I)

          else
            if (timeChoice .EQ. 1) then

              print*,"What should be the number of time steps for stage ", I,":"
              read*, num_steps(I)
              timeStep(I) = simTime(I)/ num_steps(I)

            end if

            if (timeChoice .EQ. 2) then

              print*,"What should be the time step for stage ", I," in seconds:"
              read*, timeStep(I)
              num_steps(I) = ceiling(simTime(I) / timeStep(I))

            end if

          end if

        end do


      end subroutine readTimeInfo




      subroutine readSolverSel()

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: meth,relax_m, ADI
        double precision :: tol, alpha
        common /solMethod/ meth, relax_m, ADI, tol, alpha


        print*,"Select the solution method:"
        print*,"1 for Jacobi method"
        print*,"2 for Gauss-Seidel method"
        print*,"3 for TDMA (tri-diagonal matrix algorithm)"
        print*
        if (debugmode .EQV. .TRUE.) then
          meth = 3
          print*,"meth = ", meth
        else
          read*, meth
        end if

        if (meth .EQ. 1) then
          print*,"Jacobi method was selected as a solver."
        end if
        if (meth .EQ. 2) then
          print*,"Gauss-Seidel method was selected as a solver."
        end if
        if (meth .EQ. 3) then
          print*,"TDMA (tri-diagonal matrix algorithm) was selected as a solver."
        end if


        if (meth .EQ. 3) then
          print*,"Do you want to use ADI (alternating direction implicit) method?"
          print*,"1 for Yes"
          print*,"2 for No"
          print*
          if (debugmode .EQV. .TRUE.) then
            ADI = 1
            print*,"ADI = ", ADI
          else
            read*, ADI
          end if

        end if

        print*
        print*
        print*,"What tolerance should be considered?"
        print*,"Range: 1.0D-1 - 1.0D-8"
        if (debugmode .EQV. .TRUE.) then
          tol = 1.0D-3
          print*,"tol = ", tol
        else
          read*, tol

          if (tol .GT. 1.0D-1) then
            tol = 1.0D-1
          elseif (tol .LT. 1.0D-8) then
            tol = 1.0D-8
          end if
        end if

        print*
        print*
        print*,"Do you want to use a relaxation factor?"
        print*,"1 for Yes"
        print*,"2 for No"


        if (debugmode .EQV. .TRUE.) then
          alpha = 1.0D0

          if (alpha .EQ. 1.0D0) then
            print*,"No relaxation factor will be used"
          end if
          print*,"alpha = ", alpha
        else
          read*, relax_m
          print*
          if (relax_m .EQ. 1) then
            print*,"What relaxation factor do you want to use?"
            print*,"Range: 0.01 - 2.00"
            print*,"Be cautious: the chosen value could cause divergence of solution."
            print*
            read*,alpha

            if ((alpha .LT. 0.01) .OR. (alpha .GT. 2.0)) then
              print*,"You typed a wrong value for alpha. Default value (alpha = 1.0D0) will be considered, instead."
              alpha = 1.0D0
            end if

          else
            alpha = 1.0D0
          end if
        end if
      end subroutine readSolverSel




      subroutine solveStage(stage)
        use Solver

        implicit none

        integer :: stage

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        logical :: debugmode
        common /dbgMode/ debugmode

*       Tf stores the temperatures of the neighbor fluid. If the temperature at some boundary should be kept constant,
*         just set a very high heat transfer coefficient (in variable H)
*       H stores the heat transfer coefficient at 4 sides of domain
*       Tf(1): temperature of the fluid at left boundary (West)
*       H(1): heat transfer coefficient at left face(just convection, radiation is calculated according to node temperature)

*       Tf(2): temperature of the fluid at right boundary (East)
*       H(2): heat transfer coefficient at right face(just convection, radiation is calculated according to node temperature)

*       Tf(3): temperature of the fluid at bottom boundary (South)
*       H(3): heat transfer coefficient at bottom face(just convection, radiation is calculated according to node temperature)

*       Tf(4): temperature of the fluid at top boundary (North)
*       H(4): heat transfer coefficient at top face(just convection, radiation is calculated according to node temperature)
        double precision :: Tf(1:4,1:10), H(1:4,1:10)
        common /boundCond/ Tf, H

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        double precision :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        double precision :: simTime(1:10), curTime, num_steps(1:10), timeStep(1:10)
        integer :: num_stage, timeChoice
        common /time/ simTime, curTime, num_steps, num_stage, timeStep, timeChoice

        integer :: meth,relax_m, ADI
        double precision :: tol, alpha
        common /solMethod/ meth, relax_m, ADI, tol, alpha

        double precision :: Delta_t
        logical :: timeStepChanged

        double precision :: startTime, endTime
        integer :: I

        character(LEN=50) :: name

160     format(' ', A, I6)

        startTime = 0
        endTime = 0

        do I = 1, stage, 1
          endTime = endTime + simTime(I)
        end do

        startTime = endTime - simTime(stage)

        if (timeChoice .EQ. 1) then
          Delta_t = simTime(stage) / num_steps(stage)
        elseif (timeChoice .EQ. 2) then
          Delta_t = timeStep(stage)
        end if

        curTime = startTime + Delta_t

*        if (curTime .GT. endTime) then
*          curTime = endTime
*          Delta_t = endTime - startTime
*        end if

        timeStepChanged = .FALSE.

        do while (curTime .LE. endTime)

          print*
          print*,"--------------------------------------------------------------"
          print*
          print*
          print*,"CALCULATING TEMPERATURES FOR TIME = ", curTime, " seconds."
          print*
          print*
          print*,"--------------------------------------------------------------"

*          if (curTime .GT. 1100) then
*            print*, "Current time = ", curTime
*            read*
*          end if

*         Jacobi method
          if (meth .EQ. 1) then
            name = "Jacobi Method"
            call Jacobi(Tf,H,Delta_t,tol,alpha,stage)
            if (debugmode .EQV. .TRUE.) then
*              call print_res2D(name,Tnew,Np,curTime)
            end if

            print 160,"Jacobi method has terminated - Iterations: ", iter
            print*,"=========================================================================="
            if (debugmode .EQV. .TRUE.) then
*              read*
            end if
          end if

*         Gauss-Seidel method
          if (meth .EQ. 2) then
            name = "Gauss-Seidel Method"
            call GaSe2D(Tf,H,Delta_t, tol,alpha,stage)
            if (debugmode .EQV. .TRUE.) then
*              call print_res2D(name,Tnew,Np,curTime)
            end if

            print 160,"Gauss-Seidel method has terminated - Iterations: ", iter
            print*,"=========================================================================="
            if (debugmode .EQV. .TRUE.) then
*              read*
            end if
          end if

*         TDMA method
          if (meth .EQ. 3) then
            if (ADI .EQ. 1) then
              name = "TDMA with ADI"
              call TDMA2dADI(Tf,H,Delta_t,tol,alpha,stage)
              if (debugmode .EQV. .TRUE.) then
*                call print_res2D(name,Tnew,Np,curTime)
              end if

              print 160,"TDMA method with ADI has terminated - Iterations: ", iter
              print*,"=========================================================================="
              if (debugmode .EQV. .TRUE.) then
*                read*
              end if
            elseif (ADI .EQ. 2) then
              name = "TDMA Simple (Non-ADI)"
              call TDMA2D(Tf,H,Delta_t,tol,alpha,stage)
              if (debugmode .EQV. .TRUE.) then
*                call print_res2D(name,Tnew,Np,curTime)
              end if

              print 160,"TDMA method without ADI has terminated - Iterations: ", iter
              print*,"=========================================================================="
              if (debugmode .EQV. .TRUE.) then
*                read*
              end if
            end if
          end if

          curTimeStep = curTimeStep +1
          call storeTimeSteps(curTimeStep,curTime)

          if (mod(curtimeStep, 97) .EQ. 0) then
            call SaveResults(stage,totalIter)
          end if


          if (timeStepChanged .EQV. .FALSE.) then
            if ((endTime - curTime) .GT. 1.0D0) then
              if ((endTime - curTime) .LT. Delta_t) then
                Delta_t = endTime - curTime
                timeStepChanged = .TRUE.
              end if
            end if
          end if

          curTime = curTime + Delta_t

        end do

      end subroutine solveStage




      subroutine SaveResults(stage,totalIter)

        use Properties

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: stage, totalIter

        double precision :: StartTime, EndTime

*       Lth: Lenght along each direction
*          1,2 indices refer to x and y respectively
        double precision :: Lth(1:2)
        common /domainsize/ Lth

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        double precision :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

*       Tf stores the temperatures of the neighbor fluid. If the temperature at some boundary should be kept constant,
*         just set a very high heat transfer coefficient (in variable H)
*       H stores the heat transfer coefficient at 4 sides of domain
*       Tf(1): temperature of the fluid at left boundary (West)
*       H(1): heat transfer coefficient at left face(just convection, radiation is calculated according to node temperature)

*       Tf(2): temperature of the fluid at right boundary (East)
*       H(2): heat transfer coefficient at right face(just convection, radiation is calculated according to node temperature)

*       Tf(3): temperature of the fluid at bottom boundary (South)
*       H(3): heat transfer coefficient at bottom face(just convection, radiation is calculated according to node temperature)

*       Tf(4): temperature of the fluid at top boundary (North)
*       H(4): heat transfer coefficient at top face(just convection, radiation is calculated according to node temperature)
        double precision :: Tf(1:4,1:10), H(1:4,1:10)
        common /boundCond/ Tf, H

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg

        double precision :: simTime(1:10), curTime, num_steps(1:10), timeStep(1:10)
        integer :: num_stage, timeChoice
        common /time/ simTime, curTime, num_steps, num_stage, timeStep, timeChoice

        integer :: meth,relax_m, ADI
        double precision :: tol, alpha
        common /solMethod/ meth, relax_m, ADI, tol, alpha

        double precision :: Cp, k

        character(LEN = 100) :: filename
        write(filename,"(A8,I2,A2,F8.2,A4)") "RESULT-S", stage, "-T", curTime, ".DAT"
        filename = trim(filename)

        print*,"Saving results for this time step in a text file..."

*       Writing results in RESULT.DAT file, stored in the same folder as the program
100     format(' ', 2I12,5F12.4)
200     format(' ', 7A12)
300     format(' ', A, F12.4)
400     format(' ', A, I12)
500     format(' ', A, I3, A, I3, A)
600     format(' ', A, 1F12.4)
700     format(' ', A, 1F12.4)
        open(10,file=filename,status='UNKNOWN')

        write(10,*)"===================================================================================="
        write(10,*)"              Transient state solution for a 2D rectangular plate                   "
        write(10,*)"                   for applications in steel reheating                              "
        write(10,*)"                e.g. Reheating furnaces for rolling processing                      "
        write(10,*)"===================================================================================="
        write(10,*)"Author: Dickson Alves de Souza                                                      "
        write(10,*)"                                                                                    "
        write(10,*)"Based on lectures by professor Roberto Parreiras Tavares                            "
        write(10,*)"                                                                                    "
        write(10,*)"and book Numerical Heat Transfer and Fluid Flow                                     "
        write(10,*)"by Suhas V. Patankar (1980)                                                         "
        write(10,*)"                                                                                    "
        write(10,*)"Federal University of Minas Gerais (UFMG)                                           "
        write(10,*)"October 26th, 2017                                                                  "
        write(10,*)"===================================================================================="
        write(10,*)"                             Input parameters:                                      "
        write(10,*)"                                                                                    "
        write(10,300)"Lenght in X direction (meters):     ",Lth(1)
        write(10,300)"Lenght in Y direction (meters):     ",Lth(2)
        write(10,400)"Nodes in X direction:               ",N(1)
        write(10,400)"Nodes in Y direction:               ",N(2)


        if (stage .NE. 0) then

          do I = 1, stage, 1
            startTime = 0
            endTime = 0

            do J = 1, I, 1
              endTime = endTime + simTime(J)
            end do

            startTime = endTime - simTime(I)

            write(10,*)"===================================================================================="
            write(10,*)"                                                                                    "
            write(10,500)"                      STAGE              ", I,"/",num_stage ,"                               "
            write(10,600)"                      START Time:   ", startTime
            write(10,600)"                      END Time:     ", endTime
            write(10,400)"                      Iterations:     ", totalIter
            write(10,*)"                                                                                    "
            write(10,*)"===================================================================================="
            write(10,300)"Temperature of fluid at left face (T left):          ", Tf(1,I)
            write(10,300)"Temperature of fluid at right face (T right):        ", Tf(2,I)
            write(10,300)"Temperature of fluid at bottom face (T bottom):      ", Tf(3,I)
            write(10,300)"Temperature of fluid at top face (T top):            ", Tf(4,I)
            write(10,*)"------------------------------------------------------------------------------------"
            write(10,300)"Heat transfer coefficient at left face (H left):     ", H(1,I)
            write(10,300)"Heat transfer coefficient at right face (H right):   ", H(2,I)
            write(10,300)"Heat transfer coefficient at bottom face (H bottom): ", H(3,I)
            write(10,300)"Heat transfer coefficient at top face (H top):       ", H(4,I)
            write(10,*)"===================================================================================="
            write(10,*)""
            write(10,*)""
          end do

        end if

        if (meth .EQ. 1) then
          write(10,*) "SOLUTION of LINEAR SYSTEM: Jacobi method"
        elseif (meth .EQ. 2) then
          write(10,*) "SOLUTION of LINEAR SYSTEM: Gauss-Seidel method"
        elseif (meth .EQ. 3) then

          write(10,*) "SOLUTION of LINEAR SYSTEM: TDMA (tri-diagonal matrix algorithm)"

          if (ADI .EQ. 1) then
            write(10,*) "     ADI approach employed in solution."
          elseif (ADI .EQ. 2) then
            write(10,*) "     ADI approach NOT employed in solution."
          end if

        end if

        if (relax_m .EQ. 1) then
          write(10,*)"    Relaxation factor: ", alpha
        elseif (relax_m .EQ. 2) then
          write(10,*)"    No relaxation applied to solution."
        end if

        write(10,700)"Tolerance: ", tol
        write(10,*)""

        write(10,*)"===================================================================================="
        write(10,*)"                   Calculation Results                                              "
        write(10,*)"Current time: ", curTime," seconds                                                  "
        write(10,*)"Maximum Temperature : ", Tmax
        write(10,*)"Minimum Temperature : ", Tmin
        write(10,*)"Average(algebraic, non-weigthed) Temperature : ", Tavg
        write(10,*)"===================================================================================="

        write(10,200)"I","J","X(m)","Y(m)","T(K)","Cp(J/kg.K)", "K(W/sq-m.K)"

        write(10,*)"===================================================================================="

        do I = 1, N(1), 1
          do J = 1, N(2), 1
            Cp = Cp_T(Tnew(I,J))
            k = k_Tn(Tnew(I,J))
            write(10,100)I,J,Np(I,J,1),Np(I,J,2),Tnew(I,J),Cp,k
          end do
        end do

        endfile 10
        close(10,status='KEEP')


        print*,"Results were successfully saved. Next time step..."
        print*," "
        print*," "
        print*," "
        print*," "

      end subroutine SaveResults




      subroutine storeTimeSteps(curTimeStep,curTime)
        use Properties

        implicit none

        integer :: curTimeStep
        double precision :: curTime

*       First index: 9 nodes x 9 nodes = 81 selected positions to be saved
*       Second index: up to 1000 time steps
*       Third index:
*          1: curTime
*          2: X position
*          3: Y position
*          4: Temperature
*          5: Specific heat
*          6: Thermal conductivity
        double precision :: TmRs(1:81,0:5000,1:6)
        common /transient/ TmRs

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J, Npos

        integer :: maxX, maxY, Xi(9), Yi(9)
        common /storeLimits/ maxX, maxY, Xi, Yi

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        double precision :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        Npos = 1

        do I = 1, maxX, 1
          do J = 1, maxY, 1

*            print*,Xi(I)
*            print*,Yi(J)
            TmRs(Npos,curTimeStep,1) = curTime
            TmRs(Npos,curTimeStep,2) = Np(Xi(I),Yi(J),1)
            TmRs(Npos,curTimeStep,3) = Np(Xi(I),Yi(J),2)
            TmRs(Npos,curTimeStep,4) = Tnew(Xi(I),Yi(J))
            TmRs(Npos,curTimeStep,5) = Cp_T(Tnew(Xi(I),Yi(J)))
            TmRs(Npos,curTimeStep,6) = k_Tn(Tnew(Xi(I),Yi(J)))

            Npos = Npos + 1
          end do
        end do

      end subroutine storeTimeSteps




      subroutine defineLimits()
        implicit none

        integer :: rpt, I, J
        integer :: maxX, maxY, Xi(9), Yi(9)
        common /storeLimits/ maxX, maxY, Xi, Yi

        integer :: N(1:2)
        common /gridsize/ N

        Xi(1) = 1
        Xi(2) = ceiling(1.0D0 * float(N(1)) / 8.0D0)
        Xi(3) = ceiling(2.0D0 * float(N(1)) / 8.0D0)
        Xi(4) = ceiling(3.0D0 * float(N(1)) / 8.0D0)
        Xi(5) = ceiling(4.0D0 * float(N(1)) / 8.0D0)
        Xi(6) = ceiling(5.0D0 * float(N(1)) / 8.0D0)
        Xi(7) = ceiling(6.0D0 * float(N(1)) / 8.0D0)
        Xi(8) = ceiling(7.0D0 * float(N(1)) / 8.0D0)
        Xi(9) = N(1)

*        print*,"X"

*        do I = 1, 9, 1
*          print*,I," = ",Xi(I)
*        end do

        maxX = 9
        rpt = 9

        do while (rpt .GT. 0)
          do I = 1, maxX - 2, 1

            if (I .EQ. 1) then
              rpt = 0
            end if

            if (Xi(I+1) .EQ. Xi(I)) then
              Xi(I+1) = Xi(I+2)
              rpt = rpt + 1
            end if

          end do

          if (rpt .GT. 0) then
            maxX = maxX - 1
          end if

        end do

        if (Xi(maxX - 1) .EQ. Xi(maxX)) then
          maxX = maxX - 1
        end if



*        print*
*        print*

*        do I = 1, 9, 1
*          print*,I," = ",Xi(I)
*        end do

        Yi(1) = 1
        Yi(2) = ceiling(1.0D0 * float(N(2)) / 8.0D0)
        Yi(3) = ceiling(2.0D0 * float(N(2)) / 8.0D0)
        Yi(4) = ceiling(3.0D0 * float(N(2)) / 8.0D0)
        Yi(5) = ceiling(4.0D0 * float(N(2)) / 8.0D0)
        Yi(6) = ceiling(5.0D0 * float(N(2)) / 8.0D0)
        Yi(7) = ceiling(6.0D0 * float(N(2)) / 8.0D0)
        Yi(8) = ceiling(7.0D0 * float(N(2)) / 8.0D0)
        Yi(9) = N(2)

*        print*
*        print*
*        print*,"Y"

*        do J = 1, 9, 1
*          print*,J," = ",Yi(J)
*        end do

        maxY = 9
        rpt = 9

        do while (rpt .GT. 0)
          do J = 1, maxY - 2, 1

            if (J .EQ. 1) then
              rpt = 0
            end if

            if (Yi(J+1) .EQ. Yi(J)) then
              Yi(J+1) = Yi(J+2)
              rpt = rpt + 1
            end if

          end do

          if (rpt .GT. 0) then
            maxY = maxY - 1
          end if

        end do

        if (Yi(maxY - 1) .EQ. Yi(maxY)) then
          maxY = maxY - 1
        end if



*        print*
*        print*
*
*        do J = 1, 9, 1
*          print*,J," = ",Yi(J)
*        end do



      end subroutine defineLimits




      subroutine saveTimeSteps(curTimeStep)
        implicit none

        integer :: curTimeStep, I, J, K, Npos

*       First index: 9 nodes x 9 nodes = 81 selected positions to be saved
*       Second index: up to 1000 time steps
*       Third index:
*          1: curTime
*          2: X position
*          3: Y position
*          4: Temperature
*          5: Specific heat
*          6: Thermal conductivity
        double precision :: TmRs(1:81,0:5000,1:6)
        common /transient/ TmRs

        integer :: maxX, maxY, Xi(9), Yi(9)
        common /storeLimits/ maxX, maxY, Xi, Yi

*       Lth: Lenght along each direction
*          1,2 indices refer to x and y respectively
        double precision :: Lth(1:2)
        common /domainsize/ Lth

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        double precision :: simTime(1:10), curTime, num_steps(1:10), timeStep(1:10)
        integer :: num_stage, timeChoice
        common /time/ simTime, curTime, num_steps, num_stage, timeStep, timeChoice

        double precision :: avg, total
        logical :: constVal

        character(LEN = 100) :: filename


120     format(' ', 6A12)
220     format(' ', 6F12.4)
320     format(' ', A,I2)
420     format(' ', A9,I3,A3,I3,A1,F6.2,A6,A4)

        total = 0
        do I = 1, num_stage, 1
          total = total + timeStep(I)
        end do

        avg = total / float(num_stage)

        constVal = .TRUE.
        do I = 1, num_stage - 1, 1
          if (timeStep(I) .NE. timeStep(I+1)) then
            constVal = .FALSE.
          end if
        end do

        if (constVal .EQV. .TRUE.) then
          write(filename,420) "TempEvol_", N(1), "_X_", N(2), "_", avg,"_s_CNT",".DAT"
        else
          write(filename,420) "TempEvol_", N(1), "_X_", N(2), "_", avg,"_s_AVG",".DAT"
        end if

        filename = trim(filename)

        open(20,file=filename,status='UNKNOWN')

        Npos = 1
        write(20,*)"========================================================================"
        write(20,*)"Mesh information:"
        write(20,*)"Nx =", N(1), "    -    ", "Ny = ", N(2)
        write(20,*)" "
        write(20,*)"Slab dimensions:"
        write(20,*)"X(m) = ", Lth(1), "    -    ", "Y(m) = ", Lth(2)

        write(20,*)"                                                                        "
        write(20,*)"                                                                        "
        write(20,*)"========================================================================"

        do I = 1, maxX, 1
          do J = 1, maxY, 1
            write (20,320)"Npos = ", Npos
            write(20,120)"Time(s)","X(m)","Y(m)","Temp(K)","Cp(J/kg.K)","K(W/sq-m.K)"

            do K = 0, curTimeStep, 1
              write(20,220)TmRs(Npos,K,1),TmRs(Npos,K,2),TmRs(Npos,K,3),TmRs(Npos,K,4),TmRs(Npos,K,5),TmRs(Npos,K,6)
            end do

            Npos = Npos + 1

            write(20,*)"                                                                        "
            write(20,*)"                                                                        "
            write(20,*)"========================================================================"

          end do
        end do

        endfile 20
        close(20,status='KEEP')

      end subroutine saveTimeSteps
