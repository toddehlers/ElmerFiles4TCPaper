! Get BMB as Function of Ice Thickness, GL Position and Distance.
FUNCTION GetBMB ( Model, nodenumber, x) RESULT(BMBOut)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: DepthSol,DistanceSol
   INTEGER, POINTER :: DepthPerm(:),DistancePerm(:),BMBPerm(:)
   REAL(kind=dp), POINTER :: DepthVal(:),DistanceVal(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   BMBOut, alpha, G,A,rho
   REAL(KIND=dp), ALLOCATABLE :: BMB0(:)
   LOGICAL :: FirstTime=.True., UnFoundFatal

   SAVE FirstTime
   SAVE BMB0
   !Get Depthwithout catching any error messages if fields don't exist
   DepthSol => VariableGet( Model % Variables, 'Depth',UnFoundFatal=UnFoundFatal)
   DepthPerm => DepthSol % Perm
   DepthVal => DepthSol % Values

   !Get Distance without catching any error messages if fields don't exist
   DistanceSol => VariableGet( Model % Variables, 'Distance',UnFoundFatal=UnFoundFatal)
   DistancePerm => DistanceSol % Perm
   DistanceVal => DistanceSol % Values

   alpha = 0.4
   G = 0.001
   A = 0.1
   rho = 1 - exp(-0.0001* DistanceVal(DistancePerm(nodenumber)))

   BMBOut = DepthVal(DepthPerm(nodenumber))**alpha*(rho*G+(1-rho)*A)

END FUNCTION GetBMB



