##MODULE MINE!
##
##DESCRIPTION
## A package that computes the stuff I want to
##Macros
$define EXTRA1 3
$define EXTRA2 2
$define CLOSE 4
$define CLOSE2 20
$define Increase_Digits  1/16
$define	Initial_tstep  1/8
$define TRY_CCQUAD  11 # higher accuracy in Int than this -> use _CCquad
$define	ERR1  `Newton method failed, using more Digits...`
$define mindist ComputationalGeometry[ClosestPointPair] 
$define voronoiDual ComputationalGeometry[DelaunayTriangulation]
$define searchList ListTools[Search]
$define reverseList ListTools[Reverse]
$define removeReduancy ListTools[MakeUnique]
$define cross LinearAlgebra[CrossProduct]
$define Id LinearAlgebra[IdentityMatrix]
$define zero LinearAlgebra[ZeroMatrix]
$define prod LinearAlgebra[Multiply]
$define det LinearAlgebra[Determinant]
$define sm LinearAlgebra[SubMatrix]
$define norm LinearAlgebra[Norm]
$define tp LinearAlgebra[Transpose]
$define inv LinearAlgebra[MatrixInverse]
$define linsolve LinearAlgebra[LinearSolve]
$define diag LinearAlgebra[Diagonal]
$define	delcols LinearAlgebra[DeleteColumn]
$define	coldim LinearAlgebra[ColumnDimension]
$define col LinearAlgebra[Column]
$define	row LinearAlgebra[Row]
$define	delrows LinearAlgebra[DeleteRow]
$define	rank LinearAlgebra[Rank]
$define	rowdim LinearAlgebra[RowDimension]
$define	rowop LinearAlgebra[RowOperation]
$define	colop LinearAlgebra[ColumnOperation]
$define eq LinearAlgebra[Equal]
$define	perms combinat[permute]
$define	getSubObjects combinat[choose]
$define mcoef combinat[multinomial]

KPTools := module() 
    option package; 
    #local Pmatrix, Contourintegrate,Pathintegral,Partintegral,integrate,integrand2,Acy,Permuted,Whereis,Putfirst,Analyticcontinuation, Continue, Der;
    export waveConstants,waveConstantsFromDiffs,Ksol,KPgSol, Monodromy, Homology,  PeriodMatrix, checkWCError, HBF, getDC;    
    global A1, C,p,q;
    local getPermutation, getWaveConstants, getU, thetaDet, generateEquations,processEquations, computeThetaConst,
    sortComplex,AnalyticContinuation,Monodromy_Processed,Linearmonodromy,Generalmonodromy,
    ProcessPoly,ComputeBoundary, getParameters,
    ProduceVoronoi, ProduceHurwitz,ProduceCycle,GenerateMonodromy,check_cstruct,Continue,Der,
    Distsort, Canonicalbasis,Frobeniustransform,Intersectionmatrix,Putfirst,Whereis,Tretkofftable,
    Findcycle, Listordered, Tretkofflist, Homologybasis, Makecycle, Smallest, Canbasis, Simplifycycle,
    Compresscycle, Reformcycle,Pmatrix,Contourintegrate,Pathintegral,Partintegral,integrate,integrand2,Acy;

    
    
    
$include "Monodromy.mm"
$include "AnalyticContinuation.mm"	
$include "sortComplex.mm"
$include "getPermutation.mm"
$include "Homology.mm"
$include "PeriodMatrix.mm"
$include "KPSols.mm"
$include "KPSolsTaylor.mm"
$include "KPErrorTools.mm"
	

end module:
#savelib( 'KPTools' );