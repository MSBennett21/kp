##MODULE MINE!
##
##DESCRIPTION
## A package that computes the stuff I want to
##Macros
$define EXTRA1 3
$define EXTRA2 2
$define DEFAULT_GLOBAL_LVL 12 #see david baily paper
$define DEFAULT_BASE_LVL 1
$define mindist ComputationalGeometry[ClosestPointPair] 
$define voronoiDual ComputationalGeometry[DelaunayTriangulation]
$define searchList ListTools[Search]
$define reverseList ListTools[Reverse]
$define mkUnique ListTools[MakeUnique]
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
$define input combinat[multinomial]
$define polySolve SolveTools[Polynomial]
$define polySysSolve SolveTools[PolynomialSystem]
$define getGenus algcurves[genus]
$define pProd GroupTheory[PermProduct]
$define pInv GroupTheory[PermInverse]

kptools := module() 
    option package; 
    #local Pmatrix, Contourintegrate,Pathintegral,Partintegral,integrate,integrand2,Acy,Permuted,Whereis,Putfirst,Analyticcontinuation, Continue, Der;
    export sortLex,sortArg,RealKPSolLab,intAEdge, waveConstants,lineInt, fgParameters, RCV, Monodromy, Homology,  periodMatrix, checkWCError, HBF, getDC, KPgSol; 
    local global_lvl := DEFAULT_GLOBAL_LVL, base_lvl := DEFAULT_BASE_LVL, ui_digits_origin := Digits, 
        ui_digits_processed := max(Digits, 10),   ui_digits_high := 2 * ui_digits_processed, 
        ui_digits_low := max(ui_digits_origin-3,7), s_coord, t_cord, global_abscissas_count := 0, matchUpAC, populateAbs, 
        NewtonForRational, RCV2, inf_dist_between_pts := 1, base_point := NULL,  base_preimage:= NULL, num_problem_pts := 0, 
        digits_lb := 10, plotter:={}, V:={},  E:={}, F:=table(), edgeData := table(), vertexData:=table(), pathData:=table(),
        edgeSet:={}, problemPoints, singPoints, boundaryPoints, branchPoints, diffBasis,  genusoff, num_sheets,  Pmatrix,dxx,dxxx, 
        intData, pathData, homBasis,edgeSet, 
        #procedures
        #Monodromy 
        initialize1, mono1, initialization,processPoints, linearMonodromy, generalMonodromy,ABS,intAroundBpNTimes,
        produceVoronoi, produceHurwitz, generateMonodromy, intersectionOfBLines,check_cstruct, 
        colors:=["Blue","Orange","SkyBlue","CornflowerBlue","Coral","DarkBlue","Tomato","Cyan","DarkSalmon"],
        #misc Comps
        removeVerts,minStep:=1/2^BaseLvl,
        getPermutation,permuteList,permuteFibs, wIsToRightOfz, nPtBoundary,maximalDistWRT,scrubToList,getSiteIndex,
        isNormalForm, clearCanvas, toFloat,scrub,complexToPt,complexLToPtL, ptLToComplexL,constructEdge,
        updatePlotter,intProcessing,lineIntersection,edgeProcessing,newtonY,DK,
        reverseEdgeAC,  lvlUp,intAroundBp,extractFib, getPuedoSite,getSite,intBpNTimes,
        getWaveConstants, getU, thetaDet, generateEquations,processEquations,
        computeThetaConst, analyticContinuation ,aContEdgeInt,aContEdge,updateVal1,
        processPoly, getParameters,intPrep,errorEstimate,intEdge,ContinueDK, Continue,Der,
        Canonicalbasis,Frobeniustransform,Intersectionmatrix,Putfirst,Whereis,Tretkofftable,
        Findcycle, Listordered, Tretkofflist, Homologybasis, Makecycle, Smallest, Canbasis, Simplifycycle,
        Compresscycle, Reformcycle,Acy,periodMatrixInner,intAtLvL,
        secondIntialization;

    
    
$include "analyticContinuation.mm"	  
$include "homology.mm" 
$include "integrationtools.mm"
$include "init.mm"		
$include "KPLAB.mm" 
$include "KPSolsMethod1.mm"
$include "KPSolsMethod2.mm"
$include "KPErrorTools.mm"
$include "monodromy.mm"
$include "miscComputations.mm"
$include "periodMatrix.mm"
$include "polyTools.mm" 
$include "RCV.mm" 

initialize1()



end module:
savelib( 'kptools' );