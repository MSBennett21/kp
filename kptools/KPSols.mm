###Description
###INPUT
###OUTPUT
#this assumes that 2pi I is already pulled out of the vecr
#tdo
#-add user level printing
waveConstants := proc( )
    local d, i, mm, E:=[], tempVal01, possibleSolutions, solSet, parameters := [];
    
    #print(processing input-phase 1)     
    #debug option   
    if nargs = 0 then
        error "Input is needed";
    elif type( args[1], Matrix) then   
        if not IsMatrixShape( fnormal(evalf(args[1])), symmetric) then
            error "Matrix given is not symmetric on the order of global digits";
        elif not LinearAlgebra[ IsDefinite ]( Im( fnormal(evalf(args[1])) ), query = 'positive_definite' ) then
            error "Imaginary part of the Matrix is not positive definite on the order of global digits";    
        fi;  #eventually apply siegel     
        _EnvRM := simplify( fnormal( evalf[Digits+10]( 
            args[1] 
        )), zeros);
        _EnvGenus := rank( Im( args[1] ) );  

        #Error check      
        if _EnvGenus = 1 then 
            tempVal01 := [ type( args[2], complexcons ), type( args[3], complexcons ) ];
        else
            tempVal01 := [ type( args[2], complexcons ), type( args[3], complexcons ), type( args[4], complexcons ) ];
        fi;
        _EnvUI := args[2 .. -1]; 
    
    elif type( args[1], polynom ) and type( args[2], symbol ) and type( args[3], symbol ) then
        _EnvGenus := algcurves[ genus ]( args[1], args[2], args[3] );        
        Digits := Digits + 10;
        _EnvRM := simplify( 
            fnormal( PeriodMatrix(
            args[1], args[2], args[3], 'Riemann'),
            Digits - 10), zero);           
        Digits := Digits - 10;

        #Error check  
        if _EnvGenus = 1 then 
            tempVal01 := [ type( args[4], complexcons ), type( args[5], complexcons ) ];
        else
            tempVal01 := [ type( args[4], complexcons ), type( args[5], complexcons ), type( args[6], complexcons ) ];
        fi;                   
        _EnvUI := args[4 .. -1];        
    
    else
        error "Input is incorrect: either input a Riemann Matix B and constants or an algebraic curve F,x,y and some constants";  
    fi;
    for i in tempVal01 do
        if not i then
            error "The input needs to be of the form B or a curve f,x,y and them followed by the appropriate number of constants.";
        fi;
    od;  

    #first argument in input is how do you want me to take the derivative after I get the equations
    if _EnvGenus = 1 and nops([ _EnvUI ]) >= 4  then  
        _EnvScale := _EnvUI[4];
    elif _EnvGenus > 1 and nops([_EnvUI]) >= 5 then
        _EnvScale := _EnvUI[5];            
    else
        _EnvScale := 1; 
    fi;      
    _EnvNumOfChars := 2^_EnvGenus;
    _EnvNonSingCond := _EnvGenus * ( _EnvGenus + 1 )/2 + 1;     
    _EnvIg := Id(_EnvGenus);

    _Envij := 
                [ seq( seq( [ i, j ], j = i.._EnvGenus), i = 1.._EnvGenus) ];
    _Envijk := 
                [ seq( seq( seq( [ i, j, k ], k = j.._EnvGenus), j = i.._EnvGenus), i = 1.._EnvGenus) ];
    _Envijkl := 
            [ seq( seq( seq( seq( [ i, j, k, l], l = k.._EnvGenus), k = j.._EnvGenus), j = i.._EnvGenus), i = 1.._EnvGenus) ];
    
    #PREPROCESS COMPLETE                
    possibleSolutions := getWaveConstants();
    #print(possibleSolutions);
    for solSet in possibleSolutions do      
        mm := abs( HBF( 0, 0, 0, _EnvRM, solSet[1], solSet[2], solSet[3], solSet[4],0, Digits ) );
        E := [ 
            op(E), mm 
            ];          
            parameters := [ 
            op( parameters ), simplify( fnormal( [ solSet[1], solSet[2], solSet[3], solSet[4], 0] ), zero ) 
            ];        
    end do;   
    if nops(parameters)=0 then
        return [ _EnvRM, [], nops(possibleSolutions) ];
    else 
        return [[ _EnvRM, E, nops(possibleSolutions) ],  op(parameters)  ];  
    fi;  
end proc;

getWaveConstants:= proc()
    local i,j,k, EqSetV,EqSetW, notZero, sols:=[], UU, V, W, du, Vk, Vu, Wu, USet, UVal, tempVal11,tempVal12,
        EqSetMain, DI := max(Digits, 10), Q, Qu,  M, Mu, possibleVs, 
        scrubData := X -> simplify( fnormal( evalf[DI](X)), zero);
    
         
    if _EnvGenus = 1  then #if G=1 then U1, Vk so V is at 2. G>1 U1, U2, Vk so Vk is at 3
        UU:= [ scrubData( _EnvUI[1] ) ];
        Vk := scrubData( _EnvUI[2] );                   
    else     
        Vk := scrubData( _EnvUI[3] );        
        UU := [ scrubData( _EnvUI[1] ), scrubData( _EnvUI[2] ),             
            seq( convert( cat( "u", convert( i, string ) ), symbol ), i = 3 .. _EnvGenus ) 
        ];  
    fi; 
    V := [ 
        seq( convert( cat( "v", convert( i, string ) ), symbol ), i = 1 .. _EnvGenus )
        ];   
    W := [ 
        seq( convert( cat( "w", convert( i, string ) ), symbol ), i = 1 .. _EnvGenus )
        ];  
    Q := [ 
        seq( convert( cat( "q", convert( i, string ) ), symbol ), i = 1 .. _EnvNonSingCond ) 
        ];    
    USet, M := getU(UU,DI);
     
    #print(USet);
    EqSetMain := generateEquations(_EnvGenus, _EnvNonSingCond, _Envij, _EnvIg, _EnvScale);
    if member('bug',{_EnvUI}) then
        print(USet);
        print("the generated equations:");
        print(EqSetMain);
    fi;
    for UVal in USet do
        if member('bug',{_EnvUI}) then
            print(UVal);
        fi; 
        Mu := subs( seq( UU[i]= UVal[i], i = 3.._EnvGenus ), M);
if member('bug',{_EnvUI}) then
            print(Mu);
        fi; 
        Mu := LinearAlgebra[LUDecomposition](Mu,output='U');
        if member('bug',{_EnvUI}) then
            print(Mu);
        fi;      
        Qu := LinearAlgebra[NullSpace]( Mu );   
        if  nops( {op(fnormal(UVal))} ) = 1 and fnormal(UVal[1])=0 then
            WARNING("The U vector found is zero and will be skipped");
        elif nops(Qu)>=1 then            
            Qu := scrubData( Qu[1]/Qu[1][ _EnvNonSingCond + 1 ]);           
          
        EqSetV,EqSetW, k, notZero := processEquations(_EnvGenus, _EnvNonSingCond, _Envij, _EnvIg, EqSetMain, UVal, Qu, Vk);
        possibleVs:=[]; 
        tempVal11 := getSubObjects( [ op({seq(i, i=1.._EnvGenus)} minus {k}) ] )[ 2..-1 ];      
        tempVal12 := 
        [ diag( _EnvIg ), 
            seq( diag( _EnvIg ) + add( col( _EnvIg, tempVal11[i][l] ), l = 1 .. nops(tempVal11[i])), i = 1 .. nops(tempVal11) ) 
        ];   
        for i from 1 to 2^(_EnvGenus-1) do
            tempVal11 := tempVal12[i];
            #print(tempVal01);
            Vu :=[ seq(0, j= 1.._EnvGenus ) ];    
            for j from 1 to _EnvGenus do  
                #print(j);
                #print(rhs(EqSetV[ j ]));
                if j=k then                    
                    Vu[j] := scrubData(rhs(EqSetV[ j ]));
                else 
                    #print( solve( subs( seq(V[l]=Vu[l], l=1..j-1), EqSetV[ j ])) 
                    #    );
                   Vu[j] := scrubData(solve( subs( seq(V[l]=Vu[l], l=1..j-1), EqSetV[ j ]))[ tempVal11[j]]);
                fi;             
            od;
            possibleVs := [ op(possibleVs), Vu];
        od;
        for Vu in possibleVs do
            Wu := [seq( 0, i=1.._EnvGenus )]; 
            for i from 1 to _EnvGenus do                     
                Wu[i]:= scrubData(subs(seq( V[l]=Vu[l], l=1.._EnvGenus), rhs(EqSetW[i])));
            od;
            sols :=
            [op(sols), [ UVal, Vu, Wu, Qu[ _EnvNonSingCond ]/_EnvScale^4] ];
        od;
        
        
        else
             #print("Qu is empty");
             #print(Qu);      
        fi;
           

      
    od;        
    return sols;
end proc; 

getU := proc(U,DI)
    local        
        MatrixOfConstants, CompatabilityConditions:=[],compatibleU,subIndices,subIndex,subM ,valU:=[], 
        tempVal21, tempVal22,tempVal23, i,sol,
        scrubData := X -> simplify( fnormal( evalf[DI](X), DI), zero);  
          
    
    MatrixOfConstants := scrubData(computeThetaConst(
            _EnvGenus, _EnvNumOfChars,_EnvIg, _EnvRM, U, _Envij, _Envijk,_Envijkl, _EnvScale, 2*DI
            )); 
    if member('bug',{_EnvUI}) then          
           print(MatrixOfConstants);
    fi;     
    if _EnvGenus < 3 then   
        [U], MatrixOfConstants;
    elif _EnvGenus = 3 then           
        compatibleU :=  SolveTools[Polynomial](
            thetaDet( MatrixOfConstants, 2*DI ), u3
            ); 
        #print(compatibleU);
        for sol in compatibleU do
            #print(sol);           
            valU:=[op(valU), [scrubData( U[1]), scrubData(U[2]), scrubData(sol)] ];
        od;
        return valU, MatrixOfConstants; 
    else
    #generate a system of polynomials to solve for the remaining u values
        #F=set of ints that are not rows of A. Take smallest one. Get determinate of sm formed from rows of A and row
        #bezots theorem!
        subIndices := getSubObjects( _EnvNumOfChars, _EnvNonSingCond+1 ); 
        for i from 1 to _EnvGenus-2 do
            subIndex := sort( subIndices[i]);
            subM := sm( MatrixOfConstants, subIndex, [ seq(j,j=1.._EnvNonSingCond+1)  ]);       
            CompatabilityConditions:=[
                op(CompatabilityConditions), thetaDet( subM , 2*DI)
            ];               
        od;
        compatibleU := evalf(allvalues(SolveTools[PolynomialSystem]( 
            CompatabilityConditions, [seq(U[i],i=3.._EnvGenus)], preprocess=true ))
            );
        #print(compatibleU);    
        for sol in  compatibleU do            
            valU :=[
                op(valU), subs( op(sol), scrubData~(U) ) 
            ]; 
        od;
        #print(valU);   
        valU, MatrixOfConstants;     
    fi;
end proc; 



thetaDet := proc(M :: Matrix, DI)        
    Digits:=DI;  
    return det(LinearAlgebra[LUDecomposition](M,output='U'));
end proc;    
generateEquations := proc(g, n, ij, Ig, scale)
    local u:=[seq( convert( cat( "u", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    v:=[seq( convert( cat( "v", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    w:=[seq( convert( cat( "w", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    q:=[seq( convert( cat( "q", convert( i, string ) ), symbol ), i = 1 .. n ) ], 
    IJ, j, II, J,
    eqSet:=Array(1..n-1);        
    for j from 1 to n-1 do     
        IJ := ij[j]; II := IJ[1]; J := IJ[2];          
        if II=J then
            eqSet[j] := -scale^2*u[II]*w[II]+ 3/4*scale^2*v[II]^2=q[j];            
        else
            eqSet[j] := -scale^2*(u[II]*w[J] + u[J]*w[II]) + 3/2*scale^2*v[II]*v[J]=q[j];
        fi;                                
    od;
    return convert(eqSet, list);   
end proc; 

processEquations := proc(g, n, ij, Ig, eqs,U, Qu, Vk)
    local ii, u:=[seq( convert( cat( "u", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    v:=[seq( convert( cat( "v", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    w:=[seq( convert( cat( "w", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    q:=[seq( convert( cat( "q", convert( i, string ) ), symbol ), i = 1 .. n ) ],        
    eqSet := subs( seq( u[i]=U[i], i=1.. g), seq( q[i]=Qu[i] , i=1.. n),eqs), 
    eqsV,eqsW,
    nonZeroIndex := { seq(i,i=1..g) } minus { ListTools[SearchAll]( 0,U )},
    k:= min( nonZeroIndex );
    #print("Uk,k=");
    #print(k);    
    #print("Eq set after eval");
    #print(eqSet);
    for ii from 1 to g do
        if not member(ii, nonZeroIndex) then
            eqsV[ii] := v[ii]= solve( eqSet[ searchList( [ii,ii], ij )], v[ii]);
            eqsW[ii] := w[ii]= solve( eqSet[ searchList( sort([ii,k]), ij )], w[ii]);                    
        elif ii=k then
            eqsV[ii] := v[k]=Vk;
            eqsW[ii] := w[ii]=solve( eqSet[ searchList( [ii,ii], ij )], w[ii] );
        else #k<i
        
            eqsW[ii] := w[ii]=solve( eqSet[ searchList( [ii,ii], ij )], w[ii]);
            eqsV[ii] := subs( eqsW[k], eqsW[ii], eqSet[searchList([k,ii],ij)]);                 
        fi; 
    od;

    return convert(eqsV,list),convert(eqsW,list), k, nonZeroIndex;
end proc; 



computeThetaConst := proc(g, N, Ig, BB, U, ij, ijk, ijkl, scale, DI)
    local ThetaTable := table(), theta4 := table(), thetaMatrix := Matrix( N, nops(ij) + 2 ), 
    expK, argumentVal, i, j, k, l, II, J, K, L, IJ, IJK, IJKL, ei, ej, ek, el, charVal,
    Z := [ seq( convert( cat( "z", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    scrubData := X -> simplify( fnormal( evalf[DI](X), DI/2 ), zero ), 
    tempVal21 := getSubObjects( [ seq(i, i=1..g) ] )[ 2..-1 ],
    Char := [ Vector( g ), seq( add( col( Ig/2, tempVal21[i][k] ), k = 1 .. nops(tempVal21[i])), i = 1 .. nops(tempVal21) ) ],   
    argument := 
        [ seq( convert( 2 * prod(  BB, Char[i] ), list),  i = 1..N ) ], 
    UCoef4 := 
        [ seq( seq( seq( seq( mcoef( 4, op( convert( col( _EnvIg, i ) + col( _EnvIg, j ) + col( _EnvIg, k ) + col( _EnvIg, l ), list) ) ), l = k.._EnvGenus), k = j.._EnvGenus), j = i.._EnvGenus), i = 1.._EnvGenus) ];    
    
    #_EnvUCoef2 := 
    #            [ seq( seq( mcoef( 2, op( convert( col( _EnvIg, i ) + col( _EnvIg, j ), list) ) ), j = i.._EnvGenus), i = 1.._EnvGenus) ];
    #_EnvUCoef3 := 
    #            [ seq( seq( seq( mcoef( 3, op( convert( col( _EnvIg, i ) + col( _EnvIg, j ) + col( _EnvIg, k ), list) ) ), k = j.._EnvGenus), j = i.._EnvGenus), i = 1.._EnvGenus) ];
    
    for i from 1 to N do
        argumentVal := scrubData( argument[i] );
        charVal := Vector( Char[i] ); 
        ThetaTable[ i, 0 ] := RiemannTheta( argumentVal, 2 * BB, [] );
        
        #print(prod(  BB, charVal ) + scale  *  prod(  tp( charVal ), <op(Z)>));
        expK := X -> exp( 
            (2*Pi*I) * (prod( tp(charVal), prod(  BB, charVal ))  + scale  *  prod(  tp( charVal ), <op(X)> ))
            );   #question: shound 2piI only impact the matrix here? NO: otherwise scale on derivatives must be 2piI to acount for RT function             
        for j from 1 to g do        
            ei := convert( col( Ig, j ), list ); 
            ThetaTable[i,1,j] := scale * RiemannTheta( argumentVal, 2 * BB, [ ei ] ) ;                   
        end do;
       
        for j from 1 to nops( ij ) do     
            IJ := ij[j]; II := IJ[1]; J := IJ[2];
            ei := convert( col( Ig, II ), list); ej := convert( col( Ig, J ), list);       
            ThetaTable[i,2,j] := scale^2 * RiemannTheta( argumentVal, 2 * BB, [ ei, ej ] ) ;  
        end do;
        for j from 1 to nops( ijk ) do       
            IJK := ijk[j]; II := IJK[1]; J := IJK[2]; K := IJK[3];         
            ei := convert( col( Ig, II), list); ej := convert( col( Ig, J), list); ek := convert( col( Ig, K) ,list);            
            ThetaTable[i,3,j] := scale^3 * RiemannTheta( argumentVal, 2 * BB, [ ei, ej, ek ] );
        end do;
        for j from 1 to nops( ijkl ) do     
            IJKL := ijkl[j]; II := IJKL[1]; J := IJKL[2]; K := IJKL[3]; L := IJKL[4];      
            ei := convert( col( Ig, II), list); ej := convert( col( Ig, J), list); 
            ek := convert( col( Ig, K), list); el := convert( col( Ig, L), list);           
            ThetaTable[i,4,j] := scale^4 * RiemannTheta( argumentVal, 2 * BB, [ ei, ej, ek, el ] );        
        od;        
        for j from 1 to nops(ij) do     
            IJ := ij[j]; II := IJ[1]; J := IJ[2]; 
            thetaMatrix[ i, j ] :=
                scrubData( subs( seq( Z[s]=0, s=1..g),     
                    diff( expK(Z), Z[II], Z[J] ) * ThetaTable[i,0] +
                    diff( expK(Z), Z[J] ) * ThetaTable[i,1,II] +
                    diff( expK(Z), Z[II] ) * ThetaTable[i,1,J] +
                    expK(Z) * ThetaTable[i,2,j]
                ));           
        od;
        thetaMatrix[ i, nops(ij) + 1 ] := 
            scrubData( subs( seq( Z[s]=0, s=1..g), 
            expK(Z) * ThetaTable[i,0] 
            ));            
        #Matrix construction: ith row-char, final column        
        for j from 1 to nops(ijkl) do
            IJKL := ijkl[j]; II := IJKL[1]; J := IJKL[2]; K := IJKL[3]; L := IJKL[4];    
            theta4[i,j] :=
            UCoef4[j] * U[II] * U[J] * U[K] * U[L] *  scrubData( subs(           
                seq(Z[s]=0,s=1..g),       
                diff( expK(Z), Z[II], Z[J], Z[K], Z[L]) * ThetaTable[i,0] +
                diff( expK(Z), Z[II], Z[J], Z[K]) * ThetaTable[i,1,L] +
                diff( expK(Z), Z[II], Z[J], Z[L]) * ThetaTable[i,1,K] +
                diff( expK(Z), Z[II], Z[K], Z[L]) * ThetaTable[i,1,J] +
                diff( expK(Z), Z[II], Z[J]) * ThetaTable[i,2, searchList([K,L],ij) ] +
                diff( expK(Z), Z[II], Z[K]) * ThetaTable[i,2, searchList([J,L],ij) ] +
                diff( expK(Z), Z[II], Z[L]) * ThetaTable[i,2, searchList([J,K],ij) ] +
                diff( expK(Z), Z[II]) * ThetaTable[i,3, searchList([J,K,L],ijk) ] +
                diff( expK(Z), Z[J], Z[K], Z[L]) * ThetaTable[i,1,II] +   
                diff( expK(Z), Z[J], Z[K]) * ThetaTable[i,2, searchList([II,L],ij) ] +
                diff( expK(Z), Z[J], Z[L]) * ThetaTable[i,2, searchList([II,K],ij) ] +  
                diff( expK(Z), Z[J]) * ThetaTable[i,3, searchList([II,K,L],ijk) ] +
                diff( expK(Z), Z[K], Z[L]) * ThetaTable[i,2, searchList([II,J],ij) ] +
                diff( expK(Z), Z[K]) * ThetaTable[i,3, searchList([II,J,L],ijk) ] +             
                diff( expK(Z), Z[L]) * ThetaTable[i,3, searchList([II,J,K],ijk) ] +            
                expK(Z) * ThetaTable[i,4,j]));
                thetaMatrix[i,-1] := thetaMatrix[i,-1]  +  theta4[i,j];
                #print(U[II] * U[J] * U[K] * U[L]);
                #print( thetaMatrix);
        end do;  
        thetaMatrix[i,-1] :=scrubData( thetaMatrix[i,-1] );           
        unassign(expK);
    od; 
    return thetaMatrix;
end proc;

KPgSol := proc( x , y, t,  RM::Matrix, U::list, V::list, W::list, phi, c:=0, tol:=.001)
    local Theta, ThetaX, ThetaXX, tArg := convert( <op(U)> * x + <op(V)>* y + <op(W)> * t + <op(phi)> * t,list ); 
    Theta := RiemannTheta( tArg,  RM, [], tol, output=list )[2];
    ThetaX := RiemannTheta( tArg,  RM, [U], tol, output=list )[2];
    ThetaXX := RiemannTheta( tArg,  RM, [U, U], tol, output=list )[2];
    simplify(fnormal(
        2/(2*Pi*I)^2 * ( ThetaXX * Theta - (ThetaX)^2 )/(Theta)^2 + c
    ), zero);       
end proc;
 
    
        (* These questions are my next line of investigation

    1)When is it that B, and a choice of U1 and U2 real lead to U3 being real?
    and
    2) Given B and U, real, we obtain collections of polynomials Q, P in U's components with real constant coefficients. 
    The realness of V's components depend, at minimum, on whether two of the P polynomials are negative and if we gave a real constant for normalization of V (i.e. one of the components should be set to real number a). This then leads to W being real so that we have what we need.
    Gramps Seger and Dubrovin worked on a different KP equation and did something a little bit different for their theta constants _EnvtArgs. I wonder if the fundamental region they have essentially comes from ensuring realness of u3 and on ensuing the matricies B produce polynomimals P that are in the fundmental region. If so then we should be able to produce this result here and for any other KP equation. So the question is can I obtain some conditions on B from these facts? 

    something different but if I can link what they did to dubrovins work  *)
    ####################################################################################################################
    #################################Populate/Processes Certain values relative to case data ###########################
    #Characteristics: Alpha in dubrovin. There is no order. Just combinations of the Identity Matrix multiplied by a scalar
    #Argument: Argument of theta function. Order dependendent on order of Chars.
    #UcoefK: Coeffecient of Kth order derivative in dir U,V,W: will be obtained from (x1+..xg)^k-- multinomial coeefs
    #_Envij,_Envijk,_Envijkl: An enumeration of the derivatives of order 2,3 and 4 in z1,..,zg
    

    
      
        (* These questions are my next line of investigation

    1)When is it that B, and a choice of U1 and U2 real lead to U3 being real?
    and
    2) Given B and U, real, we obtain collections of polynomials Q, P in U's components with real constant coefficients. 
    The realness of V's components depend, at minimum, on whether two of the P polynomials are negative and if we gave a real constant for normalization of V (i.e. one of the components should be set to real number a). This then leads to W being real so that we have what we need.
    Gramps Seger and Dubrovin worked on a different KP equation and did something a little bit different for their theta constants _EnvtArgs. I wonder if the fundamental region they have essentially comes from ensuring realness of u3 and on ensuing the matricies B produce polynomimals P that are in the fundmental region. If so then we should be able to produce this result here and for any other KP equation. So the question is can I obtain some conditions on B from these facts? 

    something different but if I can link what they did to dubrovins work  *)
    ####################################################################################################################
    #################################Populate/Processes Certain values relative to case data ###########################
    #Characteristics: Alpha in dubrovin. There is no order. Just combinations of the Identity Matrix multiplied by a scalar
    #Argument: Argument of theta function. Order dependendent on order of Chars.
    #UcoefK: Coeffecient of Kth order derivative in dir U,V,W: will be obtained from (x1+..xg)^k-- multinomial coeefs
    #_Envij,_Envijk,_Envijkl: An enumeration of the derivatives of order 2,3 and 4 in z1,..,zg
    

    #print(DI);
    #print DI has been adjusted
    
    ####################################################################################################################
    #################################Calculation Of Derivatives and Theta Matrix########################################
    ## need to compute derivatives of functions of the form fg, f an exponential and _EnvGenus a theta function evaluated at z=0.
    ## compute _EnvGenus's,f's separate and eval at 0. Compute f and ders symbolically call gs ders. Add to matrix.
    ## MatrixBig is 2^_EnvGenus by _EnvGenus(_EnvGenus+1)/2+2. Row:= ith char, Col:= jth derivative j in [op(_Envij),0th der,fourth order der]
    # _EnvGenus=ThetaTable[i,j,k] theta of char i, jth order derivative, number k in the order (or k dne for j=0)
    # f= expK symbolic expresion of exponential    

   #GOAL: Construct a (_EnvGenus*(_EnvGenus+1)+1)/2 square matrix that is appropriate and that can be used to obtain the remaining u value.
   #This matrix comes from our sys of equations:
   #Theta4U -  _EnvScale^2 *   ThetaUW +  _EnvScale^2*3/4 * ThetaVV +  _EnvScale^4*d * Theta = 0
   #where each computation above involves a  _EnvScale factor as well (we just divide out by a  _EnvScale^4/8 to go back to the og eq obtained from the change in coords
   #) What we do is recognize the linear dependence of the theta constant expressions in this system:
   #Theta4U + a11 theta11 +a2theta12 +a3 theta13 .... a1N theta gg a1(N+1) theta=0
   #the determinate here is clearly zero here
    