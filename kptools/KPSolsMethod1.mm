###Description
###INPUT
###OUTPUT
#this assumes that 2pi I is already pulled out of the vecr
#tdo
#-add user level printing
waveConstants := proc( )
    local d, i, mm, E:=[],
    tempVal01, possibleSolutions, solSet, parameters := [];
    #print(processing input-phase 1)     
    #debug option      
    if nargs = 0 then
        error "Input is needed";
    elif type( args[1], Matrix) then   
        #print(args[1]);
        #if not IsMatrixShape( fnormal(evalf(args[1]),Digits), symmetric) then
        #    error "Matrix given is not symmetric on the order of global digits";
        #elif not LinearAlgebra[ IsDefinite ]( Im( fnormal(evalf(args[1])) ), query = 'positive_definite' ) then
        #    error "Imaginary part of the Matrix is not positive definite on the order of global digits";    
        #fi;  #eventually apply siegel  
        UI_Digits:=max(10,Digits);
        UI_HDigits:=2*UI_Digits;
        Digits:=UI_HDigits;   
        _EnvRM := simplify( fnormal( evalf( 
            args[1] 
        )), zeros);
        g := rank( Im( args[1] ) );  
         
        #Error check      
        if g = 1 then 
            tempVal01 := [ type( args[2], complexcons ), type( args[3], complexcons ) ];
        else
            tempVal01 := [ type( args[2], complexcons ), type( args[3], complexcons ), type( args[4], complexcons ) ];
        fi;
        #clearing all input used for communication in construction and now imputing args relevant for this program
        _Env_GlobalOptions := args[2 .. -1];     
    elif type( args[1], polynom ) and type( args[2], symbol ) and type( args[3], symbol ) then
        _EnvRM := algcurves[periodmatrix](
            args[1], args[2], args[3], 'Riemann'); 
        #try
            return procname(_EnvRM, args[4..-1])
        #catch:
        #    error "You did not give curve, xcord, ycord, parameters";
        #end try;        
    else
        error "Input is incorrect";  
    fi;
    for i in tempVal01 do
        if not i then
            error "The input needs to be of the form B or a curve f,x,y and them followed by the appropriate number of constants.";
        fi;
    od;  
    #first argument in input is how do you want me to take the derivative after I get the equations
    #if g = 1 and nops([ _Env_GlobalOptions ]) >= 4  then  
    #    _EnvScale := _Env_GlobalOptions[4];
    #elif g > 1 and nops([_Env_GlobalOptions]) >= 5 then
    #    _EnvScale := _Env_GlobalOptions[5];            
    #else
    _EnvScale := 1; 
    #fi;      
    _EnvNumOfChars := 2^g;
    _EnvNonSingCond := g * ( g + 1 )/2 + 1;     
    _EnvIg := Id(g);
    _Envij := 
                [ seq( seq( [ i, j ], j = i..g), i = 1..g) ];
    _Envijk := 
                [ seq( seq( seq( [ i, j, k ], k = j..g), j = i..g), i = 1..g) ];
    _Envijkl := 
            [ seq( seq( seq( seq( [ i, j, k, l], l = k..g), k = j..g), j = i..g), i = 1..g) ];

    #PREPROCESS COMPLETE                
    possibleSolutions := getWaveConstants();

    #EXPECTED FORM
    # [s1...si   ] with
    # sj= [Uj,Vj,Wj,d] where Uj,Vj,Wj are vectors
    # d is a complex constant
    
    ErrKeys:=[];
    solDict:=table();
    for solSet in possibleSolutions do  
    print("solution obtained:");
    print(solSet);    
        mm := abs( HBF( 0, 0, 0, _EnvRM, solSet[1], solSet[2], solSet[3], solSet[4],0, Digits ) );
        print("Error");
        print(mm);
        if mm < 0.1 then
            ErrKeys := [ 
                op(ErrKeys), mm 
                ];
            solDict[mm]:=simplify( fnormal( [ solSet[1], solSet[2], solSet[3], solSet[4], 0] ), zero );
        fi;
        MM := min(MM,mm);        
    end do;   
    if nops(ErrKeys)=0 then
        return cat("no meaningfull solution found", convert(MM,string));
    else 
        ErrKeys:=sort(ErrKeys);
        return [ _EnvRM, ErrKeys, solDict ];  
    fi;  
end proc;

getWaveConstants:= proc()
    local i,j,k, EqSetV,EqSetW, notZero, sols:=[],
    UU, V, W, du, Vk, Vu, Wu, USet, UVal, tempVal11,tempVal12,
        EqSetMain, Q, Qu,  M, Mu, possibleVs;
    ##UU is the U vector in terms of knowns and unknonwns
    ## Vk is the parameter for V. We want to use it as soon as we can    
    if g = 1  then #if G=1 then U1, Vk so V is at 2. G>1 U1, U2, Vk so Vk is at 3
        UU := [ _Env_GlobalOptions[1]  ];
        Vk := _Env_GlobalOptions[2] ;                   
    else     
        Vk := _Env_GlobalOptions[3];        
        UU := [ _Env_GlobalOptions[1] , _Env_GlobalOptions[2] ,             
            seq( convert( cat( "u", convert( i, string ) ), symbol ), i = 3 .. g ) 
        ];  
    fi; 
    V := [ 
        seq( convert( cat( "v", convert( i, string ) ), symbol ), i = 1 .. g )
        ];   
    W := [ 
        seq( convert( cat( "w", convert( i, string ) ), symbol ), i = 1 .. g )
        ]; 
    # The Q vector will be the subvector of X that is nontrivial 
    Q := [ 
        seq( convert( cat( "q", convert( i, string ) ), symbol ), i = 1 .. _EnvNonSingCond ) 
        ]; 
 
    USet, M := getU(UU);
    #Expected out put
    # Uset, M
    ## Uset=[U1,... Uj] which are Ui lists
    ## M is the matrix used for the determinate condition 



    #print(USet);
    EqSetMain := generateEquations();
    
    if member('bug',{_Env_GlobalOptions}) then
        print(USet);
        print("the generated equations:");
        print(EqSetMain);
    fi;
    print("the generated equations:");
    print(EqSetMain);
    for UVal in USet do
        if member('bug',{_Env_GlobalOptions}) then
            print(UVal);
        fi; 
         print("the value of U:");
        Mu := subs( seq( UU[i]= UVal[i], i = 3..g ), M);

        Mu :=simplify( fnormal( LinearAlgebra[LUDecomposition](Mu,output='U'),UI_Digits),zeros);
        
        if member('bug',{_Env_GlobalOptions}) then
            print("matrix after rr");
            print(Mu);
        fi;      
        Qu := LinearAlgebra[NullSpace]( Mu );   
        if  nops( {op(fnormal(UVal))} ) = 1 and fnormal(UVal[1])=0 then
            WARNING("The U vector found is zero and will be skipped");
        elif nops(Qu)>=1 then            
                Qu := simplify(fnormal(evalf( Qu[1]/Qu[1][ _EnvNonSingCond + 1 ]),UI_Digits),zeros);           
            print("the X vector");
            print(Qu);
            EqSetV,EqSetW, k, notZero := processEquations(EqSetMain, UVal, Qu, Vk);
            possibleVs:=[]; 
            tempVal11 := getSubObjects( [ op({seq(i, i=1..g)} minus {k}) ] )[ 2..-1 ];      
            tempVal12 := 
            [ diag( _EnvIg ), 
                seq( diag( _EnvIg ) + add( col( _EnvIg, tempVal11[i][l] ), l = 1 .. nops(tempVal11[i])), i = 1 .. nops(tempVal11) ) 
            ];   
            for i from 1 to 2^(g-1) do
                tempVal11 := tempVal12[i];
                #print(tempVal01);
                Vu :=[ seq(0, j= 1..g ) ];    
                for j from 1 to g do  
                    #print(j);
                    #print(rhs(EqSetV[ j ]));
                    if j=k then                    
                        Vu[j] := simplify(fnormal(evalf(rhs(EqSetV[ j ])),UI_Digits),zero);
                    else 
                        #print( solve( subs( seq(V[l]=Vu[l], l=1..j-1), EqSetV[ j ])) 
                        #    );
                    Vu[j] := simplify(fnormal(evalf(solve( subs( seq(V[l]=Vu[l], l=1..j-1), EqSetV[ j ]))[ tempVal11[j]]),UI_Digits),zero);
                    fi;             
                od;
                possibleVs := [ op(possibleVs), Vu];
            od;
            for Vu in possibleVs do
                Wu := [seq( 0, i=1..g )]; 
                for i from 1 to g do                     
                    Wu[i]:= simplify(fnormal(evalf(subs(seq( V[l]=Vu[l], l=1..g), rhs(EqSetW[i]))),UI_Digits),zero);
                od;
                sols :=
                [op(sols), [ UVal, Vu, Wu, Qu[ _EnvNonSingCond ]/_EnvScale^4] ];
            od;
        
        
        else
             print("Qu is empty");
             print(Qu);      
        fi;
           

      
    od;        
    return sols;
end proc; 

getU := proc(U)
    local MatrixOfConstants, CompatabilityConditions:=[],
    compatibleU,subIndices,subIndex,subM ,valU:=[], 
        tempVal21, tempVal22,tempVal23, i,sol;  
          
    #g=1 2 by 3 square is 2 by 2
    #g=2 4 by 5 square is 4 by 4
    #g=3 8 by 8 square is both
    MatrixOfConstants := computeThetaConst(U); 
    #OUTPUT is the matrix M with right most column the one in U
    if member('bug',{_Env_GlobalOptions}) then          
        print(MatrixOfConstants);
    fi;
    
    if g < 3 then 
        print("U VALUES ALREADY OBTAINED");   
        [U], simplify(fnormal(evalf(MatrixOfConstants),Digits),zero);
    else
        print("getting U"); 
        print("Matrix");
        Mtilde:= sm( MatrixOfConstants, [seq(j,j=1.._EnvNumOfChars)], [ seq(j,j=1.._EnvNonSingCond)  ]);
        isReasonable:=false;
        DI:=Digits;
        Digits:=10;
        while not isReasonable do #just want to recover rank
            mm:=simplify(fnormal(evalf(Mtilde)),zeros);
            
            if LinearAlgebra[Rank](mm)=_EnvNonSingCond then
                mm:=simplify(fnormal(LinearAlgebra[ReducedRowEchelonForm](mm)),zeros);
                isReasonable:=true; 
                index:=[];
                indexEx:=[];
                for r from 1 to _EnvNumOfChars do 
                #RREForm piviot positions will correspond to linearly independent columns
                    for c from r to _EnvNonSingCond do
                        if mm[r,c] <> 0 then
                        #if the r,c position is non zero 
                        #it will be the first time this occurs (induction)
                            index:=[op(index),r];
                            break;
                        fi; #if this never occurs the current row is linearly dependent on the rows above
                    od;
                    if index[-1]<>r then
                        indexEx:=[op(indexEx),r];
                    fi;
                od;
            else
                if Digits > UI_HDigits then
                    error "increase digits, nonsingularity is not working"
                fi;
                Digits:=Digits+2;
                print("Digits needed to be increased");
                print(Digits);
            fi;
        od;
    fi;
    Digits:=UI_HDigits;
    index:=[op(index),indexEx[1]];
    Pu :=  solve(
            evalf(thetaDet( sm(MatrixOfConstants, index, [seq( j, j=1.._EnvNonSingCond+1)]))),
            [seq(U[j],j=3..g)]
        ); 
        #print(compatibleU   
        for sol in  Pu do            
            valU :=[
                op(valU), subs( op(sol), evalf~(U) ) 
            ]; 
        od;
        #print(valU);   
        valU, MatrixOfConstants;
        
    
end proc; 



thetaDet := proc(M :: Matrix)
    return det(LinearAlgebra[LUDecomposition](M,output='U'));
end proc;    
generateEquations := proc()
    local u:=[seq( convert( cat( "u", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    v:=[seq( convert( cat( "v", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    w:=[seq( convert( cat( "w", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    q:=[seq( convert( cat( "q", convert( i, string ) ), symbol ), i = 1 .. _EnvNonSingCond ) ], 
    IJ, j, II, J,
    eqSet:=Array(1.._EnvNonSingCond-1);        
    for j from 1 to _EnvNonSingCond-1 do     
        IJ := _Envij[j]; II := IJ[1]; J := IJ[2];          
        if II=J then
            eqSet[j] := -_EnvScale^2*u[II]*w[II]+ 3/4*_EnvScale^2*v[II]^2=q[j];            
        else
            eqSet[j] := -_EnvScale^2*(u[II]*w[J] + u[J]*w[II]) + 3/2*_EnvScale^2*v[II]*v[J]=q[j];
        fi;                                
    od;
    return convert(eqSet, list);   
end proc; 

processEquations := proc(eqs,U, Qu, Vk)
    local ii, u:=[seq( convert( cat( "u", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    v:=[seq( convert( cat( "v", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    w:=[seq( convert( cat( "w", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    q:=[seq( convert( cat( "q", convert( i, string ) ), symbol ), i = 1 .. _EnvNonSingCond) ],        
    eqSet := subs( seq( u[i]=U[i], i=1.. g), seq( q[i]=Qu[i] , i=1.. _EnvNonSingCond),eqs), 
    eqsV,eqsW,
    nonZeroIndex := { seq(i,i=1..g) } minus { ListTools[SearchAll]( 0,U )},
    k:= min( nonZeroIndex );
    #print("Uk,k=");
    #print(k);    
    #print("Eq set after eval");
    #print(eqSet);
    for ii from 1 to g do
        if not member(ii, nonZeroIndex) then
            eqsV[ii] := v[ii]= solve( eqSet[ searchList( [ii,ii], _Envij )], v[ii]);
            eqsW[ii] := w[ii]= solve( eqSet[ searchList( sort([ii,k]), _Envij )], w[ii]);                    
        elif ii=k then
            eqsV[ii] := v[k]=Vk;
            eqsW[ii] := w[ii]=solve( eqSet[ searchList( [ii,ii], _Envij )], w[ii] );
        else #k<i
        
            eqsW[ii] := w[ii]=solve( eqSet[ searchList( [ii,ii], _Envij )], w[ii]);
            eqsV[ii] := subs( eqsW[k], eqsW[ii], eqSet[searchList([k,ii],_Envij)]);                 
        fi; 
    od;

    return convert(eqsV,list),convert(eqsW,list), k, nonZeroIndex;
end proc; 


computeThetaConst := proc( U)
    local ThetaTable := table(), 
    theta4 := table(), thetaMatrix := Matrix( _EnvNumOfChars, _EnvNonSingCond + 1 ), 
    expK, argumentVal, i, j, k, l, II, J, K, L, IJ, IJK, IJKL, ei, ej, ek, el, charVal,
    Z := [ seq( convert( cat( "z", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    tempVal21 := getSubObjects( [ seq(i, i=1..g) ] )[ 2..-1 ],
    Char := [ Vector( g ), seq( add( col( _EnvIg/2, tempVal21[i][k] ), k = 1 .. nops(tempVal21[i])), i = 1 .. nops(tempVal21) ) ],   
    argument := 
        [ seq( convert( 2 * prod(  _EnvRM, Char[i] ), list),  i = 1.._EnvNumOfChars ) ], 
    UCoef4 := 
        [ seq( seq( seq( seq( mcoef( 4, op( convert( col( _EnvIg, i ) + col( _EnvIg, j ) + col( _EnvIg, k ) + col( _EnvIg, l ), list) ) ), l = k..g), k = j..g), j = i..g), i = 1..g) ];    
    
    #_EnvUCoef2 := 
    #            [ seq( seq( mcoef( 2, op( convert( col( _EnvIg, i ) + col( _EnvIg, j ), list) ) ), j = i..g), i = 1..g) ];
    #_EnvUCoef3 := 
    #            [ seq( seq( seq( mcoef( 3, op( convert( col( _EnvIg, i ) + col( _EnvIg, j ) + col( _EnvIg, k ), list) ) ), k = j..g), j = i..g), i = 1..g) ];
    #print("computing theta");
    
    for i from 1 to _EnvNumOfChars do
        argumentVal := argument[i];
        charVal := Vector( Char[i] ); 
        ThetaTable[ i, 0 ] := simplify(fnormal(evalf(RiemannTheta( argumentVal, 2 * _EnvRM, [] )),UI_Digits+EXTRA1),zero);
        #print(prod(  _EnvRM, charVal ) + _EnvScale  *  prod(  tp( charVal ), <op(Z)>));
        expK := X -> exp( 
            (2*Pi*I) * (prod( tp(charVal), prod(  _EnvRM, charVal ))  + _EnvScale  *  prod(  tp( charVal ), <op(X)> ))
            );   #question: shound 2piI only impact the matrix here? NO: otherwise _EnvScale on derivatives must be 2piI to acount for RT function             
        for j from 1 to g do        
            ei := convert( col( _EnvIg, j ), list ); 
            ThetaTable[i,1,j] := simplify(fnormal(evalf(_EnvScale * RiemannTheta( argumentVal, 2 * _EnvRM, [ ei ] )),UI_Digits+EXTRA1),zero);                   
        end do;
       
        for j from 1 to _EnvNonSingCond-1 do     
            IJ := _Envij[j]; II := IJ[1]; J := IJ[2];
            ei := convert( col( _EnvIg, II ), list); ej := convert( col( _EnvIg, J ), list);       
            ThetaTable[i,2,j] := simplify(fnormal(evalf(_EnvScale^2 * RiemannTheta( argumentVal, 2 * _EnvRM, [ ei, ej ] )),UI_Digits+EXTRA1),zero);;  
        end do;
        for j from 1 to nops( _Envijk ) do       
            IJK := _Envijk[j]; II := IJK[1]; J := IJK[2]; K := IJK[3];         
            ei := convert( col( _EnvIg, II), list); ej := convert( col( _EnvIg, J), list); ek := convert( col( _EnvIg, K) ,list);            
            ThetaTable[i,3,j] := simplify(fnormal(evalf(_EnvScale^3 * RiemannTheta( argumentVal, 2 * _EnvRM, [ ei, ej, ek ] )),UI_Digits+EXTRA1),zero);;
        end do;
        for j from 1 to nops( _Envijkl ) do     
            IJKL := _Envijkl[j]; II := IJKL[1]; J := IJKL[2]; K := IJKL[3]; L := IJKL[4];      
            ei := convert( col( _EnvIg, II), list); ej := convert( col( _EnvIg, J), list); 
            ek := convert( col( _EnvIg, K), list); el := convert( col( _EnvIg, L), list);           
            ThetaTable[i,4,j] := simplify(fnormal(evalf(_EnvScale^4 * RiemannTheta( argumentVal, 2 * _EnvRM, [ ei, ej, ek, el ] )),UI_Digits+EXTRA1),zero);        
        od;    

        for j from 1 to _EnvNonSingCond-1 do     
            IJ := _Envij[j]; II := IJ[1]; J := IJ[2]; 
            thetaMatrix[ i, j ] :=
               simplify(fnormal(evalf( subs( seq( Z[s]=0, s=1..g),     
                    diff( expK(Z), Z[II], Z[J] ) * ThetaTable[i,0] +
                    diff( expK(Z), Z[J] ) * ThetaTable[i,1,II] +
                    diff( expK(Z), Z[II] ) * ThetaTable[i,1,J] +
                    expK(Z) * ThetaTable[i,2,j])),UI_Digits+EXTRA1),zero);       
        od;
        thetaMatrix[ i, _EnvNonSingCond ] := 
              simplify(fnormal(evalf(subs( seq( Z[s]=0, s=1..g), 
            expK(Z) * ThetaTable[i,0]) 
            ),UI_Digits+EXTRA1),zero);           
        #Matrix construction: ith row-char, final column        
        for j from 1 to nops(_Envijkl) do
            IJKL := _Envijkl[j]; II := IJKL[1]; J := IJKL[2]; K := IJKL[3]; L := IJKL[4];    
            theta4[i,j] :=
             simplify(fnormal(evalf(UCoef4[j] * U[II] * U[J] * U[K] * U[L] * subs(           
                seq(Z[s]=0,s=1..g),       
                diff( expK(Z), Z[II], Z[J], Z[K], Z[L]) * ThetaTable[i,0] +
                diff( expK(Z), Z[II], Z[J], Z[K]) * ThetaTable[i,1,L] +
                diff( expK(Z), Z[II], Z[J], Z[L]) * ThetaTable[i,1,K] +
                diff( expK(Z), Z[II], Z[K], Z[L]) * ThetaTable[i,1,J] +
                diff( expK(Z), Z[II], Z[J]) * ThetaTable[i,2, searchList([K,L],_Envij) ] +
                diff( expK(Z), Z[II], Z[K]) * ThetaTable[i,2, searchList([J,L],_Envij) ] +
                diff( expK(Z), Z[II], Z[L]) * ThetaTable[i,2, searchList([J,K],_Envij) ] +
                diff( expK(Z), Z[II]) * ThetaTable[i,3, searchList([J,K,L],_Envijk) ] +
                diff( expK(Z), Z[J], Z[K], Z[L]) * ThetaTable[i,1,II] +   
                diff( expK(Z), Z[J], Z[K]) * ThetaTable[i,2, searchList([II,L],_Envij) ] +
                diff( expK(Z), Z[J], Z[L]) * ThetaTable[i,2, searchList([II,K],_Envij) ] +  
                diff( expK(Z), Z[J]) * ThetaTable[i,3, searchList([II,K,L],_Envijk) ] +
                diff( expK(Z), Z[K], Z[L]) * ThetaTable[i,2, searchList([II,J],_Envij) ] +
                diff( expK(Z), Z[K]) * ThetaTable[i,3, searchList([II,J,L],_Envijk) ] +             
                diff( expK(Z), Z[L]) * ThetaTable[i,3, searchList([II,J,K],_Envijk) ] +            
                expK(Z) * ThetaTable[i,4,j])),UI_Digits+EXTRA1),zero) ;
                thetaMatrix[i,-1] := thetaMatrix[i,-1]  +  theta4[i,j];
                #print(U[II] * U[J] * U[K] * U[L]);
                #print( thetaMatrix);
        end do;            
        unassign(expK);
    od; 
    #print("done");
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
    ## need to compute derivatives of functions of the form fg, f an exponential and g a theta function evaluated at z=0.
    ## compute g's,f's separate and eval at 0. Compute f and ders symbolically call gs ders. Add to matrix.
    ## MatrixBig is 2^g by g(g+1)/2+2. Row:= ith char, Col:= jth derivative j in [op(_Envij),0th der,fourth order der]
    # g=ThetaTable[i,j,k] theta of char i, jth order derivative, number k in the order (or k dne for j=0)
    # f= expK symbolic expresion of exponential    

   #GOAL: Construct a (g*(g+1)+1)/2 square matrix that is appropriate and that can be used to obtain the remaining u value.
   #This matrix comes from our sys of equations:
   #Theta4U -  _EnvScale^2 *   ThetaUW +  _EnvScale^2*3/4 * ThetaVV +  _EnvScale^4*d * Theta = 0
   #where each computation above involves a  _EnvScale factor as well (we just divide out by a  _EnvScale^4/8 to go back to the og eq obtained from the change in coords
   #) What we do is recognize the linear dependence of the theta constant expressions in this system:
   #Theta4U + a11 theta11 +a2theta12 +a3 theta13 .... a1N theta gg a1(N+1) theta=0
   #the determinate here is clearly zero here
    