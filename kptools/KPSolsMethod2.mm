###Description
###INPUT
###OUTPUT
#this assumes that 2pi I is already pulled out of the vecr
#tdo
#-add user level printing
fgParameters := proc(curve,x,y)
    local i,  bp, P0,Q,Qdiv,DI,
    RM, A,B,x0,coord, divisor,yprime,z0,K,AofD,place,
    utilde, vtilde, wtilde, U, V, W; 

    _Env_GlobalOptions := [ args[ 4..-1 ] ];  
    if nargs=4 or nargs>5 then
        error "Wrong input";
    fi;
    print("Computation started for");
    print(curve);
    Digits:=Digits+2*EXTRA1;
    genusCurve := algcurves[ genus ]( curve, x, y );
    Pmatrix :=  algcurves[periodmatrix]( curve,x,y ); 
    hom := algcurves[homology]( curve,x,y );
    mon := algcurves[monodromy]( curve,x,y );
    x0 := convert( hom[ basepoint ],rational); 
    Digits:=Digits-2*EXTRA1;
    numSheets := nops( hom[sheets] );    
    differentialBasis := subs(dx=1, algcurves[differentials]( curve, x, y )); 
    print(differentialBasis);
    A := sm(Pmatrix, [ 1..genusCurve ], [ 1..genusCurve ] );
    B := sm( Pmatrix, [1..genusCurve], [ genusCurve + 1..-1 ] );
    #moving bp back to an integer is fine for the purpose of integration
    if nops( _Env_GlobalOptions ) > 0 then
        divisor := _Env_GlobalOptions[2];
        coord := _Env_GlobalOptions[1];
        if type(coord, symbol) then
            error "wrong input for parameter of expansion of divisor";
        fi;
        P0 := op(
            algcurves[puiseux]( curve, x = x0, y,Digits,coord)
        ); 
        Q := rhs~( evalf( subs( coord=0, P0 ) ) );       
        try
            Qdiv:=[seq( evalf( subs( coord=0, divisor[i] ) ) , i= 1...genusCurve )];
        catch :          
            error "wrong input for divisor";
        end try; #not the best way to check if it is a divisor
        #if specials then DA(Pi)=0->DA(P0)+DA(P1)... =0
        DI:=Digits;
        Digits:=2;
        m:=max(evalf~(abs~( 
            convert( 
                add( subs( Qdiv[i], convert(differentialBasis, vector)),i=1..genusCurve)
                ,list)
                )));                
        Digits:=DI;  
        if fnormal(m)< 10^(-2) then
            WARNING("The divisor you chose is either special or very near special (WRT to this computation), this may result in failure");
        fi;              
    else        
        coord := t;
        P0 := op(
        algcurves[puiseux]( curve, x = x0, y,Digits,coord)
        );
        divisor := [seq( P0, i=1..genusCurve )];
        Q := rhs~( evalf( subs( coord=0, P0 ) ) );
    fi;    
    yprime := - normal( diff( curve, x )/diff( curve, y ) );
    
    dxx := [ 
        seq( 
            normal( diff( differentialBasis[i], x ) + diff( differentialBasis[i], y ) * yprime),
        i = 1..genusCurve ) 
    ];
    dxxx := [
            seq(
                normal( diff( dxx[i], x) + diff(dxx[i], y) * yprime),
        i = 1..genusCurve ) 
    ];
    utilde := - subs( x = Q[1], y= Q[2], differentialBasis );
    vtilde := - subs( x = Q[1], y= Q[2], dxx );
    wtilde := - subs( x = Q[1], y= Q[2], dxxx )/2;        
    U := convert( 
        simplify( fnormal( evalf( prod( inv(A), <op(utilde)>))),zero), list);
    V := convert(
         simplify( fnormal( evalf( prod(inv(A), <op(vtilde)>))),zero), list);
    W := convert(
         simplify( fnormal( evalf( prod(inv(A), <op(wtilde)>))),zero), list);      
    RM := simplify( fnormal( evalf( prod( inv(A), B ))), zero);
    if nargs=3 then
        AofD := 0;
        try

            for place in divisor do
            AofD := AofD +  
                algcurves[AbelMap]( curve, x, y, P0, place,
                coord, Digits);       
            od;
        catch:
            AofD := Vector(genusCurve); 
        end try;
    else
        AofD := 0;
        for place in divisor do
            AofD := AofD +  
                algcurves[AbelMap]( curve, x, y, P0, place,
                coord, Digits);       
        od;
    fi;  
    m:=[degree(numer(differentialBasis[1])),1];
    if genusCurve>1 then
        for i from 2 to genusCurve do
            mm:=[degree(numer(
                differentialBasis[i])),i];
            if m[1]>mm[1] then
                m:=mm;
            fi;
        od;
    fi;
    print("differential");
    print(m[1]);
    print(differentialBasis[m[2]]);
    K := RCV( curve, x, y, coord, differentialBasis[m[2]], RM, P0);   
    print("K value");
    print(K);
    z0 := - convert(AofD,Vector) - convert(K,Vector);
    print("end");
    return [ RM, U, V, W, convert(z0,list) ];      
end proc;

getDC:= proc(RM, U,V, W, DI:=Digits)
    local Theta, ThetaX, ThetaT, ThetaY, 
    ThetaXX, ThetaXT, ThetaYY, ThetaXXX, ThetaXXXX,z := convert( evalf(<op(U)> * x + <op(V)> * y + <op(W)> * t), list),
    tol := 0.1^DI;
    
    #Eval at 0 of an exven function and odd terms should be 0

    Theta := RiemannTheta( subs(x=0,y=0,t=0,z), RM, [], tol, output = list )[2]; #Even
    ThetaXX := RiemannTheta( subs(x=0,y=0,t=0,z), RM, evalf([ U, U ]), tol, output = list )[2]; 
    ThetaXT := RiemannTheta( subs(x=0,y=0,t=0,z), RM, evalf([ U, W ]), tol, output = list )[2];
    ThetaYY := RiemannTheta( subs(x=0,y=0,t=0,z), RM, evalf([ V, V ]), tol, output = list )[2];
    ThetaXXXX := RiemannTheta( subs(x=0,y=0,t=0,z), RM, evalf([ U, U, U, U ]), tol, output = list )[2];
    ThetaX := 0; #odd
    ThetaY := 0;    
    ThetaT := 0;
    ThetaXXX := 0;
    
    
    eq1 := evalf( 
        ThetaXXXX * Theta - 4 * ThetaXXX * ThetaX + 3 * ThetaXX^2 +
        4 * ThetaX * ThetaT - 4 * ThetaXT * Theta + 3 * ThetaYY * Theta
         - 3 * ThetaY^2 + 8 * d * Theta^2 + 6*c*ThetaXX*Theta -    6*c*ThetaX^2)=0; 
    Theta := RiemannTheta( subs(x=0,y=0,t=0,z), 2*RM, [], tol, output = list )[2]; #Even
    ThetaXX := RiemannTheta( subs(x=0,y=0,t=0,z), 2*RM, evalf([ U, U ]), tol, output = list )[2]; 
    ThetaXT := RiemannTheta( subs(x=0,y=0,t=0,z), 2*RM, evalf([ U, W ]), tol, output = list )[2];
    ThetaYY := RiemannTheta( subs(x=0,y=0,t=0,z), 2*RM, evalf([ V, V ]), tol, output = list )[2];
    ThetaXXXX := RiemannTheta( subs(x=0,y=0,t=0,z),2*RM, evalf([ U, U, U, U ]), tol, output = list )[2];
    ThetaX := 0; #odd
    ThetaY := 0;    
    ThetaT := 0;
    ThetaXXX := 0;
    #Must also satisfy the quartic eqs as well so     
    #arg2:= subs(x=1,y=1,t=1, GenArg);
    #Theta := RiemannTheta( arg2, RM, [], tol, output = list )[2];
    #ThetaX := RiemannTheta( arg2, RM, evalf([  Uvec  ]), tol, output = list )[2];
    #ThetaY := RiemannTheta( arg2, RM, evalf([ Vvec ]), tol, output = list )[2];
    #ThetaT := RiemannTheta( arg2, RM, evalf([ Wvec ]), tol, output = list )[2];
    #ThetaXX := RiemannTheta( arg2, RM, evalf([ Uvec, Uvec ]), tol, output = list )[2]; 
    #ThetaXT := RiemannTheta( arg2, RM, evalf([ Uvec, Wvec ]), tol, output = list )[2];
    #ThetaYY := RiemannTheta( arg2, RM, evalf([ Vvec, Vvec ]), tol, output = list )[2];
    #ThetaXXX := RiemannTheta( arg2, RM, evalf([ Uvec, Uvec, Uvec ]), tol, output = list )[2];
    #ThetaXXXX := RiemannTheta( arg2, RM, evalf([ Uvec, Uvec, Uvec, Uvec ]), tol, output = list )[2];    
    eq2 := evalf( 
        ThetaXXXX  -  ThetaXT + 3/2*c*ThetaXX +   3/4 * ThetaYY +  d * Theta)
        = 0;  
    
    sols:=solve([eq1,eq2],[d,c]);
    if nops(sols) = 1 then
        return rhs(sols[1][1]),rhs(sols[1][2]);
    else
        error cat("there is a problem:",nops(sols));
    fi;
end proc; 
# get left pt. get diffs. a canoncial divisor using P the ps series. apply alg curves. h


