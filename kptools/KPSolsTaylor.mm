###Description
###INPUT
###OUTPUT
#this assumes that 2pi I is already pulled out of the vecr
#tdo
#-add user level printing
waveConstantsFromDiffs := proc(f,x,y, P)
    local g,i,dx,dxx,dxxx,bp, RM, PM,sheetIndex, A,B,ySheets, utilde, vtilde, wtilde, U,V,W; 
        g := algcurves[ genus ]( f, x, y );        
        Digits := Digits + 10;
        PM, dx, bp :=  PeriodMatrix( f,x,y, 'getDiffs'); 
        Digits := Digits - 10;
        A := sm(PM, [1..g], [1..g]);
        B := sm(PM, [1..g], [g + 1..-1]);
        ySheets:= solve(f=0,[y]); #better to get this from Pmatrix
        mm := evalf( abs( subs(x=P[1], rhs(op( ySheets[1])) ) - P[2]) );        
        for i from 1 to nops(ySheets) do
            if evalf(abs( subs(x=P[1], rhs(op( ySheets[i])) ) - P[2])) <= mm then
                sheetIndex := i;
                mm:= evalf(abs( subs(x=P[1], rhs(op( ySheets[i])) ) - P[2]));
            fi;
        od; 
         
        dx := subs( op( ySheets[sheetIndex] ), dx);
        dxx := [ seq(diff(dx[i], x),i = 1..g ) ];
        dxxx := [ seq(diff(dx[i], x,x),i = 1..g )];
        utilde := -subs(x = P[1], dx);
        vtilde := -subs(x = P[1], dxx);
        wtilde := -subs(x = P[1], dxxx)/2;        
        U := convert(simplify(fnormal(evalf(prod(inv(A), <op(utilde)>))),zero), list);
        V := convert(simplify(fnormal(evalf(prod(inv(A), <op(vtilde)>))),zero), list);
        W := convert(simplify(fnormal(evalf(prod(inv(A), <op(wtilde)>))),zero), list);      
        RM := simplify(fnormal(evalf(prod(inv(A), B))),zero);
        if member('bug',{args}) then
            print("Pmatrix ");
            print(PM);
            print("point given");
            print(P);
            print("sheets associated");
            print(ySheets);
            print("error in sheets");
            print(mm);  
        fi; 
        if member('full',{args}) then
            return PM, A, B, RM, U, V, W;
        else     
            return [ RM, U,V,W ];  
        fi;
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
Ksol := proc(f,x,y,P )
    local P_0,ds, i,ii, mm, E:=[], tempVal01, possibleSolutions, solSet, parameters := [],abe;
    mon1:=algcurves[monodromy](f,x,y); 
    print(mon1);     
    PM :=algcurves[periodmatrix](f,x,y);
    dx1 :=  algcurves[differentials](f,x,y);
    g :=  algcurves[genus](f,x,y);
    A := sm(PM, [1..g], [1..g]);
    B := sm(PM, [1..g], [g + 1..-1]);
    RM := simplify(fnormal(evalf(prod(inv(A), B))),zero);

    
    ySheets:= solve(f=0,[y]); #better to get this from Pmatrix
    mm := evalf( abs( subs(x=P[1], rhs(op( ySheets[1])) ) - P[2]) );        
    for i from 1 to nops(ySheets) do
        if evalf(abs( subs(x=P[1], rhs(op( ySheets[i])) ) - P[2])) <= mm then
            sheetIndex := i;
            mm:= evalf(abs( subs(x=P[1], rhs(op( ySheets[i])) ) - P[2]));
        fi;
    od; 
         
    dxx := subs( op( ySheets[sheetIndex] ), dx1);
    dxxx := [ seq(diff(dxx[i], x),i = 1..g ) ];
    dxxxx := [ seq(diff(dxx[i], x,x),i = 1..g )];
    utilde := -subs(x = P[1], dxx);
    vtilde := -subs(x = P[1], dxxx);
    wtilde := -subs(x = P[1], dxxxx)/2;        
    U := convert(simplify(fnormal(evalf(prod(inv(A), <op(utilde)>))),zero), list);
    V := convert(simplify(fnormal(evalf(prod(inv(A), <op(vtilde)>))),zero), list);
    W := convert(simplify(fnormal(evalf(prod(inv(A), <op(wtilde)>))),zero), list); 
    P_0 := algcurves[puiseux](f,x=convert(mon1[1], rational),y,0,t);
    bps:={seq(mon1[3][i][1],i=1..nops(mon1[3]))};
    abe:=[seq(0,i=1..g)];
    if member('bug',{args}) then
        print("monodromy");
        print(mon1);
        print("canonical pt");
        print(P_0);
        print("branchpoints");
        print(bps);        
    fi;  
    for omegai in dx1 do
        p := subs( dx=1,numer(  omegai) );
        q := subs( dx=1, denom( omegai) );
        print(p);
        print(q);
        Rp := resultant(f, p, y);
        Rq := resultant(f, q, y);   
        possiblePts:={ op(bps), infinity, op( SolveTools[Polynomial](Rp,x)), op(SolveTools[Polynomial](Rq,x)) };
        if member('bug',{args}) then
            print("differential");
            print(omegai);
            print("possible places in divisor");
            print(possiblePts);      
        fi;  
        places:=[];
        for possibleX in possiblePts do
            possiblePt := op(algcurves[puiseux](f,x=possibleX,y,0,t));
            print(possiblePt);
            print(subs( op(possiblePt), p/q*diff(rhs(possiblePt[1]),t)));
            degX := degree( subs( op(possiblePt), p/q*diff(rhs(possiblePt[1]),t)));

            print(degX);
            if  degX <> 0 and degX <> -infinity and degX <> FAIL then
                ys:= solve( subs( x=possibleX,f),[y]);
                print(ys);
                places := [op(places),possiblePt] ;
                abe[1] :=  degX*algcurves[AbelMap](f,x,y,op(P_0),possiblePt, t, Digits);
            fi;
            if member('bug',{args}) then
                print("possible pt");
                print(possiblePt);
                print("degree of differential");
                print(degX);                
                print("record of places");
                print(places);  
                print("aMapValues");
                print(abe);  
            fi;  
        od;
        abe[i]:=abe[i]/2;
        break;
    od;
    if member('bug',{args}) then
        print("vals from vaveconst");
        print(PM, A, B, RM, U, V, W);
    fi; 
    tempVal21 := getSubObjects( [ seq(i, i=1..g) ] )[ 2..-1 ];
    Ig:=Id(g);
    HalfNs := [ Vector( g ), seq( add( col( Ig/2, tempVal21[i][k] ), k = 1 .. nops(tempVal21[i])), i = 1 .. nops(tempVal21) ) ];  
    HalfMs := [ prod(RM,Vector( g )),
    prod( RM, seq( add( col( Ig/2, tempVal21[i][k] ), k = 1 .. nops(tempVal21[i])), i = 1 .. nops(tempVal21) )) 
    ]; 
    mm:=100;
    for N in HalfNs do
        for M in HalfMs do
            h:=N+M;
            K:=convert(h-<op(A[1])>,list);
            mm1:= evalf( abs(RiemannTheta(K,RM,[])));
            if mm1<=mm then
                mm:=mm1;
                RCV:=K;
              if member('bug',{args}) then
                print("decrease in theta");
                print(mm);
                print("Rvec");
                print(RCV);               
                
                fi;  
            fi;           
        od;
    od;     

    phi:=-RCV;
    return [RM, U,V,W, phi ];  
        
end proc;

