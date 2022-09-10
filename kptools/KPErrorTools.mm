
#FACT if 0 in one direction functionis zero function
checkWCError := proc( RM::Matrix, U::list, V::list, W::list, d::complexcons, c := 0, DI:=2 )
    local mm, i, bool:=false;
   # if eq(Re(RM), zero( nops(U) ))  then
   #     for i from 1 to nops(U) do
    #        if fnormal(evalf((U[i]-conjugate(U[i]))/2), Digits-6) <> 0.0 then 
    #            bool := true;
    #            break;
    #        fi;
    #    od;        
    #fi;
    
    mm:=-1;
    if member('line',{args}) then
        mm:= abs(HBF(0, 0,0, RM, U, V, W, d, c, DI ));
        #print(mm);
        for i from 1 to 20 do
            if abs(HBF( (i-1)/19, 0, 0, RM, U, V, W, d, c, DI )) >= mm then 
                mm := abs(HBF((i-1)/19, 0,0, RM, U, V, W, d, c, DI ));
                #print(mm);
            fi;            
        od;
    else
        mm:=max([seq(abs( HBF( evalf(rand()/(10^12 - 1)), evalf(rand()/(10^12 - 1)), evalf(rand()/(10^12 - 1)), RM, U, V, W, d, c, DI )), i=1..20)] ); 
    fi;
    mm;    
end proc;

errorEstimate:=proc(e,sGrid,R)
    local Ne:=(nops(sGrid)-1)/2, AC:=R,i,j,k,x0,xs, yVec,h,sVals;	
	xs:=edgeData[e]["path"];
    h:=1/2^edgeData[e]["lvl"];
    edgeData[e]["Error"]:=<seq(0,k=1..numSheets)>;
    for i from 1 to 2*Ne+1 do #i=1-> -Ne i=N+1 0  i-N-1     
        x0:=subs(s=sGrid[i],xs);
        for j from 1 to nops(AC) do
            if scrub(UI_Digits,evalf(sGrid[i]-AC[j][1]))=0. then
                yVec:=AC[j][2];
                edgeData[e]["Error"]:=edgeData[e]["Error"]+
                    <seq(
                        evalf(subs(
                        x=x0,y=yVec[k], t=(i-Ne-1)*h,
                        gam=edgeData[e]["gam"],
                        intPrimitives["Ftt"])),                        
                    k=1..numSheets)>;
                break;
            fi;
	    od;
    od;
    Digits:=UI_Digits+2*EXTRA1;
    edgeData[e]["Error"]:= h^3/(2*Pi)* max(abs~(evalf~(
                        convert(edgeData[e]["Error"],list)
                    ))) ;
    print(edgeData[e]["Error"]);
    print(edgeData[e]["lvl"]);
	if edgeData[e]["Error"]> 10^(-UI_Digits+1) then
            if edgeData[e]["lvl"]< intPrimitives["Glvl"] then
                sVals:=intPrep(e);
                return [edgeData[e]["Error"],sVals];
            else
                error "increase digits or increase integration level"
            fi;	
	else
        return [edgeData[e]["Error"],[]];
    fi;
    
end proc;
HBF := proc(x,y,t, RM::Matrix, U::list, V::list, W::list, d::complexcons, c:=0, DI:=2 )
    local argument := convert( evalf(<op(U)> * x + <op(V)> * y + <op(W)> * t), list),
    Theta, ThetaX, ThetaT, ThetaY, 
    ThetaXX, ThetaXT, ThetaYY, 
    ThetaXXX, ThetaXXXX, tol := 0.1^DI; 
    #print(argument);
    #print(RM, [], tol, output = list);
    Theta := RiemannTheta( argument, RM, [], tol, output = list )[2];
    ThetaX := RiemannTheta( argument, RM, evalf([  U  ]), tol, output = list )[2];
    ThetaY := RiemannTheta( argument, RM, evalf([ V ]), tol, output = list )[2];
    ThetaT := RiemannTheta( argument, RM, evalf([ W ]), tol, output = list )[2];
    ThetaXX := RiemannTheta( argument, RM, evalf([ U, U ]), tol, output = list )[2]; 
    ThetaXT := RiemannTheta( argument, RM, evalf([ U, W ]), tol, output = list )[2];
    ThetaYY := RiemannTheta( argument, RM, evalf([ V, V ]), tol, output = list )[2];
    ThetaXXX := RiemannTheta( argument, RM, evalf([ U, U, U ]), tol, output = list )[2];
    ThetaXXXX := RiemannTheta( argument, RM, evalf([ U, U, U, U ]), tol, output = list )[2];

   
    evalf[DI](ThetaXXXX * Theta - 
    4 * ThetaXXX * ThetaX + 
    3 * ThetaXX^2 + 
    4 * ThetaX * ThetaT - 
    4 * ThetaXT * Theta + 
    3 * ThetaYY * Theta - 
    3 * ThetaY^2 + 
    8 * d * Theta^2 +
    6*c*ThetaXX*Theta -
    6*c*ThetaX^2);      
     
end proc;


