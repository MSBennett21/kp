###Description
###INPUT
###OUTPUT
#this assumes that 2pi I is already pulled out of the vecr
#tdo
#-add user level printing
#RVC is hardest to compute and needs
# RM <- period natrix
#Differential and all of the possible zeros <- period matrix
#a place user gives infor curve,x,y,P,D optional


#compute RCV first and pass on differentials and matricies (posibly as env vars)
#compute the the U, V, W associated to P. 

#P_o would be good env variable
#PM would be good env variable
#RM would be good env variable
#g would be good env var
#diffs would be good env var
#the left point of the monodromy <- period matrix
#a g-1 devisor

RCV := proc(curve, x, y, t, R, RM, P0)    
    local p,q,Rfdf, Rfp, possibleXs, N, AofC, i,tempVal21,Ig,HalfNs,HalfMs,K_0, parityOfVec, m,mm,
    II,J,j,vecData,KP_0, DI:=Digits, possibleChar:={}, possiblePlaces,place, px,qx,degX; 
    if g=1 then
        return  [1/2+ RM[1,1]/2];
    end if;
    print("computing RCV");
    p := numer( R );	
    q := denom( R );
    Rfdf := resultant(curve, diff(curve,y), y);
    
    Rfp := resultant(curve, p,y);

    Digits:=Digits+EXTRA1;  
    print(Rfdf);
    print(Rfp);
    possibleXs:= [ infinity,0,seq(RootOf(Rfp,x,index=i),i=1..degree(Rfp)),
        seq(RootOf(Rfdf,x,index=i),i=1..degree(Rfdf))
    ];
    possibleXs:= convert~({op(simplify(fnormal(evalf~( possibleXs,Digits)),zero))},rational) minus {0,infinity};
    possibleXs :=[0, infinity,op(possibleXs)];
    Digits:=Digits-EXTRA1;
    N := 0;
    AofC := 0;
    for i from 1 to nops(possibleXs) do
        possiblePlaces := algcurves[puiseux](curve, x=possibleXs[i], y, Digits, t );
        for place in possiblePlaces do
            px := simplify(fnormal(evalf(collect(expand(subs( op(place), p) * diff( rhs(place[1]),t )),t)),Digits),zero);        
			qx:= simplify(fnormal(evalf(collect(expand(subs(op(place), q)),t )),Digits),zero); 
            px:= convert(px,polynom);
            qx:= convert(qx,polynom);
            degX := ldegree(px)-ldegree(qx); 
            if  degX < 0 then
                error "the differential was found to be meromorphic";
            elif degX > 0 then
                N:=N+degX;
                AofC := AofC 
                    +  degX*algcurves[AbelMap](curve,x,y,P0,place, t, Digits);
            fi;                            
        od;           
    od;
    if N<>2*g-2 then
        error cat(
                    "The number of zeros counted does not match the expected value:",
                    convert(N, string)," accounted for but should be ",convert(2*g-2, string));
    fi;
    
    print("Abel map of Cannon Div");
    print(AofC);
    AofC:=-convert(evalf(AofC),Vector)/2;
    tempVal21 := getSubObjects( [ seq(i, i=1..g) ] )[ 2..-1 ];
    Ig:=Id(g);
    HalfNs := [ Vector( g ), 
        seq( add( col( Ig/2, tempVal21[i][k] ), k = 1 .. nops(tempVal21[i])), i = 1 .. nops(tempVal21) ) ];  
    HalfMs := [ seq(prod( RM, HalfNs[i]), i = 1 .. nops(HalfNs) )]; 
    
    for i from 1 to 2^g do
        for j from 1 to 2^g do
            K_0:= convert(HalfNs[i]+HalfMs[j]+AofC, list);  
            parityOfVec:=(4*prod(tp(HalfNs[i]),HalfNs[j])) mod 2;
            parityOfVec := ifelse(parityOfVec=0, "even","odd"); 

            Digits:=2;
            m:=abs(evalf(RiemannTheta(K_0,RM)));       
            if fnormal(m,Digits) = 0. then                
                possibleChar := possibleChar union { [K_0,i,j, parityOfVec] };               
            fi;  
        od;
    od; 

    if nops(possibleChar)=0 then
        print(m);
        error "increase digits."   
    elif nops(possibleChar)>1 then
        while nops(possibleChar)> 1 and Digits < min(10,DI) do
            Digits:= Digits+2;
            for vecData in possibleChar do 
                m := abs(evalf(RiemannTheta(vecData[1],RM))); 
                if fnormal(m,Digits) <> 0. then                
                    possibleChar := possibleChar minus { vecData };               
                fi;
            od;
        od;
    fi;
    Digits:=DI;
    m:=abs(evalf(RiemannTheta(possibleChar[1][1],RM)));
    KP_0:=possibleChar[1][1];
    parityOfVec :=  possibleChar[1][4];  
    II:=  possibleChar[1][2];
    J:=   possibleChar[1][3]; 
    if nops(possibleChar)>1 then
        for i from 2 to nops(possibleChar) do
            mm:=abs(evalf(RiemannTheta(possibleChar[i][1],RM)));
            if m>=mm then 
                m:=mm;
                KP_0:=possibleChar[i][1];
                parityOfVec :=  possibleChar[i][4];  
                II:=  possibleChar[i][2];
                J:=   possibleChar[i][3]; 
            fi;
        od;
    fi;   
   
    print(cat("final pass yields ",m," with parity",parityOfVec));   
    print(HalfNs[II],HalfNs[J]);
    return  KP_0;  
        
end proc;
##Description
###INPUT
###OUTPUT
#this assumes that 2pi I is already pulled out of the vecr
#tdo
#-add user level printing
#RVC is hardest to compute and needs
# RM <- period natrix
#Differential and all of the possible zeros <- period matrix
#a place user gives infor curve,x,y,P,D optional


#compute RCV first and pass on differentials and matricies (posibly as env vars)
#compute the the U, V, W associated to P. 

#P_o would be good env variable
#PM would be good env variable
#RM would be good env variable
#g would be good env var
#diffs would be good env var
#the left point of the monodromy <- period matrix
#a g-1 devisor

RCV2 := proc(curve, x, y, t, R, RM, P0)    
    local vec,Ig,MM,i,ii,j, mm, mdrmy, PM, omegaSet, A, B; 
   
        vec:=[];
        Ig:=Id(g);
        for i from 1 to g do
            MM:= LinearAlgebra[Diagonal](Ig)-col(Ig,i);
            print(MM);
            vec:=[op(vec), 1/2 + RM[i,i]/2-prod(RM,MM)[i]/2];
        od;
        return  vec;

    
        
end proc;
