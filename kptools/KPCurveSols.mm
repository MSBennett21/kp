
#For g=2, U1 is any complex number and U2,V2 are normalized. Turly you start with normalization constants and derive the relations
#on V1,W1,W2 that show that if U is chosen you have two possibilities for V1 and then W1 and W2 are completely determined

WC_Curve := proc(f, x::symbol, y::symbol)
  local RM, d, U, V, W, getPara, XList, m, g, Solutions, solSet, paramaters:=[], DI:=max(10,Digits);  
  getPara := proc(RM,g)
    local 
    ThetaTable := table(), ThetaTableEx := table(),
    MatrixOfConstantsBig := Matrix( g*( g+1 )/2+2, g*( g+1 )/2+2 ),  MatrixOfConstants,
    MatrixOfConstantsInv,  WAttained,WExpected, CompatabilityCondition, valuU, Uu, TheSign, signSys, lambda, U1, U2, V1, Wsol, Vsol,
    theta := table(), theta2 := table(), theta3 := table(), theta4 := table(), Q, Qu, P, b := Vector( g*( g+1 )/2+1 ),
    i, j, k, l, II, J, K, L, IJ, IJK, IJKL, ei, ej, ek, el, Ig := Id(g), 
    ij := [], ijk := [], ijkl := [], Ucoef2 := [], Ucoef3 := [], Ucoef4 := [],      
    expK, Char := [ Vector(g) ], argument := [ convert( 2*prod( RM, Char[-1] ), list) ], argumentEx := [], 
    Z := [ seq( convert( cat( "z", convert( i, string ) ), symbol ), i = 1 .. g ) ],
    U := [ seq( convert( cat( "u", convert( i, string) ), symbol ), i = 1 .. g ) ],
    V := [ seq( convert( cat( "v", convert( i, string) ), symbol ), i = 1 .. g ) ],
    W := [ seq( convert( cat( "w", convert( i, string) ), symbol ), i = 1 .. g ) ], 
    sol := [], alphai, betai, tempChars,tempMat, integerSet, tempInt, mm, EquationSet:=[], wExpected, charsToBeUsed,
    signs:=[< 1, 1, 1> , < -1,1,1> , < 1,-1, 1> , < 1, 1,-1>], DI:=max(10,Digits);
    #print("populating lists relative to genus");
    #########################Populate/Processes Certain relative to case data:
    #-Characteristics-there is no order. Just combinations of the Identity Matrix multiplied by a scalar
    #-Argument-Order depecendent of on order of Chars. As soon as a char is made make the argument for the theta function
    #-Ucoefk-Coeffecients of kth derivatives in U: will be obtained from (x1+..xg)^k
    # ijkl-Order of the ijkl derivatives in z1,..,zg
    for i from 1 to g do
      Char:=[ op(Char), col( Ig, i )/2 ];
      argument:=[op(argument), convert( fnormal( 2*prod( RM, Char[-1] ), DI ), list) ]; 
      for j from i to g do
        Ucoef2:=[op(Ucoef2), combinat[multinomial]( 2, op(convert( col( Ig, i) + col( Ig, j), list) ))]; 
        if j>i then
          Char:=[op(Char), col( Ig, i )/2 + col( Ig, j )/2];
          argument:=[op(argument), convert( fnormal( 2*prod( RM, Char[-1] ), DI ), list) ]; 
        fi;
        ij:=[op(ij), [i,j] ];
        for k from j to g do
          if k>j and j>i then
            Char:=[op(Char), col( Ig, i )/2 + col( Ig, j )/2 + col( Ig, k )/2 ];
            argument:=[op(argument), convert(fnormal( 2*prod( RM, Char[-1] ), DI ), list) ]; 
          fi;
          ijk:=[op(ijk), [i,j,k] ];
          Ucoef3:=[op(Ucoef3), combinat[multinomial]( 3, op(convert( col(Ig, i ) + col(Ig, j )+ col(Ig, k ), list) ))];
          for l from k to g do
            if l>k and k>j and j>i then
              Char:=[op(Char), col( Ig, i )/2 + col( Ig, j )/2 + col( Ig, k )/2 + col( Ig, l)/2 ];
              argument:=[op(argument), convert(fnormal( 2*prod( RM, Char[-1] ), DI), list) ];
            fi;
            ijkl:=[op(ijkl), [i,j,k,l] ];
            Ucoef4:=[op(Ucoef4), combinat[multinomial]( 4, op(convert( col( Ig, i) + col( Ig, j) + col( Ig, k) + col( Ig, l) , list) ))];
          end;
        end;
      end;
    end; 

    ##CALCULATION OF DERIVATIVES::
    ## IDEA:: Theta constants derivatives of functions of the form fg, f an exponential and g a theta function all evaluated at certain locations (z=0 really).
    ## compute g's separate and eval at 0. Compute f's symbolically and then derivatves will be easy to compute symbolically. Ultimately the combinations will make up the matrix
    ## of theta constants: ith row is ith char jth row is the jth derivative in ij, followed by the 0th derivative and then the associtated 4th order U derivative
    # g=ThetaTable[i,j,k] theta of char i, jth order derivative, number k in the order (or k dne for j=0)
    # f= expK symbolic expresion of exponential
    for i from 1 to nops(Char) do      
      ThetaTable[i,0] := RiemannTheta( argument[i], 2*RM, [], .1^DI ); 
      expK := Zvec -> exp( 
            (2*Pi*I) * ( prod( tp(Char[i]), prod( RM, Char[i]) ) + prod( tp( Char[i] ), <op(Zvec)> ) )
            );
      #first derivatives
      for j from 1 to g do        
        ei:=convert(col(Ig,j),list);         
        if i>1 then
          ThetaTable[i,1,j] := RiemannTheta(argument[i], 2*RM,[ei], .1^DI );         
        else
          ThetaTable[i,1,j]:=0; #placeholder only
        fi;
      end do;
      #second derivatives
      for j from 1 to nops(ij) do     
        IJ := ij[j]; II := IJ[1]; J := IJ[2];
        ei := convert( col( Ig, II ), list); ej := convert( col( Ig, J ), list);       
        ThetaTable[i,2,j] := RiemannTheta( argument[i], 2*RM, [ ei, ej ], .1^DI  ); 
      end do;
      #third derivatives
      for j from 1 to nops(ijk) do       
        IJK := ijk[j]; II := IJK[1]; J := IJK[2]; K := IJK[3];             
        if i>1 then          
          ei := convert( col( Ig, II), list); ej := convert( col( Ig, J), list); ek := convert( col( Ig, K) ,list);            
          ThetaTable[i,3,j] := RiemannTheta( argument[i], 2*RM,[ ei, ej, ek ], .1^DI );         
        else
          ThetaTable[i,3,j] := 0; #placeholder only
        fi;
      end do;  
      #fourth derivatives
      for j from 1 to nops(ijkl) do     
        IJKL := ijkl[j]; II := IJKL[1]; J := IJKL[2]; K := IJKL[3]; L := IJKL[4];      
        ei := convert( col( Ig, II), list); ej := convert( col( Ig, J), list); ek := convert( col( Ig, K), list); el := convert( col( Ig, L), list);           
        ThetaTable[i,4,j] := RiemannTheta( argument[i], 2*RM, [ ei, ej, ek, el], .1^DI  );  
      od; 
      #Matrix construction: ith row-char, jth col-der
      for j from 1 to nops(ij)+1 do       
        if j=nops(ij)+1 then 
          MatrixOfConstantsBig[i,j] :=
           eval( subs( seq(Z[s]=0,s=1..g), 
           expK(Z)*ThetaTable[i,0]
           )); 
        else   
          IJ := ij[j]; II := IJ[1]; J := IJ[2];         
          MatrixOfConstantsBig[ i, j] :=
            ( eval( subs( seq(Z[s]=0,s=1..g),     
            diff( expK(Z), Z[II], Z[J] ) * ThetaTable[i,0] +
            diff( expK(Z), Z[J] ) * ThetaTable[i,1,II] +
            diff( expK(Z), Z[II] ) * ThetaTable[i,1,J] +
            expK(Z) * ThetaTable[i,2,j]
            )) 
            )/( (2*Pi*I)^2);             
        fi;
      od;
      #Matrix construction: ith row-char, final column
      for j from 1 to nops(ijkl) do
        IJKL := ijkl[j]; II := IJKL[1]; J := IJKL[2]; K := IJKL[3]; L := IJKL[4];        
        theta4[i,j]:=
          Ucoef4[j] * U[II] * U[J] * U[K] * U[L] *
          eval( subs( seq(Z[s]=0,s=1..g),       
          diff( expK(Z), Z[II], Z[J], Z[K], Z[L]) * ThetaTable[i,0] +
          diff( expK(Z), Z[II], Z[J], Z[K]) * ThetaTable[i,1,L] +
          diff( expK(Z), Z[II], Z[J], Z[L]) * ThetaTable[i,1,K] +
          diff( expK(Z), Z[II], Z[K], Z[L]) * ThetaTable[i,1,J] +
          diff( expK(Z), Z[II], Z[J]) * ThetaTable[i,2, searchList([K,L],ij) ] +
          diff( expK(Z), Z[II], Z[K]) * ThetaTable[i,2, searchList([J,L],ij) ] +
          diff( expK(Z), Z[II], Z[L]) * ThetaTable[i,2, searchList([J,K],ij) ] +
          diff( expK(Z), Z[II]) * ThetaTable[i,3, searchList([J,K,L],ijk) ] +
          diff( expK(Z), Z[J], Z[K]) * ThetaTable[i,2, searchList([II,L],ij) ] +
          diff( expK(Z), Z[J], Z[L]) * ThetaTable[i,2, searchList([II,K],ij) ] +
          diff( expK(Z), Z[K], Z[L]) * ThetaTable[i,2, searchList([II,J],ij) ] +          
          diff( expK(Z), Z[J], Z[K], Z[L]) * ThetaTable[i,1,II] +          
          diff( expK(Z), Z[J]) * ThetaTable[i,3, searchList([II,K,L],ijk) ] +
          diff( expK(Z), Z[K]) * ThetaTable[i,3, searchList([II,J,L],ijk) ] +
          diff( expK(Z), Z[L]) * ThetaTable[i,3, searchList([II,J,K],ijk) ] +
          expK(Z) * ThetaTable[i,4,j]      
          ));
        MatrixOfConstantsBig[i,-1] := MatrixOfConstantsBig[i,-1] +  theta4[i,j];
      end do;
      MatrixOfConstantsBig[i,-1] := MatrixOfConstantsBig[i,-1]/( (2*Pi*I)^4);   
    od;
    ###VERIFICATION USING ODD CHARS: 
    #Idea: The odd chars can be used to check behavior. They can also be used to check the W value- which will have the most error.
    #Follow same procedure as above to compute values of odd char'd theta at 0 (NOT THETA CONSTANTS!). 
    #need to compute enough so that there is no problem with the calculation: we think one for every non trivial even char is enough for now.
    #g=<2^g-1 so this should be enough to get W    
    for i from 1 to nops(Char)-1 do
      #populate values
      theta[i,xxx] := 0;
      theta[i,xx] := 0;
      theta[i,yy] := 0;
      theta[i,tt] := 0;
      theta[i,x] := 0;
      theta[i,y] := 0;
      theta[i,t] := 0; 
      #odd char/theta construction
      alphai:=Char[i+1];
      for j from 1 to g do                
        if alphai[j]>0 then
          betai := col( Ig, j )/2;
          break;
        fi;
      od;
      argumentEx := [op(argumentEx), convert( prod( RM, alphai )+ betai, list) ];
      ThetaTableEx[i,0] := RiemannTheta( argumentEx[i], RM,[], .1^DI );
      if fnormal(evalf(ThetaTableEx[i,0]), DI-9) <> 0.0 then
        WARNING( cat(
            "oddness violation: Odd Theta ", convert(i,string) , " is ", convert(fnormal(evalf(ThetaTableEx[i,0]), DI),string), " at Z=0."
            ));    
      fi;
      expK := Zvec -> exp( 
        (2*Pi*I) * ( prod( tp(alphai), ( prod( RM, alphai )/2 + betai + <op(Zvec)> ) ) )
        ); 

      #first derivatives      
      for j from 1 to g do
        ei := convert( col( Ig, j), list);     
        ThetaTableEx[i,1,j] := RiemannTheta( argumentEx[i], RM,[ei], .1^DI );
        theta[i,j] := 
          eval( subs(seq(Z[s]=0,s=1..g), 
          expK(Z) * ThetaTableEx[i,1,j]  
          ));
        theta[i,x] := theta[i,x] + U[j] * theta[i,j];
        theta[i,y] := theta[i,y] + V[j] * theta[i,j];
        theta[i,t] := theta[i,t] + W[j] * theta[i,j];
      od;  
      theta[i,x] := theta[i,x]/( 2*Pi*I );  
      theta[i,y] := theta[i,y]/( 2*Pi*I );  
      theta[i,t] := theta[i,t]/( 2*Pi*I );  
      #second derivatives
      for j from 1 to nops(ij) do     
            
        ThetaTable[i,2,j] := RiemannTheta( argument[i], 2*RM, [ ei, ej ] , .1^DI ); 
      end do;
       
      for j from 1 to nops(ij) do
        IJ := ij[j]; II := IJ[1]; J := IJ[2];
        ei := convert( col( Ig, II ), list); ej := convert( col( Ig, J ), list);            
        ThetaTableEx[i,2,j] := RiemannTheta( argumentEx[i], RM, [ ei,ej ], .1^DI  );         
        theta2[i,j] :=
          eval( subs( seq(Z[s]=0,s=1..g), 
          diff( expK(Z), Z[II]) * ThetaTableEx[i,1,J] + 
          diff( expK(Z), Z[J]) * ThetaTableEx[i,1,II] +
          expK(Z) * ThetaTableEx[i,2,j]  
          ));
        theta[i,xx] := theta[i,xx] +  Ucoef2[j] * U[II] * U[J] * theta2[i,j];
        theta[i,yy] := theta[i,yy] +  Ucoef2[j] * V[II] * V[J] * theta2[i,j];
        theta[i,tt] := theta[i,tt] +  Ucoef2[j] * W[II] * W[J] * theta2[i,j]; 
      end do;      
      theta[i,xx] := theta[i,xx]/( (2*Pi*I)^2 );     
      theta[i,yy] := theta[i,yy]/( (2*Pi*I)^2 );     
      theta[i,tt] := theta[i,tt]/( (2*Pi*I)^2 );  
      #third derivatives     
      for j from 1 to nops(ijk) do
        IJK := ijk[j]; II := IJK[1]; J := IJK[2]; K := IJK[3];         
        ei := convert( col( Ig, II), list); ej := convert( col( Ig, J), list); ek := convert( col( Ig, K) ,list);          
        ThetaTableEx[i,3,j] := RiemannTheta( argumentEx[i], RM,[ ei, ej, ek ] , .1^DI );         
        theta3[i,j] :=         
          eval( subs( seq(Z[s]=0,s=1..g), 
          diff( expK(Z), Z[II], Z[J]) * ThetaTableEx[i,1,K] +
          diff( expK(Z), Z[II], Z[K]) * ThetaTableEx[i,1,J] + 
          diff( expK(Z), Z[J], Z[K]) * ThetaTableEx[i,1,II] +
          diff( expK(Z), Z[II]) * ThetaTableEx[i,2, searchList([J,K],ij) ] +
          diff( expK(Z), Z[J]) * ThetaTableEx[i,2, searchList([II,K],ij) ] + 
          diff( expK(Z), Z[K]) * ThetaTableEx[i,2, searchList([II,J],ij) ] +
          expK(Z) * ThetaTableEx[i,3,j]
          ));
          theta[i,xxx] := theta[i,xxx] +  Ucoef3[j] * U[II] * U[J] * U[K] * theta3[i,j]; 
      end do;              
      theta[i,xxx] := theta[i,xxx]/( (2*Pi*I)^3 );
    od; 
    
    #Processing of input:
    # Each g requires a certain set of parameters. There are also optional commands
    if g=1 then
      if nargs >= 4 then
        if not  ( type( args[3], complexcons) or  type( args[4], complexcons) ) then
          error "U1, V1, must be given as complex numbers"
        fi;
        U1 := args[3];
        V1 := args[4];
        MatrixOfConstantsBig := subs(u1=args[3], MatrixOfConstantsBig);
      else        
        error "For a genus 1 solution, we require 3 paramaters B,U1,V1."   
      fi;         
    elif g=2  then
      if nargs >= 3 then
        if not type( args[3], complexcons) then
          error "U1 as complex numbers"
        fi;
        U1 := args[3];
      else  
        error "For a genus 2 solution, we require 6 paramaters B (accounts for 3), U1, U2, V1";
      fi;   
      MatrixOfConstantsBig := subs( u1=U1, u2=1, MatrixOfConstantsBig);
    elif g=3 then    
      if nargs >= 4 then
        if not ( type( args[3], complexcons) or type( args[4], complexcons) ) then
          error "U1, U2 must be given as complex numbers"
        fi;
        U1 := args[3];
        U2 := args[4]; 
      else  
        error "For a genus 3 solution, we require 9 paramaters B (accounts for 6), U1, U2, V1";
      fi;   
      MatrixOfConstantsBig := subs( u1=U1, u2=U2, MatrixOfConstantsBig);
    else
      error "g must be 1,2 or 3 for now.";
    fi; 
    ## Computation of values step 1:Find the approriate submatrix.
    ## We need the system Ax=b where (A|b)(<-x,1>)=0 is the linear system involving the matricies above
    # for genus <3, there is only one choice. For g>3 just take the first non-singular one

    charsToBeUsed := [ seq(1..g*(g+1)/2+1) ];
    for i from 1 to g*(g+1)/2+1 do 
      b[i] := MatrixOfConstantsBig[ charsToBeUsed[i],-1];
    od;   
    tempMat := fnormal(evalf( sm( MatrixOfConstantsBig, charsToBeUsed, charsToBeUsed)),DI);
    mm := abs( det(tempMat) );
    if 3 > g then
      MatrixOfConstants := fnormal(evalf( sm( MatrixOfConstantsBig, charsToBeUsed, charsToBeUsed)),DI);                  
    else 
      for i from 1 to g*(g+1)/2+2 do  
        tempChars := sort( [ op({seq(1..g*(g+1)/2+2)} minus {i}) ] );         
        tempMat := fnormal( evalf( sm( MatrixOfConstantsBig, charsToBeUsed, charsToBeUsed)), DI);        
        if abs( det(tempMat) ) > 0 then      
          mm := abs( det(tempMat) );
          charsToBeUsed := tempChars;
          MatrixOfConstants := tempMat;
          break;
        fi;      
      end do;       
    fi;      
    if  fnormal(mm,DI-9) = 0.0 then
      error "Increase Digits: The submatrix needed is singular" ;   
    fi;
    #linsolve(alpha,beta,method=none,free=_t0,conjugate=true,	inplace=false,outputoptions=[]);
    if g>2 then
      Q := linsolve(MatrixOfConstants, b, method=none, free=U[-1]);
      CompatabilityCondition := collect(convert(det( MatrixOfConstantsBig),rational) ,U[-1]);
      valuU := evalf([solve( CompatabilityCondition  = 0, U[-1])]);
      print(valuU);
    else
      Q := linsolve(MatrixOfConstants, b, method=none);
      valuU := ["placeholder"];
    fi;      
    
    for i from 1 to nops(valuU) do
      if g=3 then
        Uu := subs( U[1]=U1, U[2]=U2, U[-1]=simplify( fnormal( valuU[-1], DI ), zero ), U );                 
        Qu := fnormal( subs( seq( U[i]=Uu[i], i=1..g), Q), DI);                                 
      elif g=2 then
        Uu := [U1, 1]; 
        Qu := Q;
      else 
        Uu := [U1];
        Qu := Q;
      fi;        
      P := [];
      EquationSet := [];
      TheSign := [ 1, 1, 1 ];
      if g=1 then
        Vsol := [ V1 ];
      else
        for k from 1 to g-1 do
          for j from k+1 to g do 
            #sqrts in maple need to be carefully used due to the notion of -0 vrs +0       
            P:=[ op(P),  
            simplify( fnormal( 
              Uu[k]^2 * Qu[searchList([j,j],ij)] - Uu[k] * Uu[j] * Qu[searchList([k,j],ij)]  + Uu[j]^2 * Qu[searchList([k,k],ij)],
              DI),zero)         
            ];
          od;          
        od;         
        if g>2 then     
          mm := 100;            
          for signSys in signs do
            if fnormal( evalf( abs( add((-1)^(k-1) * signSys[k] * Uu[k] * sqrt(P[g+1-k]), k=1..g) ) ), DI-9 ) < mm then
              TheSign := [ signSys[1], signSys[2], signSys[3] ];
              mm := fnormal(evalf(abs( add((-1)^(k-1) * signSys[k] * Uu[k] * sqrt(P[g+1-k]), k=1..g) ) ), DI-9 );           
            fi;
          od;     
          if mm <> 0.0 then
            WARNING(
            cat("Sign Condition is not satisfied:", convert(mm,string))
            );
          fi;
          lambda := 2*I/sqrt(3) * 1/sqrt( simplify( fnormal( Uu[1]^2+Uu[2]^2+Uu[3]^2, DI), zero) );
          Vsol:=[
              simplify( fnormal( -lambda * ( Uu[3] * TheSign[2] * sqrt(P[2]) + Uu[2] * TheSign[1] * sqrt(P[1]) ), DI), zero),
              simplify( fnormal( lambda * ( Uu[1] * TheSign[1] * sqrt(P[1]) - Uu[3] * TheSign[3] * sqrt(P[3]) ), DI), zero),
              simplify( fnormal( lambda * ( Uu[2] * TheSign[3] * sqrt(P[3]) + Uu[1] * TheSign[2] * sqrt(P[2]) ) , DI), zero)
          ];
        else
            
            Vsol:=[-simplify(fnormal( 2*I/sqrt(3)*sqrt(P[1]), DI),zero),0]
        fi;
      fi;
      #print("computation of W, both types");
      Wsol:=[];
      for j from 1 to g do #-UjWj+3/4*Vj^2=-Q -> UjWj-3/4*Vj^2=-Q
        IJ := searchList([j,j],ij);
        Wsol := [op(Wsol),
          simplify(fnormal(evalf(   
            ( Qu[IJ] + 3/4 * Vsol[j]^2 )/Uu[j]             
          ),DI),zero)
        ];
      od; 
      integerSet := { seq(1..(nops(Char)-1)) };
      mm:=1;
      while nops(integerSet)> g do
        tempInt := min(integerSet);     
        for k in integerSet do
          if evalf( abs(mm) - abs(  subs(seq( U[i]=Uu[i], i=1..g) , theta[tempInt,x] ) ) )>0 then            
            tempInt := k:
            mm := evalf( subs( seq( U[i]=Uu[i], i=1..g) , theta[tempInt,x] ) );
          fi;
        od;
        integerSet := integerSet minus { tempInt };          
      od;
      for k in integerSet do       
        if fnormal( evalf( subs( seq( U[i]=Uu[i], i=1..g), theta[k,xx])), DI-9) <> 0.0 then
          WARNING( cat(
          "oddness violation: Odd Theta ", convert(k,string) , " has xx der ", 
          convert(fnormal( evalf( subs( seq( U[i]=Uu[i], i=1..g), theta[k,xx])), DI-9),string), " at Z=0."
          ));    
        fi;         
        if fnormal( evalf( subs( seq( V[i]=Vsol[i], i=1..g), theta[k,yy])), DI-9) <> 0.0 then
          WARNING( cat(
          "oddness violation: Odd Theta ", convert(k,string) , " has yy der ", 
          convert(fnormal( evalf( subs( seq( V[i]=Vsol[i], i=1..g), theta[k,yy])), DI-9), string), " at Z=0."
          ));    
        fi;
        if fnormal( evalf( subs( seq( U[i]=Uu[i], i=1..g), theta[k,x])), DI-9) = 0.0 then
          WARNING( cat(
          "The expected W may be wrong as not enough chars are chosen to keep theta_x nonzero: ",
          convert(fnormal( evalf( subs( seq( U[i]=Uu[i], i=1..g), theta[k,x]))), string)
          ));    
        fi; 
        if fnormal( evalf( subs( seq( V[i]=Vsol[i], i=1..g), theta[k,y])), DI-9) = 0.0 then
          WARNING( cat(
          "The expected W may be wrong as not enough chars are chosen to keep theta_y nonzero: ",
          convert(fnormal( evalf( subs( seq( V[i]=Vsol[i], i=1..g), theta[k,y]))), string)
          ));    
        fi;           
      od; 
      EquationSet := [seq(         
        subs( seq( U[k]=Uu[k], k=1..g), seq( V[k]=Vsol[k], k=1..g),
        -4*theta[integerSet[i],xxx]*theta[integerSet[i],x]-3*(theta[integerSet[i],y])^2+
        4*theta[integerSet[i],t]*theta[integerSet[i],x])=0,
        i=1..g)];   
      wExpected := solve({op(EquationSet)},[seq(W[j],j=1..g)])[1];
      wExpected := [seq(
        simplify(fnormal(evalf(rhs(wExpected[j])),DI),zero),
        j=1..g)           
      ];
      sol:=
        [op(sol),
        [Uu,Vsol,<seq(Wsol[j],j=1..g)>,simplify(fnormal(Qu[-1],DI),zero),<seq(wExpected[j],j=1..g)>]
        ];
    od;
    return sol;  
  end proc;
  if algcurves[genus](f,x,y) > 3 then
    error `genus is not 3 or less`;
  else
    g:=algcurves[genus](f,x,y);
    #print(cat("genus is ",convert(g, string)));
  fi;
  Digits := Digits+9;
	RM := fnormal(algcurves[periodmatrix](f,x,y,'Riemann'),DI);
  Digits := Digits-9;  
  #print("Riemann Matrix:"  );
  #print(RM);
  Solutions:= getPara(RM,g,args[4..-1]);
  #print(Solutions);
  for solSet in Solutions do
    U:=solSet[1];    
    V:=solSet[2];
    W:=convert(solSet[3], list);
    d:=-solSet[4];
    if fnormal(LinearAlgebra[Norm](solSet[3]-solSet[5],infinity),DI-9) >0.0 then
      WARNING(cat(
        "W is not within the desired limits of the expected value: ",
      convert(fnormal(LinearAlgebra[Norm](solSet[3]-solSet[5],infinity),DI-9),string)
      ));
    fi;
    if (nargs > 5 ) then
      checkWCError(g,RM,U,V,W, d,0,args[6]);       
    else
      checkWCError(g,RM,U,V,W, d,0);  
    fi;
          
      paramaters:=[op(paramaters),[U,V,W, d,0]]; 
  od;   
  [RM,op(paramaters)];    
end proc;

checkWCError := proc( g, RM::Matrix, U::list, V::list, W::list, d::complexcons, c:=0, sampleSize:=1 )
  local ThetaTable1:=table(), ThetaTable2:=table(), tval1:=table(), tval2:=table(), argument1:=table(), argument2:=table(),
    i, j, k, l, II, J, K, L, IJ, IJK, IJKL, ei, ej, ek, el, ij:=[], ijk:=[], ijkl:=[], Ig:=Id(g), sampleVals:=[], 
    Ucoef2:=[], Ucoef3:=[], Ucoef4:=[], dir1:=Vector(g), dir2:=Vector(g), Mat:=Matrix( [ convert( dir1, list ), convert( dir2, list ) ] ),
    Theta1, Theta2,Theta1X:=0, Theta1T:=0, Theta1Y:=0, Theta1XX:=0, Theta1XT:=0, Theta1YY:=0, Theta1XXX:=0, Theta1XXXX:=0,
    Theta2X:=0, Theta2T:=0, Theta2Y:=0, Theta2XX:=0, Theta2XT:=0, Theta2YY:=0, Theta2XXX:=0, Theta2XXXX:=0,
    DI:=max(10, Digits), tol:= 0.1^DI; 
    for i from 1 to g do     
        for j from i to g do
            Ucoef2 := [ op(Ucoef2), combinat[multinomial]( 2, op( convert(col(Ig, i) + col(Ig, j), list ) ) ) ];  
            ij := [op(ij), [i,j] ]; 
            for k from j to g do        
                ijk := [ op(ijk), [i,j,k] ];
                Ucoef3 := [ op(Ucoef3), combinat[multinomial]( 3, op( convert(col(Ig, i) + col(Ig, j) + col(Ig, k), list ) ) ) ];
                for l from k to g do           
                    ijkl := [ op(ijkl), [i,j,k,l] ];
                    Ucoef4 := [ op(Ucoef4), combinat[multinomial]( 4, op( convert( col(Ig, i) + col(Ig, j) + col(Ig, k) + col(Ig, l), list ) ) ) ];
                end;
            end;
        end;
    end;   
    for i from 1 to sampleSize + 1 do
        Theta1X := 0;Theta1T := 0; Theta1Y := 0;
        Theta2X := 0;Theta2T := 0; Theta2Y := 0; 
        Theta1XX := 0; Theta1XT := 0; Theta1YY := 0; 
        Theta2XX := 0; Theta2XT := 0; Theta2YY := 0; 
        Theta1XXX := 0; Theta1XXXX := 0;
        Theta2XXX := 0; Theta2XXXX := 0;
        if i=1 then
            tval1[i] := 0;
            tval2[i] := 0;
        else
            tval1[i] := (-1)^rand()*evalf( rand()/( 10^12 - 1 ));
            tval2[i] := (-1)^rand()*evalf( rand()/( 10^12 - 1 ));
        fi;    
        argument1 := convert( <op(U)> * tval1[i] + <op(V)> * 0 + <op(W)> * tval2[i], list );
        argument2 := convert( <op(U)> * 0 + <op(V)> * tval2[i] + <op(W)> * tval1[i], list );
        #print("KP1Errorvals");
        #g1KPError(tval1[i],0, tval2[i], RM, U[1], V[1], W[1], d, 0 );
        #g1KPError(0,tval2[i], tval1[i], RM, U[1], V[1], W[1], d, 0 );
        ThetaTable1[i,0] := RiemannTheta( argument1, RM, [], tol, output=list )[2];
        ThetaTable2[i,0] := RiemannTheta( argument2, RM, [], tol, output=list )[2];
        Theta1 := ThetaTable1[i,0];
        Theta2 := ThetaTable2[i,0];        
        for j from 1 to g do 
            ei := convert( col(Ig,j),list );  
            ThetaTable1[ i, 1, j ] := RiemannTheta( argument1, RM, [ei], tol, output=list )[2];
            ThetaTable2[ i, 1, j ] := RiemannTheta( argument2, RM, [ei], tol, output=list )[2];
            Theta1X := Theta1X + U[j] * ThetaTable1[i,1,j];
            Theta2X := Theta2X + U[j] * ThetaTable2[i,1,j];
            Theta1T := Theta1T + W[j] * ThetaTable1[i,1,j];
            Theta2T := Theta2T + W[j] * ThetaTable2[i,1,j];
            Theta1Y := Theta1Y + V[j] * ThetaTable1[i,1,j];
            Theta2Y := Theta2Y + V[j] * ThetaTable2[i,1,j];
        od; 
        for j from 1 to nops(ij) do  
            IJ := ij[j]; II := IJ[1]; J := IJ[2]; 
            ei := convert( col(Ig, II), list ); ej := convert( col(Ig, J), list );
            ThetaTable1[ i, 2, j ] := RiemannTheta( argument1, RM, [ ei, ej ], tol, output=list )[2];
            ThetaTable2[ i, 2, j ] := RiemannTheta( argument2, RM, [ ei, ej ], tol, output=list )[2];
            Theta1XX := Theta1XX + Ucoef2[j]*U[II]*U[J]*ThetaTable1[ i, 2, j ];
            Theta2XX := Theta2XX + Ucoef2[j]*U[II]*U[J]*ThetaTable2[ i, 2, j ];
            Theta1YY := Theta1YY + Ucoef2[j]*V[II]*V[J]*ThetaTable1[ i, 2, j ];
            Theta2YY := Theta2YY + Ucoef2[j]*V[II]*V[J]*ThetaTable2[ i, 2, j ];
            if II<J then
                Theta1XT := Theta1XT + ( U[II]*W[J] + U[J]*W[II] )*ThetaTable1[i, 2, j ];
                Theta2XT := Theta2XT + ( U[II]*W[J] + U[J]*W[II] )*ThetaTable2[i, 2, j ];
            else
                Theta1XT := Theta1XT + U[II]*W[J]*ThetaTable1[ i, 2, j ];
                Theta2XT := Theta2XT + U[II]*W[J]*ThetaTable2[ i, 2, j ];
            fi;      
        od;
        for j from 1 to nops(ijk) do  
            IJK := ijk[j]; II := IJK[1];  J := IJK[2]; K := IJK[3];
            ei := convert( col(Ig, II) , list ); ej := convert( col(Ig, J), list ); ek := convert( col(Ig, K),list );
            ThetaTable1[ i, 3, j ] := RiemannTheta( argument1, RM, [ ei, ej, ek ], tol, output=list )[2];
            ThetaTable2[ i, 3, j ] := RiemannTheta( argument2, RM, [ ei, ej, ek ], tol, output=list )[2];
            Theta1XXX := Theta1XXX + Ucoef3[j]*U[II]*U[J]*U[K]*ThetaTable1[ i, 3, j ];
            Theta2XXX := Theta2XXX + Ucoef3[j]*U[II]*U[J]*U[K]*ThetaTable2[ i, 3, j ];
        od;
        for j from 1 to nops(ijkl) do
            IJKL := ijkl[j]; II := IJKL[1]; J := IJKL[2]; K := IJKL[3]; L := IJKL[4];
            ei := convert( col(Ig, II), list); ej := convert( col(Ig, J), list );
            ek := convert( col(Ig, K), list ); el := convert( col(Ig, L), list );
            ThetaTable1[ i, 4, j ] := RiemannTheta( argument1, RM, [ ei, ej, ek, el ], tol, output=list )[2];
            ThetaTable2[ i, 4, j ] := RiemannTheta( argument2, RM, [ ei, ej, ek, el ], tol, output=list )[2];
            Theta1XXXX := Theta1XXXX + Ucoef4[j]*U[II]*U[J]*U[K]*U[L]*ThetaTable1[ i, 4, j ];
            Theta2XXXX := Theta2XXXX + Ucoef4[j]*U[II]*U[J]*U[K]*U[L]*ThetaTable2[ i, 4, j ];      
        end do;  
        Theta1X := Theta1X/( 2*Pi*I );
        Theta2X := Theta2X/( 2*Pi*I );  
        Theta1Y := Theta1Y/( 2*Pi*I );
        Theta2Y := Theta2Y/( 2*Pi*I );  
        Theta1T := Theta1T/( 2*Pi*I );
        Theta2T := Theta2T/( 2*Pi*I ); 
        Theta1XX := Theta1XX/( (2*Pi*I)^2 ); 
        Theta2XX := Theta2XX/( (2*Pi*I)^2 ); 
        Theta1YY := Theta1YY/( (2*Pi*I)^2 );
        Theta2YY := Theta2YY/( (2*Pi*I)^2 );  
        Theta1XT := Theta1XT/( (2*Pi*I)^2 ); 
        Theta2XT := Theta2XT/( (2*Pi*I)^2 );   
        Theta1XXX := Theta1XXX/( (2*Pi*I)^3 ); 
        Theta2XXX := Theta2XXX/( (2*Pi*I)^3 );
        Theta1XXXX := Theta1XXXX/( (2*Pi*I)^4 ); 
        Theta2XXXX := Theta2XXXX/( (2*Pi*I)^4 );         
        sampleVals := [
            op(sampleVals),
            max( 
                abs(evalf(
                    Theta1XXXX*Theta1 - 4*Theta1XXX*Theta1X + 3*Theta1XX^2 + 4*Theta1X*Theta1T - 4*Theta1XT*Theta1 + 
                    3*Theta1YY*Theta1 - 3*Theta1Y^2 + 8*d*Theta1^2
                )),
                abs(evalf(
                    Theta2XXXX*Theta2 - 4*Theta2XXX*Theta2X + 3*Theta2XX^2 + 4*Theta2X*Theta2T - 4*Theta2XT*Theta2 + 
                    3*Theta2YY*Theta2 - 3*Theta2Y^2 + 8*d*Theta2^2
                ))
                )
            ];         
    od;
  if max( op(sampleVals) ) > 10^( -(DI-9) ) then
    WARNING(
            cat(
            "Bilinear form is not as small as desired for given number of digits:",
            convert( simplify(
                max( op(sampleVals)), 
                zero ),
                string )
            )
        );
  fi;  
end proc;

KPgSol := proc(x , y, t, RM::Matrix, U::list, V::list, W::list, c:=0)
    local ThetaTable:=table(), g:=nops(U), argument:=convert( <op(U)>*x + <op(V)>*y + <op(W)>*t ,list ), Ig:=Id(g),
        i, j, II, J, IJ, ei, ej, ij:=[], Ucoef2:=[],
        Theta, ThetaX:=0, ThetaXX:=0,DI:=max(10, Digits), tol:= 0.1^DI;   
    for i from 1 to g do     
        for j from i to g do
            Ucoef2 := [ 
                op(Ucoef2), combinat[multinomial]( 2, op( convert( col(Ig, i) + col(Ig, j), list ) ) ) 
                ];  
            ij := [ op(ij), [i,j] ];       
        end;
    end;  
    ThetaTable[0] := RiemannTheta( argument, RM, [], tol, output=list )[2];
    Theta := ThetaTable[0];
    for j from 1 to g do 
        ei := convert( col(Ig, j), list );  
        ThetaTable[ 1, j ] := RiemannTheta( argument, RM, [ei], tol, output=list )[2];
        ThetaX := ThetaX + U[j]*ThetaTable[ 1, j ];
    od; 
    for j from 1 to nops(ij) do  
        IJ := ij[j]; II := IJ[1]; J := IJ[2]; 
        ei := convert( col(Ig,II), list ); ej := convert( col(Ig,J), list );
        ThetaTable[ 2,j ] := RiemannTheta( argument, RM, [ei,ej], tol, output=list )[2];
        ThetaXX := ThetaXX + Ucoef2[j]*U[II]*U[J]*ThetaTable[ 2, j ];   
    od;    
    ThetaX := ThetaX/( 2*Pi*I );   
    ThetaXX := ThetaXX/( (2*Pi*I)^2 );    
    2 * ( ThetaXX*Theta - (ThetaX)^2 )/( (Theta)^2 ) + c;
end proc;
