
#For g=2, U1 is any complex number and U2,V2 are normalized. Turly you start with normalization constants and derive the relations
#on V1,W1,W2 that show that if U is chosen you have two possibilities for V1 and then W1 and W2 are completely determined

WC_Mat := proc(RM::Matrix)
  local d, U, V, W, getPara, m, g, Solutions, solSet, paramaters:=[],DI:=max(10, Digits);  
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
      argument:=[op(argument), convert( simplify(fnormal( 2*prod( RM, Char[-1] ), DI ),zero), list) ]; 
      for j from i to g do
        Ucoef2:=[op(Ucoef2), combinat[multinomial]( 2, op(convert( col( Ig, i) + col( Ig, j), list) ))]; 
        if j>i then
          Char:=[op(Char), col( Ig, i )/2 + col( Ig, j )/2];
          argument:=[op(argument), convert( simplify(fnormal( 2*prod( RM, Char[-1] ), DI ),zero), list) ];  
        fi;
        ij:=[op(ij), [i,j] ];
        for k from j to g do
          if k>j and j>i then
            Char:=[op(Char), col( Ig, i )/2 + col( Ig, j )/2 + col( Ig, k )/2 ];
            argument:=[op(argument), convert( simplify(fnormal( 2*prod( RM, Char[-1] ), DI ),zero), list) ]; 
          fi;
          ijk:=[op(ijk), [i,j,k] ];
          Ucoef3:=[op(Ucoef3), combinat[multinomial]( 3, op(convert( col(Ig, i ) + col(Ig, j )+ col(Ig, k ), list) ))];
          for l from k to g do
            if l>k and k>j and j>i then
              Char:=[op(Char), col( Ig, i )/2 + col( Ig, j )/2 + col( Ig, k )/2 + col( Ig, l)/2 ];
              argument:=[op(argument), convert( simplify(fnormal( 2*prod( RM, Char[-1] ), DI ),zero), list) ]; 
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
          ThetaTable[i,1,j] := simplify(fnormal(RiemannTheta(argument[i], 2*RM,[ei], .1^DI )),zero);         
        else
          ThetaTable[i,1,j] := 0; #placeholder only
        fi;
      end do;
      #second derivatives
      for j from 1 to nops(ij) do     
        IJ := ij[j]; II := IJ[1]; J := IJ[2];
        ei := convert( col( Ig, II ), list); ej := convert( col( Ig, J ), list);       
        ThetaTable[i,2,j] := simplify(fnormal(RiemannTheta( argument[i], 2*RM, [ ei, ej ], .1^DI  )),zero); 
      end do;
      #third derivatives
      for j from 1 to nops(ijk) do       
        IJK := ijk[j]; II := IJK[1]; J := IJK[2]; K := IJK[3];             
        if i>1 then          
          ei := convert( col( Ig, II), list); ej := convert( col( Ig, J), list); ek := convert( col( Ig, K) ,list);            
          ThetaTable[i,3,j] := simplify(fnormal(RiemannTheta( argument[i], 2*RM,[ ei, ej, ek ], .1^DI )),zero);         
        else
          ThetaTable[i,3,j] := 0; #placeholder only
        fi;
      end do;  
      #fourth derivatives
      for j from 1 to nops(ijkl) do     
        IJKL := ijkl[j]; II := IJKL[1]; J := IJKL[2]; K := IJKL[3]; L := IJKL[4];      
        ei := convert( col( Ig, II), list); ej := convert( col( Ig, J), list); ek := convert( col( Ig, K), list); el := convert( col( Ig, L), list);           
        ThetaTable[i,4,j] := simplify(fnormal( RiemannTheta( argument[i], 2*RM, [ ei, ej, ek, el], .1^DI  )),zero);  
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
      ThetaTableEx[i,0] := simplify(fnormal(RiemannTheta( argumentEx[i], RM,[], .1^DI )),zero);
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
        ThetaTableEx[i,1,j] := simplify(fnormal(RiemannTheta( argumentEx[i], RM,[ei], .1^DI )),zero);
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
            
        ThetaTable[i,2,j] := simplify(fnormal(RiemannTheta( argument[i], 2*RM, [ ei, ej ] , .1^DI )),zero);
      end do;
       
      for j from 1 to nops(ij) do
        IJ := ij[j]; II := IJ[1]; J := IJ[2];
        ei := convert( col( Ig, II ), list); ej := convert( col( Ig, J ), list);            
        ThetaTableEx[i,2,j] := simplify(fnormal(RiemannTheta( argumentEx[i], RM, [ ei,ej ], .1^DI  )),zero);         
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
        ThetaTableEx[i,3,j] := simplify(fnormal(RiemannTheta( argumentEx[i], RM,[ ei, ej, ek ] , .1^DI )),zero);         
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
    if  fnormal(mm,DI-6) = 0.0 then
      error "Increase Digits: The submatrix needed is singular" ;   
    fi;
    #linsolve(alpha,beta,method=none,free=_t0,conjugate=true,	inplace=false,outputoptions=[]);
    if g>2 then
      Q := simplify(fnormal( linsolve(MatrixOfConstants, b, method=none, free=U[-1]) ),zero);
      CompatabilityCondition := collect(convert(det( MatrixOfConstantsBig),rational) ,U[-1]);
      valuU := evalf([solve( CompatabilityCondition  = 0, U[-1])]);
      #print(valuU);
    else
      Q := simplify(fnormal(linsolve(MatrixOfConstants, b, method=none) ),zero);
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
  g := rank(simplify(fnormal(evalf(RM),DI-6),zero));
  if g > 3 then
    error "either the genus is not correct or the given matrix is not of full rank. Genus must be 3 or less";  
  elif not IsMatrixShape(simplify(fnormal(evalf(RM),DI-6),zero), symmetric) then
    error "Matrix must be symmetric";
  elif not LinearAlgebra[IsDefinite](Im(simplify(fnormal(evalf(RM),DI-6),zero)), query = 'positive_definite') then
   error "Matrix must have positive definite imaginary part";    
  fi;
  if nargs>1 then
    Solutions := getPara(simplify(fnormal(evalf(RM),DI-3),zero),g,args[2..-1]);
    else
    Solutions:= getPara(simplify(fnormal(evalf(RM),DI-3),zero),g);
  fi;
  for solSet in Solutions do
    U := solSet[1];    
    V := solSet[2];
    W := convert(solSet[3], list);
    d := -solSet[4];
     if fnormal(LinearAlgebra[Norm](solSet[3]-solSet[5],infinity),DI-3) >0.0 then
      WARNING(cat(
        "W is not within 10^(-20) of the expected value: ",
      convert(fnormal(LinearAlgebra[Norm](solSet[3]-solSet[5],infinity),DI-3),string)
      ));
    fi;
    if nargs > 3 then
      checkWCError(g,RM,U,V,W, d,0,args[4]);     
    else
      checkWCError(g,RM,U,V,W, d,0);  
    fi;          
    paramaters:=[op(paramaters),[U,V,W, d,0]]; 
  od;     
  [RM,op(paramaters)];    
end proc;
