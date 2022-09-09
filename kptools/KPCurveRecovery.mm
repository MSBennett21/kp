GetCurve := proc(RM::Matrix,z) 
    local g:=rank(RM), ThetaTable := table(), 
    MatrixOfConstantsBig := Matrix( g*( g+1 )/2+2, g*( g+1 )/2+2 ),   
    theta := table(), theta2 := table(), theta4 := table(), 
    i, j, k, l, II, J, K, L, IJ, IJK, IJKL, 
    ei, ej, ek, el, ij := [], ijk := [], ijkl := [], Ig := Id(g),   
    expK, Z := [ seq( convert( cat( "z", convert( i, string ) ), symbol ), i = 1 .. g ) ],     
    Char := [ Vector(g) ], argument := [ convert( 2*prod( RM, Char[-1] ), list) ],
    Ucoef2 := [], Ucoef3 := [], Ucoef4 := [],
    U := [ seq( convert( cat( "x", convert( i, string) ), symbol ), i = 1 .. g-1 ), z ];
    
    #populate certain variables
    for i from 1 to g do
        Char:=[ op(Char), col(Ig,i)/2]; #standard basis vectors are char elements
        argument:=[op(argument), convert(simplify(fnormal(2*prod(RM,Char[-1])),zero),list)]; # need to make argument to go in theta for newly formed char
        for j from i to g do
            Ucoef2:=[op(Ucoef2),combinat[multinomial](2, op(convert(col(Ig, i) + col(Ig, j), list)))]; #the coeeficients will be coeef of xixj in:(x1+..xg)^2 
            if j>i then
                Char:=[op(Char),col(Ig,i)/2+col(Ig,j)/2]; #if j>i then this is componetent with two nonzero components::half char
                argument:=[op(argument), convert(simplify(fnormal(2*prod(RM,Char[-1])),zero),list)];#need to add this 
            fi;
            ij:=[op(ij), [i,j]]; #construct ij in canonical order
            for k from j to g do
                if k>j and j>i then
                    Char:=[op(Char), col(Ig,i)/2+col(Ig,j)/2+col(Ig,k)/2];
                    argument:=[op(argument), convert(simplify(fnormal(2*prod(RM,Char[-1])),zero),list)];
                fi;
                ijk:=[op(ijk), [i,j,k]];
                Ucoef3:=[op(Ucoef3),combinat[multinomial](3, op(convert(col(Ig, i) + col(Ig, j)+ col(Ig, k), list)))];
                for l from k to g do
                    if l>k and k>j and j>i then
                        Char:=[op(Char), col(Ig,i)/2+col(Ig,j)/2+col(Ig,k)/2+col(Ig,l)/2];
                        argument:=[op(argument),convert(simplify(fnormal(2*prod(RM,Char[-1])),zero),list)];
                    fi;
                    ijkl:=[op(ijkl), [i,j,k,l]];
                    Ucoef4:=[op(Ucoef4), combinat[multinomial](4, op(convert(col(Ig, i) + col(Ig, j)+ col(Ig, k)+col(Ig, l), list)))];
                end;
            end;
        end;
    end; 
    #THETA CONSTANTS   
    for i from 1 to nops(Char) do      
        ThetaTable[i,0] := fnormal(evalf(RiemannTheta(argument[i], 2*RM,[]))); #need theta portion of theta constant
        expK:=Zvec-> exp( 
            (2*Pi*I) * ( prod(tp(Char[i]),prod(RM,Char[i])) + prod(tp(Char[i]), <op(Zvec)>))
            );      
        for j from 1 to g do        
            ei:=convert(col(Ig,j),list);         
            if i>1 then
                ThetaTable[i,1,j]:=fnormal(evalf(RiemannTheta(argument[i], 2*RM,[ei])));         
            else
                ThetaTable[i,1,j]:=0; 
            fi;
        end do;
        for j from 1 to nops(ij) do     
            IJ:=ij[j]; II:=IJ[1]; J:=IJ[2];
            ei:=convert(col(Ig,II),list); ej:=convert(col(Ig,J),list);       
            ThetaTable[i,2,j]:=fnormal(evalf(RiemannTheta(argument[i], 2*RM,[ei,ej]))); 
        end do;      
        for j from 1 to nops(ijk) do       
            IJK:=ijk[j];II:=IJK[1];J:=IJK[2]; K:=IJK[3];             
            if i>1 then          
                ei:=convert(col(Ig,II),list); ej:=convert(col(Ig,J),list); ek:=convert(col(Ig,K),list);            
                ThetaTable[i,3,j]:=fnormal(evalf(RiemannTheta(argument[i], 2*RM,[ei,ej,ek])));         
            else
                ThetaTable[i,3,j]:=0;
            fi;
        end do;        
        for j from 1 to nops(ijkl) do     
            IJKL:=ijkl[j];II:=IJKL[1];J:=IJKL[2];K:=IJKL[3];L:=IJKL[4];      
            ei:=convert(col(Ig,II),list); ej:=convert(col(Ig,J),list); ek:=convert(col(Ig,K),list); el:=convert(col(Ig,L),list);           
            ThetaTable[i,4,j]:=fnormal(evalf(RiemannTheta(argument[i], 2*RM,[ei,ej,ek,el])));  
        od;      
        for j from 1 to nops(ij)+1 do        
            if j=nops(ij)+1 then 
                MatrixOfConstantsBig[i,j]:=simplify(fnormal(
                evalf(eval(subs(seq(Z[s]=0,s=1..g), 
                expK(Z)*ThetaTable[i,0])))
                ),zero); 
            else   
                IJ:=ij[j];II:=IJ[1];J:=IJ[2];          
                MatrixOfConstantsBig[i,j]:=simplify(fnormal(
                evalf(eval(subs(seq(Z[s]=0,s=1..g),     
                diff(expK(Z),Z[II],Z[J]) * ThetaTable[i,0]
                + diff(expK(Z),Z[J]) * ThetaTable[i,1,II] #i
                + diff(expK(Z),Z[II]) * ThetaTable[i,1,J] #j
                + expK(Z) * ThetaTable[i,2,j] #ij
                ))/((2*Pi*I)^2))
                ),zero);     
            fi;
        od;
        for j from 1 to nops(ijkl) do
            IJKL:=ijkl[j];II:=IJKL[1];J:=IJKL[2];K:=IJKL[3];L:=IJKL[4];       
            theta4[i,j]:=
            Ucoef4[j]*U[II] * U[J] * U[K] * U[L]*
            eval(subs(seq(Z[s]=0,s=1..g),       
            diff(expK(Z),Z[II],Z[J],Z[K],Z[L]) * ThetaTable[i,0] +
            diff(expK(Z),Z[II],Z[J],Z[K]) * ThetaTable[i,1,L] +
            diff(expK(Z),Z[II],Z[J],Z[L]) * ThetaTable[i,1,K] +
            diff(expK(Z),Z[II],Z[K],Z[L]) * ThetaTable[i,1,J] +
            diff(expK(Z),Z[J],Z[K],Z[L]) * ThetaTable[i,1,II] +
            diff(expK(Z),Z[II],Z[J]) * ThetaTable[i,2,searchList([K,L],ij)] +
            diff(expK(Z),Z[II],Z[K]) * ThetaTable[i,2,searchList([J,L],ij)] +
            diff(expK(Z),Z[II],Z[L]) * ThetaTable[i,2,searchList([J,K],ij)] +
            diff(expK(Z),Z[J],Z[K]) * ThetaTable[i,2,searchList([II,L],ij)] +
            diff(expK(Z),Z[J],Z[L]) * ThetaTable[i,2,searchList([II,K],ij)] +
            diff(expK(Z),Z[K],Z[L]) * ThetaTable[i,2,searchList([II,J],ij)] +
            diff(expK(Z),Z[II]) * ThetaTable[i,3,searchList([J,K,L],ijk)] +
            diff(expK(Z),Z[J]) * ThetaTable[i,3,searchList([II,K,L],ijk)] +
            diff(expK(Z),Z[K]) * ThetaTable[i,3,searchList([II,J,L],ijk)] +
            diff(expK(Z),Z[L]) * ThetaTable[i,3,searchList([II,J,K],ijk)] +
            expK(Z)* ThetaTable[i,4,j]      
            ));
            MatrixOfConstantsBig[i,-1]:= MatrixOfConstantsBig[i,-1] +  theta4[i,j];
        end do;
        MatrixOfConstantsBig[i,-1]:=simplify(fnormal(evalf(
        MatrixOfConstantsBig[i,-1]/((2*Pi*I)^4)
        )),zero);   
    od;
    if g=3 then
        return det(MatrixOfConstantsBig): 
    else
        return "nathan jason";
    fi;         
end proc;
