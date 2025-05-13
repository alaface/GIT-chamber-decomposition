// Input: an ideal I or an integer 
// Output: The F-faces indices of I.
//
// Description:
// - If I is an integer, it returns all the non-empty subsets of {1, ..., I}.
// - If I is an ideal, it computes the F-faces indices based on Remark 3.1.1.11.

Ffaces := function(I)
    if Type(I) eq RngIntElt then
        return Subsets({1..I}) diff {{}};

    elif IsPrincipal(I) then
	g := Basis(I)[1];
	R := Parent(I.1);
	n := Rank(R);
	faces := {};
	for S in Subsets({1..n}) diff {{}} do
		gS := Evaluate(g, [(i in S) select R.i else 0 : i in [1..n]]);
        	if #Monomials(gS) ne 1 then
            		Include(~faces, S);
        	end if;
    end for;
    else
	B := Basis(I);
	R := Parent(I.1);
	n := Rank(R);
	faces := {};

	for S in Subsets({1..n}) diff {{}} do
        	BS := [Evaluate(g, [(i in S) select R.i else 0 : i in [1..n]]) : g in B];
        	if &*[R.i : i in S] notin Radical(Ideal(BS)) then
        	    Include(~faces, S);
        	end if;
    	end for;
    end if;
    return faces;
end function;


// Input: A grading matrix Q
// Output: The effective cone of the grading matrix.
//
// Description:
// The effective cone is computed using the transpose of Q to create a list 
// of lattice vectors and constructing the corresponding cone.

Eff := function(Q)
    n := Ncols(Q);
    K := ToricLattice(Nrows(Q));
    return Cone([K!Eltseq(r) : r in Rows(Transpose(Q))]);
end function;


// Input: A grading matrix Q
// Output: The moving cone of the grading matrix.
//
// Description:
// The moving cone is obtained as the intersection of cones obtained by 
// removing one generator at a time.

Mov := function(Q)
    n := Ncols(Q);
    K := ToricLattice(Nrows(Q));
    L := [K!Eltseq(r) : r in Rows(Transpose(Q))];
    return &meet([Cone([L[j] : j in Remove([1..#L], i)]) : i in [1..#L]]);
end function;


// Input: F-faces and a grading matrix Q
// Output: The orbit cones associated with the F-faces.
//
// Description:
// For each subset in F, constructs a cone using the corresponding lattice vectors.

OrbitCones := function(F, Q)
    n := Ncols(Q);
    K := ToricLattice(Nrows(Q));
    w := [K!Eltseq(r) : r in Rows(Transpose(Q))];

    if #F eq 0 then 
        F := Subsets({1..n}) diff {{}};
    end if;

    return {Cone([w[i] : i in S]) : S in F};
end function;


// Input: Orbit cones and a class w
// Output: The GIT chamber containing the class w.
//
// Description:
// Computes the intersection of all cones in the orbit cones that contain w.

GitChamber := function(orb, w)
    K := Ambient(Random(orb));
    w := K!Eltseq(w);
    return &meet{C : C in orb | w in C};
end function;


// Input: Orbit cones and a class w
// Output: The bunch of cones containing the class w.
//
// Description:
// Returns all cones from the orbit cones that contain w.

BunchCones := function(orb, w)
    return {C : C in orb | w in C};
end function;


// Input: A bunch of cones and two classes w1, w2
// Output: A boolean value.
//
// Description:
// Returns true if the two classes w1 and w2 have the same stable base locus.

SameSbl := function(bun, w1, w2)
    K := Ambient(Random(bun));
    w1 := K!Eltseq(w1);
    w2 := K!Eltseq(w2);
    return {C : C in bun | w1 in C} eq {C : C in bun | w2 in C};
end function;


// Input: Orbit cones
// Output: The GIT fan.
//
// Description:
// Constructs the GIT fan by iteratively finding chambers and updating 
// the fan until all chambers are found.

GitFan := function(orb)
    Eff := Cone(&cat[Rays(C) : C in orb]);
    SH := {SupportingHyperplane(Eff, C) : C in Facets(Eff)};

    repeat
        W := Random(orb);
        K := Ambient(W);
        w := &+Rays(W);
        la := &meet{C : C in orb | w in C};
    until Dimension(la) eq Dimension(K);

    L := {la};
    F := {C : C in Facets(la) | SupportingHyperplane(la, C) notin SH};
    if #F eq 0 then return L; end if;

    repeat
        ff := Random(F);
        la := Random([C : C in L | IsFace(C, ff)]); 
        H := K!SupportingHyperplane(la, ff);
        w := &+Rays(ff);

        if Dimension(Cone([w, w+H]) meet la) eq 0 then 
            e := 1; 
        else 
            e := -1;
        end if;

        n := 1;
        repeat
            u := w + e/10^n * H;
            n := n + 1;
        until u in Eff;

        repeat
            lb := &meet{C : C in orb | u in C};
            u := w + e/10^n * H;
            n := n + 1;
        until Dimension(lb) eq Dimension(K) and lb meet la eq ff;

        L := L join {lb};
        Fb := {C : C in Facets(lb) | SupportingHyperplane(lb, C) notin SH};
        F := (F join Fb) diff (F meet Fb);
    until IsEmpty(F);

    return L;
end function;
