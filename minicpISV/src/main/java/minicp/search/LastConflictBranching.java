/*
 * mini-cp is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License  v3
 * as published by the Free Software Foundation.
 *
 * mini-cp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY.
 * See the GNU Lesser General Public License  for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with mini-cp. If not, see http://www.gnu.org/licenses/lgpl-3.0.en.html
 *
 * Copyright (c)  2018. by Laurent Michel, Pierre Schaus, Pascal Van Hentenryck
 */

package minicp.search;


import minicp.cp.BranchingScheme;

import static minicp.cp.Factory.equal;
import static minicp.cp.Factory.notEqual;
import minicp.engine.core.IntVar;
import minicp.util.Procedure;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.function.Supplier;

public class LastConflictBranching implements Supplier<Procedure[]> {

    private final Supplier<IntVar> vs;
    private final Function<IntVar, Integer> valueS;
    private IntVar lastConflictVar = null;

    public LastConflictBranching(Supplier<IntVar> variableSelector,
                                 Function<IntVar, Integer> valueSelector){
        this.vs = variableSelector;
        this.valueS = valueSelector;
    }
    @Override
    public Procedure[] get() {
        IntVar X = nextVarToBranchOn();
        lastConflictVar = X;
        return BranchingScheme.branch(
                ()-> equal(X, valueS.apply(X)),    // X = v
                ()-> notEqual(X, valueS.apply(X))  // X != v
                );
    }

    private IntVar nextVarToBranchOn() {
        if (lastConflictVar == null) {
            return vs.get();
        } else {
            if(lastConflictVar.isBound()){
                return vs.get();
            }
            return lastConflictVar;
        }
    }
}
