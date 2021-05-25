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

import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;

public class ConflictOrderingSearch implements Supplier<Procedure[]> {

    private final Supplier<IntVar> vs;
    private final Function<IntVar, Integer> valueS;
    private int timestamp=0;
    private PriorityQueue<Pair> nextVars;

    public ConflictOrderingSearch(Supplier<IntVar> variableSelector,
                                 Function<IntVar, Integer> valueSelector){
        this.vs = variableSelector;
        this.valueS = valueSelector;
        this.nextVars=new PriorityQueue<Pair>();
    }
    @Override
    public Procedure[] get() {
        IntVar X = nextVarToBranchOn();
        timestamp++;
        nextVars.add(new Pair(X,timestamp));
        return BranchingScheme.branch(
                ()-> equal(X, valueS.apply(X)),    // X = v
                ()-> notEqual(X, valueS.apply(X))  // X != v
        );
    }

    private IntVar nextVarToBranchOn() {
        if (timestamp == 0) {
            return vs.get();
        } else {
            IntVar nextvar = nextVars.peek().key;
            if(nextvar.isBound()){
                timestamp--;
                ArrayList<Pair> putBack = new ArrayList<Pair>();
                putBack.add(nextVars.poll());
                if(timestamp==0) // bah en fait pas besoin de putback
                    return vs.get();
                // find next unbound
                while (true){
                    Pair tmp = nextVars.poll();
                    if(tmp==null){
                        nextVars.addAll(putBack);
                        return vs.get();
                    }
                    putBack.add(tmp);
                    if(!tmp.key.isBound()) {
                        nextVars.addAll(putBack);
                        return tmp.key;
                    }
                }
            }
            return nextvar;
        }
    }
    private class Pair implements Comparable{
        IntVar key = null;
        int value = 0;

        Pair(IntVar newVar,int newVal){
            this.key = newVar;
            this.value = newVal;
        }

        @Override
        public int compareTo(Object o) {
            return this.value-((Pair)o).value;
        }
    }
}

