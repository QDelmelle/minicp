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


package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.IntVar;
import minicp.state.StateInt;
import minicp.util.GraphUtil;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;

import static minicp.cp.Factory.allDifferent;

/**
 * Hamiltonian Circuit Constraint with a successor model
 */
public class Circuit extends AbstractConstraint {

    private final IntVar[] x;
    private final StateInt[] dest;
    private final StateInt[] orig;
    private final StateInt[] lengthToDest;


    private GraphUtil.Graph g = new GraphUtil.Graph() {
        @Override
        public int n() {
            return dest.length;
        }
        @Override
        public Iterable<Integer> in(int idx) {
            ArrayList ret = new ArrayList<>();
            for(int i=0; i<n(); i++){
                if(x[i].contains(idx))
                    ret.add(i);
            }
            return ret;
        }
        @Override
        public Iterable<Integer> out(int idx) {
            ArrayList ret = new ArrayList<>();
            for(int i=0; i<n(); i++){
                if(x[idx].contains(i))
                    ret.add(i);
            }
            return ret;
        }
    };

    /**
     * Creates an Hamiltonian Circuit Constraint
     * with a successor model.
     *
     * @param x the variables representing the successor array that is
     *          {@code x[i]} is the city visited after city i
     */
    public Circuit(IntVar[] x) {
        super(x[0].getSolver());
        assert (x.length > 0);
        this.x = x;
        dest = new StateInt[x.length];
        orig = new StateInt[x.length];
        lengthToDest = new StateInt[x.length];
        for (int i = 0; i < x.length; i++) {
            dest[i] = getSolver().getStateManager().makeStateInt(i);
            orig[i] = getSolver().getStateManager().makeStateInt(i);
            lengthToDest[i] = getSolver().getStateManager().makeStateInt(0);
        }
    }

    @Override
    public void post() {
        getSolver().post(allDifferent(x));
        // TODO
        for (int i=0; i<x.length; i++){
            final int w = i;
            x[i].whenBind(()->bind(w));
        }
        SCC();
    }


    private void bind(int i) {
        // TODO
        int j = x[i].min();
        dest[orig[i].value()] = dest[j];
        orig[dest[j].value()] = orig[i];
        lengthToDest[orig[i].value()].setValue(
                lengthToDest[orig[i].value()].value() + lengthToDest[j].value() + 1);
        if(lengthToDest[orig[i].value()].value() < x.length){
            x[dest[j].value()].remove(orig[i].value());
        }

        SCC();
    }

    private void SCC(){
        int[] scc = GraphUtil.stronglyConnectedComponents(g);
        for(int k:scc){
            if(k!=scc[0]) throw new InconsistencyException();
        }
    }

}

