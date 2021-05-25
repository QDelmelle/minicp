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

import minicp.cp.Factory;
import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.Constraint;
import minicp.engine.core.IntVar;
import minicp.state.StateInt;
import minicp.state.StateManager;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * Element Constraint modeling {@code array[y] = z}
 *
 */
public class Element1D extends AbstractConstraint {

    protected final int[] t;
    protected final IntVar y;
    protected final IntVar z;
    protected final StateInt low;
    protected final StateInt up;
    protected final ArrayList<Element1D.Pair> yz;

    protected final static class Pair implements Comparable<Element1D.Pair> {
        private final int y, z;

        private Pair(int y, int z) {
            this.y = y;
            this.z = z;
        }

        @Override
        public int compareTo(Element1D.Pair p) {
            return z - p.z;
        }
    }


    /**
     * Creates an element constraint {@code array[y] = z}
     *
     * @param array the array to index
     * @param y the index variable
     * @param z the result variable
     */
    public Element1D(int[] array, IntVar y, IntVar z) {
        super(y.getSolver());
        this.t = array;
        yz = new ArrayList<Pair>();
        for (int j = 0; j < t.length; j++)
            yz.add(new Pair(j, t[j]));
        Collections.sort(yz);
        this.y = y;
        this.z = z;
        StateManager sm = getSolver().getStateManager();
        low = sm.makeStateInt(0);
        up = sm.makeStateInt(t.length - 1);
    }

    @Override
    public void post() {
        y.removeBelow(0);
        y.removeAbove(t.length-1);
        y.propagateOnDomainChange(this);
        z.propagateOnBoundChange(this);
        propagate();
    }

    public void propagate(){
        int l = low.value(), u = up.value();
        int zMin = z.min(), zMax = z.max();

        while (yz.get(l).z < zMin || !y.contains(yz.get(l).y)) {
            y.remove(yz.get(l).y);
            l++;
            if (l > u) throw new InconsistencyException();
        }
        while (yz.get(u).z > zMax || !y.contains(yz.get(u).y)) {
            y.remove(yz.get(u).y);
            u--;
            if (l > u) throw new InconsistencyException();
        }
        z.removeBelow(yz.get(l).z);
        z.removeAbove(yz.get(u).z);
        low.setValue(l);
        up.setValue(u);


    }
}
