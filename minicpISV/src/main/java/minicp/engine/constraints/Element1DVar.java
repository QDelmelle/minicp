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
import minicp.util.exception.NotImplementedException;

import java.util.HashSet;

public class Element1DVar extends AbstractConstraint {

    private final IntVar[] array;
    private final IntVar y;
    private final IntVar z;

    

    public Element1DVar(IntVar[] array, IntVar y, IntVar z) {
        super(y.getSolver());
        this.array = array;
        this.y = y;
        this.z = z;
    }

    @Override
    public void post() {
        y.removeBelow(0);
        y.removeAbove(array.length-1);
        y.propagateOnDomainChange(this);
        z.propagateOnBoundChange(this);
        for(int i=0; i< array.length; i++){
            final int w = i;
            array[i].propagateOnBind(this);
        }
        propagate();
    }

    @Override
    public void propagate() {
        if(y.isBound()) {
            z.removeAbove(array[y.min()].max());
            z.removeBelow(array[y.min()].min());
            array[y.min()].removeAbove(z.max());
            array[y.min()].removeBelow(z.min());
            return;
        }
        int[] domy = new int[y.size()];
        y.fillArray(domy);
        for (int i = 0; i < domy.length; i++) {
            if(array[domy[i]].max() < z.min()) y.remove(domy[i]);
            if(array[domy[i]].min() > z.max()) y.remove(domy[i]);
        }
        y.fillArray(domy);
        int Max = Integer.MIN_VALUE, Min = Integer.MAX_VALUE;
        for (int i = 0; i < domy.length; i++) {
            Max = Math.max(Max, array[domy[i]].max());
            Min = Math.min(Min, array[domy[i]].min());
        }
        z.removeAbove(Max);
        z.removeBelow(Min);
        if(y.isBound()) {
            array[y.min()].removeAbove(z.max());
            array[y.min()].removeBelow(z.min());
        }
    }
}
