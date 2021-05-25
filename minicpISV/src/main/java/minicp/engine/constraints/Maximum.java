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

/**
 * Maximum Constraint
 */
public class Maximum extends AbstractConstraint {

    private final IntVar[] x;
    private final IntVar y;

    /**
     * Creates the maximum constraint y = maximum(x[0],x[1],...,x[n])?
     *
     * @param x the variable on which the maximum is to be found
     * @param y the variable that is equal to the maximum on x
     */
    public Maximum(IntVar[] x, IntVar y) {
        super(x[0].getSolver());
        assert (x.length > 0);
        this.x = x;
        this.y = y;
    }


    @Override
    public void post() {
        for(int i=0;i<x.length;i++) {
            x[i].propagateOnBoundChange(this);
        }
        y.propagateOnBoundChange(this);
        propagate();
    }


    @Override
    public void propagate() {
        int maxmax = Integer.MIN_VALUE;
        int maxmin = Integer.MIN_VALUE;
        int nvalidx = 0;
        int lastvalidx = -1;
        for(int i = 0;i < x.length;i++) {
            x[i].removeAbove(y.max());
            maxmax = Math.max(maxmax, x[i].max());
            maxmin = Math.max(maxmin, x[i].min());
            if(x[i].max()>=y.min()){
                nvalidx++;
                lastvalidx = i;
            }
        }
        if(nvalidx == 1){   // if only one xi can provide a solution, its domain should be the same as y
            x[lastvalidx].removeBelow(y.min());
        }
        y.removeAbove(maxmax);
        y.removeBelow(maxmin);
    }

}
