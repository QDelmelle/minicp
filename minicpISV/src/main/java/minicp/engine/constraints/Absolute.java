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

import java.util.Arrays;

/**
 * Absolute value constraint
 */
public class Absolute extends AbstractConstraint {

    private final IntVar x;
    private final IntVar y;

    /**
     * Creates the absolute value constraint {@code y = |x|}.
     *
     * @param x the input variable such that its absolut value is equal to y
     * @param y the variable that represents the absolute value of x
     */
    public Absolute(IntVar x, IntVar y) {
        super(x.getSolver());
        this.x = x;
        this.y = y;
    }

    public void post() {
        x.propagateOnBoundChange(this);
        y.propagateOnBoundChange(this);
        y.removeBelow(0);
        propagate();
    }

    @Override
    public void propagate() {
        x.removeBelow(-y.max());
        x.removeAbove(y.max());
        int maxX = Math.max(Math.abs(x.max()), Math.abs(x.min()));
        y.removeAbove(maxX);
        if(x.max()<=0){
            y.removeBelow(-x.max());
        } else if(x.min()>=0){
            y.removeBelow(x.min());
        } else{
            int[] a= new int[x.size()];
            x.fillArray(a);
            Arrays.sort(a);
            int holemin = 0;    // repérage des bornes du trou
            int holemax = 0;
            for(int i = 0; i<x.size(); i++){
                if(a[i]>0) {
                    holemin = a[i - 1];
                    holemax = a[i];
                    break;
                }
            }
            y.removeBelow(Math.min(holemax, -holemin));
        }
        for(int i = -y.min()+1; i < y.min(); i++){
            x.remove(i);
        }
    }

}
