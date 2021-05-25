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
import minicp.util.GraphUtil;
import minicp.util.GraphUtil.Graph;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Arc Consistent AllDifferent Constraint
 *
 * Algorithm described in
 * "A filtering algorithm for constraints of difference in CSPs" J-C. Régin, AAAI-94
 */
public class AllDifferentAC extends AbstractConstraint {

    private IntVar[] x;

    private final MaximumMatching maximumMatching;

    private final int nVar;
    private int nVal;

    // residual graph
    private ArrayList<Integer>[] in;
    private ArrayList<Integer>[] out;
    private int nNodes;
    private Graph g = new Graph() {
        @Override
        public int n() {
            return nNodes;
        }

        @Override
        public Iterable<Integer> in(int idx) {
            return in[idx];
        }

        @Override
        public Iterable<Integer> out(int idx) {
            return out[idx];
        }
    };

    private int[] match;
    private boolean[] matched;

    private int minVal;
    private int maxVal;

    public AllDifferentAC(IntVar... x) {
        super(x[0].getSolver());
        this.x = x;
        maximumMatching = new MaximumMatching(x);
        match = new int[x.length];
        this.nVar = x.length;
    }

    @Override
    public void post() {
        for (int i = 0; i < nVar; i++) {
            x[i].propagateOnDomainChange(this);
        }
        updateRange();

        matched = new boolean[nVal];
        nNodes = nVar + nVal + 1;
        in = new ArrayList[nNodes];
        out = new ArrayList[nNodes];
        for (int i = 0; i < nNodes; i++) {
            in[i] = new ArrayList<>();
            out[i] = new ArrayList<>();
        }
        propagate();
    }

    private void updateRange() {
        minVal = Integer.MAX_VALUE;
        maxVal = Integer.MIN_VALUE;
        for (int i = 0; i < nVar; i++) {
            minVal = Math.min(minVal, x[i].min());
            maxVal = Math.max(maxVal, x[i].max());
        }
        nVal = maxVal - minVal + 1;
    }

    private void updateGraph() {
        nNodes = nVar + nVal + 1;
        int sink = nNodes - 1;
        for (int i = 0; i < nNodes; i++) {
            in[i].clear();
            out[i].clear();
        }
        for (int i = 0; i<nVar; i++){ // directed var/val edges
            if(match[i] != MaximumMatching.NONE){
                in[i].add(match[i]-minVal+nVar);
                in[match[i]-minVal+nVar].add(sink);
                out[match[i]-minVal+nVar].add(i);
                out[sink].add(match[i]-minVal+nVar);

                matched[match[i]-minVal] = true;
            }
            int[] arr = new int[x[i].size()];
            x[i].fillArray(arr);
            for(int k = 0; k< arr.length; k++){
                int valIndex = arr[k] - minVal + nVar;
                if(match[i] != arr[k]){
                    out[i].add(valIndex);
                    in[valIndex].add(i);
                }
            }
        }
        for (int i = nVar; i<sink; i++){ // sink edges
            if(out[i].isEmpty()) { // if M-free
                out[i].add(sink);
                in[sink].add(i);
            } else{
                in[i].add(sink);
                out[sink].add(i);
            }
        }
    }


    @Override
    public void propagate() {
        // hint: use maximumMatching.compute(match) to update the maximum matching
        maximumMatching.compute(match);
        //       use updateRange() to update the range of values
        updateRange();
        //       use updateGraph() to update the residual graph
        updateGraph();
        //       use  GraphUtil.stronglyConnectedComponents to compute SCC's
        int[] comp = GraphUtil.stronglyConnectedComponents(g);

        for(int i = 0; i<nVal; i++){
            for(int k:in[i+nVar]){
                if(k<nVar)
                    if(comp[i+nVar]!= comp[k])
                        x[k].remove(i+minVal);
            }
        }
    }
}
