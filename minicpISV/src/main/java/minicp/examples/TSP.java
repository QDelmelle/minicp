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

package minicp.examples;

import minicp.engine.constraints.Circuit;
import minicp.engine.constraints.Element1D;
import minicp.engine.core.IntVar;
import minicp.engine.core.Solver;
import minicp.search.DFSearch;
import minicp.search.Objective;
import minicp.search.SearchStatistics;
import minicp.util.io.InputReader;

import java.util.Random;

import java.util.Random;
import java.util.stream.IntStream;

import static minicp.cp.BranchingScheme.*;
import static minicp.cp.Factory.*;

import static minicp.cp.Factory.*;

/**
 * Traveling salesman problem.
 * <a href="https://en.wikipedia.org/wiki/Travelling_salesman_problem">Wikipedia</a>.
 */
public class TSP {
    public static void main(String[] args) {


        // instance gr17 https://people.sc.fsu.edu/~jburkardt/datasets/tsp/gr17_d.txt
        InputReader reader = new InputReader("data/tsp.txt");

        int n = reader.getInt();

        int[][] distanceMatrix = reader.getMatrix(n, n);


        Solver cp = makeSolver(false);
        IntVar[] succ = makeIntVarArray(cp, n, n);
        IntVar[] distSucc = makeIntVarArray(cp, n, 1000);


        cp.post(new Circuit(succ));

        for (int i = 0; i < n; i++) {
            cp.post(new Element1D(distanceMatrix[i], succ[i], distSucc[i]));
        }

        IntVar totalDist = sum(distSucc);

        Objective obj = cp.minimize(totalDist);

       // do first fail variable heuristique
        DFSearch dfs = makeDfs(cp, () -> {
            IntVar xs=null;
            int p = -1;
            int min = Integer.MAX_VALUE;
            for(int i=0;i<n;i++) {
                if (!succ[i].isBound()) {
                    int [] tmp = new int[succ[i].size()];
                    succ[i].fillArray(tmp);
                    int[] next = minRClosest(i, tmp, distanceMatrix);
                    int v = distanceMatrix[i][next[0]]-distanceMatrix[i][next[1]];
                    if(v<min) { // replace succi.size to v for closest // doesn't work
                        xs = succ[i];
                        p = next[0];
                        min = v;//succ[i].size();
                    }
                }
            }
            if (xs == null)
                return EMPTY;
            else {
                int v = p;
                final IntVar xt = xs;
                return branch(() -> equal(xt,v),
                        () -> notEqual(xt, v));
            }
        });


        //dfs = makeDfs(cp, firstFail(succ));

        /*
        dfs.onSolution(() ->
                System.out.println(totalDist)
        );

        // take a while (optimum = 291)
        SearchStatistics stats = dfs.solve();
        System.out.println(stats);

        */
        // --- Large Neighborhood Search ---

        // Current best solution
        int[] succBest = IntStream.range(0, n).toArray();
        int[] totBest = new int[1];
        dfs.onSolution(() -> {
            // Update the current best solution
            for (int i = 0; i < n; i++) {
                succBest[i] = succ[i].min();
            }
            System.out.println("objective:" + totalDist.min());
            totBest[0]=totalDist.min();
        });

        //dfs.optimize(obj);

        int nRestarts = 10000;
        int failureLimit = 100;
        Random rand = new java.util.Random(0);
        int[] percentage = new int[2];
        percentage[0] = 5; // for simulated annealing % = 50->5
        percentage[1] = 50;

        for (int i = 0; i < nRestarts; i++) {
            int w = i;
            if (i%100==0)
                System.out.println("restart number #"+w);
            // Record the state such that the fragment constraints can be cancelled
            dfs.optimizeSubjectTo(obj,statistics -> statistics.numberOfFailures() >= failureLimit,
                    () -> {
                        // Assign the fragment 10% of the variables randomly chosen
                        for (int j = 0; j < n; j++) {
                            if (rand.nextInt(100) < ((percentage[1]-percentage[0])/((double)(nRestarts))*w+percentage[0])) { // kinda smulated annealing
                                equal(succ[j],succBest[j]);
                            }
                        }
                    });
        }
        System.out.println(totBest[0]);
    }
    private static int closest(int p, int[] dom,int[][] distanceMatrix){
        int closest = 0;
        int min = Integer.MAX_VALUE;
        for(int j=0;j<dom.length;j++){
            if(p==dom[j])
                continue;
            if(distanceMatrix[p][dom[j]]<min){
                min = distanceMatrix[p][dom[j]];
                closest = dom[j];
            }
        }
        return closest;
    }

    private static int[] minRClosest(int p, int[] dom,int[][] distanceMatrix){ // variable needs to be unbound
        int[] closests = new int[2]; // en 0 : plus proche // en 1 suivant
        closests[0] = -1;
        closests[1] = -1;
        int min = Integer.MAX_VALUE;
        int secondMin = Integer.MAX_VALUE;
        for(int j=0;j<dom.length;j++){
            if(p==dom[j])
                continue;
            if(distanceMatrix[p][dom[j]]<min){
                secondMin = min;
                min = distanceMatrix[p][dom[j]];
                closests[1] = closests[0];
                closests[0] = dom[j];
            }
            if(distanceMatrix[p][dom[j]]<secondMin) {
                secondMin = distanceMatrix[p][dom[j]];
                closests[1] = dom[j];
            }
        }
        return closests;
    }
}
