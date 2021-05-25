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
import minicp.engine.core.BoolVar;
import minicp.engine.core.IntVar;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;

import java.util.Arrays;
import java.util.Comparator;

import static minicp.cp.Factory.*;

/**
 * Disjunctive Scheduling Constraint:
 * Any two pairs of activities cannot overlap in time.
 */
public class Disjunctive extends AbstractConstraint {

    private final IntVar[] start;
    private final int[] duration;
    private final IntVar[] end;

    

    /**
     * Creates a disjunctive constraint that enforces
     * that for any two pair i,j of activities we have
     * {@code start[i]+duration[i] <= start[j] or start[j]+duration[j] <= start[i]}.
     *
     * @param start the start times of the activities
     * @param duration the durations of the activities
     */
    public Disjunctive(IntVar[] start, int[] duration) {
        this(start, duration, true);
    }


    private Disjunctive(IntVar[] start, int[] duration, boolean postMirror) {
        super(start[0].getSolver());
        this.start = start;
        this.duration = duration;
        this.end = Factory.makeIntVarArray(start.length, i -> plus(start[i], duration[i]));

        
    }


    @Override
    public void post() {

        int[] demands = new int[start.length];
        for (int i = 0; i < start.length; i++) {
            demands[i] = 1;
        }
        getSolver().post(new Cumulative(start, duration, demands, 1), false);


        // TODO 1: replace by  posting  binary decomposition using IsLessOrEqualVar
        for(int i = 0; i<start.length; i++){
            for(int j = i+1; j<start.length; j++){
                BoolVar bij = makeBoolVar(this.getSolver());
                BoolVar bji = makeBoolVar(this.getSolver());
                getSolver().post(new IsLessOrEqualVar(bij, end[i], start[j]));
                getSolver().post(new IsLessOrEqualVar(bji, end[j], start[i])); //lol
                getSolver().post(notEqual(bij,bji));
            }
        }
        // TODO 2: add the mirror filtering as done in the Cumulative Constraint
        IntVar[] startMirror = Factory.makeIntVarArray(start.length, i -> minus(end[i]));
        IntVar[] endMirror = Factory.makeIntVarArray(end.length, i -> minus(start[i]));
        for(int i = 0; i<startMirror.length; i++){
            for(int j = i+1; j<startMirror.length; j++){
                BoolVar bij = makeBoolVar(this.getSolver());
                BoolVar bji = makeBoolVar(this.getSolver());
                getSolver().post(new IsLessOrEqualVar(bij, endMirror[i], startMirror[j]));
                getSolver().post(new IsLessOrEqualVar(bji, endMirror[j], startMirror[i]));
                getSolver().post(notEqual(bij,bji));
            }
        }
    }

    @Override
    public void propagate() {

        // HINT: for the TODO 1-4 you'll need the ThetaTree data-structure

        // TODO 3: add the OverLoadCheck algorithms

        // TODO 4: add the Detectable Precedences algorithm

        // TODO 5: add the Not-Last algorithm

        // TODO 6 (optional, for a bonus): implement the Lambda-Theta tree and implement the Edge-Finding        overLoadChecker();

        boolean fixed = false;
        while (!fixed) {
            fixed = true;
            overLoadChecker();
            fixed =  fixed && !detectablePrecedence();
            fixed =  fixed && !notLast();
        }

    }
    

    private void overLoadChecker() {
        int startmin = Integer.MAX_VALUE;
        int endmax = Integer.MIN_VALUE;
        int sumL = 0;
        for(int i = 0; i<start.length; i++){
            startmin = Math.min(startmin, start[i].min());
            endmax = Math.max(endmax, end[i].max());
            sumL += duration[i];
        }
        if(endmax-startmin < sumL){
            throw new InconsistencyException();
        }
    }

    /**
     * @return true if one domain was changed by the detectable precedence algo
     */
    private boolean detectablePrecedence() {
         throw new NotImplementedException("Disjunctive");
    }

    /**
     * @return true if one domain was changed by the not-last algo
     */
    private boolean notLast() {
         throw new NotImplementedException("Disjunctive");
    }


}
