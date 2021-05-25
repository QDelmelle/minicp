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
import minicp.engine.constraints.Profile.Rectangle;
import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.IntVar;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;

import java.util.ArrayList;

import static minicp.cp.Factory.*;

/**
 * Cumulative constraint with time-table filtering
 */
public class Cumulative extends AbstractConstraint {

    private final IntVar[] start;
    private final IntVar[] end;
    private final int[] demand;
    private final int capa;
    private final boolean postMirror;


    /**
     * Creates a cumulative constraint with a time-table filtering.
     * At any time-point t, the sum of the demands
     * of the activities overlapping t do not overlap the capacity.
     *
     * @param start the start time of each activity
     * @param duration the duration of each activity (non negative, constant)
     * @param demand the demand of each activity, non negative
     * @param capa the capacity of the constraint
     */
    public Cumulative(IntVar[] start, int[] duration, int[] demand, int capa) {
        this(start, duration, demand, capa, true);
    }

    private Cumulative(IntVar[] start, int[] duration, int[] demand, int capa, boolean postMirror) {
        super(start[0].getSolver());
        this.start = start;
        this.end = Factory.makeIntVarArray(start.length, i -> plus(start[i], duration[i]));
        this.demand = demand;
        this.capa = capa;
        this.postMirror = postMirror;
    }

    /**
     * Creates a cumulative constraint with a time-table filtering.
     * At any time-point t, the sum of the demands
     * of the activities overlapping t do not overlap the capacity.
     *
     * @param start the start time of each activity
     * @param end the end time of each activity (for variable durations)
     * @param demand the demand of each activity, non negative
     * @param capa the capacity of the constraint
     */
    public Cumulative(IntVar[] start, IntVar[] end, int[] demand, int capa) {
        this(start, end, demand, capa, true);
    }

    private Cumulative(IntVar[] start, IntVar[] end, int[] demand, int capa, boolean postMirror) {
        super(start[0].getSolver());
        this.start = start;
        this.end = end;
        this.demand = demand;
        this.capa = capa;
        this.postMirror = postMirror;
    }


    @Override
    public void post() {
        for (int i = 0; i < start.length; i++) {
            start[i].propagateOnBoundChange(this);
            lessOrEqual(plus(start[i], 1), end[i]); // durations must be > 0
        }

        if (postMirror) {
            IntVar[] startMirror = Factory.makeIntVarArray(start.length, i -> minus(end[i]));
            IntVar[] endMirror = Factory.makeIntVarArray(start.length, i -> minus(start[i]));
            getSolver().post(new Cumulative(startMirror, endMirror, demand, capa, false), false);
        }

        propagate();
    }

    @Override
    public void propagate() {
        Profile profile = buildProfile();
        // TODO 2: check that the profile is not exceeding the capa otherwise throw an INCONSISTENCY
        for (int i = 0; i < profile.size(); i++) {
            if(profile.get(i).height()>capa)
                throw new InconsistencyException();
        }

        for (int i = 0; i < start.length; i++) {
            if (!start[i].isBound()) {
                // j is the index of the profile rectangle overlapping t
                //-------int j = profile.rectangleIndex(start[i].min());
                // TODO 3: push i to the right
                // hint:
                // Check that at every-point on the interval
                // [start[i].getMin() ... start[i].getMin()+duration[i]-1] there is enough space.
                // You may have to look-ahead on the next profile rectangle(s)
                // Be careful that the activity you are currently pushing may have contributed to the profile.
                //throw new NotImplementedException("todo 3");
                for(int time=start[i].min();time<end[i].min();time++) {
                    if(start[i].max()<end[i].min() && time >= start[i].max() && time <= end[i].min()) // if(dans mandatory)
                        continue;
                    int j = profile.rectangleIndex(time);
                    int actualC = profile.get(j).height();
                    if(demand[i]+actualC>capa)
                        start[i].removeBelow(time+1);
                }
            }
        }
    }

    public Profile buildProfile() {
        ArrayList<Rectangle> mandatoryParts = new ArrayList<Rectangle>();
        for (int i = 0; i < start.length; i++) {
            // TODO 1: add mandatory part of activity i if any
            if(start[i].max()<end[i].min())
                mandatoryParts.add(new Rectangle(start[i].max(),end[i].min(),demand[i]));
        }
        return new Profile(mandatoryParts.toArray(new Profile.Rectangle[0]));
    }

}
