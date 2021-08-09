package minicp.engine.constraints;

import minicp.cp.Factory;
import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.InsertionSequenceVar;
import minicp.engine.core.IntVar;
import minicp.util.exception.InconsistencyException;

import java.util.List;

public class TransitionTimes extends AbstractConstraint {
    private InsertionSequenceVar V;
    private IntVar[] start;
    private int[] dur;
    private int[][] tt;
    private int n;

    /**
     * @param V              an ISV defined on a set of size n.
     * @param start          the time windows of each activity. Length must be n.
     * @param duration       size n.
     * @param transitionTime size n x n.
     *                       This constraint ensures that start[b] >= start[a] + duration[a] + trans[a][b] for all pairs a,b of consecutive elements in V.
     */
    public TransitionTimes(InsertionSequenceVar V, IntVar[] start, int[] duration, int[][] transitionTime) {
        super(V.getSolver());
        this.V = V;
        n = V.domainSize();
        this.start = start;
        dur = duration;
        tt = transitionTime;
    }

    public void post() {
        V.propagateOnInsert(this);
        for (IntVar s : start) s.propagateOnBoundChange(this);
        propagate();
    }

    public void propagate() {
        timeWindowUpdate();
        insertionUpdate();
    }

    // removes all inserts that would violate a transition time.
    private void insertionUpdate() {
        if (V.isEmpty()) return;
        for (int e = 0; e < n; e++) {
            if (!V.isMember(e)) {
                List<Integer> inserts = V.getInserts(e);
                for (int p : inserts) {
                    if (V.isMember(p)) {
                        if (p < n && start[p].min() + dur[p] + tt[p][e] > start[e].max()) V.remInsert(e, p);
                        int s = V.nextMember(p);
                        if (s < n && start[e].min() + dur[e] + tt[e][s] > start[s].max()) V.remInsert(e, p);
                    }
                }
            }
        }
    }

    // iterates the sequence and fails if a transition time is violated.
    private void timeWindowUpdate() {
        if (V.size() < 2) return;
        this.setActive(false); // the bound update goes over the whole sequence, so we don't need this to propagate itself.
        int last = V.nextMember(n);
        int current = V.nextMember(last);
        while (current != n) {
            cp.post(Factory.lessOrEqual(Factory.plus(start[last], dur[last] + tt[last][current]), start[current]));
            last = current;
            current = V.nextMember(current);
        }
        this.setActive(true);
    }

    // for any element r that is required in V but not inserted yet, checks if it still has time window feasible inserts.
    private void feasiblePathChecking() {
        //TODO
    }
}
