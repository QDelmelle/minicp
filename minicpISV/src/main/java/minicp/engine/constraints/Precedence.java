package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.InsertionSequenceVar;
import minicp.engine.core.Solver;
import minicp.util.exception.InconsistencyException;

public class Precedence extends AbstractConstraint {
    int[] O;
    int[] reverseO;
    InsertionSequenceVar V;

    public Precedence(InsertionSequenceVar V, int[] O) {
        super(V.getSolver());
        this.V = V;
        this.O = O;
        reverseO = new int[O.length];
        for (int i=0;i<O.length;i++) reverseO[i] = O[O.length - i -1];
    }

    public void post(){
        propagate();
        V.propagateOnInsert(this);
    }

    public void propagate() {
        intersectionConsistency();
        precFiltering(false);
        precFiltering(true);
    }

    private void intersectionConsistency() {
        int begin = V.domainSize();
        int e = begin;
        int i = 0;
        while (V.nextMember(e) != begin) {
            e = V.nextMember(e);
            while (i < O.length && (e == O[i] || !V.isMember(O[i]))) i++;
        }
        for (;i<O.length;i++) {
            if (V.isMember(O[i])) throw new InconsistencyException();
        }
    }

    private void precFiltering(boolean reversed) {
        int begin = V.domainSize();
        int prec = -1;
        int[] tab = reversed ? reverseO : O;
        for (int e: tab) {
            if (V.isMember(e)) {
                prec = e;
            } else {
                if (prec != -1) {
                    if (reversed) {
                        int p = prec;
                        while (p != begin) {
                            V.remInsert(e, p);
                            p = V.nextMember(p);
                        }
                    } else {
                        int p = prec;
                        do {
                            p = V.prevMember(p);
                            V.remInsert(e, p);
                        } while (p != begin);
                    }
                }
            }
        }
    }
}
