package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.InsertionSequenceVar;
import minicp.util.exception.InconsistencyException;

/**
 * Dependency constraint on Insertion Sequence Variable:
 * ensures that either all the elements in U are in the sequence,
 * or none of them are.
 */

public class Dependency extends AbstractConstraint {
    int[] U;
    InsertionSequenceVar V;

    public Dependency(InsertionSequenceVar V, int[] U){
        super(V.getSolver());
        this.V = V;
        this.U = U;
    }

    public void post(){
        propagate();
        V.propagateOnExclude(this);
        V.propagateOnRequire(this);
        V.propagateOnInsert(this);
    }

    public void propagate(){
        if (V.isBound()) return;

        int[] status = new int[U.length];
        int[] counts = new int[3];
        for (int i=0; i<U.length; i++){
            int s = V.getStatus(U[i]);
            counts[s]++;
            status[i] = s;
        }
        if (counts[2] == U.length) return;

        if (counts[0] > 0){
            if (counts[2] > 0) throw new InconsistencyException();
            for (int i=0; i<U.length; i++){
                if (status[i] == 1){
                    V.require(U[i]);
                }
            }
        } else if (counts[2] > 0){
            for (int i=0; i<U.length; i++){
                if (status[i] == 1){
                    V.exclude(U[i]);
                }
            }
        }
        setActive(false);
    }
}
