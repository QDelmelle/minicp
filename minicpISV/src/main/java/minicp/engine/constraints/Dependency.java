package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.InsertionSequenceVar;
import minicp.util.exception.InconsistencyException;

/**
 * Dependency constraint on Insertion Sequence Variable:
 * ensures that either all the elements in U are in the sequence,
 * or none of them are.
 *
 * @author Quentin Delmelle qdelmelle@gmail.com
 */

public class Dependency extends AbstractConstraint {
    int[] U;
    InsertionSequenceVar V;

    public Dependency(InsertionSequenceVar V, int[] U) {
        super(V.getSolver());
        this.V = V;
        this.U = U;
    }

    public void post() {
        propagate();
        V.propagateOnExclude(this);
        //V.propagateOnRequire(this);
        V.propagateOnInsert(this);
    }

    public void propagate() {
        if (V.isBound()) return;

        //TODO
    }
}
