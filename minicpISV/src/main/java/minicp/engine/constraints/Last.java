package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.InsertionSequenceVar;
import minicp.util.exception.InconsistencyException;

/**
 * @author Quentin Delmelle qdelmelle@gmail.com
 */

public class Last extends AbstractConstraint {
    private InsertionSequenceVar V;
    private int l;

    /**
     * @param V
     * @param l ensures that l is the las element of V.
     */
    public Last(InsertionSequenceVar V, int l) {
        super(V.getSolver());
        this.V = V;
        this.l = l;
    }

    public void post() {
        propagate();
    }

    public void propagate() {
        if (V.isMember(l)) {
            if (V.nextMember(l) != V.domainSize()) throw new InconsistencyException();
        } else {
            V.insert(l, V.prevMember(V.domainSize()));
        }
        for (int i = 0; i < V.domainSize(); i++) {
            V.remInsert(i, l);
        }
    }
}
