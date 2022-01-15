package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.InsertionSequenceVar;
import minicp.engine.core.Solver;
import minicp.util.exception.InconsistencyException;

/**
 * @author Quentin Delmelle qdelmelle@gmail.com
 */

public class First extends AbstractConstraint {
    private InsertionSequenceVar V;
    private int f;

    /**
     * @param V
     * @param f ensures that f is the first element of V.
     */
    public First(InsertionSequenceVar V, int f) {
        super(V.getSolver());
        this.V = V;
        this.f = f;
    }

    public void post() {
        propagate();
    }

    public void propagate() {
        if (V.isMember(f)) {
            if (V.nextMember(V.domainSize()) != f) throw new InconsistencyException();
        } else {
            V.insert(f, -1);
        }
        for (int i = 0; i < V.domainSize(); i++) {
            V.remInsert(i, V.domainSize());
        }
    }
}
