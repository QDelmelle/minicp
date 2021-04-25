package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.InsertionSequenceVar;
import minicp.engine.core.Solver;
import minicp.util.exception.InconsistencyException;

public class First extends AbstractConstraint {
    private InsertionSequenceVar V;
    private int f;

    public First(InsertionSequenceVar V, int f) {
        super(V.getSolver());
        this.V = V;
        this.f = f;
    }

    public void post(){
        propagate();
    }

    public void propagate(){
        if (V.isMember(f)) {
            if (V.nextMember(V.size()) != f) throw new InconsistencyException();
        } else {
            V.insert(f, -1);
        }
        for (int i=0; i<V.size(); i++){
            V.remInsert(i, V.size());
        }
    }
}
