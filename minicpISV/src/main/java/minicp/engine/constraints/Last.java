package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.InsertionSequenceVar;
import minicp.util.exception.InconsistencyException;

public class Last extends AbstractConstraint {
    private InsertionSequenceVar V;
    private int l;

    public Last(InsertionSequenceVar V, int l) {
        super(V.getSolver());
        this.V = V;
        this.l = l;
    }

    public void post(){
        propagate();
    }

    public void propagate(){
        if (V.isMember(l)) {
            if (V.nextMember(l) != V.size()) throw new InconsistencyException();
        } else {
            int x = V.size();
            while(V.nextMember(x) != V.size()){
                x = V.nextMember(x);
            }
            V.insert(l, x);
        }
        for (int i=0; i<V.size(); i++){
            V.remInsert(i, l);
        }
    }
}
