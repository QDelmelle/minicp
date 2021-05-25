package minicp.engine.constraints;

import minicp.engine.core.AbstractConstraint;
import minicp.engine.core.IntVar;
import minicp.engine.core.Solver;
import minicp.state.StateInt;

import java.util.stream.IntStream;

public class AllDifferentFW extends AbstractConstraint {
    private IntVar[] x;
    private StateInt nUnBounds;
    private int[] unBounds;


    public AllDifferentFW(IntVar... x) {
        super(x[0].getSolver());
        this.x = x;
        int n = x.length;
        nUnBounds = getSolver().getStateManager().makeStateInt(n);
        unBounds = IntStream.range(0, n).toArray();
    }

    @Override
    public void post() {
        for (IntVar var : x)
            var.propagateOnBind(this);
        propagate();
    }

    @Override
    public void propagate() {
        int nU = nUnBounds.value();
        for(int i= nU-1; i>=0; i--){
            int idx = unBounds[i];
            IntVar y = x[idx];
            if(y.isBound()){
                // filter
                for(int unb : unBounds){
                    if (unb != idx)
                        x[unb].remove(y.max());
                }
                // update unbound variables
                unBounds[i] = unBounds[nU-1];
                unBounds[nU-1] = idx;
                nU--;
            }
        }
        nUnBounds.setValue(nU);
    }
}
