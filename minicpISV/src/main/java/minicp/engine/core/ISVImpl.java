package minicp.engine.core;


import minicp.state.StateStack;
import minicp.util.exception.NotImplementedException;
import minicp.cp.Factory;

public class ISVImpl implements InsertionSequenceVar {

    private Solver cp;
    private int n; // size of the domain
    private ISVDomain dom;

    private StateStack<Constraint> onExclude;
    private StateStack<Constraint> onInsert;
    private StateStack<Constraint> onRequire;

    public ISVImpl(Solver cp, int n){
        this.cp = cp; this.n = n;
        dom = new ISVDomain(cp.getStateManager(), n);
        onRequire = new StateStack<>(cp.getStateManager());
        onInsert = new StateStack<>(cp.getStateManager());
        onExclude = new StateStack<>(cp.getStateManager());
    }

    @Override
    public void propagateOnExclude(Constraint c) { onExclude.push(c);}

    @Override
    public void propagateOnInsert(Constraint c) { onInsert.push(c);}

    @Override
    public void propagateOnRequire(Constraint c) { onRequire.push(c);}

    protected void scheduleAll(StateStack<Constraint> constraints) {
        for (int i = 0; i < constraints.size(); i++)
            cp.schedule(constraints.get(i));
    }

    protected void notifyExclude() { scheduleAll(onExclude); }
    protected void notifyInsert() { scheduleAll(onInsert); }
    protected void notifyRequire() { scheduleAll(onRequire); }

    @Override
    public Solver getSolver() { return cp; }

    @Override
    public int size() {
        return n;
    }

    @Override
    public boolean isBound() { return dom.isBound(); }

    @Override
    public boolean isMember(int e)  { return dom.isMember(e); }

    @Override
    public String allMembers()  { return dom.allMembers(); }

    @Override
    public String allCurrentInserts()  { return dom.allCurrentInserts(); }

    @Override
    public String allInserts()  { return dom.allInserts(); }

    /**
     * return the current successor of e in S.
     */
    @Override
    public int nextMember(int e)  { return dom.nextMember(e); }

    /**
     * returns true iif (e, p) is in I.
     */
    @Override
    public boolean canInsert(int e, int p)  { return dom.canInsert(e, p); }

    /**
     * remove (e, p) from I.
     */
    @Override
    public void remInsert(int e, int p)  { dom.remInsert(e, p); }

    @Override
    public void insert(int e, int p) {
        dom.insert(e, p);
        notifyInsert();
    }

    @Override
    public void exclude(int e) {
        dom.exclude(e);
        notifyExclude();
    }

    @Override
    public void require(int e) {
        dom.require(e);
        notifyRequire();
    }

    @Override
    public int getStatus(int e) { return dom.getStatus(e); }

    @Override
    public boolean isEmpty() { return dom.isEmpty(); }
}
