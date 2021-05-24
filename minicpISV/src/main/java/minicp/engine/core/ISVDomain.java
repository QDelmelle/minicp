package minicp.engine.core;

import minicp.state.StateInt;
import minicp.state.StateManager;
import minicp.state.StateSparseSet;
import minicp.state.StateStack;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;

import java.util.ArrayList;
import java.util.List;

/**
 * an object representing the domain of the ISV.
 */
public class ISVDomain {
    private int n;
    private StateInt[] succ;
    private StateInt[] pred;
    private StateSparseSet elems;
    private StateInt[] elemPos;
    private StateInt r, p;
    private StateSparseSet[] posPreds;

    public ISVDomain(StateManager sm, int n){
        this.n = n;
        succ = new StateInt[n+1];
        pred = new StateInt[n+1];
        for (int i = 0; i <= n; i++){
            succ[i] = sm.makeStateInt(i);
            pred[i] = sm.makeStateInt(i);
        }
        elemPos = new StateInt[n];
        for (int i = 0; i < n; i++){
            elemPos[i] = sm.makeStateInt(i);
        }
        elems = new StateSparseSet(sm, n, 0);
        r = sm.makeStateInt(0);
        p = sm.makeStateInt(n);
        posPreds = new StateSparseSet[n];
        for (int i = 0; i < n; i++){
            posPreds[i] = new StateSparseSet(sm, n+1, 0);
            posPreds[i].remove(i); // can't insert i after i
        }
    }

    public boolean isBound() { return r.value()==p.value(); }

    public boolean isEmpty() { return !isMember(n); }

    public boolean isMember(int e) { return nextMember(e) != e; }

    public String allMembers() {
        StringBuilder b = new StringBuilder();
        b.append("[");
        int current = n;
        boolean ok = true;
        if (!isMember(current)) ok = false;
        while(ok){
            if (current == n) b.append("$");
            else b.append(" " + current);
            current = succ[current].value();
            if (current == n) ok = false;
        }
        b.append("]");
        return b.toString();
    }

    public String allCurrentInserts()  {
        StringBuilder b = new StringBuilder();
        b.append("[");
        for (int e = 0; e < n; e++){
            for(int p : posPreds[e].toArray()){
                if (isMember(p) || p == n) { // we can always insert after n = $
                    if (p == n) b.append("(" + e + ", $), ");
                    else b.append("(" + e + ", " + p + "), ");
                }
            }
        }
        b.append("]");
        return b.toString();
    }

    public String allInserts()  {
        StringBuilder b = new StringBuilder();
        b.append("[");
        for (int e = 0; e < n; e++){
            for(int p : posPreds[e].toArray()){
                if (p == n) b.append( "("+e+", $), ");
                else b.append( "("+e+", "+p+"), ");
            }
        }
        b.append("]");
        return b.toString();
    }

    public int nextMember(int e)  { return succ[e].value(); }

    public int prevMember(int e) { return pred[e].value(); }

    public boolean canInsert(int e, int p)  { return posPreds[e].contains(p); }

    public void remInsert(int e, int p)  { posPreds[e].remove(p); }

    public void insert(int e, int p) {
        if (p < 0) p = n;
        if (canInsert(e, p)) {
            // remove all inserts of e
            posPreds[e].removeAll();

            succ[e].setValue(succ[p].value());
            pred[succ[e].value()].setValue(e);
            pred[e].setValue(p);
            succ[p].setValue(e);
            int x = elems.get(r.value());
            if (x != e) {
                elems.exchangeFromIndex(r.value(), elemPos[e].value());
                elemPos[x].setValue(elemPos[e].value());
                elemPos[e].setValue(r.value());
            }
            r.increment();
        } else throw new InconsistencyException();
    }

    public void require(int e) { // ?
        if (getStatus(e) == 2) throw new InconsistencyException();
        if (getStatus(e) == 0) return;
        int x = elems.get(r.value());
        if (x != e) {
            elems.exchangeFromIndex(r.value(), elemPos[e].value());
            elemPos[x].setValue(elemPos[e].value());
            elemPos[e].setValue(r.value());
        }
        r.increment();
    }

    public void exclude(int e) {
        if (getStatus(e) == 0) throw new InconsistencyException();
        if (getStatus(e) == 2) return;

        int x = elems.get(p.value()-1);
        int pos = elemPos[e].value();
        elems.remove(e);
        elemPos[x].setValue(pos);
        elemPos[e].setValue(p.value()-1);
        p.decrement();

        posPreds[e].removeAll();
        for (int i=0; i<n; i++){
            remInsert(i, e);
        }
    }

    public int getStatus(int e) {
        if (elemPos[e].value() < r.value()){
            return 0;
        } else if (elemPos[e].value() >= p.value()){
            return 2;
        } else return 1;
    }

    public List<Integer> getInserts(int e) {
        List<Integer> ret = new ArrayList<Integer>();
        for (int p: posPreds[e].toArray()) {
            if (isMember(p) || p==n) ret.add(p);
        }
        return ret;
    }

    public int size() {
        return r.value();
    }
}
