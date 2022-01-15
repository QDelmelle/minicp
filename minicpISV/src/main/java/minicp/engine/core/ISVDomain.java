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
 *
 *  @author Quentin Delmelle qdelmelle@gmail.com
 */
public class ISVDomain {
    private int n;
    private StateInt[] succ;
    private StateInt[] pred;
    private int[] elems;
    private int[] elemPos;
    private StateInt r, p;
    private StateSparseSet[] posPreds;

    public ISVDomain(StateManager sm, int n) {
        this.n = n;
        succ = new StateInt[n + 1];
        pred = new StateInt[n + 1];
        for (int i = 0; i <= n; i++) {
            succ[i] = sm.makeStateInt(i);
            pred[i] = sm.makeStateInt(i);
        }
        elems = new int[n];
        for (int i = 0; i < n; i++) {
            elems[i] = i;
        }
        elemPos = new int[n];
        for (int i = 0; i < n; i++) {
            elemPos[i] = i;
        }
        r = sm.makeStateInt(0);
        p = sm.makeStateInt(n);
        posPreds = new StateSparseSet[n];
        for (int i = 0; i < n; i++) {
            posPreds[i] = new StateSparseSet(sm, n + 1, 0);
            posPreds[i].remove(i); // can't insert i after i
        }
    }

    public boolean isBound() {
        return r.value() == p.value();
    }

    public boolean isEmpty() {
        return !isMember(n);
    }

    public boolean isMember(int e) {
        return nextMember(e) != e;
    }

    public String allMembers() {
        StringBuilder b = new StringBuilder();
        b.append("[");
        int current = n;
        boolean ok = true;
        if (!isMember(current)) ok = false;
        while (ok) {
            if (current == n) b.append("$");
            else b.append(" " + current);
            current = succ[current].value();
            if (current == n) ok = false;
        }
        b.append("]");
        return b.toString();
    }

    public String allCurrentInserts() {
        StringBuilder b = new StringBuilder();
        b.append("[");
        for (int e = 0; e < n; e++) {
            for (int p : posPreds[e].toArray()) {
                if (isMember(p) || p == n) { // we can always insert after n = $
                    if (p == n) b.append("(" + e + ", $), ");
                    else b.append("(" + e + ", " + p + "), ");
                }
            }
        }
        b.append("]");
        return b.toString();
    }

    public String allInserts() {
        StringBuilder b = new StringBuilder();
        b.append("[");
        for (int e = 0; e < n; e++) {
            for (int p : posPreds[e].toArray()) {
                if (p == n) b.append("(" + e + ", $), ");
                else b.append("(" + e + ", " + p + "), ");
            }
        }
        b.append("]");
        return b.toString();
    }

    public int nextMember(int e) {
        return succ[e].value();
    }

    public int prevMember(int e) {
        return pred[e].value();
    }

    public boolean canInsert(int e, int p) {
        if (p == -1) return posPreds[e].contains(n);
        else return posPreds[e].contains(p);
    }

    public void remInsert(int e, int p) {
        posPreds[e].remove(p);
    }

    public void insert(int e, int p) {
        if (p < 0) p = n;
        if (canInsert(e, p)) {
            // remove all inserts of e
            posPreds[e].removeAll();

            succ[e].setValue(succ[p].value());
            pred[succ[e].value()].setValue(e);
            pred[e].setValue(p);
            succ[p].setValue(e);

            require(e);
        } else throw new InconsistencyException();
    }

    public void require(int e) { // ?
        if (isExcluded(e)) throw new InconsistencyException();
        if (isRequired(e)) return;

        swap(e, elems[r.value()]);
        r.increment();
        if (!isRequired(e)) {
            //System.out.println("require foire");
            //printDomain();
            throw new InconsistencyException();
        }
    }

    public void exclude(int e) {
        if (isRequired(e)) throw new InconsistencyException();
        if (isExcluded(e)) return;

        swap(e, elems[p.value() - 1]);
        p.decrement();

        posPreds[e].removeAll();
        for (int i = 0; i < n; i++) {
            remInsert(i, e);
        }
    }

    // swap the positions of a and b in elems, update elemPos accordingly.
    private void swap(int a, int b) {
        if (a == b) return;
        int posa = elemPos[a];
        int posb = elemPos[b];
        elems[posa] = b;
        elems[posb] = a;
        elemPos[a] = posb;
        elemPos[b] = posa;
    }

    // return true e is in R
    public boolean isRequired(int e) {
        return elemPos[e] < r.value();
    }

    // return true e is in E
    public boolean isExcluded(int e) {
        return elemPos[e] >= p.value();
    }

    public List<Integer> getInserts(int e) {
        List<Integer> ret = new ArrayList<Integer>();
        for (int p : posPreds[e].toArray()) {
            if (isMember(p) || p == n) ret.add(p);
        }
        return ret;
    }

    private void printDomain() {
        System.out.println("r,p = " + r.value() + ", " + p.value());
        System.out.println("elems = ");
        for (int i = 0; i < n; i++) {
            System.out.print(elems[i] + ", ");
        }
        System.out.println("");
        System.out.println("elemPos = ");
        for (int i = 0; i < n; i++) {
            System.out.print(elemPos[i] + ", ");
        }
        System.out.println("");
    }

    public int size() {
        return r.value();
    }
}
