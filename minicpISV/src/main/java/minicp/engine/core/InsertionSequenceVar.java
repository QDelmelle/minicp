package minicp.engine.core;

import minicp.util.Procedure;
import minicp.util.exception.InconsistencyException;

import java.util.List;

/**
 * A variable representing a sequence of integers.
 * Its domain is represented by <S, I, P, E, R>
 *     S is the current sequence.
 *     I is the set of possible insertions (e, p) == (insert e after p in S).
 *     P is the set of optional members.
 *     E is the set of excluded members.
 *     R is the set of mandatory members.
 */

public interface InsertionSequenceVar {

    /**
     * Returns the solver in which this variable was created.
     *
     * @return the solver in which this variable was created
     */
    Solver getSolver();

    /**
     * return the max size of S.
     */
    int domainSize();

    /**
     * return the number of members in S.
     */
    int size();

    /**
     * returns true iif the variable is bound.
     */
    boolean isBound();

    /**
     * returns true iif the element e is part of S.
     */
    boolean isMember(int e);

    /**
     * return a string enumerating all the elements in S.
     */
    String allMembers();

    /**
     * return a string enumerating all inserts (e, p) in I where p is in S.
     */
    String allCurrentInserts();

    /**
     * return a string enumerating all inserts (e, p) in I.
     */
    String allInserts();

    /**
     * @param e an element not yet in S.
     * @return a list containing all possible insertions of e in the current sequence.
     */
    List<Integer> getInserts(int e);

    /**
     * return the current successor of e in S. if e == n, e is considered
     * as the start symbol $
     */
    int nextMember(int e);

    /**
     * return the current predecessor of e in S. if e == n, e is considered
     * as the start symbol $
     */
    int prevMember(int e);

    /**
     * returns true iif (e, p) is in I.
     */
    boolean canInsert(int e, int p);

    /**
     * remove (e, p) from I.
     */
    void remInsert(int e, int p);

    /**
     * insert e after p in S. If p < 0 or p == n, insert e after $ (as first elem).
     */
    void insert(int e, int p);

    /**
     * exclude e from the domain.
     */
    void exclude(int e);

    /**
     * if e is not in the sequence yet, put e in R.
     */
    void require(int e);

    /**
     * return 0 if e is required (R), 2 if e is excluded (E) and 1 otherwise (P)
     */
    int getStatus(int e);

    /**
     * return true iif the sequence is empty.
     */
    boolean isEmpty();

    /**
     * Asks that {@link Constraint#propagate()} is called whenever an element is excluded.
     * We say that an <i>exclude</i> event occurs.
     *
     * @param c the constraint for which the {@link Constraint#propagate()}
     *          method should be called on exclude events of this variable.
     */
    void propagateOnExclude(Constraint c);

    /**
     * Asks that {@link Constraint#propagate()} is called whenever the insertion set I is modified.
     * We say that a <i>IChange</i> event occurs.
     *
     * @param c the constraint for which the {@link Constraint#propagate()}
     *          method should be called on IChange events of this variable.
     */
    void propagateOnIChange(Constraint c);

    /**
     * Asks that {@link Constraint#propagate()} is called whenever an element is inserted.
     * We say that an <i>insert</i> event occurs.
     *
     * @param c the constraint for which the {@link Constraint#propagate()}
     *          method should be called on insert events of this variable.
     */
    void propagateOnInsert(Constraint c);

    /**
     * creates a new empty domain for this variable
     */
    void resetDomain();
}
