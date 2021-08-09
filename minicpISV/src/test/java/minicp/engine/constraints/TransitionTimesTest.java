package minicp.engine.constraints;

import minicp.engine.SolverTest;
import minicp.engine.core.*;
import minicp.search.DFSearch;
import minicp.search.SearchStatistics;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;
import minicp.util.NotImplementedExceptionAssume;
import org.junit.Assert;
import org.junit.Test;

import static minicp.cp.BranchingScheme.firstFail;
import static minicp.cp.Factory.*;
import static org.junit.Assert.*;

public class TransitionTimesTest extends SolverTest {

    public int[][] constDistMatrix(int n, int val) {
        int[][] ret = new int[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) ret[i][j] = 0;
                else ret[i][j] = val;
            }
        }
        return ret;
    }

    @Test
    public void tt1() {
        try {
            Solver cp = solverFactory.get();
            int n = 3;
            InsertionSequenceVar V = new ISVImpl(cp, n);
            IntVar[] start = makeIntVarArray(cp, n, 10);
            int[] d = new int[n];
            for (int i = 0; i < n; i++) d[i] = 1;
            int[][] tt = constDistMatrix(n, 1);
            cp.post(new TransitionTimes(V, start, d, tt));
            V.insert(0, -1);
            cp.fixPoint();
            assertEquals(0, start[0].min());
            assertEquals(9, start[0].max());
            V.insert(1, -1);
            cp.fixPoint();
            assertEquals(2, start[0].min());
            assertEquals(7, start[1].max());

        } catch (InconsistencyException e) {
            fail("should not fail here");
        }
    }

    @Test
    public void tt2() {
        Solver cp = solverFactory.get();
        int n = 3;
        InsertionSequenceVar V = new ISVImpl(cp, n);
        IntVar[] start = makeIntVarArray(cp, n, 10);
        int[] d = new int[n];
        for (int i = 0; i < n; i++) d[i] = 5;
        int[][] tt = constDistMatrix(n, 5);
        cp.post(new TransitionTimes(V, start, d, tt));
        try {
            V.insert(0, -1);
            cp.fixPoint();
        } catch (InconsistencyException e) {
            fail("should not fail here");
        }
        try {
            V.insert(1, -1);
            cp.fixPoint();
            fail("should fail here");
        } catch (InconsistencyException e) {}
    }

    @Test
    public void tt3() {
        Solver cp = solverFactory.get();
        int n = 3;
        InsertionSequenceVar V = new ISVImpl(cp, n);
        IntVar[] start = makeIntVarArray(cp, n, 10);
        int[] d = new int[n];
        for (int i = 0; i < n; i++) d[i] = 1;
        int[][] tt = constDistMatrix(n, 1);
        cp.post(new TransitionTimes(V, start, d, tt));
        try {
            V.insert(0, -1);
            V.insert(1, -1);
            start[0].removeAbove(7);
            cp.fixPoint();
            assertEquals(5, start[1].max());
            start[1].removeBelow(2);
            cp.fixPoint();
            assertEquals(4, start[0].min());
        } catch (InconsistencyException e) {
            fail("should not fail here");
        }
    }

    @Test
    public void tt4() {
        Solver cp = solverFactory.get();
        int n = 3;
        InsertionSequenceVar V = new ISVImpl(cp, n);
        IntVar[] start = makeIntVarArray(cp, n, 10);
        int[] d = new int[n];
        for (int i = 0; i < n; i++) d[i] = 1;
        int[][] tt = constDistMatrix(n, 1);
        cp.post(new TransitionTimes(V, start, d, tt));
        try {
            V.insert(0, -1);
            start[0].removeBelow(5);
            start[1].removeAbove(6);
            cp.fixPoint();
            assertFalse(V.canInsert(1, 0));
            assertTrue(V.canInsert(1, -1));
        } catch (InconsistencyException e) {
            fail("should not fail here");
        }
    }

    @Test
    public void tt5() {
        Solver cp = solverFactory.get();
        int n = 3;
        InsertionSequenceVar V = new ISVImpl(cp, n);
        IntVar[] start = makeIntVarArray(cp, n, 10);
        int[] d = new int[n];
        for (int i = 0; i < n; i++) d[i] = 1;
        int[][] tt = constDistMatrix(n, 1);
        cp.post(new TransitionTimes(V, start, d, tt));
        try {
            V.insert(0, -1);
            V.insert(1, -1);
            V.insert(2, -1);
            cp.fixPoint();
            assertEquals(4, start[0].min());
            assertEquals(5, start[2].max());
            start[0].removeAbove(7);
            cp.fixPoint();
            assertEquals(3, start[2].max());
        } catch (InconsistencyException e) {
            fail("should not fail here");
        }
    }
}
