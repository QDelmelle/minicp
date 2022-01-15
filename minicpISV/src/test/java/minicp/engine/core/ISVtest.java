package minicp.engine.core;
import minicp.engine.SolverTest;
import minicp.util.exception.InconsistencyException;
import org.junit.Test;

import static minicp.cp.Factory.*;
import static org.junit.Assert.*;

public class ISVtest extends SolverTest {

    @Test
    public void testISV0(){
        Solver cp = solverFactory.get();
        int n = 3;

        InsertionSequenceVar V = new ISVImpl(cp, n);

        assertFalse(V.isBound());
        assertTrue(V.isEmpty());
        assertEquals("[(0, $), (0, 1), (0, 2), (1, 0), (1, $), (1, 2), (2, 0), (2, 1), (2, $), ]"
                , V.allInserts());
        assertEquals("[]", V.allMembers());
        assertEquals("[(0, $), (1, $), (2, $), ]", V.allCurrentInserts());

        V.insert(0, -1);
        assertFalse(V.isBound());
        assertEquals("[(1, 0), (1, $), (1, 2), (2, 0), (2, 1), (2, $), ]", V.allInserts());
        assertEquals("[$ 0]", V.allMembers());
        assertEquals("[(1, 0), (1, $), (2, 0), (2, $), ]", V.allCurrentInserts());

        try{
            V.insert(0, 2);
            fail("should fail here");
        } catch (InconsistencyException e) {
        }

        V.insert( 1, -1);
        V.insert(2, 1);
        assertTrue(V.isBound());
        assertFalse(V.isEmpty());
        assertEquals("[]", V.allInserts());
        assertEquals("[$ 1 2 0]", V.allMembers());
        assertEquals("[]", V.allCurrentInserts());
    }

    @Test
    public void testISV1(){
        Solver cp = solverFactory.get();
        int n = 3;
        InsertionSequenceVar V = new ISVImpl(cp, n);
        assertEquals(0, V.size());
        for (int i=0;i<10;i++) {
            cp.getStateManager().saveState();

            V.insert(0, -1);
            assertTrue(V.isRequired(0));
            assertEquals(1, V.size());
            V.insert(1, 0);
            assertTrue(V.isRequired(1));
            assertTrue(V.isRequired(0));
            assertEquals(2, V.size());
            V.insert(2, 1);
            assertTrue(V.isRequired(2));
            assertTrue(V.isRequired(1));
            assertTrue(V.isRequired(0));
            assertEquals(3, V.size());

            cp.getStateManager().restoreState();

            assertEquals(0, V.size());
        }
    }

    @Test
    public void testISV2() {
        Solver cp = solverFactory.get();
        int n = 3;
        InsertionSequenceVar V = new ISVImpl(cp, n);
        assertFalse(V.isRequired(1));
        V.require(1);
        assertTrue(V.isRequired(1));
    }
}
