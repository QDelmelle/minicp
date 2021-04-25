package minicp.engine.constraints;

import minicp.engine.SolverTest;
import minicp.engine.core.ISVImpl;
import minicp.engine.core.InsertionSequenceVar;
import minicp.engine.core.Solver;
import minicp.util.exception.InconsistencyException;
import org.junit.Test;

import static minicp.cp.Factory.*;
import static org.junit.Assert.*;

public class LastTest extends SolverTest {
    @Test
    public void lastTest1(){

        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            assertFalse(V.isBound());
            cp.post(new Last(V, 1));
            assertEquals("[(0, $), (0, 2), (2, 0), (2, $), ]", V.allInserts());
            assertEquals("[$ 1]", V.allMembers());
        } catch (InconsistencyException e){
            fail("should not fail here");
        }

    }

    @Test
    public void lastTest2(){

        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            V.insert(2, -1);
            cp.post(new Last(V, 2));
            assertEquals("[(0, $), (0, 1), (1, 0), (1, $), ]", V.allInserts());
            assertEquals("[$ 2]", V.allMembers());
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }

    @Test
    public void lastTest3(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            V.insert(2, -1);
            cp.post(new Last(V, 0));
            assertEquals("[(1, 2), (1, $), ]", V.allInserts());
            assertEquals("[$ 2 0]", V.allMembers());
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }

    @Test
    public void lastTest4(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            cp.post(new Last(V, 1));
            cp.post(new Last(V, 1));
            assertEquals("[(0, $), (0, 2), (2, 0), (2, $), ]", V.allInserts());
            assertEquals("[$ 1]", V.allMembers());
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }

    @Test
    public void lastTest5(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            cp.post(new Last(V, 1));
            cp.post(new Last(V, 2));
            fail("should fail here");
        } catch (InconsistencyException e){
        }
    }

    @Test
    public void lastTest6(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            V.insert(2, -1);
            V.remInsert(1, 2);
            cp.post(new Last(V, 1));
            fail("should fail here");
        } catch (InconsistencyException e){
        }
    }

    @Test
    public void lastTest7(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            V.insert(0, -1);
            V.insert(1, -1);
            cp.post(new Last(V, 1));
            fail("should fail here");
        } catch (InconsistencyException e){
        }
    }
}
