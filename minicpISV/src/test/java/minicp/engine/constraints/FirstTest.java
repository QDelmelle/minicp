package minicp.engine.constraints;

import minicp.engine.SolverTest;
import minicp.engine.core.ISVImpl;
import minicp.engine.core.InsertionSequenceVar;
import minicp.engine.core.Solver;
import minicp.util.exception.InconsistencyException;
import org.junit.Test;

import static minicp.cp.Factory.*;
import static org.junit.Assert.*;

public class FirstTest extends SolverTest {
    @Test
    public void firstTest1(){

        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            assertFalse(V.isBound());
            cp.post(new First(V, 1));
            assertEquals("[(0, 2), (0, 1), (2, 0), (2, 1), ]", V.allInserts());
            assertEquals("[$ 1]", V.allMembers());
        } catch (InconsistencyException e){
            fail("should not fail here");
        }

    }

    @Test
    public void firstTest2(){

        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            V.insert(2, -1);
            cp.post(new First(V, 2));
            assertEquals("[(0, 2), (0, 1), (1, 0), (1, 2), ]", V.allInserts());
            assertEquals("[$ 2]", V.allMembers());
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }

    @Test
    public void firstTest3(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            V.insert(2, -1);
            cp.post(new First(V, 0));
            assertEquals("[(1, 0), (1, 2), ]", V.allInserts());
            assertEquals("[$ 0 2]", V.allMembers());
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }

    @Test
    public void firstTest4(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            cp.post(new First(V, 1));
            cp.post(new First(V, 1));
            assertEquals("[(0, 2), (0, 1), (2, 0), (2, 1), ]", V.allInserts());
            assertEquals("[$ 1]", V.allMembers());
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }

    @Test
    public void firstTest5(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            cp.post(new First(V, 1));
            cp.post(new First(V, 2));
            fail("should fail here");
        } catch (InconsistencyException e){
        }
    }

    @Test
    public void firstTest6(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            V.remInsert(1, 3);
            cp.post(new First(V, 1));
            fail("should fail here");
        } catch (InconsistencyException e){
        }
    }

    @Test
    public void firstTest7(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            V.insert(0, -1);
            V.insert(1, 0);
            cp.post(new First(V, 1));
            fail("should fail here");
        } catch (InconsistencyException e){
        }
    }
}
