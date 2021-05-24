package minicp.engine.constraints;

import minicp.engine.SolverTest;
import minicp.engine.core.ISVImpl;
import minicp.engine.core.InsertionSequenceVar;
import minicp.engine.core.Solver;
import minicp.util.exception.InconsistencyException;
import org.junit.Test;

import static org.junit.Assert.*;

public class PrecedenceTest extends SolverTest {
    @Test
    public void precedenceTest0(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 5);
            cp.post(new Precedence(V, new int[]{2,3,4}));
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }
    @Test
    public void precedenceTest1(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 5);
            V.insert(0, -1); V.insert(1, 0); V.insert(2, 1);
            V.insert(3, 2); V.insert(4, 3);
            cp.post(new Precedence(V, new int[]{2,3,4}));
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }

    @Test
    public void precedenceTest2(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 5);
            V.insert(0, -1); V.insert(1, 0); V.insert(2, 1);
            V.insert(3, 2); V.insert(4, 3);
            cp.post(new Precedence(V, new int[]{2,4,3}));
            fail("should fail here");
        } catch (InconsistencyException e){
        }
    }
    @Test
    public void precedenceTest3(){
        try {
            Solver cp = solverFactory.get();
            InsertionSequenceVar V = new ISVImpl(cp, 3);
            cp.post(new Precedence(V, new int[]{0,1}));
            assertTrue(V.canInsert(1, 3));
            V.insert(0, -1);
            cp.fixPoint();
            assertFalse(V.canInsert(1, 3));
        } catch (InconsistencyException e){
            fail("should not fail here");
        }
    }
}
