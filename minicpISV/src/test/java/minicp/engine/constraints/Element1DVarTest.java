/*
 * mini-cp is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License  v3
 * as published by the Free Software Foundation.
 *
 * mini-cp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY.
 * See the GNU Lesser General Public License  for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with mini-cp. If not, see http://www.gnu.org/licenses/lgpl-3.0.en.html
 *
 * Copyright (c)  2018. by Laurent Michel, Pierre Schaus, Pascal Van Hentenryck
 */

package minicp.engine.constraints;

import com.github.guillaumederval.javagrading.GradeClass;
import minicp.engine.SolverTest;
import minicp.engine.core.IntVar;
import minicp.engine.core.Solver;
import minicp.search.DFSearch;
import minicp.search.SearchStatistics;
import minicp.util.exception.InconsistencyException;
import minicp.util.exception.NotImplementedException;
import minicp.util.NotImplementedExceptionAssume;
import org.junit.Test;

import static minicp.cp.BranchingScheme.firstFail;
import static minicp.cp.Factory.*;
import static minicp.cp.Factory.makeIntVarArray;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

@GradeClass(totalValue = 1, defaultCpuTimeout = 1000)
public class Element1DVarTest extends SolverTest {

    @Test
    public void element1dVarTest1() {

        try {

            Solver cp = solverFactory.get();
            IntVar y = makeIntVar(cp, -3, 10);
            IntVar z = makeIntVar(cp, 2, 40);

            IntVar[] T = new IntVar[]{makeIntVar(cp, 9, 9), makeIntVar(cp, 8, 8), makeIntVar(cp, 7, 7), makeIntVar(cp, 5, 5), makeIntVar(cp, 6, 6)};

            cp.post(new Element1DVar(T, y, z));

            assertEquals(0, y.min());
            assertEquals(4, y.max());


            assertEquals(5, z.min());
            assertEquals(9, z.max());

            z.removeAbove(7);
            cp.fixPoint();

            assertEquals(2, y.min());


            y.remove(3);
            cp.fixPoint();

            assertEquals(7, z.max());
            assertEquals(6, z.min());


        } catch (InconsistencyException e) {
            fail("should not fail");
        } catch (NotImplementedException e) {
            NotImplementedExceptionAssume.fail(e);
        }
    }

    @Test
    public void element1dVarTest2() {

        try {

            Solver cp = solverFactory.get();
            IntVar y = makeIntVar(cp, -3, 10);
            IntVar z = makeIntVar(cp, -4, 40);

            IntVar[] T = new IntVar[]{makeIntVar(cp, 1, 2),
                    makeIntVar(cp, 3, 4),
                    makeIntVar(cp, 5, 6),
                    makeIntVar(cp, 7, 8),
                    makeIntVar(cp, 9, 10)};

            cp.post(new Element1DVar(T, y, z));

            assertEquals(0, y.min());
            assertEquals(4, y.max());

            assertEquals(1, z.min());
            assertEquals(10, z.max());

            y.removeAbove(2);
            cp.fixPoint();

            assertEquals(6, z.max());

            y.assign(2);
            cp.fixPoint();

            assertEquals(5, z.min());
            assertEquals(6, z.max());


        } catch (InconsistencyException e) {
            fail("should not fail");
        } catch (NotImplementedException e) {
            NotImplementedExceptionAssume.fail(e);
        }
    }


    @Test
    public void element1dVarTest3() {

        try {

            Solver cp = solverFactory.get();
            IntVar y = makeIntVar(cp, -3, 10);
            IntVar z = makeIntVar(cp, -20, 40);

            IntVar[] T = new IntVar[]{makeIntVar(cp, 9, 9), makeIntVar(cp, 8, 8), makeIntVar(cp, 7, 7), makeIntVar(cp, 5, 5), makeIntVar(cp, 6, 6)};

            cp.post(new Element1DVar(T, y, z));

            DFSearch dfs = makeDfs(cp, firstFail(y, z));
            dfs.onSolution(() ->
                    assertEquals(T[y.min()].min(), z.min())
            );
            SearchStatistics stats = dfs.solve();

            assertEquals(5, stats.numberOfSolutions());


        } catch (InconsistencyException e) {
            fail("should not fail");
        } catch (NotImplementedException e) {
            NotImplementedExceptionAssume.fail(e);
        }
    }

    @Test
    public void miniTSPtest() {

        try {
            Solver cp = solverFactory.get();
            int n = 5;
            IntVar[] next = makeIntVarArray(cp, n, n); // next[x] = successor of node x
            IntVar[] prev = makeIntVarArray(cp, n, n); // prev[x] = predecessor of node x
            cp.post(new Circuit(next));
            // i = next[prev[i]] = prev[next[i]]
            for(int i=0; i<n; i++) {
                IntVar uselessVar = makeIntVar(cp, i, i);
                cp.post(new Element1DVar(next, prev[i], uselessVar));
                cp.post(new Element1DVar(prev, next[i], uselessVar));
            }

            next[0].assign(4);
            cp.fixPoint();
            assertEquals(prev[4].min(), 0);
            assertEquals(prev[4].max(), 0);
            /*
            for(int i=1; i<n; i++) {
                lessOrEqual(next[i], i);
            }*/


        } catch (InconsistencyException e) {
            fail("should not fail");
        } catch (NotImplementedException e) {
            NotImplementedExceptionAssume.fail(e);
        }
    }
}
