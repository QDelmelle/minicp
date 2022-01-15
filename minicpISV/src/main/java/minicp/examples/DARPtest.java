package minicp.examples;

import minicp.examples.DARPDataModel.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Quentin Delmelle qdelmelle@gmail.com
 */

public class DARPtest {

    public static void main(String[] args) {
        int runtime = 120;
        int k = 15;
        DARPInstance instance = DARPParser.parseInstance("data/DARP/Cordeau/a2-16.txt");
        //DARPInstance instance = DARPParser.parseInstance("data/DARP/sample2.txt");
        //optimalityTest(instance);
        //slackTest(instance, 0, 10000, 100);
        //perfTest(instance, runtime, k);
        DARPSolution sol = DARPModelVH.run(instance, runtime, 40, 1, false);
    }

    // warning: very slow
    private static void optimalityTest(DARPInstance instance) {
        System.out.println("Running optimality test for instance " + instance.name);
        long t = System.currentTimeMillis();
        //DARPSolution sol = DARPModelVH.run(instance, true);
        long t1 = (System.currentTimeMillis() - t) / 1000;
        //System.out.println(isSolValid(instance, sol)+", "+sol);
        System.out.println("brute force done in " + t1 + " seconds.");
        t = System.currentTimeMillis();
        DARPSolution sol2 = DARPModelVH.run(instance, 900);
        long t2 = (System.currentTimeMillis() - t) / 1000;
        //System.out.println(isSolValid(instance, sol2)+", "+sol2);
        System.out.println("no ISV done in " + t2 + " seconds.");
        t = System.currentTimeMillis();
        DARPSolution sol3 = DARP_LNSFFPA_Seqvar.run(instance, 900);
        long t3 = (System.currentTimeMillis() - t) / 1000;
        //System.out.println(isSolValid(instance, sol3)+", "+sol3);
        System.out.println("ISV done in " + t3 + " seconds.");
        System.out.println("all done.");
    }

    private static void perfTest(DARPInstance instance, int runtime, int k) {
        System.out.println("Running perf test for instance " + instance.name + ", mean of " + k + " runs.");
        double meanCostSeqvar = 0;
        double failsISV = 0;
        double meanCostNoSeqvar = 0;
        double failsNoISV = 0;
        for (int i = 0; i < k; i++) {
            System.out.println("run " + (i + 1) + ", seqvar:");
            DARPSolution sol = DARP_LNSFFPA_Seqvar.run(instance, runtime);
            if (!isSolValid(instance, sol)) {
                System.out.println("illegal solution!");
                System.out.println(sol);
                return;
            }
            meanCostSeqvar += sol.cost;
            failsISV += sol.fails;
            System.out.println("run " + (i + 1) + ", no seqvar:");
            sol = DARPModelVH.run(instance, runtime);
            if (!isSolValid(instance, sol)) {
                System.out.println("illegal solution!");
                System.out.println(sol);
                return;
            }
            meanCostNoSeqvar += sol.cost;
            failsNoISV += sol.fails;
        }
        meanCostSeqvar /= k;
        meanCostNoSeqvar /= k;
        failsISV /= k;
        failsNoISV /= k;
        System.out.println("Perf test done for instance " + instance.name + ", mean of " + k + " runs.");
        System.out.println("mean cost with seqvar =" + meanCostSeqvar + ", fails: " + failsISV);
        System.out.println("mean cost without seqvar =" + meanCostNoSeqvar + ", fails: " + failsNoISV);
    }

    /**
     * @param instance
     * @param start
     * @param end
     * @param step     solves for the first solution of instance for every value of aplha/beta between start and end, with step.
     *                 writes the plot data in a file.
     *                 /!\ path should be replaced
     */
    private static void slackTest(DARPInstance instance, int start, int end, int step) {
        String path = "D:/Documents/UNIF2020-21/THESE/PAPER/slackTest/" + instance.name; //+"_"+start+"_"+end+"_"+step;
        writeToFile("", path, false);
        System.out.println("Running slack test for instance " + instance.name + ", alpha = " + start + "-->" + end + ", step = " + step);
        int beta = 1;
        double bestObj = Integer.MAX_VALUE;
        int best = start;
        for (int alpha = start; alpha <= end; alpha += step) {
            System.out.println("alpha=" + alpha + ", beta=" + beta);
            DARPSolution sol;
            sol = DARP_LNSFFPA_Seqvar.run(instance, 30, alpha, beta, true);
            if (!isSolValid(instance, sol)) {
                System.out.println("illegal solution!");
                System.out.println(sol);
            }
            writeToFile(alpha + " " + sol.cost + "\n", path, true);
            if (sol.cost <= bestObj) {
                bestObj = sol.cost;
                best = alpha;
            }
        }
        DARPSolution sol;
        sol = DARP_LNSFFPA_Seqvar.run(instance, 30, 1, 0, true);
        if (!isSolValid(instance, sol)) {
            System.out.println("illegal solution!");
            System.out.println(sol);
        }
        if (sol.cost <= bestObj) {
            bestObj = sol.cost;
            best = Integer.MAX_VALUE;
        }
        System.out.println("Slack test done for instance " + instance.name + ", alpha = " + start + "-->" + end + ", step = " + step);
        System.out.println("best value for alpha = " + best);
    }

    /**
     * @return true iif sol is a valid solution to the instance and satisfies all the constraints.
     */
    public static boolean isSolValid(DARPInstance instance, DARPSolution sol) {
        int n = instance.nSites();
        DARPStop[] sites = instance.getSites();
        int[] g = new int[n];
        int[] vehicle = new int[n];
        for (int j = 0; j < n; j++) {
            vehicle[j] = -1;
        }
        DARPStep[] stops = new DARPStep[n];
        int i = 0;
        for (int v = 0; v < instance.nVehicles; v++) {
            int lastStop = -1;
            int lastStart = -1;
            int lastEnd = -1;
            int load = 0;
            for (DARPStep step : sol.paths[v].steps) {
                stops[step.stop] = step;
                if (step.starttime < lastStart || step.endtime < lastEnd) {
                    System.out.println("Serving time clash at stops "
                            + lastStop + "[" + lastStart + ", " + lastEnd + "]/" + step.stop + "[" + step.starttime + ", " + step.endtime + "]");
                    return false;
                }
                lastStart = step.starttime;
                lastEnd = step.endtime;
                lastStop = step.stop;
                if (step.starttime < sites[step.stop].winStart || step.endtime > sites[step.stop].winEnd) {
                    System.out.println("time window violated for site " + step.stop
                            + ", " + step.starttime + " < " + sites[step.stop].winStart + ", " + step.endtime + " > " + sites[step.stop].winEnd);
                    return false;
                }
                load += sites[step.stop].load;
                if (load < 0) {
                    System.out.println("negative load in vehicle " + v + " at stop " + step.stop);
                    return false;
                } else if (load > instance.vCapacity) {
                    System.out.println("overload in vehicle " + v + " at stop " + step.stop);
                    return false;
                }
                g[step.stop] = i;
                vehicle[step.stop] = v;
                i++;
            }
            if (sol.paths[v].steps.get(sol.paths[v].len - 1).starttime - sol.paths[v].steps.get(0).endtime > instance.maxDuration) {
                System.out.println("max route duration violated for vehicle " + v);
                return false;
            }
        }
        if (i != n) {
            List<Integer> missingSites = new ArrayList<>();
            for (int j = 0; j < n; j++) {
                if (vehicle[j] == -1) missingSites.add(j);
            }
            System.out.println("some sites are not visited: " + missingSites);
            return false;
        }
        for (int r = 0; r < instance.nRequests(); r++) {
            if (vehicle[r] != vehicle[r + instance.nRequests()]) {
                System.out.println("request " + r + " is served by 2 vehicles!");
                return false;
            }
            if (g[r] > g[r + instance.nRequests()]) {
                System.out.println("request " + r + " is not visited in the right order!");
                return false;
            }
            if (stops[r + instance.nRequests()].starttime - stops[r].endtime - sites[r].service > instance.maxRideTime) {
                System.out.println("max ride time violated for request " + r);
                return false;
            }
        }
        return true;
    }

    public static void writeToFile(String s, String path, boolean append) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(path, append));
            writer.append(s);
            writer.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    // WORK IN PROGRESS
    private static void dynamicDARPTest() {
        int runtime = 5;
        int nRuns = 10;
        DynamicDARPInstance instance = new DynamicDARPInstance("empty", 3, 3, 60, 240, 0);
        instance.timeHorizon = 60;
        System.out.println("Solving the DARP...");
        //instance = DARPParser.parseInstance("data/DARP/Cordeau/a3-24.txt");
        DARPSolution sol = DynamicDARP.run(instance, runtime);
        System.out.println("Solution: \n" + sol); // display solution ?

        for (int i = 0; i < nRuns; i++) {
            System.out.println("Making new instance...");
            DynamicDARPInstance instance2 = makeDynamicInstance(instance, sol);
            System.out.println("Solving the DARP...");
            DARPSolution sol2 = DynamicDARP.run(instance2, runtime);
            if (sol2.isEmpty()) {
                System.out.println("no Solution found!");
                //i--;
            } else {
                sol = sol2;
                instance = instance2;
                if (!isSolValid(instance, sol)) {
                    System.out.println("invalid solution!");
                    System.out.println("Solution: \n" + sol); // display solution ?
                    return;
                }
            }
            System.out.println("Solution: \n" + sol); // display solution ?
        }
    }

    /**
     * @param instance
     * @param sol
     * @return modify instance to simulate vehicles moving along their predicted routes, and add some random requests.
     */
    private static DynamicDARPInstance makeDynamicInstance(DynamicDARPInstance instance, DARPSolution sol) {
        DynamicDARPInstance ret = instance.copy();
        Random rn = new Random();
        int a = instance.nRequests();
        int timePassed = randomBetween(30, 60);
        ret.startTime += timePassed;
        ret.timeHorizon = instance.timeHorizon + timePassed;
        System.out.println("time horizon: " + ret.timeHorizon);
        // add request(s)
        ret.addRequest(makeRandomRequest(ret, 5, 3));
        int b = ret.nRequests();
        // compute new vehicle positions
        for (int v = 0; v < instance.nVehicles; v++) {
            for (DARPStep s : sol.paths[v].steps) {
                if (s.endtime <= instance.startTime) ret.addStep(v, renumber(s, a, b));
                else {
                    ret.addStep(v, renumber(s, a, b));
                    break;
                }
            }
            System.out.println("Vehicle " + v + " moved " + ret.paths[v].len + " steps.");
        }
        return ret;
    }

    /**
     * @param s
     * @param a the number of requests in the last instance
     * @param b the number of requests in the new instance
     * @return a new renumbered DARPStep adapted to the new instance.
     */
    private static DARPStep renumber(DARPStep s, int a, int b) {
        int s2 = 0;
        if (s.stop < a) s2 = s.stop;
        else if (s.stop >= a && s.stop < 2 * a) s2 = s.stop + b - a;
        else s2 = s.stop + 2 * (b - a);
        return new DARPStep(s2, s.starttime, s.endtime);
    }

    /**
     * @param instance
     * @param size     random coordinates will be generated between -size and size.
     * @return a random request similar to the ones in instance.
     */
    private static DARPRequest makeRandomRequest(DARPInstance instance, double size, int service) {
        int ofs = (int) size * 2;
        int start = instance.startTime;
        int end = instance.timeHorizon - ofs;
        boolean inbound = Math.random() < 0.5;
        int startP = inbound ? randomBetween(start, end) : start;
        int endP = inbound ? randomBetween(startP, end) : end;
        int startD = inbound ? start : randomBetween(Math.min(startP + ofs, end - 1), end);
        int endD = inbound ? end : randomBetween(startD, end);
        DARPStop pickup = new DARPStop(Math.random() * 2 * size - size, Math.random() * 2 * size - size, service, 1, startP, endP);
        DARPStop drop = new DARPStop(Math.random() * 2 * size - size, Math.random() * 2 * size - size, service, -1, startD, endD);
        //pickup =
        DARPRequest ret = new DARPRequest(pickup, drop);
        System.out.println("new request: " + ret);
        return ret;
    }

    static int randomBetween(int a, int b) {
        if (b <= a) {
            System.out.println("wrong bounds!");
            return 0;
        }
        Random rn = new Random();
        return rn.nextInt(b - a) + a;
    }
}
