package minicp.examples;

import com.sun.net.httpserver.Authenticator;
import minicp.engine.constraints.Equal;
import minicp.engine.constraints.LessOrEqual;
import minicp.engine.core.*;

import minicp.examples.DARPParser.*;
import minicp.examples.DARPDataModel.*;
import minicp.search.DFSearch;
import minicp.search.SearchStatistics;
import minicp.state.*;
import minicp.util.Procedure;
import minicp.util.exception.InconsistencyException;

import java.util.*;

import static minicp.cp.Factory.*;
import static minicp.cp.BranchingScheme.*;

/**
 * @author Quentin Delmelle qdelmelle@gmail.com
 * transtlated from the scala work of
 * Roger Kameugne rkameugne@gmail.com
 * and Charles Thomas cftmthomas@gmail.com
 */

class DARPModelVH {
    static DARPInstance instance;
    static Solver cp;
    static DFSearch search;

    // meta-parameters
    static int alpha = 1;
    static int beta = 0;
    static int gamma = 200;
    static int tau = 1000;
    static int maxSize;
    static int range = 4;
    static int numIter = 300;
    static double d = 0.07;
    static int maxTime = 30;
    static long remainTime = maxTime * 1000L;
    static int SCALING = 100;

//    solPath: Option[String],
//    logPath: Option[String],

    // utility parameters
    static boolean firstSolOnly = false;
    static boolean debug = false;

    // solution stats
    static DARPSolution bestSolution = null;
    static int bestSolutionObjective = Integer.MAX_VALUE;
    static DARPSolution currentSolution = null;
    static int currentSolutionObjective = Integer.MAX_VALUE;
    static int totalNumFails = 0;

    // Data

    // Parameters of the problem
    static int numVars;
    static int numRequests;
    static int vehicleCapacity;
    static int numVehicles;
    static int maxRideTime;
    static int maxRouteDuration;
    static int timeHorizon = 0;

    // Variables necessary for posting constraints
    static int[][] dist;
    static int[] load;
    static int[] timeWindowStarts;
    static int[] timeWindowEnds;
    static int[] servingDuration;

    // Variables and structures related to modelling
    static IntVar[] servingTime;
    static IntVar[] servingVehicle;
    static StateInt[] succ;
    static StateInt[] pred;
    static StateInt[] capacityLeftInRoute;
    static HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>[] insertionObjChange;
    static StateSparseSet customersLeft;

    public static DARPSolution run(DARPInstance I, int runtime, int a, int b) {
        bestSolution = null;
        bestSolutionObjective = Integer.MAX_VALUE;
        currentSolution = null;
        currentSolutionObjective = Integer.MAX_VALUE;
        totalNumFails = 0;
        maxTime = runtime;
        remainTime = maxTime * 1000L;
        alpha = a;
        beta = b;
        instance = I;
        cp = makeSolver();
        System.out.println("firstSolOnly= " + firstSolOnly);

        // Data

        // Parameters of the problem
        numVars = instance.nSites();
        numRequests = instance.nRequests();
        vehicleCapacity = instance.vCapacity;
        numVehicles = instance.nVehicles;
        maxRideTime = instance.maxRideTime * SCALING;
        maxRouteDuration = instance.maxDuration * SCALING;
        timeHorizon = instance.timeHorizon * SCALING;
        maxSize = Math.max(numRequests / 2, range * 2);

        System.out.println("numVars = " + numVars);
        System.out.println("numRequests = " + numRequests);
        System.out.println("numVehicles = " + numVehicles);

        // Variables necessary for posting constraints
        DARPStop[] sites = instance.getSites();
        dist = instance.getDistances(sites, SCALING);

        load = new int[numVars];
        for (int i = 0; i < numVars; i++) load[i] = sites[i].load;
        timeWindowStarts = new int[numVars];
        for (int i = 0; i < numVars; i++) timeWindowStarts[i] = sites[i].winStart * SCALING;
        timeWindowEnds = new int[numVars];
        for (int i = 0; i < numVars; i++) timeWindowEnds[i] = sites[i].winEnd * SCALING;
        servingDuration = new int[numVars];
        for (int i = 0; i < numVars; i++) servingDuration[i] = sites[i].service * SCALING;
        insertionObjChange = new HashMap[numVars];

        // Variables and structures related to modelling
        servingTime = new IntVar[numVars];
        servingVehicle = makeIntVarArray(cp, numVars, 0, numVehicles - 1);
        succ = new StateInt[numVars];
        pred = new StateInt[numVars];
        capacityLeftInRoute = new StateInt[numVars];
        for (int i = 0; i < numVars; i++) capacityLeftInRoute[i] = cp.getStateManager().makeStateInt(vehicleCapacity);
        customersLeft = new StateSparseSet(cp.getStateManager(), numRequests, 0);

        //Start solving the problem
        // initialise variables
        initCpVars();

        // post constraints
        postConstraints();

        // search algorithm
        search = makeDfs(cp, () -> {
            printDebug("" + exportSol(0));
            if (customersLeft.isEmpty()) return EMPTY;
            else {
                int request = getUnassignedRequest();
                int[][] points = getInsertionPoints(request);
                //printInsertionPoints(points);
                printDebug("branching on request " + request + " with " + points.length + " insertion points...");
                Procedure[] branches = new Procedure[points.length];
                for (int i = 0; i < points.length; i++) {
                    int[] p = points[i];
                    branches[i] = () -> branchRequestPoint(request, p);
                }
                return branch(branches);
            }
        });

        // On solution
        search.onSolution(() -> {
            double rand = Math.random();
            int obj = getDistanceObjective();
            if (obj < currentSolutionObjective || rand < d) {
                currentSolution = exportSol(obj);
                currentSolutionObjective = obj;
                if (currentSolutionObjective < bestSolutionObjective) {
                    bestSolution = currentSolution;
                    bestSolutionObjective = currentSolutionObjective;
                    System.out.println("new best solution found:= " + bestSolutionObjective / (double) SCALING);
                    System.out.println("remainingTime := " + remainTime);
                }
            }
        });
        boolean solFound = false;

        //search.solve(); System.out.println("best obj = "+bestSolutionObjective);
        // find initial solution
        while (!solFound && remainTime > 0) {
            long t1 = System.currentTimeMillis();
            SearchStatistics stats = search.solve(statistics -> statistics.numberOfSolutions() == 1);
            // + (int)Math.round(remainTime / 1000.0)); ??
            if (stats.numberOfSolutions() == 1) solFound = true;
            else System.out.println("restart! " + remainTime);
            //System.out.println(stats);
            //solFound = true; // if we don't want restarts
            remainTime -= System.currentTimeMillis() - t1;
            totalNumFails += stats.numberOfFailures();
        }

        if (!firstSolOnly) lns();
        System.out.println("total failures: " + totalNumFails);
        return bestSolution;

        // write sol in a file
        //  val pw = new PrintWriter(new File("data/DARP/Cordeau2003/" + name + "Solution.txt"))
        //  pw.write(name + " " + bestSolutionObjective/100.0 +"\n")
        //  for(i <- numVarsRange){
        //    pw.write(bestSolution.get.succ(i) + " ")
        //    pw.write(bestSolution.get.pred(i) + " ")
        //    pw.write(bestSolution.get.minServingTime(i) + " ")
        //    pw.write(bestSolution.get.maxServingTime(i) + " ")
        //    pw.write(bestSolution.get.servingVehicle(i)+"\n")
        //  }
        //  pw.close
    }

    // lns
    static void lns() {
        int i = 2;
        while (remainTime > 0 && i <= maxSize - range) {
            int j = 0;
            int finalI = i;
            int finalJ = j;
            while (remainTime > 0 && j <= range) {
                int k = 1;
                while (remainTime > 0 && k <= numIter) {
                    long t1 = System.currentTimeMillis();
                    SearchStatistics stats = search.solveSubjectTo(
                            statistics -> statistics.numberOfSolutions() == 1,
                            () -> {
                                relax(finalI + finalJ);
                            });
                    k++;
                    remainTime -= System.currentTimeMillis() - t1;
                    totalNumFails += stats.numberOfFailures();
                }
                j++;
            }
            i++;
        }
    }

    static void relax(int nRelax) {
        Set<Integer> relaxedCustomers = selectRelaxedCustomers(nRelax);
        for (int r = 0; r < numRequests; r++) {
            if (!relaxedCustomers.contains(r)) customersLeft.remove(r);
        }
        for (int v = 0; v < numVehicles; v++) {
            int prev = getBeginDepot(v);
            for (DARPStep s : currentSolution.paths[v].steps) {
                int current = s.stop;
                int r = getCorrespondingRequest(current);
                if (current >= numRequests * 2 || !relaxedCustomers.contains(r)) {
                    succ[prev].setValue(current);
                    pred[current].setValue(prev);
                    cp.post(lessOrEqual(getArrivalTime(prev, current), servingTime[current]));
                    equal(servingVehicle[current], v);
                    prev = current;
                }
            }
            updateCapacityLeftInRoute(v, -1);
        }
    }

    static int getUnassignedRequest() {
        int[] cl = customersLeft.toArray();
        for (int r : cl) {
            insertionObjChange[r] = new HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>();
            for (int v = 0; v < numVehicles; v++) {
                setInsertionCost(r, v);
                if (!insertionObjChange[r].containsKey(v)) servingVehicle[r].remove(v);
            }
        }
        int minVehicles = numVehicles + 1;

        // step 1: minimize number of routes
        for (int r : cl) {
            minVehicles = Math.min(minVehicles, servingVehicle[r].size());
        }
        List<Integer> minRoutesRequests = new LinkedList<Integer>();
        for (int r : cl) {
            if (servingVehicle[r].size() == minVehicles) minRoutesRequests.add(r);
        }
        if (minRoutesRequests.size() == 1) return minRoutesRequests.get(0); // shortcut

        // step 2: minimize number of insertion points
        int minChoices = Integer.MAX_VALUE;
        int[] bestInsertionCost = new int[numRequests];
        int[] numInsertions = new int[numRequests];
        for (int r : minRoutesRequests) {
            int[][] iP = getInsertionPoints(r);
            if (iP.length == 0) System.out.println(minVehicles);
            bestInsertionCost[r] = iP[0][3];
            numInsertions[r] = iP.length;
            minChoices = Math.min(minChoices, iP.length);
        }
        List<Integer> minInsPointsRequests = new LinkedList<Integer>();
        for (int r : minRoutesRequests) {
            if (numInsertions[r] == minChoices) minInsPointsRequests.add(r);
        }
        if (minInsPointsRequests.size() == 1) return minInsPointsRequests.get(0); // shortcut

        // step 3: minimize cost of best insertion
        int bestChange = Integer.MAX_VALUE;
        for (int r : minInsPointsRequests) {
            bestChange = Math.min(bestChange, bestInsertionCost[r]);
        }
        List<Integer> minBestCostRequests = new LinkedList<Integer>();
        for (int r : minInsPointsRequests) {
            if (bestInsertionCost[r] == bestChange) minBestCostRequests.add(r);
        }

        Random rn = new Random();
        return minBestCostRequests.get(rn.nextInt(minBestCostRequests.size()));
    }

    /**
     * @param request
     * @return all the insertion points of request, sorted in increasing order of insertion cost.
     */
    static int[][] getInsertionPoints(int request) {
        List<int[]> vehiclePointsChangeBuffer = new ArrayList<int[]>();
        for (int v : insertionObjChange[request].keySet()) {
            for (int pickup : insertionObjChange[request].get(v).keySet()) {
                for (int drop : insertionObjChange[request].get(v).get(pickup).keySet()) {
                    vehiclePointsChangeBuffer.add(new int[]
                            {v, pickup, drop, insertionObjChange[request].get(v).get(pickup).get(drop)});
                }
            }
        }
        int[][] vehiclePointsChange = vehiclePointsChangeBuffer.toArray(new int[0][0]);
        Arrays.sort(vehiclePointsChange, new Comparator<int[]>() {
            @Override
            public int compare(int[] a, int[] b) {
                return a[3] - b[3];
            }
        });
        return vehiclePointsChange;
    }

    static void branchRequestPoint(int request, int[] point) {
        int vehicle = point[0], pPrev = point[1], dPrev = point[2], e = point[3];
        printDebug("insertion point: " + vehicle + " " + pPrev + " " + dPrev + " " + e);
        int pickup = request, drop = request + numRequests;
        equal(servingVehicle[pickup], vehicle);
        equal(servingVehicle[drop], vehicle);
        // /!\ insert pickup first
        if (!insertVertexIntoRoute(pickup, pPrev)) {
            printDebug("fail1");
            throw new InconsistencyException();
        }
        if (!insertVertexIntoRoute(drop, dPrev)) {
            printDebug("fail2");
            throw new InconsistencyException();
        }

        updateCapacityLeftInRoute(vehicle, pickup);
        customersLeft.remove(request);
        // recompute the insertions of all remaining unassigned requests for this vehicle
        // and check consistency.
        for (int r : customersLeft.toArray()) {
            if (insertionObjChange[r].containsKey(vehicle)) {
                insertionObjChange[r].remove(vehicle);
                setInsertionCost(r, vehicle);
                if (!insertionObjChange[r].containsKey(vehicle)) servingVehicle[r].remove(vehicle);
            }
        }
    }

    // insert i after j in route v.
    // return false if the insertion fails.
    static boolean insertVertexIntoRoute(int i, int j) {
        int s = succ[j].value();
        if (!tryPost(new LessOrEqual(getArrivalTime(j, i), servingTime[i])) ||
                !tryPost(new LessOrEqual(getArrivalTime(i, s), servingTime[s])))
            return false;
        succ[j].setValue(i);
        pred[s].setValue(i);
        succ[i].setValue(s);
        pred[i].setValue(j);
        return true;
    }

    // finds all possible insertions of request in route v and compute their cost.
    static void setInsertionCost(int request, int v) {
        int pickup = request, drop = request + numRequests;
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int pPred = begin;
        while (pPred != end) {
            int pSucc = succ[pPred].value();
            int pMinServingTime = Math.max(getArrivalTime(pPred, pickup).min(), servingTime[pickup].min());
            int pMaxServingTime = Math.min(servingTime[pSucc].max() - dist[pickup][pSucc] - servingDuration[pickup], servingTime[pickup].max());
            if (pMaxServingTime >= pMinServingTime && capacityLeftInRoute[pPred].value() >= load[pickup]) {
                // drop inserted just after pickup
                int dMinServingTime = Math.max(pMinServingTime + servingDuration[pickup] + dist[pickup][drop], servingTime[drop].min());
                int dMaxServingTime = Math.min(servingTime[pSucc].max() - dist[pSucc][drop] - servingDuration[drop], servingTime[drop].max());
                if (dMaxServingTime >= dMinServingTime) {
                    int slack = dMaxServingTime - servingTime[pPred].min() - dist[pPred][pickup] - servingDuration[pPred] - dist[pickup][pSucc] - servingDuration[pickup]
                            + servingTime[pSucc].max() - pMinServingTime - dist[pickup][drop] - servingDuration[pickup] - dist[drop][pSucc] - servingDuration[drop];
                    int costIncrease = dist[pPred][pickup] + dist[pickup][drop] + dist[drop][pSucc] - dist[pPred][pSucc];
                    addToInsertionObjChange(request, v, pPred, pickup, alpha * costIncrease - beta * slack);
                }

                int dPred = pSucc; // drop inserted further
                while (dPred != end && capacityLeftInRoute[dPred].value() >= load[pickup]) {
                    int dSucc = succ[dPred].value();
                    dMinServingTime = Math.max(getArrivalTime(dPred, drop).min(), servingTime[drop].min());
                    dMaxServingTime = Math.min(servingTime[dSucc].max() - dist[dSucc][drop] - servingDuration[drop], servingTime[drop].max());
                    if (dMaxServingTime >= dMinServingTime) {
                        int slack = servingTime[pSucc].max() - servingTime[pPred].min() - dist[pPred][pickup] - servingDuration[pPred] - dist[pickup][pSucc] - servingDuration[pickup]
                                + servingTime[dSucc].max() - servingTime[dPred].min() - dist[dPred][drop] - servingDuration[dPred] - dist[drop][dSucc] - servingDuration[drop];
                        int costIncrease = dist[pPred][pickup] + dist[pickup][pSucc] - dist[pPred][pSucc]
                                + dist[dPred][drop] + dist[drop][dSucc] - dist[dPred][dSucc];
                        addToInsertionObjChange(request, v, pPred, dPred, alpha * costIncrease - beta * slack);
                    }
                    dPred = dSucc;
                }
            }
            pPred = pSucc;
        }
    }

    static void addToInsertionObjChange(int request, int v, int pickup, int drop, int cost) {
        printDebug("insertion point found: " + request + ", " + v + ", " + pickup + ", " + drop + ", " + cost);
        if (!insertionObjChange[request].containsKey(v)) {
            insertionObjChange[request].put(v, new HashMap<Integer, HashMap<Integer, Integer>>());
        }
        if (!insertionObjChange[request].get(v).containsKey(pickup)) {
            insertionObjChange[request].get(v).put(pickup, new HashMap<Integer, Integer>());
        }
        insertionObjChange[request].get(v).get(pickup).put(drop, cost);
    }

    // update capacity for vehicule v starting from node <start> in route.
    // if start == -1, start from begin depot.
    static void updateCapacityLeftInRoute(int v, int start) {
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int index = start;

        if (index == -1) {
            index = succ[begin].value();
        }
        int capacity = capacityLeftInRoute[pred[index].value()].value();
        while (index != end) {
            //System.out.println(index);
            capacity -= load[index];
            capacityLeftInRoute[index].setValue(capacity);
            index = succ[index].value();
        }
    }

    static void printCapacityLeftInRoute(int v) {
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int index = succ[begin].value();
        System.out.println("capa left in route " + v + ":");
        while (index != end) {
            System.out.print(capacityLeftInRoute[index].value() + ", ");
            index = succ[index].value();
        }
        System.out.print("\n");
    }

    /**
     * @param numCustomersToRelax
     * @return a set of [numCustomersToRelax] randomly selected customers.
     */
    static Set<Integer> selectRelaxedCustomers(int numCustomersToRelax) {
        int[] customers = new int[numRequests];
        Arrays.setAll(customers, i -> i);
        Random rn = new Random();
        int relaxEnd = 0;
        while (relaxEnd < numCustomersToRelax && relaxEnd < customers.length) {
            int toRelax = relaxEnd + rn.nextInt(numRequests - relaxEnd);
            int cRelaxed = customers[toRelax];
            customers[toRelax] = customers[relaxEnd];
            customers[relaxEnd] = cRelaxed;
            relaxEnd++;
        }
        Set<Integer> ret = new HashSet<Integer>();
        for (int i = 0; i < relaxEnd; i++) ret.add(customers[i]);
        return ret;
    }

    static boolean tryPost(Constraint c) {
        try {
            cp.post(c);
        } catch (InconsistencyException e) {
            return false;
        }
        return true;
    }

    // variables && constraints initialization

    static void initCpVars() {
        for (int i = 0; i < numVars; i++) {
            if (i < 2 * numRequests) { // i is a site
                succ[i] = cp.getStateManager().makeStateInt(i);
                pred[i] = cp.getStateManager().makeStateInt(i);
            } else {
                if (isBeginDepot(i)) { // i is a start depot
                    succ[i] = cp.getStateManager().makeStateInt(getEndDepot(getVehicleOfDepot(i)));
                    pred[i] = cp.getStateManager().makeStateInt(-1);
                    servingVehicle[i].assign(getVehicleOfDepot(i));
                } else { // i is an end depot
                    succ[i] = cp.getStateManager().makeStateInt(-1);
                    pred[i] = cp.getStateManager().makeStateInt(getBeginDepot(getVehicleOfDepot(i)));
                    servingVehicle[i].assign(getVehicleOfDepot(i));
                }
            }
            servingTime[i] = makeIntVar(cp, timeWindowStarts[i], timeWindowEnds[i]);
        }
    }

    static void postConstraints() {
        // dependency
        /*for (int i = 0; i < numRequests; i++) {
            cp.post(equal(servingVehicle[i], servingVehicle[i + numRequests]));
        }*/
        // max ride time
        for (int i = 0; i < numRequests; i++) {
            cp.post(lessOrEqual(servingTime[numRequests + i], plus(servingTime[i], servingDuration[i] + maxRideTime)));
        }
        // max route duration
        for (int i = 0; i < numVehicles; i++) {
            cp.post(lessOrEqual(servingTime[getEndDepot(i)], plus(servingTime[getBeginDepot(i)], maxRouteDuration)));
        }
    }

    // helper methods
    // return the total distance traveled by all vehicles
    static int getDistanceObjective() {
        int length = 0;
        for (int v = 0; v < numVehicles; v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int i = begin;
            while (i != end) {
                length += dist[i][succ[i].value()];
                i = succ[i].value();
            }
        }
        return length;
    }

    static boolean isPositive(StateInt[] array) {
        for (int i = 0; i < array.length; i++) {
            if (array[i].value() < 0) {
                return false;
            }
        }
        return true;
    }

    static int getVehicleOfDepot(int i) {
        if (isBeginDepot(i))
            return i - 2 * numRequests;
        else return i - (2 * numRequests + numVehicles);
    }

    static int getBeginDepot(int i) {
        return 2 * numRequests + i;
    }

    static int getEndDepot(int i) {
        return numVars - numVehicles + i;
    }

    static boolean isBeginDepot(int i) {
        return (i >= 2 * numRequests && i < numVars - numVehicles);
    }

    boolean isEndDepot(int i) {
        return i >= numVars - numVehicles && i < numVars;
    }

    static IntVar getArrivalTime(int vertex, int successor) {
        return plus(servingTime[vertex], servingDuration[vertex] + dist[vertex][successor]);
    }

    static int getCorrespondingRequest(int i) {
        return (i < numRequests) ? i : i - numRequests;
    }

    static int getCorrespondingPickup(int i) {
        return i - numRequests;
    }

    static int getCorrespondingDelivery(int i) {
        return i + numRequests;
    }

    static boolean isInbound(int r) {
        return timeWindowStarts[r] > 0 || timeWindowEnds[r] < timeHorizon;
    }

    static int getCriticalVertex(int request) {
        return (isInbound(request)) ?
                request : request + numRequests;
    }

    static int getCorrespondingVertex(int i) {
        return (isPickup(i)) ?
                getCorrespondingDelivery(i) : getCorrespondingPickup(i);
    }

    static boolean isPickup(int i) {
        return i < numRequests;
    }

    boolean isDelivery(int i) {
        return i >= numRequests && i < 2 * numRequests;
    }

    boolean isCustomerVertex(int i) {
        return i < 2 * numRequests;
    }

    boolean isCriticalVertex(int vertex) {
        return timeWindowStarts[vertex] > 0 || timeWindowEnds[vertex] < timeHorizon;
    }
    /*
    // returns the minimum duration of route v.
    static int getRouteDuration(int v) {
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int p = begin;
        int i = bestSolution.succ[begin];
        int rd = 0;
        while (p != end) {
            rd += servingDuration[i] + dist[p][i];
            p = i;
            i = bestSolution.succ[i];
        }
        return rd;
    }*/

    /**
     * wraps the current state into an object representing the current solution.
     */
    static DARPSolution exportSol(double cost) {
        DARPPath[] paths = new DARPPath[numVehicles];
        for (int v = 0; v < numVehicles; v++) {
            paths[v] = new DARPPath(v);
            int current = getBeginDepot(v);
            int end = getEndDepot(v);
            while (current != end) {
                int st = (int) ((double) servingTime[current].min() / SCALING);
                int et = (int) ((double) servingTime[current].max() / SCALING);
                paths[v].addStep(new DARPStep(current, st, et));
                current = succ[current].value();
            }
            paths[v].addStep(new DARPStep(end, (int) ((double) servingTime[end].min() / SCALING), (int) ((double) servingTime[end].max() / SCALING)));
        }
        return new DARPSolution(paths, cost / SCALING, totalNumFails);
    }

    static void printsequence(int v) {
        System.out.println("Current sequence for vehicle v: ");
        for (int i = 0; i < numVars; i++) {
            if (servingVehicle[i].isBound() && servingVehicle[i].min() == v)
                System.out.println("(" + pred[i] + ", " + i + ", " + succ[i] + ")");
        }
    }

    private static void printInsertionPoints(int[][] points) {
        System.out.println("points: ");
        System.out.println("");
        for (int[] p : points) {
            System.out.print(p[0] + " " + p[1] + " " + p[2] + " " + p[3]);
            System.out.println("");
        }
    }

    // prints s only if debug is on.
    private static void printDebug(String s) {
        if (debug) System.out.println(s);
    }

}
