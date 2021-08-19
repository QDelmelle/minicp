package minicp.examples;

import minicp.engine.constraints.*;
import minicp.engine.core.*;

import minicp.examples.DARPParser.*;
import minicp.examples.DARPDataModel.*;
import minicp.search.DFSearch;
import minicp.search.SearchStatistics;
import minicp.search.StopSearchException;
import minicp.state.*;
import minicp.util.Procedure;
import minicp.util.exception.InconsistencyException;

import java.util.*;

import static minicp.cp.Factory.*;
import static minicp.cp.BranchingScheme.*;

/**
 * @author Quentin Delmelle qdelmelle@gmail.com
 */

class DynamicDARP {
    static DynamicDARPInstance instance;
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
    static int maxTime;
    static long remainTime;
    static int SCALING = 100;

    // utility parameters
    static boolean firstSolOnly = true;
    static boolean debug = false;

    // solution stats
    static DARPSolution bestSolution;
    static int bestSolutionObjective;
    static DARPSolution currentSolution;
    static int currentSolutionObjective;
    static int totalNumFails;

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
    static int[] sitesVisited; // for each request, the number of sites already visited (0,1 or 2).
    static int[] incompleteRequests; // the set of all requests that are completed yet.
    static int[] currentPosition; // the current position of each vehicle = the next site it will visit.

    // Variables and structures related to modelling
    static IntVar[] servingTime;
    static IntVar[] servingVehicle;
    static InsertionSequenceVar[] route;
    static StateInt[] capacityLeftInRoute;
    static StateInt[] routeLength;
    static HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>[] insertionObjChange;
    static StateSparseSet customersLeft;

    public static DARPSolution run(DynamicDARPInstance I, int runtime) {
        bestSolution = null;
        bestSolutionObjective = Integer.MAX_VALUE;
        currentSolution = null;
        currentSolutionObjective = Integer.MAX_VALUE;
        totalNumFails = 0;
        cp = makeSolver();
        instance = I;
        maxTime = runtime;
        remainTime = maxTime * 1000L;
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
        currentPosition = new int[numVehicles];

        //Start solving the problem
        try {
            // initialise variables
            initCpVars();
            // post constraints
            postConstraints();
        } catch (InconsistencyException e) {
            System.out.println("bad instance!");
            return emptySol();
        }

        // search algorithm
        search = makeDfs(cp, () -> {
            printDebug("cL: " + customersLeft.toString());
            printRoutes();
            if (customersLeft.isEmpty()) return EMPTY;
            else {
                int request = getUnassignedRequest();
                int[][] points = getInsertionPoints(request);
                if (points.length == 0) {
                    printDebug("aborting! request "+request+" has no insertion points!");
                    throw new StopSearchException();
                }
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
        SearchStatistics stats = search.solve(statistics -> statistics.numberOfSolutions() == 1);
        if (stats.numberOfSolutions() == 0) return emptySol();

        /*
        while (!solFound && remainTime > 0) {
            long t1 = System.currentTimeMillis();
            SearchStatistics stats = search.solve(statistics -> statistics.numberOfSolutions() == 1
                    && statistics.numberOfFailures() <= Math.max(gamma * numVehicles, tau));
            if (stats.numberOfSolutions() == 1) solFound = true;
            else System.out.println("restart! " + remainTime);
            //System.out.println(stats);
            //solFound = true; // if we don't want restarts
            remainTime -= System.currentTimeMillis() - t1;
            totalNumFails += stats.numberOfFailures();
        }*/

        if (!firstSolOnly) lns();
        //printBestRoutesSolution();
        //System.out.println(maxRouteDuration);
        //for (int v = 0; v < numVehicles; v++) System.out.println(getSolutionRouteLength(v));
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
            int prev = -1;
            for (DARPStep s : currentSolution.paths[v].steps) {
                int current = s.stop;
                int r = getCorrespondingRequest(current);
                if (current >= numRequests * 2 || !relaxedCustomers.contains(r) || (sitesVisited[r] == 1 && isPickup(current))) {
                    route[v].insert(current, prev);
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
            if (iP.length == 0) {
                return r;
            }
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
        int vehicle = point[0], pPrev = point[1], dPrev = point[2], costIncrease = point[3];
        printDebug("insertion point: " + vehicle + " " + pPrev + " " + dPrev + " " + costIncrease);
        int pickup = request, drop = request + numRequests;
        // /!\ insert pickup first
        if (sitesVisited[request] == 0) {
            if (!insertVertexIntoRoute(pickup, pPrev, vehicle)) {
                throw new InconsistencyException();
            }
            equal(servingVehicle[pickup], vehicle);
        }
        if (!insertVertexIntoRoute(drop, dPrev, vehicle)) {
            printDebug("fail2");
            throw new InconsistencyException();
        }
        equal(servingVehicle[drop], vehicle);

        routeLength[vehicle].setValue(getRouteLength(vehicle));
        updateCapacityLeftInRoute(vehicle, pickup);
        customersLeft.remove(request);
        // recompute the insertions of all remaining unassigned requests for this vehicle
        // and check consistency.

        /*
        for (int r : customersLeft.toArray()) {
            if (insertionObjChange[r].containsKey(vehicle)) {
                insertionObjChange[r].remove(vehicle);
                setInsertionCost(r, vehicle);
                int[][] iP = getInsertionPoints(r);
                if (iP.length == 0) {
                    System.out.println("fail4"); // never fails here, why?
                    throw new InconsistencyException();
                }
            }
        }*/
    }

    // insert i after j in route v.
    // return false if the insertion fails.
    static boolean insertVertexIntoRoute(int i, int j, int v) {
        int s = route[v].nextMember(j);
        if (!tryPost(new LessOrEqual(getArrivalTime(j, i), servingTime[i])) ||
                !tryPost(new LessOrEqual(getArrivalTime(i, s), servingTime[s])))
            return false;
        try {
            route[v].insert(i, j);
        } catch (InconsistencyException e) {
            printDebug("insertion fail");
            return false;
        }
        return true;
    }

    /**
     * finds all possible insertions of request in route v and compute their cost.
     */
    static void setInsertionCost(int request, int v) {
        int pickup = request;
        if (sitesVisited[request] == 0) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int pPred = begin;
            while (pPred != end) {
                int pSucc = route[v].nextMember(pPred);
                if (capacityLeftInRoute[pPred].value() >= load[pickup] && route[v].canInsert(pickup, pPred)) {
                    int pMinServingTime = Math.max(getArrivalTime(pPred, pickup).min(), servingTime[pickup].min());
                    int pMaxServingTime = Math.min(servingTime[pSucc].max() - dist[pickup][pSucc] - servingDuration[pickup], servingTime[pickup].max());
                    if (pMaxServingTime < pMinServingTime) {
                        pPred = pSucc;
                        continue;
                    }
                    findDropInsert(request, v, pPred, pSucc, pPred);
                }
                pPred = pSucc;
            }
        } else if (sitesVisited[request] == 1) { // if pickup was already visited
            if (servingVehicle[pickup].isBound() && servingVehicle[pickup].min() == v)
                findDropInsert(request, v, route[v].prevMember(pickup), route[v].nextMember(pickup), currentPosition[v]);
        }
    }

    /**
     * computes all possible insertion points for request, assuming that the pickup vertex
     * has already been inserted between pPred and pSucc. Insertions for the drop vertex are
     * tried at every position from start to the end of the route.
     */
    static void findDropInsert(int request, int v, int pPred, int pSucc, int start) {
        int pickup = request, drop = request + numRequests;
        int end = getEndDepot(v);
        int dPred = start;
        if (start == pPred) { // drop inserted just after pickup
            int dMinServingTime = Math.max(getArrivalTime(pickup, drop).min(), servingTime[drop].min());
            int dMaxServingTime = Math.min(servingTime[pSucc].max() - dist[pSucc][drop] - servingDuration[drop], servingTime[drop].max());
            printDebug("dmin: " + dMinServingTime + ", dmax: " + dMaxServingTime);

            if (dMaxServingTime >= dMinServingTime && route[v].canInsert(drop, pickup)) {
                int slack = Math.min(servingTime[pickup].max() - servingTime[pickup].min(), servingTime[drop].max() - servingTime[drop].min());
                int costIncrease = dist[pPred][pickup] + dist[pickup][drop] + dist[drop][pSucc] - dist[pPred][pSucc];
                // check max route duration
                if (costIncrease + servingDuration[pickup] + servingDuration[drop] + routeLength[v].value() <= maxRouteDuration) {
                    addToInsertionObjChange(request, v, pPred, pickup, alpha * costIncrease - beta * slack);
                }
            }
            dPred = route[v].nextMember(dPred);
        }

        // drop inserted further
        boolean done = false;
        while (dPred != end && (capacityLeftInRoute[dPred].value() >= load[pickup] || sitesVisited[request] == 1) && !done) {
            int dSucc = route[v].nextMember(dPred);
            if (route[v].canInsert(drop, dPred)) {
                int dMinServingTime = Math.max(getArrivalTime(dPred, drop).min(), servingTime[drop].min());
                int dMaxServingTime = Math.min(servingTime[dSucc].max() - dist[dSucc][drop] - servingDuration[drop], servingTime[drop].max());
                // check max ride time
                //if (dMinServingTime - (pMaxServingTime + servingDuration[pickup]) > maxRideTime) done = true;
                if (dMaxServingTime >= dMinServingTime) {
                    int slack = Math.min(servingTime[pickup].max() - servingTime[pickup].min(), servingTime[drop].max() - servingTime[drop].min());
                    int costIncrease = dist[pPred][pickup] + dist[pickup][pSucc] - dist[pPred][pSucc]
                            + dist[dPred][drop] + dist[drop][dSucc] - dist[dPred][dSucc];
                    // check max route duration
                    if (costIncrease + servingDuration[pickup] + servingDuration[drop] + routeLength[v].value() <= maxRouteDuration)
                        addToInsertionObjChange(request, v, pPred, dPred, alpha * costIncrease - beta * slack);
                }
            }
            dPred = dSucc;
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

    // variables && constraints initialization

    static void initCpVars() {
        servingTime = new IntVar[numVars];
        servingVehicle = makeIntVarArray(cp, numVars, 0, numVehicles);
        for (int i = 0; i < numVars; i++) {
            if (i >= 2 * numRequests) { // i is a depot
                servingVehicle[i].assign(getVehicleOfDepot(i));
            }
            servingTime[i] = makeIntVar(cp, timeWindowStarts[i], timeWindowEnds[i]);
        }
        sitesVisited = new int[numRequests];
        for (int i = 0; i < numRequests; i++) sitesVisited[i] = 0;
        route = new InsertionSequenceVar[numVehicles];
        for (int v = 0; v < numVehicles; v++) { // rebuild routes in progress from instance
            route[v] = makeISV(cp, numVars);
            int last = -1;
            for (DARPStep s : instance.paths[v].steps) {
                route[v].insert(s.stop, last);
                last = s.stop;
                if (last < 2 * numRequests) { // if not depot
                    servingTime[s.stop].removeBelow(s.starttime * SCALING);
                    servingTime[s.stop].removeAbove(s.endtime * SCALING);
                    servingVehicle[s.stop].assign(v);
                    sitesVisited[getCorrespondingRequest(s.stop)]++;
                }
            }
            currentPosition[v] = last;
        }
        customersLeft = new StateSparseSet(cp.getStateManager(), numRequests, 0);
        ArrayList<Integer> IR = new ArrayList<Integer>();
        for (int i = 0; i < numRequests; i++) {
            if (sitesVisited[i] < 2) IR.add(i);
            else customersLeft.remove(i);
        }
        incompleteRequests = new int[IR.size()];
        for (int i = 0; i < IR.size(); i++) {
            incompleteRequests[i] = IR.get(i);
        }

        capacityLeftInRoute = new StateInt[numVars];
        for (int i = 0; i < numVars; i++)
            capacityLeftInRoute[i] = cp.getStateManager().makeStateInt(vehicleCapacity);
    }

    static void postConstraints() {
        for (int i = 0; i < numVehicles; i++) {
            cp.post(new First(route[i], getBeginDepot(i)));
            cp.post(new Last(route[i], getEndDepot(i)));
        }
        routeLength = new StateInt[numVars];
        for (int v = 0; v < numVehicles; v++) {
            routeLength[v] = cp.getStateManager().makeStateInt(getRouteLength(v));
            updateCapacityLeftInRoute(v, -1);
        }
        // dependency
        /*for (int i = 0; i < numRequests; i++)
            cp.post(equal(servingVehicle[i], servingVehicle[i + numRequests]));
        */
        // max ride time
        for (int i = 0; i < numRequests; i++) {
            cp.post(lessOrEqual(servingTime[numRequests + i], plus(servingTime[i], servingDuration[i] + maxRideTime)));
        }
        // max route duration
        /*for (int i = 0; i < numVehicles; i++) {
            cp.post(lessOrEqual(servingTime[getEndDepot(i)], plus(servingTime[getBeginDepot(i)], maxRouteDuration)));
        }*/
    }


    // auxiliary methods

    /**
     * update capacity for vehicule v starting from node <start> in route.
     * start MUST be part of route[v].
     * if start == -1, start from begin depot.
     */
    static void updateCapacityLeftInRoute(int v, int start) {
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int index = start;

        if (index == -1) {
            index = route[v].nextMember(begin);
        }
        if (!route[v].isMember(index)) return;
        int capacity = capacityLeftInRoute[route[v].prevMember(index)].value();
        while (index != end) {
            //System.out.println(index);
            capacity -= load[index];
            capacityLeftInRoute[index].setValue(capacity);
            index = route[v].nextMember(index);
        }
    }

    /**
     * @param numCustomersToRelax
     * @return a set of [numCustomersToRelax] randomly selected customers.
     */
    static Set<Integer> selectRelaxedCustomers(int numCustomersToRelax) {
        int[] customers = new int[incompleteRequests.length];
        Arrays.setAll(customers, i -> incompleteRequests[i]);
        Random rn = new Random();
        int relaxEnd = 0;
        while (relaxEnd < numCustomersToRelax && relaxEnd < customers.length) {
            int toRelax = relaxEnd + rn.nextInt(customers.length - relaxEnd);
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

    // return the total distance traveled by all vehicles
    static int getDistanceObjective() {
        int length = 0;
        for (int v = 0; v < numVehicles; v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int i = begin;
            while (i != end) {
                length += dist[i][route[v].nextMember(i)];
                i = route[v].nextMember(i);
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

    // helper methods
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
        int distance = dist[vertex][successor];
        if (isBeginDepot(vertex)) return new IntVarViewOffset(servingTime[vertex], distance);
        else return new IntVarViewOffset(servingTime[vertex], servingDuration[vertex] + distance);
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

    static int getArrivalTimeValue(int vertex, int successor, boolean getMin) {
        if (isBeginDepot(vertex)) {
            if (getMin) return servingTime[vertex].min() + dist[vertex][successor];
            else return servingTime[vertex].max() + dist[vertex][successor];
        }
        if (getMin) return servingTime[vertex].min() + servingDuration[vertex]
                + dist[vertex][successor];
        else return servingTime[vertex].max() + servingDuration[vertex]
                + dist[vertex][successor];
    }

    // returns the current duration of route v.
    static int getRouteLength(int v) {
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int p = begin;
        int i = route[v].nextMember(begin);
        int rd = 0;
        while (p != end && p != route[v].domainSize()) {
            rd += servingDuration[i] + dist[p][i];
            p = i;
            i = route[v].nextMember(i);
        }
        return rd;
    }

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
                current = route[v].nextMember(current);
            }
            paths[v].addStep(new DARPStep(end, (int) ((double) servingTime[end].min() / SCALING), (int) ((double) servingTime[end].max() / SCALING)));
        }
        return new DARPSolution(paths, cost / SCALING, totalNumFails);
    }

    static DARPSolution emptySol() {
        DARPPath[] paths = new DARPPath[numVehicles];
        for (int v = 0; v < numVehicles; v++) {
            paths[v] = new DARPPath(v);
        }
        return new DARPSolution(paths, Double.MAX_VALUE, totalNumFails);
    }

    private static void printInsertionPoints(int[][] points) {
        System.out.println("points: ");
        System.out.println("");
        for (int[] p : points) {
            System.out.print(p[0] + " " + p[1] + " " + p[2] + " " + p[3]);
            System.out.println("");
        }
    }

    private static void printRoutes() {
        for (int v = 0; v < numVehicles; v++) {
            printDebug("v: " + v + ", route: " + route[v].allMembers());
        }
    }

    static void printDistances() {
        for (int i = 0; i < numVars; i++) {
            System.out.println("");
            for (int j = 0; j <= i; j++) {
                System.out.print(dist[i][j] + " ");
            }
        }
    }

    static void printCustomersLeft() {
        System.out.println("customers left: ");
        for (int r : customersLeft.toArray()) System.out.print(r + " ");
    }

    // prints s only if debug is on.
    private static void printDebug(String s) {
        if (debug) System.out.println(s);
    }

    private static class NoSolutionException extends Exception {
    }
}
