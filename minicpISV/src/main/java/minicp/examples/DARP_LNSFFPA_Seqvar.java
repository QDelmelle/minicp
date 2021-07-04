package minicp.examples;

import minicp.engine.constraints.*;
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

class DARP_LNSFFPA_Seqvar {
    static DARPInstance instance;
    static Solver cp;
    static DFSearch search;

    // meta-parameters
    static int searchTime;
    static int alpha = 1;
    static int beta = 0; // slack doesn't work for now
    static int gamma = 200;
    static int tau = 1000;
    static int maxSize;
    static int range = 4;
    static int numIter = 300;
    static double d = 0.07;
    static int maxTime = 30;
    static long remainTime = maxTime * 1000L;

    // utility parameters
    static boolean firstSolOnly = false;
    static boolean debug = false;

    // solution stats
    static DarpSol bestSolution = null;
    static int bestSolutionObjective = Integer.MAX_VALUE;
    static DarpSol currentSolution = null;
    static int currentSolutionObjective = Integer.MAX_VALUE;
    static int totalNumFails = 0;
    static boolean constraintsPosted = false;

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
    static InsertionSequenceVar[] route;
    static StateInt[] capacityLeftInRoute;
    static StateInt[] routeLength;
    static HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>[] insertionObjChange;
    static StateSparseSet customersLeft;

    public static void main(String[] args) {
        cp = makeSolver();
        String path = "data/DARP/Cordeau/b6-60.txt";
        //path = "data/DARP/Cordeau/a2-5.txt";
        //path = "data/DARP/sample.txt";
        instance = DARPParser.parseInstance(path);
        System.out.println("firstSolOnly= " + firstSolOnly);

        // Data

        // Parameters of the problem
        numVars = instance.nSites;
        numRequests = instance.nRequests;
        vehicleCapacity = instance.vehicles[0].capacity;
        numVehicles = instance.nVehicles;
        maxRideTime = instance.maxRideTime;
        maxRouteDuration = instance.maxDuration;
        for (DARPStop site : instance.sites) timeHorizon = Math.max(timeHorizon, site.winEnd);
        maxSize = Math.max(numRequests / 2, range * 2);

        System.out.println("numVars = " + numVars);
        System.out.println("numRequests = " + numRequests);
        System.out.println("numVehicles = " + numVehicles);

        // Variables necessary for posting constraints
        dist = instance.distances;
        //printDistances();
        load = new int[instance.nSites];
        for (int i = 0; i < numVars; i++) load[i] = instance.sites[i].load;
        timeWindowStarts = new int[instance.nSites];
        for (int i = 0; i < numVars; i++) timeWindowStarts[i] = instance.sites[i].winStart;
        timeWindowEnds = new int[instance.nSites];
        for (int i = 0; i < numVars; i++) timeWindowEnds[i] = instance.sites[i].winEnd;

        // Variables and structures related to modelling
        servingTime = new IntVar[numVars];
        servingVehicle = makeIntVarArray(cp, numVars, 0, numVehicles);
        route = new InsertionSequenceVar[numVehicles];
        for (int i = 0; i < numVehicles; i++) route[i] = makeISV(cp, numVars);
        capacityLeftInRoute = new StateInt[numVars];
        for (int i = 0; i < numVars; i++) capacityLeftInRoute[i] = cp.getStateManager().makeStateInt(vehicleCapacity);
        routeLength = new StateInt[numVars];
        for (int i = 0; i < numVehicles; i++) routeLength[i] = cp.getStateManager().makeStateInt(0);
        servingDuration = new int[numVars];
        for (int i = 0; i < numVars; i++) servingDuration[i] = instance.sites[i].service;
        insertionObjChange = new HashMap[numVars];
        customersLeft = new StateSparseSet(cp.getStateManager(), numRequests, 0);

        //Start solving the problem
        // initialise variables
        initCpVars();

        // post constraints
        postConstraints();

        // search algorithm
        search = makeDfs(cp, () -> {
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
            int[] servingVP = new int[numVars];
            Arrays.setAll(servingVP, i -> servingVehicle[i].min());
            int[] succP = new int[numVars];
            Arrays.setAll(succP, i -> route[servingVP[i]].nextMember(i));
            int[] predP = new int[numVars];
            Arrays.setAll(predP, i -> route[servingVP[i]].prevMember(i));
            double[] minServingTime = new double[numVars];
            Arrays.setAll(minServingTime, (i) -> servingTime[i].min() / 100);
            double[] maxServingTime = new double[numVars];
            Arrays.setAll(maxServingTime, (i) -> servingTime[i].max() / 100);
            double rand = Math.random();
            int obj = getDistanceObjective();
            if (obj < currentSolutionObjective || rand < d) {
                currentSolution = new DarpSol(succP, predP, servingVP, obj / 100.0, minServingTime, maxServingTime);
                currentSolutionObjective = obj;
                if (currentSolutionObjective < bestSolutionObjective) {
                    bestSolution = currentSolution;
                    bestSolutionObjective = currentSolutionObjective;
                    System.out.println("new best solution found:= " + bestSolutionObjective / 100.0);
                    System.out.println("remainingTime := " + remainTime);
                }
            }
        });
        boolean solFound = false;

        //search.solve(); System.out.println("best obj = "+bestSolutionObjective);
        // find initial solution
        while (!solFound && remainTime > 0) {
            long t1 = System.currentTimeMillis();
            SearchStatistics stats = search.solve(statistics -> statistics.numberOfSolutions() == 1
                    && statistics.numberOfFailures() <= Math.max(gamma * numVehicles, tau));
            // + (int)Math.round(remainTime / 1000.0)); ??
            if (stats.numberOfSolutions() == 1) solFound = true;
            else System.out.println("restart! " + remainTime);
            //System.out.println(stats);
            //solFound = true; // if we don't want restarts
            remainTime -= System.currentTimeMillis() - t1;
            totalNumFails += stats.numberOfFailures();
        }

        if (!firstSolOnly) lns();
        printBestRoutesSolution();
        System.out.println(maxRouteDuration);
        for (int v = 0; v < numVehicles; v++) System.out.println(getRouteDuration(v));
        System.out.println("valid solution: "+exportSol().isValid());
        System.out.println("total failures: " + totalNumFails);

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
        Set<Integer> relaxedCustomers = selectRelaxedCustomers(currentSolution, nRelax);
        int[] solSucc = new int[numVars];
        Arrays.setAll(solSucc, i -> currentSolution.succ[i]);
        for (int r = 0; r < numRequests; r++) {
            if (!relaxedCustomers.contains(r)) customersLeft.remove(r);
        }
        for (int v = 0; v < numVehicles; v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int current = solSucc[begin];
            int prev = begin;
            while (current != end) {
                int r = getCorrespondingRequest(current);
                if (!relaxedCustomers.contains(r)) {
                    route[v].insert(current, prev);
                    lessOrEqual(getArrivalTime(prev, current), servingTime[current]);
                    equal(servingVehicle[current], v);
                    routeLength[v].setValue(routeLength[v].value() + dist[prev][current] + servingDuration[prev]);
                    prev = current;
                }
                current = solSucc[current];
            }
            lessOrEqual(getArrivalTime(prev, end), servingTime[end]);
            routeLength[v].setValue(routeLength[v].value() + dist[prev][end] + servingDuration[prev]);
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
        int bestChange = Integer.MIN_VALUE;
        for (int r : minInsPointsRequests) {
            bestChange = Math.max(bestChange, bestInsertionCost[r]);
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

        equal(servingVehicle[pickup], vehicle);
        // /!\ insert pickup first
        if (!insertVertexIntoRoute(pickup, pPrev, vehicle)) {
            printDebug("fail1");
            throw new InconsistencyException();
        }
        if (!insertVertexIntoRoute(drop, dPrev, vehicle)) {
            printDebug("fail2");
            throw new InconsistencyException();
        }

        routeLength[vehicle].setValue(routeLength[vehicle].value() + costIncrease + servingDuration[pickup] + servingDuration[drop]);
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

    // finds all possible insertions of request in route v and compute their cost.
    static void setInsertionCost(int request, int v) {
        int pickup = request, drop = request + numRequests;
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
                // drop inserted just after pickup
                int dMinServingTime = Math.max(getArrivalTime(pickup, drop).min(), servingTime[drop].min());
                int dMaxServingTime = Math.min(servingTime[pSucc].max() - dist[pSucc][drop] - servingDuration[drop], servingTime[drop].max());
                if (dMaxServingTime >= dMinServingTime && route[v].canInsert(drop, pickup)) {
                    int slack = servingTime[pSucc].max() - servingTime[pickup].min() + servingTime[drop].max() - servingTime[pPred].min();
                    int costIncrease = dist[pPred][pickup] + dist[pickup][drop] + dist[drop][pSucc] - dist[pPred][pSucc];
                    // check max route duration
                    if (costIncrease + servingDuration[pickup] + servingDuration[drop] + routeLength[v].value() <= maxRouteDuration)
                        addToInsertionObjChange(request, v, pPred, pickup, alpha * costIncrease - beta * slack);
                }

                int dPred = pSucc; // drop inserted further
                boolean done = false;
                while (dPred != end && capacityLeftInRoute[dPred].value() >= load[pickup] && !done) {
                    int dSucc = route[v].nextMember(dPred);
                    if (route[v].canInsert(drop, dPred)) {
                        dMinServingTime = Math.max(getArrivalTime(dPred, drop).min(), servingTime[drop].min());
                        dMaxServingTime = Math.min(servingTime[dSucc].max() - dist[dSucc][drop] - servingDuration[drop], servingTime[drop].max());
                        // check max ride time
                        if (dMinServingTime - (pMaxServingTime + servingDuration[pickup]) > maxRideTime) done = true;
                        else if (dMaxServingTime >= dMinServingTime) {
                            int slack = servingTime[dSucc].max() - servingTime[dPred].min() + servingTime[pSucc].max() - servingTime[pPred].min();
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
            pPred = pSucc;
        }
    }

    static void addToInsertionObjChange(int request, int v, int pickup, int drop, int cost) {
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
        for (int i = 0; i < numVars; i++) {
            if (i >= 2 * numRequests) { // i is a depot
                servingVehicle[i].assign(getVehicleOfDepot(i));
            }
            servingTime[i] = makeIntVar(cp, timeWindowStarts[i], timeWindowEnds[i]);
        }
    }

    static void postDependency() { // pas besoin?
        for (int i = 0; i < numRequests; i++) {
            for (int s = 0; s < numVehicles; s++) {
                cp.post(new Dependency(route[s], new int[]{i, i + numRequests}));
            }
        }
    }

    static void postConstraints() {
        if (constraintsPosted) {
            return;
        }
        constraintsPosted = true;
        for (int i = 0; i < numVehicles; i++) {
            cp.post(new First(route[i], getBeginDepot(i)));
            cp.post(new Last(route[i], getEndDepot(i)));
            cp.post(lessOrEqual(servingTime[getEndDepot(i)],
                    plus(servingTime[getBeginDepot(i)], timeHorizon)));
        }
        // same pickup & drop vehicle
        for (int i = 0; i < numRequests; i++)
            cp.post(equal(servingVehicle[i], servingVehicle[i + numRequests]));
        // max ride time
        for (int i = 0; i < numRequests; i++) {
            cp.post(lessOrEqual(servingTime[numRequests + i], plus(servingTime[i], servingDuration[i] + maxRideTime)));
        }
        // max route duration
        for (int i = 0; i < numVehicles; i++) {
            cp.post(lessOrEqual(servingTime[getEndDepot(i)], plus(servingTime[getBeginDepot(i)], maxRouteDuration)));
        }
    }

    // update capacity for vehicule v starting from node <start> in route.
    // start MUST be part of route[v].
    // if start == -1, start from begin depot.
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

    // auxiliary methods

    /**
     * @param s
     * @param numCustomersToRelax
     * @return a set of [numCustomersToRelax] randomly selected customers.
     */
    static Set<Integer> selectRelaxedCustomers(DarpSol s, int numCustomersToRelax) {
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

    // returns the duration of route v in the solution.
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
    }

    static DARPSolution exportSol() {
        DARPPath[] paths = new DARPPath[numVehicles];
        for(int v = 0; v < numVehicles; v++){
            paths[v] = new DARPPath(v);
            int current = getBeginDepot(v);
            int end = getEndDepot(v);
            while (current != end) {
                paths[v].addStep(new DARPStep(current, (int)(bestSolution.minServingTime[current]*100), (int)(bestSolution.maxServingTime[current] * 100)));
                current = bestSolution.succ[current];
            }
            paths[v].addStep(new DARPStep(end, (int)(bestSolution.minServingTime[end]*100), (int)(bestSolution.maxServingTime[end] * 100)));
        }
        return new DARPSolution(instance, paths);
    }

    private static void printInsertionPoints(int[][] points) {
        System.out.println("points: ");
        System.out.println("");
        for (int[] p : points) {
            System.out.print(p[0] + " " + p[1] + " " + p[2] + " " + p[3]);
            System.out.println("");
        }
    }

    static void printBestRoutesSolution() {
        for (int v = 0; v < numVehicles; v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int i = begin;
            System.out.print("vehicle " + v + ": ");
            while (i != end && i != bestSolution.succ[i]) {
                System.out.print(i + ", ");
                i = bestSolution.succ[i];
            }
            System.out.println(end);
        }
    }

    private static void printRoutes() {
        for (int v = 0; v < numVehicles; v++) {
            System.out.println("v: " + v + ", route: " + route[v].allMembers());
        }
    }

    static void printCurrentSolution() {
        System.out.println("Darp Current Solution");
        for (int v = 0; v < numVehicles; v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int i = begin;
            System.out.print("vehicle: " + v + ": ");
            while (i != end && i != currentSolution.succ[i]) {
                System.out.print(i + ", ");
                i = currentSolution.succ[i];
            }
            System.out.println(end);
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
}

    /*

    def printRoutes(): Unit = {
        for (int v=0;v<numVehicles;v++) {
            val begin = getBeginDepot(v)
            val end = getEndDepot(v)
            var i = begin
            print("vehicle: " + v + ": ")
            while (i != end && i != succ(i).value) {
                print("(i : " + i + " v: " + servingVehicle(i) + ") -> ")
                i = succ(i).value
            }
            println("(i : " + end + " v: " + servingVehicle(end) + ")")
        }
    }
    def printRoutes(solSucc: Array[Int], servingV : Array[Int]): Unit = {
        for (int v=0;v<numVehicles;v++) {
            val begin = getBeginDepot(v)
            val end = getEndDepot(v)
            var i = begin
            print("vehicle: " + v + ": ")
            while (i != end && i != solSucc(i)) {
                print("(i : " + i + " v: " + servingV(i) + ") -> ")
                i = solSucc(i)
            }
            println("(i : " + end + " v: " + servingV(end) + ")")
        }
    }

    def exportBestSol: DARPSolution = exportSol(bestSolution.get)

}
*/
