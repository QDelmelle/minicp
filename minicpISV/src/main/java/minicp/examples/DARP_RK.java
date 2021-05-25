package minicp.examples;

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
    static int searchTime;
    static int alpha = 80;
    static int beta = 1;
    static int gamma = 200;
    static int tau = 1000;
    static int maxSize;
    static int range = 4;
    static int numIter = 300;
    static double d = 0.07;
    static int maxTime = 300;
    static long remainTime = maxTime * 1000L;

//    solPath: Option[String],
//    logPath: Option[String],

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
    static int[] vertexLoadChange;
    static int[] timeWindowStarts;
    static int[] timeWindowEnds;

    // Variables and structures related to modelling
    static IntVar[] servingTime;
    static IntVar[] servingVehicle;
    static StateInt[] succ;
    static StateInt[] pred;
    static StateInt[] capacityLeftInRoute;
    static int[] servingDuration;
    static HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>[] insertionObjChange;
    static StateBool[] customersLeft;

    public static void main(String[] args) {
        cp = makeSolver();
        String path = "C:\\Users\\Utilisateur\\Documents\\UNIF2020\\TFE\\DARP RK\\DARP\\Cordeau\\a3-24.txt";
        //path = "C:\\Users\\Utilisateur\\Documents\\UNIF2020\\TFE\\DARP RK\\DARP\\sample.txt";
        instance = DARPParser.parseInstance(path);
        System.out.println("firstSolOnly= " + firstSolOnly);

        // Data

        // Parameters of the problem
        numVars = instance.nSites;
        numRequests = instance.nRequests;
        vehicleCapacity = instance.vehicles[0].capacity;
        numVehicles = instance.nVehicles;
        maxRideTime = instance.requests[0].maxRideTime;
        maxRouteDuration = instance.vehicles[0].maxDuration;
        for (DARPStop site : instance.sites) timeHorizon = Math.max(timeHorizon, site.winEnd);
        maxSize = numRequests / 2;

        System.out.println("numVars = " + numVars);
        System.out.println("numRequests = " + numRequests);
        System.out.println("numVehicles = " + numVehicles);

        // Variables necessary for posting constraints
        dist = instance.distances;
        vertexLoadChange = new int[instance.nSites];
        for (int i = 0; i < numVars; i++) vertexLoadChange[i] = instance.sites[i].load;
        timeWindowStarts = new int[instance.nSites];
        for (int i = 0; i < numVars; i++) timeWindowStarts[i] = instance.sites[i].winStart;
        timeWindowEnds = new int[instance.nSites];
        for (int i = 0; i < numVars; i++) timeWindowEnds[i] = instance.sites[i].winEnd;

        // Variables and structures related to modelling
        servingTime = new IntVar[numVars];
        servingVehicle = makeIntVarArray(cp, numVars, 0, numVehicles);
        succ = new TrailInt[numVars];
        pred = new TrailInt[numVars];
        capacityLeftInRoute = new TrailInt[numVars];
        for (int i = 0; i < numVars; i++)
            capacityLeftInRoute[i] =
                    new TrailInt((Trailer) cp.getStateManager(), vehicleCapacity);
        servingDuration = new int[numVars];
        for (int i = 0; i < numVars; i++) servingDuration[i] = instance.sites[i].service;
        insertionObjChange = new HashMap[numVars];
        customersLeft = new StateBool[numRequests];
        for (int i = 0; i < numRequests; i++) customersLeft[i] = new TrailBool((Trailer) cp.getStateManager(), false);

        //Start solving the problem
        // initialise variables
        initCpVars();
        //Add request to set of unassigned request "customersLeft"
        for (int i = 0; i < numRequests; i++) {
            customersLeft[i].setValue(true);
        }

        // post constraints
        postConstraints();

        // search algorithm
        search = makeDfs(cp, () -> {
            boolean CL = false;
            for (StateBool C : customersLeft) if (C.value()) CL = true;
            if (!CL) {
                return EMPTY;
            } else {
                int request = getUnassignedRequest();
                int[][] points = getInsertionPoints(request);
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
            //System.out.println("Total route cost: " + getDistanceObjective()/100.0);
            int[] succP = new int[numVars];
            Arrays.setAll(succP, (i) -> succ[i].value());
            int[] predP = new int[numVars];
            Arrays.setAll(predP, (i) -> pred[i].value());
            int[] servingVP = new int[numVars];
            Arrays.setAll(servingVP, (i) -> servingVehicle[i].min());
            double[] minServingTime = new double[numVars];
            Arrays.setAll(minServingTime, (i) -> servingTime[i].min() / 100);
            double[] maxServingTime = new double[numVars];
            Arrays.setAll(maxServingTime, (i) -> servingTime[i].max() / 100);
//    println("maxRideTime: " + maxRideTime/100.0)
//    println("serving times:")
//    println(numVarsRange.map(i => i + " (" + minServintgTime(i) + ":" + maxServintgTime(i) + ")").mkString("\n"))
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
            //System.out.println("-------------------------");
        });
        boolean solFound = false;

        //search.solve(); System.out.println("best obj = "+bestSolutionObjective);
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
        }

        if (!firstSolOnly) lns();
        printBestRoutesSolution();

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
            while (remainTime > 0 && j <= range) {
                int k = 1;
                while (remainTime > 0 && k <= numIter) {
                    int finalI = i;
                    int finalJ = j;
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
        int[] solServingVehicle = new int[numVars];
        Arrays.setAll(solServingVehicle, i -> currentSolution.servingVehicle[i]);
        clearCustomerLeft();
        for (int r : relaxedCustomers) {
            succ[r].setValue(r);
            pred[r].setValue(r);
            succ[r + numRequests].setValue(r + numRequests);
            pred[r + numRequests].setValue(r + numRequests);
            customersLeft[r].setValue(true);
        }
        for (int v = 0; v < numVehicles; v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int current = begin;
            int prev = -1;
            while (current != end) {
                int currentRequest = getCorrespondingRequest(current);
                if (!relaxedCustomers.contains(currentRequest)) {
                    if (prev != -1) {
                        succ[prev].setValue(current);
                        pred[current].setValue(prev);
                        lessOrEqual(getArrivalTime(prev, current), servingTime[current]);
                        equal(servingVehicle[current], v);
                    }
                    prev = current;
                }
                current = solSucc[current];
            }
            succ[prev].setValue(end);
            pred[end].setValue(prev);
            lessOrEqual(getArrivalTime(prev, end), servingTime[end]);
            updateCapacityLeftInRoute(v, -1);
        }
    }

    // original version of algo 2
    static int getUnassignedMinVehicleMinInsertionPointsRequest2() {
        int minChoices = Integer.MAX_VALUE;
        int minVehicles = numVehicles + 1;
        List<BranchingChoice> bQueueBuffer = new ArrayList<BranchingChoice>();
        int bestChange = Integer.MIN_VALUE;
        for (int i = 0; i < numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            Mut numChoices = new Mut(0);
            List<BranchingChoice> branchingQueue = getInsertionCost(i, numChoices);
            int tempChange = branchingQueue.get(0).change;
            if (servingVehicle[i].size() < minVehicles) {
                minVehicles = servingVehicle[i].size();
                bQueueBuffer.clear();
                bQueueBuffer.addAll(branchingQueue);
            } else if (servingVehicle[i].size() == minVehicles && numChoices.value < minChoices) {
                minChoices = numChoices.value;
                bQueueBuffer.clear();
                bQueueBuffer.addAll(branchingQueue);
            } else if (servingVehicle[i].size() == minVehicles && numChoices.value == minChoices && tempChange > bestChange) {
                bestChange = tempChange;
                bQueueBuffer.clear();
                bQueueBuffer.addAll(branchingQueue);
            } else if (servingVehicle[i].size() == minVehicles && numChoices.value == minChoices && tempChange == bestChange) {
                bQueueBuffer.addAll(branchingQueue);
            }
        }
        Random rn = new Random();
        return bQueueBuffer.get(rn.nextInt(bQueueBuffer.size())).request;
    }

    //
    static List<BranchingChoice> getInsertionCost(int request, Mut numChoices) {
        ArrayList<BranchingChoice> branchingQueue = new ArrayList<BranchingChoice>();
        numChoices.value = 0;
        int bestChange = Integer.MIN_VALUE;

        for (int v = 0; v < numVehicles; v++) {
            if (insertionObjChange[request].containsKey(v) && servingVehicle[request].contains(v)) {
                for (int cvi : insertionObjChange[request].get(v).keySet()) {
                    for (int ncvi : insertionObjChange[request].get(v).get(cvi).keySet()) {
                        numChoices.value++;
                        if (insertionObjChange[request].get(v).get(cvi).get(ncvi) > bestChange) {
                            bestChange = insertionObjChange[request].get(v).get(cvi).get(ncvi);
                            branchingQueue.clear();
                            branchingQueue.add(new BranchingChoice(request, cvi, ncvi, bestChange,
                                    servingVehicle[cvi].min()));
                        } else if (insertionObjChange[request].get(v).get(cvi).get(ncvi) == bestChange) {
                            branchingQueue.add(new BranchingChoice(request, cvi, ncvi, bestChange,
                                    servingVehicle[cvi].min()));
                        }
                    }
                }
            }
        }
        return branchingQueue;
    }

    // my version of algo 2, works just as good
    static int getUnassignedMinVehicleMinInsertionPointsRequest() {
        int minVehicles = numVehicles + 1;
        List<Integer> unassignedRequests = new LinkedList<Integer>();

        // step 1: minimize number of routes
        for (int i = 0; i < numRequests; i++) {
            if (customersLeft[i].value()) {
                unassignedRequests.add(i);
                minVehicles = Math.min(minVehicles, servingVehicle[i].size());
            }
        }
        List<Integer> minRoutesRequests = new LinkedList<Integer>();
        for (int i : unassignedRequests) {
            if (servingVehicle[i].size() == minVehicles) minRoutesRequests.add(i);
        }

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

    static int getUnassignedRequest() {
        for (int i = 0; i < numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            insertionObjChange[i] = new HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>();
            for (int v = 0; v < numVehicles; v++) {
                setInsertionCost(i, v);
            }
        }
        return getUnassignedMinVehicleMinInsertionPointsRequest();
    }

    /**
     * @param request
     * @return all the insertion points of request, sorted by increasing order of insertion cost.
     */
    static int[][] getInsertionPoints(int request) {
        List<int[]> vehiclePointsChangeBuffer = new ArrayList<int[]>();
        for (int v : insertionObjChange[request].keySet()) {
            for (int cv : insertionObjChange[request].get(v).keySet()) {
                for (int ncv : insertionObjChange[request].get(v).get(cv).keySet()) {
                    vehiclePointsChangeBuffer.add(new int[]
                            {v, cv, ncv, insertionObjChange[request].get(v).get(cv).get(ncv)});
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
        int vehicle = point[0], cvSucc = point[1], ncvSucc = point[2], change = point[3];
        printDebug("branching point: " + request + " " + cvSucc + " " + ncvSucc + " " + change);
        int pickup = request;
        int cv = getCriticalVertex(request), ncv = getCorrespondingVertex(cv);

        equal(servingVehicle[pickup], vehicle);
        // /!\ insert cv first
        if (!insertVertexIntoRoute(cv, cvSucc)) {
            printDebug("fail1");
            throw new InconsistencyException();
        }
        if (!insertVertexIntoRoute(ncv, ncvSucc)) {
            printDebug("fail2");
            throw new InconsistencyException();
        }
        updateCapacityLeftInRoute(vehicle, pickup);
        if (!isPositive(capacityLeftInRoute)) {
            printDebug("fail3");
            throw new InconsistencyException();
        }
        customersLeft[request].setValue(false);

        // recompute the insertions of all remaining requests for this vehicle
        // and check consistency.
        for (int i = 0; i < numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            if (insertionObjChange[i].containsKey(vehicle)) {
                insertionObjChange[i].remove(vehicle);
            }
            setInsertionCost(i, vehicle);
            int[][] insertionPoint = getInsertionPoints(i);
            if (insertionPoint.length == 0) {
                printDebug("fail4");
                throw new InconsistencyException();
            }
        }
    }

    // prints s only if debug is on.
    private static void printDebug(String s) {
        if (debug) System.out.println(s);
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
            capacity -= vertexLoadChange[index];
            capacityLeftInRoute[index].setValue(capacity);
            index = succ[index].value();
        }
    }

    static boolean insertVertexIntoRoute(int i, int j) {
        //System.out.println(i+" "+j);
        if (!tryPost(new LessOrEqual(getArrivalTime(i, j), servingTime[j])) ||
                !tryPost(new LessOrEqual(getArrivalTime(pred[j].value(), i), servingTime[i])))
            return false;
        succ[i].setValue(j);
        pred[i].setValue(pred[j].value());
        succ[pred[i].value()].setValue(i);
        pred[j].setValue(i);
        return true;
    }

    static void setInsertionCost(int request, int v) {
        int cvv = getCriticalVertex(request);
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int start = succ[begin].value();
        boolean done = false;
        while (!done) {
            if (start == end) done = true;
            if (getArrivalTimeValue(pred[start].value(), cvv, true) <= servingTime[cvv].max()
                    && getArrivalTimeValue(cvv, start, true) <= servingTime[start].max()) {
                //if (!isInbound(request) ||
                //      capacityLeftInRoute[pred[start].value()].value() >= vertexLoadChange[cvv])
                setBestServingTimeFail2(request, v, start);
            }
            start = succ[start].value();
        }
    }

    // finds all possible insertions of request in route v, with critical vertex before cvSucc
    // and compute their cost.
    static void setBestServingTimeFail(int request, int v, int cvSucc) {
        int begin = getBeginDepot(v);
        int cv = getCriticalVertex(request);
        int ncv = getCorrespondingVertex(cv);
        int cvPred = pred[cvSucc].value();
        // simulate insertion of cv between cvPred and cvSucc
        succ[cv].setValue(cvSucc);
        pred[cv].setValue(cvPred);
        succ[cvPred].setValue(cv);
        pred[cvSucc].setValue(cv);

        // we want to insert ncv before s
        if (isInbound(request)) {
            int s = cvSucc;
            int end = getEndDepot(v);
            boolean done = false;
            while (!done) {
                // simulate insertion of ncv
                succ[ncv].setValue(s);
                pred[ncv].setValue(pred[s].value());
                succ[pred[ncv].value()].setValue(ncv);
                pred[s].setValue(ncv);
                done = checkInsertionInBound(request, v, cv, ncv, s);
                // /!\ undo insertion of ncv
                succ[pred[ncv].value()].setValue(succ[ncv].value());
                pred[succ[ncv].value()].setValue(pred[ncv].value());
                succ[ncv].setValue(ncv);
                pred[ncv].setValue(ncv);
                //
                if (s == end) done = true;
                s = succ[s].value();
            }
        } else {
            int s = cv;
            int start = getBeginDepot(v);
            boolean done = false;
            while (s != start && !done) {
                // simulate insertion of ncv
                succ[ncv].setValue(s);
                pred[ncv].setValue(pred[s].value());
                succ[pred[ncv].value()].setValue(ncv);
                pred[s].setValue(ncv);
                done = checkInsertionOutBound(request, v, cv, ncv, s);
                // /!\ undo insertion of ncv
                succ[pred[ncv].value()].setValue(succ[ncv].value());
                pred[succ[ncv].value()].setValue(pred[ncv].value());
                succ[ncv].setValue(ncv);
                pred[ncv].setValue(ncv);
                //
                s = pred[s].value();
            }
        }
        // /!\ undo insertion of cv
        succ[pred[cv].value()].setValue(succ[cv].value());
        pred[succ[cv].value()].setValue(pred[cv].value());
        succ[cv].setValue(cv);
        pred[cv].setValue(cv);
    }

    // checks if insertion is feasible. if it is, compute its cost. the request must be inbound.
    // return true if the search must be stopped, false otherwise.
    private static boolean checkInsertionInBound(int r, int v, int cv, int ncv, int ncvSucc) {
        int cvPred = pred[cv].value();
        int ncvPred = pred[ncv].value();
        int cvSucc = succ[cv].value();
        int cvMinServingTime = Math.max(servingTime[cv].min(),
                getArrivalTimeValue(cvPred, cv, true));
        int ncvMaxServingTime = Math.min(servingTime[ncv].max(),
                servingTime[ncvSucc].max() - dist[ncv][ncvSucc] - servingDuration[ncv]);
        int cvMaxServingTime, ncvMinServingTime;
        if (cvSucc == ncv) { // edge case: we insert drop right after pickup
            cvMaxServingTime = Math.min(servingTime[cv].max(),
                    ncvMaxServingTime - servingDuration[cv] - dist[cv][ncv]);
            ncvMinServingTime = Math.max(servingTime[ncv].min(),
                    cvMinServingTime + servingDuration[cv] + dist[cv][ncv]);
        } else {
            cvMaxServingTime = Math.min(servingTime[cv].max(),
                    servingTime[cvSucc].max() - servingDuration[cv] - dist[cv][cvSucc]);
            ncvMinServingTime = Math.max(servingTime[ncv].min(),
                    getArrivalTimeValue(ncvPred, ncv, true));
        }
        // check feasibility
        if (ncvMaxServingTime < cvMinServingTime) return false;
        // check max ride time
        int minRideTime = Math.max(dist[cv][ncv],
                ncvMinServingTime - cvMaxServingTime - servingDuration[cv]);
        if (minRideTime > maxRideTime) return true;
        // check capacity --> no need, already done in setInsertionCost
        // compute insertion cost e
        int costIncrease = dist[cvPred][cv] + dist[cv][cvSucc] - dist[cvPred][cvSucc]
                + dist[ncvPred][ncv] + dist[ncv][ncvSucc] - dist[ncvPred][ncvSucc];
        int slackAfterInsert;
        if (cvSucc == ncv) {
            slackAfterInsert = servingTime[ncvSucc].max() - cvMinServingTime;
            //-servingDuration[ncv]-dist[ncv][ncvSucc]-servingDuration[ncvPred]-dist[ncvPred][ncv];
            slackAfterInsert += ncvMaxServingTime - servingTime[cvPred].min();
            //-servingDuration[cv]-dist[cv][cvSucc]-servingDuration[cvPred]-dist[cvPred][cv];
        } else {
            slackAfterInsert = servingTime[ncvSucc].max() - servingTime[ncvPred].min();
            //-servingDuration[ncv]-dist[ncv][ncvSucc]-servingDuration[ncvPred]-dist[ncvPred][ncv];
            slackAfterInsert += servingTime[cvSucc].max() - servingTime[cvPred].min();
            //-servingDuration[cv]-dist[cv][cvSucc]-servingDuration[cvPred]-dist[cvPred][cv];
        }
        if (cvSucc == ncv) cvSucc = ncvSucc; // we want to be able to insert cv first
        addToInsertionObjChange(r, v, cvSucc, ncvSucc, alpha * costIncrease - beta * slackAfterInsert);
        return false;
    }

    // checks if insertion is feasible. if it is, compute its cost. the request must be outbound.
    // return true if the search must be stopped, false otherwise.
    private static boolean checkInsertionOutBound(int r, int v, int cv, int ncv, int ncvSucc) {
        int cvPred = pred[cv].value();
        int ncvPred = pred[ncv].value();
        int cvSucc = succ[cv].value();
        int cvMaxServingTime = Math.min(servingTime[cv].max(),
                servingTime[cvSucc].max() - servingDuration[cv] - dist[cv][cvSucc]);
        int ncvMinServingTime = Math.max(servingTime[ncv].min(),
                getArrivalTimeValue(ncvPred, ncv, true));
        int ncvMaxServingTime, cvMinServingTime;
        if (ncvSucc == cv) { // edge case: we insert drop right after pickup
            cvMinServingTime = Math.max(servingTime[cv].min(),
                    ncvMinServingTime + servingDuration[ncv] + dist[ncv][cv]);
            ncvMaxServingTime = Math.min(servingTime[ncv].max(),
                    cvMaxServingTime - dist[ncv][cv] - servingDuration[ncv]);
        } else {
            cvMinServingTime = Math.max(servingTime[cv].min(),
                    getArrivalTimeValue(cvPred, cv, true));
            ncvMaxServingTime = Math.min(servingTime[ncv].max(),
                    servingTime[ncvSucc].max() - dist[ncv][ncvSucc] - servingDuration[ncv]);
        }
        // check feasibility
        if (cvMaxServingTime < ncvMinServingTime) return false;
        // check max ride time
        int minRideTime = Math.max(dist[cv][ncv],
                cvMinServingTime - ncvMaxServingTime - servingDuration[ncv]);
        if (minRideTime > maxRideTime) return true;
        // check capacity
        if (capacityLeftInRoute[ncvPred].value() < vertexLoadChange[ncv]) return false;
        // compute insertion cost e
        int costIncrease = dist[cvPred][cv] + dist[cv][cvSucc] - dist[cvPred][cvSucc]
                + dist[ncvPred][ncv] + dist[ncv][ncvSucc] - dist[ncvPred][ncvSucc];
        int slackAfterInsert;
        if (ncvSucc == cv) {
            slackAfterInsert = cvMaxServingTime - servingTime[ncvPred].min();
            //-servingDuration[ncv]-dist[ncv][ncvSucc]-servingDuration[ncvPred]-dist[ncvPred][ncv];
            slackAfterInsert += servingTime[cvSucc].max() - ncvMinServingTime;
            //-servingDuration[cv]-dist[cv][cvSucc]-servingDuration[cvPred]-dist[cvPred][cv];
        } else {
            slackAfterInsert = servingTime[ncvSucc].max() - servingTime[ncvPred].min();
            //-servingDuration[ncv]-dist[ncv][ncvSucc]-servingDuration[ncvPred]-dist[ncvPred][ncv];
            slackAfterInsert += servingTime[cvSucc].max() - servingTime[cvPred].min();
            //-servingDuration[cv]-dist[cv][cvSucc]-servingDuration[cvPred]-dist[cvPred][cv];
        }
        addToInsertionObjChange(r, v, cvSucc, ncvSucc, alpha * costIncrease - beta * slackAfterInsert);
        return false;
    }


    static void setBestServingTimeFail2(int request, int v, int start) {
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int cvv = getCriticalVertex(request);
        int ncv = getCorrespondingVertex(cvv);
        int cvvMinServingTime = Math.max(
                getArrivalTimeValue(pred[start].value(), cvv, true),
                timeWindowStarts[cvv]);
        int cvvMaxServingTime = Math.min(
                servingTime[start].max() - dist[cvv][start] - servingDuration[cvv],
                timeWindowEnds[cvv]);
        if (cvvMaxServingTime < cvvMinServingTime) {
            return;
        }
        int changeCvv = servingTime[start].max()
                - (getArrivalTimeValue(pred[start].value(), cvv, true)
                + servingDuration[cvv] + dist[cvv][start]);
        if (changeCvv < 0) {
            return;
        }
        changeCvv -= 80 * (dist[cvv][start] + dist[pred[start].value()][cvv]
                - dist[pred[start].value()][start]);
        int changeNcv = 0;
        succ[cvv].setValue(start);
        pred[cvv].setValue(pred[start].value());
        succ[pred[cvv].value()].setValue(cvv);
        pred[start].setValue(cvv);
        if (isPickup(cvv)) {
            int index = start;
            int p = cvv;
            int minRideTime = maxRideTime;
            boolean done = false;
            while (index != begin && !done) {
                if (p == cvv) {
                    minRideTime = dist[cvv][ncv];
                } else {
                    minRideTime = servingTime[p].min() + servingDuration[p]
                            + dist[p][ncv] - (cvvMaxServingTime + servingDuration[cvv]);
                }
                if (minRideTime > maxRideTime) {
                    done = true;
                }
                int ncvMinServingTime = Math.max(getArrivalTimeValue(p, ncv, true),
                        timeWindowStarts[ncv]);
                int ncvMaxServingTime = Math.min(servingTime[index].max() - dist[index][ncv]
                        - servingDuration[ncv], timeWindowEnds[ncv]);

                changeNcv = servingTime[index].max()
                        - (getArrivalTimeValue(pred[index].value(), ncv, true)
                        + servingDuration[ncv] + dist[ncv][index]);

                if (ncvMaxServingTime >= ncvMinServingTime && changeNcv >= 0) {
                    changeNcv -= 80 * (dist[ncv][index] + dist[pred[index].value()][ncv]
                            - dist[pred[index].value()][index]);
                    addToInsertionObjChange(request, v, start, index, -changeCvv - changeNcv);
                    if (capacityLeftInRoute[index].value() < vertexLoadChange[cvv]) {
                        done = true;
                    }
                }
                p = index;
                index = succ[index].value();
            }
        } else {
            int p = pred[cvv].value();
            int index = cvv;
            int minRideTime = maxRideTime;
            boolean done = false;
            while (index != begin && capacityLeftInRoute[index].value() >= vertexLoadChange[ncv] && !done) {
                if (index == cvv) {
                    minRideTime = dist[ncv][cvv];
                } else {
                    minRideTime = cvvMinServingTime - (getArrivalTimeValue(p, ncv, false)
                            + servingDuration[ncv]);
                }
                if (minRideTime > maxRideTime) {
                    done = true;
                }
                int ncvMinServingTime = Math.max(getArrivalTimeValue(p, ncv, true),
                        timeWindowStarts[ncv]);
                int ncvMaxServingTime = Math.min(servingTime[index].max() - dist[index][ncv]
                        - servingDuration[ncv], timeWindowEnds[ncv]);
                changeNcv = servingTime[index].max() - (getArrivalTimeValue(pred[index].value(), ncv, true)
                        + servingDuration[ncv] + dist[ncv][index]);
                if (ncvMaxServingTime >= ncvMinServingTime && changeNcv >= 0) {
                    changeNcv -= 80 * (dist[ncv][index]
                            + dist[pred[index].value()][ncv]
                            - dist[pred[index].value()][index]);
                    addToInsertionObjChange(request, v, start, index, -changeCvv - changeNcv);
                }
                index = p;
                p = pred[p].value();
            }
        }
        succ[pred[cvv].value()].setValue(succ[cvv].value());
        pred[succ[cvv].value()].setValue(pred[cvv].value());
        succ[cvv].setValue(cvv);
        pred[cvv].setValue(cvv);
    }

    static void addToInsertionObjChange(int request, int v, int cvi, int ncvi, int change) {
        if (!insertionObjChange[request].containsKey(v)) {
            insertionObjChange[request].put(v, new HashMap<Integer, HashMap<Integer, Integer>>());
        }
        if (!insertionObjChange[request].get(v).containsKey(cvi)) {
            insertionObjChange[request].get(v).put(cvi, new HashMap<Integer, Integer>());
        }
        insertionObjChange[request].get(v).get(cvi).put(ncvi, change);
    }

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

    static void clearCustomerLeft() {
        for (StateBool cl : customersLeft) cl.setValue(false);
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
        int routeLength = 0;
        for (int v = 0; v < numVehicles; v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int i = begin;
            while (i != end) {
                routeLength += dist[i][succ[i].value()];
                i = succ[i].value();
            }
        }
        return routeLength;
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

    // variables && constraints initialization

    static void initCpVars() {
        for (int i = 0; i < numVars; i++) {
            if (i < 2 * numRequests) { // i is a site
                succ[i] = new TrailInt((Trailer) cp.getStateManager(), i);
                pred[i] = new TrailInt((Trailer) cp.getStateManager(), i);
            } else {
                if (isBeginDepot(i)) { // i is a start depot
                    succ[i] = new TrailInt((Trailer) cp.getStateManager(), getEndDepot(getVehicleOfDepot(i)));
                    pred[i] = new TrailInt((Trailer) cp.getStateManager(), succ[i].value());
                    servingVehicle[i].assign(getVehicleOfDepot(i));
                } else { // i is an end depot
                    succ[i] = new TrailInt((Trailer) cp.getStateManager(), getBeginDepot(getVehicleOfDepot(i)));
                    pred[i] = new TrailInt((Trailer) cp.getStateManager(), succ[i].value());
                    servingVehicle[i].assign(getVehicleOfDepot(i));
                }
            }
            servingTime[i] = makeIntVar(cp, timeWindowStarts[i], timeWindowEnds[i]);
        }
    }

    static void postPrecedence() {
        for (int i = 0; i < numRequests; i++) {
            cp.post(lessOrEqual(servingTime[i], minus(servingTime[numRequests + i],
                    dist[i][i + numRequests] + servingDuration[i])));
        }
        for (int i = 0; i < numVehicles; i++) {
            cp.post(lessOrEqual(servingTime[getEndDepot(i)],
                    plus(servingTime[getBeginDepot(i)], timeHorizon)));
        }
    }

    static void postRideTime() {
        for (int i = 0; i < numRequests; i++) {
            cp.post(lessOrEqual(servingTime[numRequests + i],
                    plus(servingTime[i], servingDuration[i] + maxRideTime)));
        }
    }

    static void postMaxRouteDuration() {
        for (int i = 0; i < numVehicles; i++) {
            cp.post(lessOrEqual(servingTime[getEndDepot(i)],
                    plus(servingTime[getBeginDepot(i)], maxRouteDuration)));
        }
    }

    static void postConstraints() {
        if (constraintsPosted) {
            return;
        }
        constraintsPosted = true;
        for (int i = 0; i < numRequests; i++) {
            cp.post(equal(servingVehicle[i], servingVehicle[i + numRequests]));
        }
        postPrecedence();
        postRideTime();
        postMaxRouteDuration();
    }

    static void printBestRoutesSolution() {
        for (int v = 0; v < numVehicles; v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int i = begin;
            System.out.print("vehicle: " + v + ": ");
            while (i != end && i != bestSolution.succ[i]) {
                System.out.print(i + ", ");
                i = bestSolution.succ[i];
            }
            System.out.println(end);
        }
    }

    static void printsequence(int v) {
        System.out.println("Current sequence for vehicle v: ");
        for (int i = 0; i < numVars; i++) {
            if (servingVehicle[i].isBound() && servingVehicle[i].min() == v)
                System.out.println("(" + pred[i] + ", " + i + ", " + succ[i] + ")");
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
                System.out.print("(i : " + i + " v: " + currentSolution.servingVehicle[i] + ") -> ");
                i = currentSolution.succ[i];
            }
            System.out.println("(i : " + end + " v: " + currentSolution.servingVehicle[end] + ")");
        }
    }

}

class Mut {
    int value;

    public Mut(int i) {
        value = i;
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

    def printBestRoutesSolution(): Unit = {
        for (int v=0;v<numVehicles;v++) {
            val begin = getBeginDepot(v)
            val end = getEndDepot(v)
            var i = begin
            print("vehicle: " + v + ": ")
            while (i != end && i != bestSolution.get.succ(i)) {
                print("(i : " + i + " v: " + bestSolution.get.servingVehicle(i) + ") -> ")
                i = bestSolution.get.succ(i)
            }
            println("(i : " + end + " v: " + bestSolution.get.servingVehicle(end) + ")")
        }
    }

    def printCustomerLeft(): Unit = {
        val custLeft: Array[Int] = new Array[Int](customersLeft.size)
                var h = 0
        for (i <- customersLeft) {
            custLeft(h) = i
            h += 1
        }
        println("customers Left := " + custLeft.mkString(","))
    }

    def exportSol(sol: DarpSol): DARPSolution = {
        val pathBuffer = Array.fill(numVehicles)(mutable.ArrayBuffer[DARPStep]())
        for(i <- sol.servingVehicle.indices){
            pathBuffer(sol.servingVehicle(i)) += DARPStep(i, ((sol.minServingTime(i) * 100).round.toInt, (sol.maxServingTime(i) * 100).round.toInt))
        }
        val darpSol = DARPSolution(instance, pathBuffer.indices.map(v => DARPPath(v, pathBuffer(v).toArray.sortBy(_.time._1))).toArray)
        println(darpSol.paths.map(_.steps.mkString(" -> ")).mkString("\n"))
        println()
        darpSol
    }

    def exportBestSol: DARPSolution = exportSol(bestSolution.get)


    static BranchingChoice getBestBranchingDecision(int request) {
        ArrayList<BranchingChoice> branchingQueue = new ArrayList<BranchingChoice>();
        int cvv = getCriticalVertex(request);
        int ncv = getCorrespondingVertex(cvv);
        int bestCvi = -1;
        int bestNcvi = -1;
        int bestChange = Integer.MIN_VALUE;
        for (int v=0;v<numVehicles;v++) {
            if (insertionObjChange[request].containsKey(v) && servingVehicle[request].contains(v)) {
                for (int cvi : insertionObjChange[request].get(v).keySet()) {
                    for (int ncvi : insertionObjChange[request].get(v).get(cvi).keySet()) {
                        if (insertionObjChange[request].get(v).get(cvi).get(ncvi) > bestChange) {
                            bestChange = insertionObjChange[request].get(v).get(cvi).get(ncvi);
                            bestCvi = cvi;
                            bestNcvi = ncvi;
                            branchingQueue.clear();
                            branchingQueue.add(new BranchingChoice(request, cvi, ncvi, bestChange,
                                    servingVehicle[cvi].min()));
                        }
                        else if (insertionObjChange[request].get(v).get(cvi).get(ncvi) == bestChange) {
                            branchingQueue.add(new BranchingChoice(request, cvi, ncvi, bestChange,
                                    servingVehicle[cvi].min()));
                        }
                    }
                }
            }
        }
        BranchingChoice[] bQueue = branchingQueue.toArray(new BranchingChoice[0]);
        Random rn = new Random();
        return bQueue[rn.nextInt(bQueue.length)];
    }

}
*/
