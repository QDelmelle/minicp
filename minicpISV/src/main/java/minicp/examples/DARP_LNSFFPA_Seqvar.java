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
    static int alpha = 80;
    static int beta = 1;
    static int gamma = 200;
    static int tau = 1000;
    static int maxSize;
    static int range = 4;
    static int numIter = 300;
    static double d = 0.07;
    static int maxTime = 30;
    static long remainTime = maxTime * 1000L;

//    solPath: Option[String],
//    logPath: Option[String],

    // utility parameters
    static boolean firstSolOnly = false;
    boolean silent = true;

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
    static HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>[] insertionObjChange;
    static StateBool[] customersLeft;

    public static void main(String[] args) {
        cp = makeSolver();
        String path = "C:\\Users\\Utilisateur\\Documents\\UNIF2020\\TFE\\DARP RK\\DARP\\Cordeau\\a3-24.txt";
        //path = "C:\\Users\\Utilisateur\\Documents\\UNIF2020\\TFE\\DARP RK\\DARP\\sample.txt";
        instance = DARPParser.parseInstance(path);
        System.out.println("firstSolOnly= "+firstSolOnly);

        // Data

        // Parameters of the problem
        numVars = instance.nSites;
        numRequests = instance.nRequests;
        vehicleCapacity = instance.vehicles[0].capacity;
        numVehicles = instance.nVehicles;
        maxRideTime = instance.requests[0].maxRideTime;
        maxRouteDuration = instance.vehicles[0].maxDuration;
        for (DARPStop site : instance.sites) timeHorizon = Math.max(timeHorizon, site.winEnd);
        maxSize = numRequests/2;

        System.out.println("numVars = "+numVars);
        System.out.println("numRequests = "+numRequests);
        System.out.println("numVehicles = "+numVehicles);

        // Variables necessary for posting constraints
        dist = instance.distances;
        load = new int[instance.nSites];
        for (int i=0;i<numVars;i++) load[i] = instance.sites[i].load;
        timeWindowStarts = new int[instance.nSites];
        for (int i=0;i<numVars;i++) timeWindowStarts[i] = instance.sites[i].winStart;
        timeWindowEnds = new int[instance.nSites];
        for (int i=0;i<numVars;i++) timeWindowEnds[i] = instance.sites[i].winEnd;

        // Variables and structures related to modelling
        servingTime = new IntVar[numVars];
        servingVehicle = makeIntVarArray(cp, numVars, 0, numVehicles);
        route = new InsertionSequenceVar[numVehicles];
        for (int i=0;i<numVehicles;i++) route[i] = makeISV(cp, numVars);
        capacityLeftInRoute = new StateInt[numVars];
        for (int i=0;i<numVars;i++)
            capacityLeftInRoute[i] = cp.getStateManager().makeStateInt(vehicleCapacity);
        servingDuration = new int[numVars];
        for (int i=0;i<numVars;i++) servingDuration[i] = instance.sites[i].service;
        insertionObjChange = new HashMap[numVars];
        customersLeft = new StateBool[numRequests];
        for (int i=0;i<numRequests;i++)
            customersLeft[i] = cp.getStateManager().makeStateBool(false);

        //Start solving the problem
        // initialise variables
        initCpVars();
        //Add request to the set of unassigned requests "customersLeft"
        for (int i=0;i<numRequests;i++) {
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
                //System.out.println("branching on request "+request+" with "+points.length+" insertion points...");
                Procedure[] branches = new Procedure[points.length];
                for (int i=0;i<points.length;i++) {
                    int[] p = points[i];
                    branches[i] = () -> branchRequestPoint(request, p);
                }
                return branch(branches);
            }
        });

        // On solution
        search.onSolution(() -> {
            //System.out.println("Total route cost: " + getDistanceObjective()/100.0);
            //printRoutes();
            //System.out.println("remainingTime := "+remainTime);
            int[] servingVP = new int[numVars];
            Arrays.setAll(servingVP, i -> servingVehicle[i].min());
            int[] succP = new int[numVars];
            Arrays.setAll(succP, i -> route[servingVP[i]].nextMember(i));
            int[] predP = new int[numVars];
            Arrays.setAll(predP, i -> route[servingVP[i]].prevMember(i));
            double[] minServingTime = new double[numVars];
            Arrays.setAll(minServingTime, (i) -> servingTime[i].min()/100);
            double[] maxServingTime = new double[numVars];
            Arrays.setAll(maxServingTime, (i) -> servingTime[i].max()/100);
//    println("maxRideTime: " + maxRideTime/100.0)
//    println("serving times:")
//    println(numVarsRange.map(i => i + " (" + minServintgTime(i) + ":" + maxServintgTime(i) + ")").mkString("\n"))
            double rand = Math.random();
            int obj = getDistanceObjective();
            if (obj < currentSolutionObjective || rand < d) {
                currentSolution = new DarpSol(succP, predP, servingVP,obj/100.0, minServingTime, maxServingTime);
                currentSolutionObjective = obj;
                if (currentSolutionObjective < bestSolutionObjective) {
                    bestSolution = new DarpSol(succP, predP, servingVP,
                            obj/100.0, minServingTime, maxServingTime);
                    bestSolutionObjective = currentSolutionObjective;
                    System.out.println("new best solution found:= " + bestSolutionObjective/100.0);
                    System.out.println("remainingTime := "+remainTime);
                }
            }
            //System.out.println("-------------------------");
        });
        boolean solFound = false;

        //search.solve(); System.out.println("best obj = "+bestSolutionObjective);
        while(!solFound && remainTime > 0){
            long t1 = System.currentTimeMillis();
            SearchStatistics stats = search.solve(statistics -> statistics.numberOfSolutions() == 1
                    && statistics.numberOfFailures() <= Math.max(gamma * numVehicles, tau));
            // + (int)Math.round(remainTime / 1000.0)); ??
            if(stats.numberOfSolutions() == 1) solFound = true;
            else System.out.println("restart! "+ remainTime);
            //System.out.println(stats);
            //solFound = true; // if we don't want restarts
            remainTime -= System.currentTimeMillis() - t1;
        }

        if(!firstSolOnly) lns();
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
        while(remainTime > 0 && i <= maxSize - range) {
            int j = 0;
            while(remainTime > 0 && j <= range) {
                int k = 1;
                while(remainTime > 0 && k <= numIter) {
                    int finalI = i;
                    int finalJ = j;
                    long t1 = System.currentTimeMillis();
                    SearchStatistics stats = search.solveSubjectTo(
                            statistics -> statistics.numberOfSolutions() == 1,
                            () -> {
                                relax(finalI + finalJ);
                            });
                    k++;
                    remainTime -= System.currentTimeMillis()-t1;
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
        clearCustomerLeft(); // needed?
        for (int r : relaxedCustomers) customersLeft[r].setValue(true);
        for (int v=0;v<numVehicles;v++) {
            route[v].resetDomain();
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            // these are only propagated on post, so we must re-post them
            cp.post(new First(route[v], begin));
            cp.post(new Last(route[v], end));
            int current = solSucc[begin];
            int prev = begin;
            while (current != end) {
                int r = getCorrespondingRequest(current);
                if(!relaxedCustomers.contains(r)){
                    if (!insertVertexIntoRoute(current, prev, v)) System.out.println("relax fail!!!");
                    equal(servingVehicle[current], v);
                    prev = current;
                }
                current = solSucc[current];
            }
            updateCapacityLeftInRoute(v, -1);
        }
    }

    static int getUnassignedMinVehicleMinInsertionPointsRequest() {
        int minVehicles = numVehicles + 1;
        List<Integer> unassignedRequests = new LinkedList<Integer>();

        // step 1: minimize number of routes
        for (int i=0 ; i<numRequests; i++) {
            if (customersLeft[i].value()) {
                unassignedRequests.add(i);
                minVehicles = Math.min(minVehicles, servingVehicle[i].size());
            }
        }
        List<Integer> minRoutesRequests = new LinkedList<Integer>();
        for (int i: unassignedRequests) {
            if (servingVehicle[i].size() == minVehicles) minRoutesRequests.add(i);
        }

        // step 2: minimize number of insertion points
        int minChoices = Integer.MAX_VALUE;
        int[] bestInsertionCost = new int[numRequests];
        int[] numInsertions = new int[numRequests];
        for (int r: minRoutesRequests) {
            int[][] iP = getInsertionPoints(r);
            bestInsertionCost[r] = iP[0][3];
            numInsertions[r] = iP.length;
            minChoices = Math.min(minChoices, iP.length);
        }
        List<Integer> minInsPointsRequests = new LinkedList<Integer>();
        for (int r: minRoutesRequests) {
            if (numInsertions[r] == minChoices) minInsPointsRequests.add(r);
        }

        // step 3: minimize cost of best insertion
        int bestChange = Integer.MIN_VALUE;
        for (int r: minInsPointsRequests) {
            bestChange = Math.max(bestChange, bestInsertionCost[r]);
        }
        List<Integer> minBestCostRequests = new LinkedList<Integer>();
        for (int r: minInsPointsRequests) {
            if (bestInsertionCost[r] == bestChange) minBestCostRequests.add(r);
        }

        Random rn = new Random();
        return minBestCostRequests.get(rn.nextInt(minBestCostRequests.size()));
    }

    static int getUnassignedRequest() {
        for (int i=0 ; i<numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            insertionObjChange[i] = new HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>();
            for (int v=0;v<numVehicles;v++) {
                setInsertionCost(i, v);
            }
        }
        return getUnassignedMinVehicleMinInsertionPointsRequest();
    }

    /**
     *
     * @param request
     * @return all the insertion points of request, sorted in increasing order of insertion cost.
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
                return a[3]-b[3];
            }
        });
        return vehiclePointsChange;
    }

    static void branchRequestPoint(int request, int[] point) {
        int vehicle = point[0], cvPrev = point[1], ncvPrev = point[2], change = point[3];
        //System.out.println(request+" "+cvPrev+" "+ncvPrev+" "+change);
        int pickup = request, drop = request + numRequests;
        int pPrev = isInbound(request) ? cvPrev : ncvPrev;
        int dPrev = isInbound(request) ? ncvPrev : cvPrev;

        equal(servingVehicle[pickup], vehicle);
        // /!\ insert pickup first
        if (!insertVertexIntoRoute(pickup, pPrev, vehicle)) {
            throw new InconsistencyException();
        }
        if (!insertVertexIntoRoute(drop, dPrev, vehicle)) {
            throw new InconsistencyException();
        }
        updateCapacityLeftInRoute(vehicle, pickup);
        if (!isPositive(capacityLeftInRoute)) { // not sure this is useful
            throw new InconsistencyException();
        }
        customersLeft[request].setValue(false);

        // recompute the insertions of all remaining unassigned requests for this vehicule
        // and check consistency.
        for (int i=0 ; i<numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            if (insertionObjChange[i].containsKey(vehicle)) {
                insertionObjChange[i].remove(vehicle);
            }
            setInsertionCost(i, vehicle);
            int[][] iP = getInsertionPoints(i);
            if(iP.length == 0){
                throw new InconsistencyException();
            }
        }
    }

    // update capacity for vehicule v starting from node <start> in route.
    // start MUST be part of the route of v.
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

    // insert i after j in route v.
    // return false if the insertion fails.
    static boolean insertVertexIntoRoute(int i, int j, int v) {
        //System.out.println(i+" "+j);
        int s = route[v].nextMember(j);
        if (!tryPost(new LessOrEqual(getArrivalTime(j, i), servingTime[i])) ||
                !tryPost(new LessOrEqual(getArrivalTime(i, s), servingTime[s])))
            return false;
        route[v].insert(i, j);
        return true;
    }

    // finds all possible insertions of request in route v and compute their cost.
    static void setInsertionCost(int request, int v) {
        int cv = getCriticalVertex(request);
        int ncv = getCorrespondingVertex(cv);
        for (int cvPred : route[v].getInserts(cv)) {
            if (!isInbound(request) || capacityLeftInRoute[cvPred].value() >= load[cv]) {
                // simulate insertion of cv and look for insertions of ncv
                cp.getStateManager().saveState();
                if (insertVertexIntoRoute(cv, cvPred, v)) {
                    updateCapacityLeftInRoute(v, cv);
                    for (int ncvPred : route[v].getInserts(ncv)) {
                        if (isInbound(request) || capacityLeftInRoute[ncvPred].value() >= load[ncv]) {
                            // simulate insertion of ncv
                            cp.getStateManager().saveState();
                            if (insertVertexIntoRoute(ncv, ncvPred, v)) {
                                updateCapacityLeftInRoute(v, ncv);
                                computeInsertionCost(request, v, cv, ncv);
                            }
                            cp.getStateManager().restoreState();
                        }
                    }
                }
                cp.getStateManager().restoreState();
            }
        }
    }

    private static void computeInsertionCost(int r, int v, int cv, int ncv) {
        int cvPred = route[v].prevMember(cv);
        int cvSucc = route[v].nextMember(cv);
        int ncvPred = route[v].prevMember(ncv);
        int ncvSucc = route[v].nextMember(ncv);
        // compute insertion cost e
        int costIncrease = dist[cvPred][cv] + dist[cv][cvSucc] - dist[cvPred][cvSucc]
                + dist[ncvPred][ncv] + dist[ncv][ncvSucc] - dist[ncvPred][ncvSucc];
        int slackAfterInsert;
        slackAfterInsert = servingTime[ncvSucc].max() - servingTime[ncvPred].min();
        //-servingDuration[ncv]-dist[ncv][ncvSucc]-servingDuration[ncvPred]-dist[ncvPred][ncv];
        slackAfterInsert += servingTime[cvSucc].max() - servingTime[cvPred].min();
        //-servingDuration[cv]-dist[cv][cvSucc]-servingDuration[cvPred]-dist[cvPred][cv];
        addToInsertionObjChange(r, v, cvPred, ncvPred, alpha*costIncrease - beta*slackAfterInsert);
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
        while(relaxEnd < numCustomersToRelax && relaxEnd < customers.length){
            int toRelax = relaxEnd + rn.nextInt(numRequests-relaxEnd);
            int cRelaxed = customers[toRelax];
            customers[toRelax] = customers[relaxEnd];
            customers[relaxEnd] = cRelaxed;
            relaxEnd++;
        }
        Set<Integer> ret = new HashSet<Integer>();
        for (int i=0;i<relaxEnd;i++) ret.add(customers[i]);
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
        for (int v=0;v<numVehicles;v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int i = begin;
            while (i != end) {
                routeLength += dist[i][route[v].nextMember(i)];
                i = route[v].nextMember(i);
            }
        }
        return routeLength;
    }

    static boolean isPositive(StateInt[] array) {
        for (int i=0; i<array.length; i++) {
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

    static int getBeginDepot(int i) { return  2 * numRequests + i; }

    static int getEndDepot(int i) { return numVars - numVehicles + i; }

    static boolean isBeginDepot(int i) { return (i >= 2 * numRequests && i < numVars - numVehicles); }

    boolean isEndDepot(int i) { return i >= numVars - numVehicles && i < numVars; }

    static IntVar getArrivalTime(int vertex, int successor) {
        int distance = dist[vertex][successor];
        if (isBeginDepot(vertex)) return new IntVarViewOffset(servingTime[vertex], distance);
        else return new IntVarViewOffset(servingTime[vertex], servingDuration[vertex] + distance);
    }

    static int getCorrespondingRequest(int i) { return (i < numRequests) ? i : i - numRequests; }

    static int getCorrespondingPickup(int i) { return i - numRequests; }

    static int getCorrespondingDelivery(int i) { return i + numRequests; }

    static boolean isInbound(int r) { return timeWindowStarts[r] > 0 || timeWindowEnds[r] < timeHorizon; }

    static int getCriticalVertex(int request) {
        return (isInbound(request)) ?
                request : request + numRequests; }

    static int getCorrespondingVertex(int i) { return (isPickup(i)) ?
            getCorrespondingDelivery(i) : getCorrespondingPickup(i); }

    static boolean isPickup(int i) { return i < numRequests; }

    boolean isDelivery(int i) { return i >= numRequests && i < 2 * numRequests; }

    boolean isCustomerVertex(int i) { return i < 2 * numRequests; }

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
        for (int i=0;i<numVars;i++) {
            if (i>=2*numRequests) { // i is a depot
                servingVehicle[i].assign(getVehicleOfDepot(i));
            }
            servingTime[i] = makeIntVar(cp, timeWindowStarts[i], timeWindowEnds[i]);
        }
    }

    static void postDependency() { // pas besoin?
        for (int i=0;i<numRequests;i++) {
            equal(servingVehicle[i], servingVehicle[i+numRequests]);
            for (int s=0;s<numVehicles;s++) {
                cp.post(new Dependency(route[s], new int[]{i, i+numRequests}));
            }
        }
    }

    static void postPrecedence() {
        for (int i=0;i<numRequests;i++) {
            cp.post(lessOrEqual(servingTime[i], minus(servingTime[numRequests + i],
                    dist[i][i + numRequests] + servingDuration[i])));
            for (int v=0;v<numVehicles;v++)
                cp.post(new Precedence(route[v], new int[]{i, i+numRequests}));
        }
    }

    static void postRideTime() {
        for (int i=0;i<numRequests;i++) {
            cp.post(lessOrEqual(servingTime[numRequests + i],
                    plus(servingTime[i], servingDuration[i] + maxRideTime)));
        }
    }

    static void postMaxRouteDuration() {
        for (int i=0;i<numVehicles;i++) {
            cp.post(lessOrEqual(servingTime[getEndDepot(i)],
                    plus(servingTime[getBeginDepot(i)], maxRouteDuration)));
        }
    }

    static void postTransitionTimes() {} // not needed?

    static void postConstraints() {
        if (constraintsPosted) {
            return;
        }
        constraintsPosted = true;
        for (int i=0;i<numVehicles;i++) {
            cp.post(new First(route[i], getBeginDepot(i)));
            cp.post(new Last(route[i], getEndDepot(i)));
            cp.post(lessOrEqual(servingTime[getEndDepot(i)],
                    plus(servingTime[getBeginDepot(i)], timeHorizon)));
        }
        for (int i=0;i<numRequests;i++)
            cp.post(equal(servingVehicle[i], servingVehicle[i+numRequests]));
        //postDependency();
        postPrecedence();
        postTransitionTimes();
        postRideTime();
        postMaxRouteDuration();
    }

    static void printBestRoutesSolution() {
        for (int v=0;v<numVehicles;v++) {
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

    private static void printRoutes() {
        for (int v=0;v<numVehicles;v++){
            System.out.println("v: "+v+", route: "+route[v].allMembers());
        }
    }

    /*
    static void printsequence(int v) {
        System.out.println("Current sequence for vehicle v: ");
        for (int i=0;i<numVars;i++){
            if (servingVehicle[i].isBound() && servingVehicle[i].min() == v)
                System.out.println("("+pred[i]+", "+i+", "+succ[i]+")");
        }
    }*/

    static void printCurrentSolution() {
        System.out.println("Darp Current Solution");
        for (int v=0;v<numVehicles;v++) {
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
