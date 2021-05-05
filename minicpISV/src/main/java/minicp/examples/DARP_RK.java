package minicp.examples;

import minicp.engine.constraints.Equal;
import minicp.engine.constraints.LessOrEqual;
import minicp.engine.core.*;

import minicp.examples.DARPParser.*;
import minicp.examples.DARPDataModel.*;
import minicp.search.DFSearch;
import minicp.search.SearchStatistics;
import minicp.state.CopyBool;
import minicp.state.CopyInt;
import minicp.state.StateBool;
import minicp.state.StateInt;
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
    static boolean firstSolOnly = true;
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
    static int timeHorizon = 0;

    // Variables necessary for posting constraints
    static int[][] travelTimeMatrix;
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

    // Ranges
    /*
    val numVarsRange = 0 until numVars
    val numRequestsRange = 0 until numRequests
    val beginDepotRange = 2 * numRequests until numVars - numVehicles
    val endDepotRange = numVars - numVehicles until numVars
    val numVehiclesRange = 0 until numVehicles
    val customerVertexRange = 0 until 2 * numRequests
    val succRange = 0 until numVars - numVehicles
    val predRange = customerVertexRange ++ endDepotRange
    */

    public static void main(String[] args) {
        cp = makeSolver();
        instance = DARPParser.parseInstance("C:\\Users\\Utilisateur\\Documents\\UNIF2020\\TFE\\DARP RK\\DARP\\Cordeau\\a2-16.txt");
        System.out.println("firstSolOnly= "+firstSolOnly);

        // Data

        // Parameters of the problem
        numVars = instance.nSites;
        numRequests = instance.nRequests;
        vehicleCapacity = instance.vehicles[0].capacity;
        numVehicles = instance.nVehicles;
        maxRideTime = instance.requests[0].maxRideTime;
        for (DARPStop site : instance.sites) timeHorizon = Math.max(timeHorizon, site.winEnd);
        maxSize = numRequests/2;

        // Variables necessary for posting constraints
        travelTimeMatrix = instance.distances;
        vertexLoadChange = new int[instance.nSites];
        for (int i=0;i<numVars;i++) vertexLoadChange[i] = instance.sites[i].load;
        timeWindowStarts = new int[instance.nSites];
        for (int i=0;i<numVars;i++) timeWindowStarts[i] = instance.sites[i].winStart;
        timeWindowEnds = new int[instance.nSites];
        for (int i=0;i<numVars;i++) timeWindowEnds[i] = instance.sites[i].winEnd;

        // Variables and structures related to modelling
        servingTime = new IntVar[numVars];
        servingVehicle = makeIntVarArray(cp, numVars, 0, numVehicles);
        succ = new CopyInt[numVars];
        pred = new CopyInt[numVars];
        capacityLeftInRoute = new CopyInt[numVars];
        for (int i=0;i<numVars;i++) capacityLeftInRoute[i].setValue(vehicleCapacity);
        servingDuration = new int[numVars];
        for (int i=0;i<numVars;i++) servingDuration[i] = instance.sites[i].service;
        insertionObjChange = new HashMap[numVars];
        customersLeft = new StateBool[numRequests];
        for (int i=0;i<numRequests;i++) customersLeft[i] = new CopyBool(false);

        //Start solving the problem
        // initialise variables
        initCpVars();
        //Add request to set of unassigned request "customersLeft"
        for (int i=0;i<numRequests;i++) {
            customersLeft[i].setValue(true);
        }

        // post constraints
        postConstraints();
        postUnaryConstraint();
        postCumulativeConstraint();

        // search algorithm
        search = makeDfs(cp, () -> {
            boolean CL = false;
            for (StateBool C : customersLeft) if (C.value()) CL = true;
            if (!CL) {
                return EMPTY;
            } else {
                int request = getUnassignedRequest();
                int[][] points = getInsertionPoints(request);
                Procedure[] branches = new Procedure[points.length];
                for (int i=0;i<points.length;i++) {
                    int[] p = points[i];
                    branches[i] = () -> branchRequestPoint(request, new int[]{p[1], p[2], p[3], p[4]});
                }
                return branch(branches);
            }
        });

        // On solution
        search.onSolution(() -> {
            System.out.println("Total route cost: " + getDistanceObjective()/100.0);
            System.out.println("best solution found:= " + bestSolutionObjective/100.0);
            System.out.println("remainingTime := "+remainTime);
            int[] succP = new int[numVars];
            Arrays.setAll(succP, (i) -> succ[i].value());
            int[] predP = new int[numVars];
            Arrays.setAll(predP, (i) -> pred[i].value());
            int[] servingVP = new int[numVars];
            Arrays.setAll(servingVP, (i) -> servingVehicle[i].min());
            double[] minServingTime = new double[numVars];
            Arrays.setAll(minServingTime, (i) -> servingTime[i].min()/100);
            double[] maxServingTime = new double[numVars];
            Arrays.setAll(maxServingTime, (i) -> servingTime[i].max()/100);
//    println("maxRideTime: " + maxRideTime/100.0)
//    println("serving times:")
//    println(numVarsRange.map(i => i + " (" + minServintgTime(i) + ":" + maxServintgTime(i) + ")").mkString("\n"))
            double rand = Math.random();
            if (getDistanceObjective() < currentSolutionObjective || rand < d) {
                currentSolution = new DarpSol(succP, predP, servingVP, getDistanceObjective()/100.0, minServingTime, maxServingTime);
                currentSolutionObjective = getDistanceObjective();
                if (currentSolutionObjective < bestSolutionObjective) {
                    bestSolution = new DarpSol(succP, predP, servingVP, getDistanceObjective()/100.0, minServingTime, maxServingTime);
                    bestSolutionObjective = currentSolutionObjective;
                }
            }
            System.out.println("-------------------------");
        });
        boolean solFound = false;

//  println(start(nSols = Int.MaxValue, failureLimit = Int.MaxValue, timeLimit = (remainTime/1000.0).round.toInt))

        while(!solFound && remainTime > 0){
            SearchStatistics stats = search.solve();
                    //(1, Math.max(gamma * numVehicles, tau), (int)Math.round(remainTime / 1000.0)); ??
            if(stats.numberOfSolutions() == 1) solFound = true;
            System.out.println(stats);
            //remainTime -= (int)Math.round(stats.time/1000.0); ?? no time in stats?
        }

        if(!firstSolOnly) lns();

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
        int dist = travelTimeMatrix[vertex][successor];
        if (isBeginDepot(vertex)) return new IntVarViewOffset(servingTime[vertex], dist);
        else return new IntVarViewOffset(servingTime[vertex], servingDuration[vertex] + dist);
    }

    static int getCorrespondingRequest(int i) { return (i < numRequests) ? i : i - numRequests; }

    static int getCorrespondingPickup(int i) { return i - numRequests; }

    static int getCorrespondingDelivery(int i) { return i + numRequests; }

    static int getCriticalVertex(int request) {
        return (timeWindowStarts[request] > 0 || timeWindowEnds[request] < timeHorizon) ?
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
            if (getMin) return servingTime[vertex].min() + travelTimeMatrix[vertex][successor];
            else return servingTime[vertex].max() + travelTimeMatrix[vertex][successor];
        }
        if (getMin) return servingTime[vertex].min() + servingDuration[vertex] + travelTimeMatrix[vertex][successor];
        else return servingTime[vertex].max() + servingDuration[vertex] + travelTimeMatrix[vertex][successor];
    }

    // variables && constraints initialization

    static void initCpVars() {
        for (int i=0;i<numVars;i++) {
            if (i<2*numRequests) { // i is a site
                succ[i] = new CopyInt(i);
                pred[i] = new CopyInt(i);
            } else {
                if (isBeginDepot(i)) { // i is a start depot
                    succ[i] = new CopyInt(getEndDepot(getVehicleOfDepot(i)));
                    pred[i] = new CopyInt(succ[i].value());
                    servingVehicle[i].assign(getVehicleOfDepot(i));
                } else { // i is an end depot
                    succ[i] = new CopyInt(getBeginDepot(getVehicleOfDepot(i)));
                    pred[i] = new CopyInt(succ[i].value());
                    servingVehicle[i].assign(getVehicleOfDepot(i));
                }
            }
            servingTime[i] = makeIntVar(cp, timeWindowStarts[i], timeWindowEnds[i]);
        }
    }

    static void postPrecedence() {
        for (int i=0;i<numRequests;i++) {
            cp.post(lessOrEqual(servingTime[i], minus(servingTime[numRequests + i],
                    travelTimeMatrix[i][i + numRequests] + servingDuration[i])));
        }
        for (int i=0;i<numVehicles;i++) {
            cp.post(lessOrEqual(servingTime[getEndDepot(i)],
                    plus(servingTime[getBeginDepot(i)], timeHorizon)));
        }
    }

    static void postRideTime() {
        for (int i=0;i<numRequests;i++) {
            cp.post(lessOrEqual(servingTime[numRequests + i],
                    plus(servingTime[i], servingDuration[i] + maxRideTime)));
        }
    }

    static void postConstraints() {
        if (constraintsPosted) {
            return;
        }
        constraintsPosted = true;
        for (int i=0;i<numRequests;i++) {
            cp.post(equal(servingVehicle[i], servingVehicle[i + numRequests]));
        }
        postPrecedence();
        postRideTime();
    }

    static void postCumulativeConstraint() {
        /* ???
        val travelStart: Array[CPIntVar] = Array.tabulate(numRequests)(i => servingTime(i))
        val travelEnd: Array[CPIntVar] = Array.tabulate(numRequests)(i => servingTime(i+numRequests))
        val travelDuration: Array[CPIntVar] = Array.tabulate(numRequests)(i => travelEnd(i) - travelStart(i))
        val travelLoad = Array.tabulate(numRequests)(i => CPIntVar(vertexLoadChange(i)))
        val travelVehicle = Array.tabulate(numRequests)(i => servingVehicle(i))
        for(int i=0;i<numRequests;i++) add(travelDuration(i) > 0)
        for(int v=0;v<numVehicles;v++) {
            //Adding cumulative constraint:
            val capVar = CPIntVar(vehicleCapacity)
            add(maxCumulativeResource(travelStart, travelDuration, travelEnd,
                    travelLoad, travelVehicle, capVar, v))
        }
        */
    }

    static void postUnaryConstraint() {
        /* ???
        val start: Array[CPIntVar] = Array.tabulate(2*numRequests)(i => servingTime(i))
        val duration: Array[CPIntVar] = Array.tabulate(2*numRequests)(i => CPIntVar(servingDuration(i)))
        val end: Array[CPIntVar] = Array.tabulate(2*numRequests)(i => start(i) + duration(i))

        val resource = Array.tabulate(2*numRequests)(i => servingVehicle(i))
        for(int v=0;v<numVehicles;v++){
            //Adding unary constraint:
            add(unaryResource(start, duration, end, resource, v))
        }
        */
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
                    SearchStatistics stats = search.solveSubjectTo(
                            statistics -> statistics.numberOfSolutions() == 1,
                            () -> {
                                relax(finalI + finalJ);
                            });

                    //val stats = startSubjectTo(1, Int.MaxValue, (remainTime/1000.0).round.toInt){ relax(i + j) }
                    k++;
                    //remainTime -= stats.time // stats have no time??
                    totalNumFails += stats.numberOfFailures();
                    System.out.println(stats);
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
        for (int v=0;v<numVehicles;v++) {
            int begin = getBeginDepot(v);
            int end = getEndDepot(v);
            int current = begin;
            int prev = -1;
            while (current != end) {
                int currentRequest = getCorrespondingRequest(current);
                if(!relaxedCustomers.contains(currentRequest)){
                    if(prev != -1){
                        succ[prev].setValue(current);
                        pred[current].setValue(prev);
                        lessOrEqual(getArrivalTime(prev, current), servingTime[current]);
                        equal(servingVehicle[current], solServingVehicle[current]);
                    }
                    prev = current;
                }
                current = solSucc[current];
            }
            if(prev != -1){
                succ[prev].setValue(end);
                pred[end].setValue(prev);
                lessOrEqual(getArrivalTime(prev, end), servingTime[end]);
            }
            updateCapacityLeftInRoute(v, -1);
        }
    }

    static int getUnassignedMinVehicleMinInsertionPointsRequest() {
        int minChoices = Integer.MAX_VALUE;
        int minVehicles = numVehicles + 1;
        List<BranchingChoice> bQueueBuffer = new ArrayList<BranchingChoice>();
        int bestChange = Integer.MIN_VALUE;
        for (int i=0 ; i<numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            int tempChange = Integer.MIN_VALUE;
            Mut numChoices = new Mut(0);
            Iterable<BranchingChoice> branchingQueue = getInsertionCost(i, numChoices);
            if (servingVehicle[i].size() < minVehicles) {
                minVehicles = servingVehicle[i].size();
                bQueueBuffer.clear();
                for (BranchingChoice j : branchingQueue) {
                    bQueueBuffer.add(j);
                }
            }
            else if (servingVehicle[i].size() == minVehicles && numChoices.value < minChoices) {
                minChoices = numChoices.value;
                bQueueBuffer.clear();
                for (BranchingChoice j : branchingQueue) {
                    bQueueBuffer.add(j);
                }
            }
            else if (servingVehicle[i].size() == minVehicles && numChoices.value == minChoices && tempChange >= bestChange) {
                for (BranchingChoice j : branchingQueue) {
                    bQueueBuffer.add(j);
                }
            }
        }
        BranchingChoice[] bQueue = bQueueBuffer.toArray(new BranchingChoice[0]);
        Random rn = new Random();
        int currentRequest = bQueue[rn.nextInt(bQueue.length)].request;
        BranchingChoice bc = getBestBranchingDecision(currentRequest);
        return bc.request;
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
        int vehicle = point[0], cvSucc = point[1], ncvSucc = point[2], change = point[3];
        int cvv = getCriticalVertex(request);
        int ncv = getCorrespondingVertex(cvv);

        if (!tryPost(new Equal(servingVehicle[request], servingVehicle[cvSucc]))) {
            throw new InconsistencyException();
        }
        if (!insertVertexIntoRoute(cvv, cvSucc)) {
            throw new InconsistencyException();
        }
        if (!insertVertexIntoRoute(ncv, ncvSucc)) {
            throw new InconsistencyException();
        }
        updateCapacityLeftInRoute(servingVehicle[request].min(), request);
        if (!isPositive(capacityLeftInRoute)) {
            throw new InconsistencyException();
        }
        customersLeft[request].setValue(false);
        int v = servingVehicle[cvSucc].min();
        for (int i=0 ; i<numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            if (insertionObjChange[i].containsKey(v)) {
                insertionObjChange[i].remove(v);
            }
            setInsertionCost(i, v);
            int[][] insertionPoint = getInsertionPoints(i);
            if(insertionPoint.length == 0){
                throw new InconsistencyException();
            }
        }
    }

    static void updateCapacityLeftInRoute(int v, int start) {
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int star = start;
        if (star == -1) {
            star = succ[begin].value();
        }
        int capacity = capacityLeftInRoute[pred[star].value()].value();
        while (star != end) {
            capacity -= vertexLoadChange[star];
            capacityLeftInRoute[star].setValue(capacity);
            star = succ[star].value();
        }
    }

    static boolean insertVertexIntoRoute(int i, int j) {
        if (!tryPost(new LessOrEqual(getArrivalTime(i, j), servingTime[j])) ||
                !tryPost(new LessOrEqual(getArrivalTime(pred[j].value(), i), servingTime[i])))
            return false;
        succ[i].setValue(j);
        pred[i].setValue(pred[j].value());
        succ[pred[i].value()].setValue(i);
        pred[j].setValue(i);
        return true;
    }

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

    static Iterable<BranchingChoice> getInsertionCost(int request, Mut numChoices) {
        ArrayList<BranchingChoice> branchingQueue = new ArrayList<BranchingChoice>();
        numChoices.value = 0;
        int cvv = getCriticalVertex(request);
        int ncv = getCorrespondingVertex(cvv);
        int bestCvi = -1;
        int bestNcvi = -1;
        int bestChange = Integer.MIN_VALUE;

        for (int v=0;v<numVehicles;v++) {
            if (insertionObjChange[request].containsKey(v) && servingVehicle[request].contains(v)) {
                for (int cvi : insertionObjChange[request].get(v).keySet()) {
                    for (int ncvi : insertionObjChange[request].get(v).get(cvi).keySet()) {
                        numChoices.value++;
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
        return branchingQueue;
    }

    static void setInsertionCost(int request, int v) {
        int cvv = getCriticalVertex(request);
        int ncv = getCorrespondingVertex(cvv);
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int start = succ[begin].value();
        while(start != end){
            if(getArrivalTimeValue(pred[start].value(), cvv, true) <= servingTime[cvv].max()
                    && getArrivalTimeValue(cvv, start, true) <= servingTime[start].max()){
                setBestServingTimeFail(request, v, start);
            }
            start = succ[start].value();
        }
        if(getArrivalTimeValue(pred[start].value(), cvv, true) <= servingTime[cvv].max()
                && getArrivalTimeValue(cvv, start, true) <= servingTime[start].max()){
            setBestServingTimeFail(request, v, start);
        }
    }

    static void setBestServingTimeFail(int request, int v, int start) {
        int begin = getBeginDepot(v);
        int end = getEndDepot(v);
        int cvv = getCriticalVertex(request);
        int ncv = getCorrespondingVertex(cvv);
        int cvvMinServingTime = Math.max(
                getArrivalTimeValue(pred[start].value(), cvv, true),
                timeWindowStarts[cvv]);
        int cvvMaxServingTime = Math.min(
                servingTime[start].max() - travelTimeMatrix[start][cvv] - servingDuration[cvv],
                timeWindowEnds[cvv]);
        if (cvvMaxServingTime < cvvMinServingTime) {
            return;
        }
        int changeCvv = servingTime[start].max()
                - (getArrivalTimeValue(pred[start].value(), cvv, true)
                + servingDuration[cvv] + travelTimeMatrix[cvv][start]);
        if (changeCvv < 0) {
            return;
        }
        changeCvv -= 80 * (travelTimeMatrix[cvv][start] + travelTimeMatrix[pred[start].value()][cvv]
                - travelTimeMatrix[pred[start].value()][start]);
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
                    minRideTime = travelTimeMatrix[cvv][ncv];
                }
                else {
                    minRideTime = servingTime[p].min() + servingDuration[p]
                            + travelTimeMatrix[p][ncv] - (cvvMaxServingTime + servingDuration[cvv]);
                }
                if (minRideTime > maxRideTime) {
                    done = true;
                }
                int ncvMinServingTime = Math.max(getArrivalTimeValue(p, ncv, true),
                        timeWindowStarts[ncv]);
                int ncvMaxServingTime = Math.min(servingTime[index].max() - travelTimeMatrix[index][ncv]
                        - servingDuration[ncv], timeWindowEnds[ncv]);

                changeNcv = servingTime[index].max()
                        - (getArrivalTimeValue(pred[index].value(), ncv, true)
                        + servingDuration[ncv] + travelTimeMatrix[ncv][index]);

                if (ncvMaxServingTime >= ncvMinServingTime && changeNcv >= 0) {
                    changeNcv -= 80 * (travelTimeMatrix[ncv][index] + travelTimeMatrix[pred[index].value()][ncv]
                            - travelTimeMatrix[pred[index].value()][index]);
                    addToInsertionObjChange(request, v, start, index, changeCvv + changeNcv);
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
                    minRideTime = travelTimeMatrix[ncv][cvv];
                }
                else {
                    minRideTime = cvvMinServingTime - (getArrivalTimeValue(p, ncv, false)
                            + servingDuration[ncv]);
                }
                if (minRideTime > maxRideTime) {
                    done = true;
                }
                int ncvMinServingTime = Math.max(getArrivalTimeValue(p, ncv, true),
                        timeWindowStarts[ncv]);
                int ncvMaxServingTime = Math.min(servingTime[index].max() - travelTimeMatrix[index][ncv]
                        - servingDuration[ncv], timeWindowEnds[ncv]);
                changeNcv = servingTime[index].max() - (getArrivalTimeValue(pred[index].value(), ncv, true)
                        + servingDuration[ncv] + travelTimeMatrix[ncv][index]);
                if(ncvMaxServingTime >= ncvMinServingTime && changeNcv >= 0){
                    changeNcv -= 80 * (travelTimeMatrix[ncv][index] + travelTimeMatrix[pred[index].value()][ncv]
                            - travelTimeMatrix[pred[index].value()][index]);
                    addToInsertionObjChange(request, v, start, index, changeCvv + changeNcv);
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
        int relaxEnd = 0;
        while(relaxEnd < numCustomersToRelax && relaxEnd < customers.length){
            Random rn = new Random();
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
                routeLength += travelTimeMatrix[i][succ[i].value()];
                i = succ[i].value();
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

    def printCurrentSolution(): Unit = {
        println("Darp Current Solution")
        for (int v=0;v<numVehicles;v++) {
            val begin = getBeginDepot(v)
            val end = getEndDepot(v)
            var i = begin
            print("vehicle: " + v + ": ")
            while (i != end && i != currentSolution.get.succ(i)) {
                print("(i : " + i + " v: " + currentSolution.get.servingVehicle(i) + ") -> ")
                i = currentSolution.get.succ(i)
            }
            println("(i : " + end + " v: " + currentSolution.get.servingVehicle(end) + ")")
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

}
*/
