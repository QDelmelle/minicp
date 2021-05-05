package minicp.examples;

import minicp.engine.constraints.Equal;
import minicp.engine.constraints.LessOrEqual;
import minicp.engine.core.BoolVar;
import minicp.engine.core.IntVar;
import minicp.engine.core.IntVarViewOffset;
import minicp.engine.core.Solver;

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

        // On solution
        /*
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
                    branches[i] = () -> branchRequestPoint(request, (p[1], p[2], p[3], p[4]);
                }
                return branch(branches);
            }
        });

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
        */
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

    int getCorrespondingPickup(int i) { return i - numRequests; }

    int getCorrespondingDelivery(int i) { return i + numRequests; }

    int getCriticalVertex(int request) {
        return (timeWindowStarts[request] > 0 || timeWindowEnds[request] < timeHorizon) ?
                request : request + numRequests; }

    int getCorrespondingVertex(int i) { return (isPickup(i)) ?
            getCorrespondingDelivery(i) : getCorrespondingPickup(i); }

    boolean isPickup(int i) { return i < numRequests; }

    boolean isDelivery(int i) { return i >= numRequests && i < 2 * numRequests; }

    boolean isCustomerVertex(int i) { return i < 2 * numRequests; }

    boolean isCriticalVertex(int vertex) {
        return timeWindowStarts[vertex] > 0 || timeWindowEnds[vertex] < timeHorizon;
    }

    int getArrivalTimeValue(int vertex, int successor, boolean getMin) {
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

    void relax(int nRelax) {
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

    int getUnassignedMinVehicleMinInsertionPointsRequest() {
        int minChoices = Integer.MAX_VALUE;
        int minVehicles = numVehicles + 1;
        List<BranchingChoice> bQueueBuffer = new ArrayList<BranchingChoice>();
        int bestChange = Integer.MIN_VALUE;
        for (int i=0 ; i<numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            int tempChange = Integer.MIN_VALUE;
            Mut numChoices = new Mut(0);
            Queue<BranchingChoice> branchingQueue = getInsertionCost(i, numChoices);
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

    int getUnassignedRequest() {
        for (int i=0 ; i<numRequests; i++) {
            if (!customersLeft[i].value()) continue;
            insertionObjChange[i] = new HashMap<Integer, HashMap<Integer, HashMap<Integer, Integer>>>();
            for (int v=0;v<numVehicles;v++) {
                setInsertionCost(i, v);
            }
        }
        return getUnassignedMinVehicleMinInsertionPointsRequest();
    }

    int[][] getInsertionPoints(int request) {
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

    void branchRequestPoint(int request, int[] point) {
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

    void updateCapacityLeftInRoute(int v, int start) {
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

    boolean insertVertexIntoRoute(int i, int j) {
        if (!tryPost(new LessOrEqual(getArrivalTime(i, j), servingTime[j])) ||
                !tryPost(new LessOrEqual(getArrivalTime(pred[j].value(), i), servingTime[i])))
            return false;
        succ[i].setValue(j);
        pred[i].setValue(pred[j].value());
        succ[pred[i].value()].setValue(i);
        pred[j].setValue(i);
        return true;
    }

    BranchingChoice getBestBranchingDecision(int request) {
        var branchingQueue: ArrayBuffer[BranchingChoice] = ArrayBuffer[BranchingChoice]()
        val cvv = getCriticalVertex(request)
        val ncv = getCorrespondingVertex(cvv)
        var bestCvi = -1
        var bestNcvi = -1
        var bestChange = Int.MinValue
        for (v <- 0 until numVehicles) {
            if (insertionObjChange(request).contains(v) && servingVehicle(request).hasValue(v)) {
                for (cvi <- insertionObjChange(request)(v).keys) {
                    for (ncvi <- insertionObjChange(request)(v)(cvi).keys) {
                        if (insertionObjChange(request)(v)(cvi)(ncvi) > bestChange) {
                            bestChange = insertionObjChange(request)(v)(cvi)(ncvi)
                                    bestCvi = cvi
                            bestNcvi = ncvi
                            branchingQueue.clear()
                            branchingQueue += BranchingChoice(request, cvi, ncvi, bestChange, servingVehicle(cvi).value)
                        }
            else if (insertionObjChange(request)(v)(cvi)(ncvi) == bestChange) {
                            branchingQueue += BranchingChoice(request, cvi, ncvi, bestChange, servingVehicle(cvi).value)
                        }
                    }
                }
            }
        }
        val branchingQ: Array[BranchingChoice] = branchingQueue.toArray
        branchingQ(scala.util.Random.nextInt(branchingQ.length))
    }
}

class Mut {
    int value;

    public Mut(int i) {
        value = i;
    }
}

    /*

    def getInsertionCost(request: Int, numChoices: Mut[Int]): scala.collection.mutable.Queue[BranchingChoice] = {
        val branchingQueue: scala.collection.mutable.Queue[BranchingChoice] = scala.collection.mutable.Queue.empty[BranchingChoice]
        numChoices.value = 0
        val cvv = getCriticalVertex(request)
        val ncv = getCorrespondingVertex(cvv)
        var bestCvi = -1
        var bestNcvi = -1
        var bestChange = Int.MinValue

        for (v <- 0 until numVehicles) {
            if (insertionObjChange(request).contains(v) && servingVehicle(request).hasValue(v)) {
                for (cvi <- insertionObjChange(request)(v).keys) {
                    for (ncvi <- insertionObjChange(request)(v)(cvi).keys) {
                        numChoices.value += 1
                        if (insertionObjChange(request)(v)(cvi)(ncvi) > bestChange) {
                            bestChange = insertionObjChange(request)(v)(cvi)(ncvi)
                                    bestCvi = cvi
                            bestNcvi = ncvi
                            branchingQueue.clear()
                            branchingQueue.enqueue(BranchingChoice(request, cvi, ncvi, bestChange, servingVehicle(cvi).value))
                        }
            else if (insertionObjChange(request)(v)(cvi)(ncvi) == bestChange) {
                            branchingQueue.enqueue(BranchingChoice(request, cvi, ncvi, bestChange, servingVehicle(cvi).value))
                        }
                    }
                }
            }
        }
        branchingQueue
    }

    def setInsertionCost(request: Int, v: Int): Unit = {
        val cvv = getCriticalVertex(request)
        val ncv = getCorrespondingVertex(cvv)
        val begin = getBeginDepot(v)
        val end = getEndDepot(v)
        var start = succ(begin).value
        while(start != end){
            if(getArrivalTimeValue(pred(start).value, cvv, true) <= servingTime(cvv).max && getArrivalTimeValue(cvv, start, true) <= servingTime(start).max){
                setBestServingTimeFail(request, v, start)
            }
            start = succ(start).value
        }
        if(getArrivalTimeValue(pred(start).value, cvv, true) <= servingTime(cvv).max && getArrivalTimeValue(cvv, start, true) <= servingTime(start).max){
            setBestServingTimeFail(request, v, start)
        }
    }

    def setBestServingTimeFail(request: Int, v: Int, start: Int): Unit = {
        val begin = getBeginDepot(v)
        val end = getEndDepot(v)
        val cvv = getCriticalVertex(request)
        val ncv = getCorrespondingVertex(cvv)
        val cvvMinServingTime = math.max(getArrivalTimeValue(pred(start).value, cvv, true), timeWindowStarts(cvv))
        val cvvMaxServingTime = math.min(servingTime(start).max - getTravelTime(start, cvv) - servingDuration(cvv), timeWindowEnds(cvv))
        if (cvvMaxServingTime < cvvMinServingTime) {
            return
        }
        var changeCvv = servingTime(start).max - (getArrivalTimeValue(pred(start).value, cvv, true) + servingDuration(cvv) + getTravelTime(cvv, start))
        if (changeCvv < 0) {
            return
        }
        changeCvv -= 80 * (getTravelTime(cvv, start) + getTravelTime(pred(start).value, cvv) - getTravelTime(pred(start).value, start))
        var changeNcv = 0
        succ(cvv).setValue(start)
        pred(cvv).setValue(pred(start).value)
        succ(pred(cvv).value).setValue(cvv)
        pred(start).setValue(cvv)
        if (isPickup(cvv)) {
            var index = start
            var p = cvv
            var minRideTime = maxRideTime
            var done = false
            while (index != begin && !done) {
                if (p == cvv) {
                    minRideTime = getTravelTime(cvv, ncv)
                }
                else {
                    minRideTime = servingTime(p).min + servingDuration(p) + getTravelTime(p, ncv) - (cvvMaxServingTime + servingDuration(cvv))
                }
                if (minRideTime > maxRideTime) {
                    done = true
                }
                val ncvMinServingTime = math.max(getArrivalTimeValue(p, ncv, true), timeWindowStarts(ncv))
                val ncvMaxServingTime = math.min(servingTime(index).max - getTravelTime(index, ncv) - servingDuration(ncv), timeWindowEnds(ncv))

                changeNcv = servingTime(index).max - (getArrivalTimeValue(pred(index).value, ncv, true) + servingDuration(ncv) + getTravelTime(ncv, index))



                if (ncvMaxServingTime >= ncvMinServingTime && changeNcv >= 0) {
                    changeNcv -= 80 * (getTravelTime(ncv, index) + getTravelTime(pred(index).value, ncv) - getTravelTime(pred(index).value, index))
                    addToInsertionObjChange(request, v, start, index, changeCvv + changeNcv)
                    if (capacityLeftInRoute(index).value < vertexLoadChange(cvv)) {
                        done = true
                    }
                }
                p = index
                index = succ(index).value
            }
        }
        else {
            var p = pred(cvv).value
            var index = cvv
            var minRideTime = maxRideTime
            var done = false
            while (index != begin && capacityLeftInRoute(index).value >= vertexLoadChange(ncv) && !done) {
                if (index == cvv) {
                    minRideTime = getTravelTime(ncv, cvv)
                }
                else {
                    minRideTime = cvvMinServingTime - (getArrivalTimeValue(p, ncv, false) + servingDuration(ncv))
                }
                if (minRideTime > maxRideTime) {
                    done = true
                }
                val ncvMinServingTime = math.max(getArrivalTimeValue(p, ncv, true), timeWindowStarts(ncv))
                val ncvMaxServingTime = math.min(servingTime(index).max - getTravelTime(index, ncv) - servingDuration(ncv), timeWindowEnds(ncv))
                changeNcv = servingTime(index).max - (getArrivalTimeValue(pred(index).value, ncv, true) + servingDuration(ncv) + getTravelTime(ncv, index))
                if(ncvMaxServingTime >= ncvMinServingTime && changeNcv >= 0){
                    changeNcv -= 80 * (getTravelTime(ncv, index) + getTravelTime(pred(index).value, ncv) - getTravelTime(pred(index).value, index))
                    addToInsertionObjChange(request, v, start, index, changeCvv + changeNcv)
                }
                index = p
                p = pred(p).value
            }
        }

        succ(pred(cvv).value).setValue(succ(cvv).value)
        pred(succ(cvv).value).setValue(pred(cvv).value)
        succ(cvv).setValue(cvv)
        pred(cvv).setValue(cvv)
    }

    def addToInsertionObjChange(request: Int, v: Int, cvi: Int, ncvi: Int, change: Int): Unit = {
        if (!insertionObjChange(request).contains(v)) {
            insertionObjChange(request)(v) = HashMap[Int, HashMap[Int, Int]]()
        }
        if (!insertionObjChange(request)(v).contains(cvi)) {
            insertionObjChange(request)(v)(cvi) = HashMap[Int, Int]()
        }
        insertionObjChange(request)(v)(cvi)(ncvi) = change
    }


    def selectRelaxedCustomers(s: DarpSol, numCustomersToRelax: Int): Set[Int] = {
        val customers = Array.tabulate(numRequests)(i => i)
        var relaxEnd = 0
        while(relaxEnd < numCustomersToRelax && relaxEnd < customers.length){
            val toRelax = relaxEnd + Random.nextInt(numRequests-relaxEnd)
            val cRelaxed = customers(toRelax)
            customers(toRelax) = customers(relaxEnd)
            customers(relaxEnd) = cRelaxed
            relaxEnd += 1
        }
        customers.take(relaxEnd).toSet
    }

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

    def clearCustomerLeft(): Unit = {
        if (customersLeft.nonEmpty)
            for (i <- customersLeft) customersLeft.remove(i)
    }

    // Getter
    def getVehicleOfDepot(int i) { return if (beginDepotRange.contains(i)) i - 2 * numRequests else i - (2 * numRequests + numVehicles)

    def getBeginDepot(int i) { return 2 * numRequests + i

    def getEndDepot(int i) { return numVars - numVehicles + i

    def isBeginDepot(i: Int): Boolean = i >= 2 * numRequests && i < numVars - numVehicles

    def isEndDepot(i: Int): Boolean = i >= numVars - numVehicles && i < numVars

    def getArrivalTime(vertex: Int, successor: Int): CPIntVar = {
        val dist = travelTimeMatrix(vertex)(successor)
        if (isBeginDepot(vertex)) servingTime(vertex) + dist else servingTime(vertex) + servingDuration(vertex) + dist
    }

    def getCorrespondingRequest(int i) { return if (i < numRequests) i else i - numRequests

    def getCorrespondingPickup(int i) { return i - numRequests

    def getCorrespondingDelivery(int i) { return i + numRequests

    def getCriticalVertex(request: Int): Int = if (timeWindowStarts(request) > 0 || timeWindowEnds(request) < timeHorizon) request else request + numRequests

    def getCorrespondingVertex(int i) { return if (isPickup(i)) getCorrespondingDelivery(i) else getCorrespondingPickup(i)

    def isPickup(i: Int): Boolean = i < numRequests

    def isDelivery(i: Int): Boolean = i >= numRequests && i < 2 * numRequests

    def isCustomerVertex(i: Int): Boolean = i < 2 * numRequests

    def isCriticalVertex(vertex: Int): Boolean = timeWindowStarts(vertex) > 0 || timeWindowEnds(vertex) < timeHorizon

    def getTravelTimeMatrix: Array[Array[Int]] = travelTimeMatrix

    def getTravelTime(i: Int, j: Int): Int = travelTimeMatrix(i)(j)

    def getNumVehicles: Int = numVehicles

    def getBestSolution: DarpSol = bestSolution.get

    def getArrivalTimeValue(vertex: Int, successor: Int, getMin: Boolean): Int = {
        if (isBeginDepot(vertex)) {
            if (getMin) return servingTime(vertex).min + travelTimeMatrix(vertex)(successor) else return servingTime(vertex).max + travelTimeMatrix(vertex)(successor)
        }
        if (getMin)  servingTime(vertex).min + servingDuration(vertex) + travelTimeMatrix(vertex)(successor) else  servingTime(vertex).max + servingDuration(vertex) + travelTimeMatrix(vertex)(successor)
    }

    def tryPost(c: Constraint): Boolean = {
        try {
            solver.post(c)
        } catch {
            case _: Inconsistency => {
                return false
            }
        }
        true
    }

    def getDistanceObjective: Int = {
        var routeLength = 0
        for (int v=0;v<numVehicles;v++) {
            val begin = getBeginDepot(v)
            val end = getEndDepot(v)
            var i = begin
            while (i != end) {
                routeLength += travelTimeMatrix(i)(succ(i).value)
                i = succ(i).value
            }
        }
        routeLength
    }

    def getBestSolutionObjective: Int = bestSolutionObjective

    def isPositive(array: Array[ReversibleInt]): Boolean = {
        for (i <- array.indices) {
            if (array(i).value < 0) {
                return false
            }
        }
        true
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
