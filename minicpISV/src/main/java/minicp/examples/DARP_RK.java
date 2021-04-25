package minicp.examples;

import minicp.engine.core.Solver;

import static minicp.cp.Factory.makeSolver;
import minicp.examples.DARPParser.*;
import minicp.examples.DARPDataModel.*;

/**
 * @author Quentin Delmelle qdelmelle@gmail.com
 * transtlated from the scala work of
 * Roger Kameugne rkameugne@gmail.com
 * and Charles Thomas cftmthomas@gmail.com
 */

public class DARPModelVH {
    DARPInstance instance;
    int searchTime;
//    solPath: Option[String],
//    logPath: Option[String],
    boolean firstSolOnly = true;

    public void main(){
        Solver cp = makeSolver();
        instance =;
    }
}

    solver.silent = true

    var bestSolutionObjective = Int.MaxValue
    var bestSolution: Option[DarpSol] = None
    var currentSolution: Option[DarpSol] = None
    var currentSolutionObjective = Int.MaxValue
    var totalNumFails = 0
    var constraintsPosted: Boolean = false

    println(firstSolOnly)


    /* Data */

    // Parameters of the problem
    val numVars = instance.nSites
    val numRequests = instance.nRequests
    val vehicleCapacity = instance.vehicles.head.capacity
    val numVehicles = instance.nVehicles
    val maxRideTime = instance.requests.head.maxRideTime
    val timeHorizon = instance.sites.map(_.winEnd).max

    // Ranges
    val numVarsRange = 0 until numVars
    val numRequestsRange = 0 until numRequests
    val beginDepotRange = 2 * numRequests until numVars - numVehicles
    val endDepotRange = numVars - numVehicles until numVars
    val numVehiclesRange = 0 until numVehicles
    val customerVertexRange = 0 until 2 * numRequests
    val succRange = 0 until numVars - numVehicles
    val predRange = customerVertexRange ++ endDepotRange

    // Variables necessary for posting constraints
    val travelTimeMatrix = instance.distances
    val vertexLoadChange = instance.sites.map(_.load)
    val timeWindowStarts = instance.sites.map(_.winStart)
    val timeWindowEnds = instance.sites.map(_.winEnd)

    // Variables related to modelling
    var servingTime: Array[CPIntVar] = new Array[CPIntVar](numVars)
    val servingVehicle: Array[CPIntVar] = new Array[CPIntVar](numVars)
    val succ: Array[ReversibleInt] = new Array[ReversibleInt](numVars)
    val pred: Array[ReversibleInt] = new Array[ReversibleInt](numVars)

    val capacityLeftInRoute: Array[ReversibleInt] = Array.tabulate(numVars)(_ => new ReversibleInt(solver, vehicleCapacity))
    val servingDuration = instance.sites.map(_.service)

    val insertionObjChange: Array[scala.collection.mutable.HashMap[Int, scala.collection.mutable.HashMap[Int, scala.collection.mutable.HashMap[Int, Int]]]] = new Array[scala.collection.mutable.HashMap[Int, scala.collection.mutable.HashMap[Int, scala.collection.mutable.HashMap[Int, Int]]]](numRequests)
    val customersLeft: ReversibleSet = new ReversibleSet(solver)

    def initCpVars(): Unit = {
        for (i <- numVarsRange) {
            if (customerVertexRange.contains(i)) {
                succ(i) = new ReversibleInt(solver, i, numVars)
                pred(i) = new ReversibleInt(solver, i, numVars)
                servingVehicle(i) = CPIntVar(0, numVehicles)
            } else {
                if (beginDepotRange.contains(i)) {
                    succ(i) = new ReversibleInt(solver, getEndDepot(getVehicleOfDepot(i)))
                    pred(i) = new ReversibleInt(solver, succ(i).value)
                    servingVehicle(i) = CPIntVar(getVehicleOfDepot(i))
                } else {
                    succ(i) = new ReversibleInt(solver, getBeginDepot(getVehicleOfDepot(i)))
                    pred(i) = new ReversibleInt(solver, succ(i).value)
                    servingVehicle(i) = CPIntVar(getVehicleOfDepot(i))
                }
            }
            servingTime(i) = CPIntVar(timeWindowStarts(i), timeWindowEnds(i))
        }
    }

    def postPrecedence(): Unit = {
        for (i <- numRequestsRange) {
            add(servingTime(i) <= (servingTime(numRequests + i) - travelTimeMatrix(i)(i + numRequests) - servingDuration(i)))
        }
        for (i <- numVehiclesRange) {
            add(servingTime(getEndDepot(i)) - servingTime(getBeginDepot(i)) <= timeHorizon)
        }
    }

    def postRideTime(): Unit = {
        for (i <- numRequestsRange) {
            add(servingTime(numRequests + i) - (servingTime(i) + servingDuration(i)) <= maxRideTime)
        }
    }


    def postConstraints(): Unit = {
        if (constraintsPosted) {
            return
        }
        constraintsPosted = true
        for (i <- numRequestsRange) {
            add(servingVehicle(i) === servingVehicle(i + numRequests))
        }
        postPrecedence()
        postRideTime()
    }

    def postCumulativeConstraint(): Unit = {
        val travelStart: Array[CPIntVar] = Array.tabulate(numRequests)(i => servingTime(i))
        val travelEnd: Array[CPIntVar] = Array.tabulate(numRequests)(i => servingTime(i+numRequests))
        val travelDuration: Array[CPIntVar] = Array.tabulate(numRequests)(i => travelEnd(i) - travelStart(i))
        val travelLoad = Array.tabulate(numRequests)(i => CPIntVar(vertexLoadChange(i)))
        val travelVehicle = Array.tabulate(numRequests)(i => servingVehicle(i))
        for(i <- numRequestsRange) add(travelDuration(i) > 0)
        for(v <- numVehiclesRange){
            //Adding cumulative constraint:
            val capVar = CPIntVar(vehicleCapacity)
            add(maxCumulativeResource(travelStart, travelDuration, travelEnd, travelLoad, travelVehicle, capVar, v))
        }
    }

    def postUnaryConstraint(): Unit = {
        val start: Array[CPIntVar] = Array.tabulate(2*numRequests)(i => servingTime(i))
        val duration: Array[CPIntVar] = Array.tabulate(2*numRequests)(i => CPIntVar(servingDuration(i)))
        val end: Array[CPIntVar] = Array.tabulate(2*numRequests)(i => start(i) + duration(i))

        val resource = Array.tabulate(2*numRequests)(i => servingVehicle(i))
        for(v <- numVehiclesRange){
            //Adding unary constraint:
            add(unaryResource(start, duration, end, resource, v))
        }
    }



    // parameters
    val gamma = 200
    val tau = 1000
    val maxSize = numRequests / 2
    val range = 4
    val numIter = 300
    val d = 0.07.toFloat
    val maxTime = 300
    var remainTime = maxTime * 1000L


    /**
     * Start Solve the problem
     */
    // initialise variable
    initCpVars()

    //Add request to set of unassigned request "customersLeft"
  for (i <- numRequestsRange) {
        customersLeft.add(i)
    }

    // post constraints
    postConstraints()
    //postUnaryConstraint()
    //postCumulativeConstraint()

    // On solution
    onSolution {
        println("Total route cost: " + getDistanceObjective/100.0)
        println("best solution found:= " + bestSolutionObjective/100.0)
        println("remainingTime := "+remainTime)
        val succP = Array.tabulate(numVars)(i => succ(i).value)
        val predP = Array.tabulate(numVars)(i => pred(i).value)
        val servingVP = Array.tabulate(numVars)(i => servingVehicle(i).value)
        val minServintgTime = Array.tabulate(numVars)(i => servingTime(i).min/100.0)
        val maxServintgTime = Array.tabulate(numVars)(i => servingTime(i).max/100.0)
//    println("maxRideTime: " + maxRideTime/100.0)
//    println("serving times:")
//    println(numVarsRange.map(i => i + " (" + minServintgTime(i) + ":" + maxServintgTime(i) + ")").mkString("\n"))
        val rand = scala.util.Random.nextFloat()
        if (getDistanceObjective < currentSolutionObjective || rand < d) {
            currentSolution = Some(new DarpSol(succP, predP, servingVP, getDistanceObjective/100.0, minServintgTime, maxServintgTime))
            currentSolutionObjective = getDistanceObjective
            if (currentSolutionObjective < bestSolutionObjective) {
                bestSolution = Some(new DarpSol(succP, predP, servingVP, getDistanceObjective/100.0, minServintgTime, maxServintgTime))
                bestSolutionObjective = currentSolutionObjective
            }
        }
        println("-------------------------")
    }

    def relax(nRelax: Int): Unit = {
        val relaxedCustomers = selectRelaxedCustomers(currentSolution.get, nRelax)
        val solSucc = Array.tabulate(numVars)(i => currentSolution.get.succ(i))
        val solServingVehicle = Array.tabulate(numVars)(i => currentSolution.get.servingVehicle(i))
        clearCustomerLeft()
        for (r <- relaxedCustomers) {
            succ(r).setValue(r)
            pred(r).setValue(r)
            succ(r + numRequests).setValue(r + numRequests)
            pred(r + numRequests).setValue(r + numRequests)
            customersLeft.add(r)
        }
        for (v <- numVehiclesRange) {
            val begin = getBeginDepot(v)
            val end = getEndDepot(v)
            var current = begin
            var prev = -1
            while (current != end) {
                val currentRequest = getCorrespondingRequest(current)
                if(!relaxedCustomers.contains(currentRequest)){
                    if(prev != -1){
                        succ(prev).setValue(current)
                        pred(current).setValue(prev)
                        add(getArrivalTime(prev, current) <= servingTime(current))
                        add(servingVehicle(current) === solServingVehicle(current))
                    }
                    prev = current
                }
                current = solSucc(current)
            }
            if(prev != -1){
                succ(prev).setValue(end)
                pred(end).setValue(prev)
                add(getArrivalTime(prev, end) <= servingTime(end))
            }
            updateCapacityLeftInRoute(v, -1)
        }
    }

    // lns
    def lns(): Unit = {
        var i = 2
        while(remainTime > 0 && i <= maxSize - range) {
            var j = 0
            while(remainTime > 0 && j <= range) {
                var k = 1
                while(remainTime > 0 && k <= numIter) {
                    val stats = startSubjectTo(1, Int.MaxValue, (remainTime/1000.0).round.toInt){
                        relax(i + j)
                    }
                    k += 1
                    remainTime -= stats.time
                    totalNumFails += stats.nFails
                    println(stats)
                }
                j += 1
            }
            i += 1
        }
    }

    var solFound = false

    search {
        if (customersLeft.isEmpty) {
            noAlternative
        } else {
            val request = getUnassignedRequest
            val point= getInsertionPoints(request)
            branchAll(point)(p => branchRequestPoint(request, (p._1, p._2, p._3, p._4)))
        }
    }

//  println(start(nSols = Int.MaxValue, failureLimit = Int.MaxValue, timeLimit = (remainTime/1000.0).round.toInt))

  while(!solFound && remainTime > 0){
        val stats = start(1, Math.max(gamma * numVehicles, tau), (remainTime/1000.0).round.toInt)
        if(stats.nSols == 1) solFound = true
        println(stats)
        remainTime -= (stats.time/1000.0).round.toInt
    }

  if(!firstSolOnly) lns()

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


    //-----------------------------------------


    def getUnassignedMinVehicleMinInsertionPointsRequest: Int = {
        var minChoices = Int.MaxValue
        var minVehicles = numVehicles + 1
        val bQueueBuffer: ArrayBuffer[BranchingChoice] = ArrayBuffer[BranchingChoice]()
        val bestChange = Int.MinValue
        for (i <- customersLeft) {
            val tempChange = Int.MinValue
            val numChoices: Mut[Int] = new Mut[Int](0)
            val branchingQueue: scala.collection.mutable.Queue[BranchingChoice] = getInsertionCost(i, numChoices)
            if (servingVehicle(i).size < minVehicles) {
                minVehicles = servingVehicle(i).size
                bQueueBuffer.clear()
                for (j <- branchingQueue) {
                    bQueueBuffer += j
                }
            }
            else if (servingVehicle(i).size == minVehicles && numChoices.value < minChoices) {
                minChoices = numChoices.value
                bQueueBuffer.clear()
                for (j <- branchingQueue) {
                    bQueueBuffer += j
                }
            }
            else if (servingVehicle(i).size == minVehicles && numChoices.value == minChoices && tempChange >= bestChange) {
                for (j <- branchingQueue) {
                    bQueueBuffer += j
                }
            }
        }
        val bQueue: Array[BranchingChoice] = bQueueBuffer.toArray
        val currentRequest = bQueue(scala.util.Random.nextInt(bQueue.length)).request
        val bc: BranchingChoice = getBestBranchingDecision(currentRequest)
        bc.request
    }


    def getUnassignedRequest: Int = {
        for (i <- customersLeft) {
            insertionObjChange(i) = new scala.collection.mutable.HashMap[Int, scala.collection.mutable.HashMap[Int, scala.collection.mutable.HashMap[Int, Int]]]()
            for (v <- numVehiclesRange) {
                setInsertionCost(i, v)
            }
        }
        getUnassignedMinVehicleMinInsertionPointsRequest
    }

    def getInsertionPoints(request: Int): Array[(Int, Int, Int, Int)] = {
        val vehiclePointsChangeBuffer: ArrayBuffer[(Int, Int, Int, Int)] = ArrayBuffer[(Int, Int, Int, Int)]()
        for (v <- insertionObjChange(request).keys) {
            for (cv <- insertionObjChange(request)(v).keys) {
                for (ncv <- insertionObjChange(request)(v)(cv).keys) {
                    vehiclePointsChangeBuffer += ((v, cv, ncv, insertionObjChange(request)(v)(cv)(ncv)))
                }
            }
        }
        val vehiclePointsChange: Array[(Int, Int, Int, Int)] = vehiclePointsChangeBuffer.toArray.sortBy(-_._4)
        vehiclePointsChange
    }

    def branchRequestPoint(request: Int, point: (Int, Int, Int, Int)): Unit = {
        val (vehicle, cvSucc, ncvSucc, change): (Int, Int, Int, Int) = point
        val cvv = getCriticalVertex(request)
        val ncv = getCorrespondingVertex(cvv)


        if (!tryPost(servingVehicle(request) === servingVehicle(cvSucc))) {
            throw Inconsistency
        }
        if (!insertVertexIntoRoute(cvv, cvSucc)) {
            throw Inconsistency
        }
        if (!insertVertexIntoRoute(ncv, ncvSucc)) {
            throw Inconsistency
        }
        updateCapacityLeftInRoute(servingVehicle(request).value, request)
        if (!isPositive(capacityLeftInRoute)) {
            throw Inconsistency
        }
        customersLeft.remove(request)
        val v = servingVehicle(cvSucc).value
        for (i <- customersLeft) {
            if (insertionObjChange(i).contains(v)) {
                insertionObjChange(i).remove(v)
            }
            setInsertionCost(i, v)
            val insertionPoint = getInsertionPoints(i)
            if(insertionPoint.isEmpty){
                throw Inconsistency
            }
        }

    }




    def updateCapacityLeftInRoute(v: Int, start: Int): Unit = {
        val begin = getBeginDepot(v)
        val end = getEndDepot(v)
        var star = start
        if (star == -1) {
            star = succ(begin).value
        }
        var capacity = capacityLeftInRoute(pred(star).value).value
        while (star != end) {
            capacity -= vertexLoadChange(star)
            capacityLeftInRoute(star).setValue(capacity)
            star = succ(star).value
        }
    }

    def insertVertexIntoRoute(i: Int, j: Int): Boolean = {
        if (!tryPost(getArrivalTime(i, j) <= servingTime(j)) || !tryPost(getArrivalTime(pred(j).value, i) <= servingTime(i))) return false
        succ(i).setValue(j)
        pred(i).setValue(pred(j).value)
        succ(pred(i).value).setValue(i)
        pred(j).setValue(i)
        true
    }

    def getBestBranchingDecision(request: Int): BranchingChoice = {
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
            insertionObjChange(request)(v) = scala.collection.mutable.HashMap[Int, scala.collection.mutable.HashMap[Int, Int]]()
        }
        if (!insertionObjChange(request)(v).contains(cvi)) {
            insertionObjChange(request)(v)(cvi) = scala.collection.mutable.HashMap[Int, Int]()
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
        for (v <- numVehiclesRange) {
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
        for (v <- numVehiclesRange) {
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
        for (v <- numVehiclesRange) {
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
        for (v <- numVehiclesRange) {
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
    def getVehicleOfDepot(i: Int): Int = if (beginDepotRange.contains(i)) i - 2 * numRequests else i - (2 * numRequests + numVehicles)

    def getBeginDepot(i: Int): Int = 2 * numRequests + i

    def getEndDepot(i: Int): Int = numVars - numVehicles + i

    def isBeginDepot(i: Int): Boolean = i >= 2 * numRequests && i < numVars - numVehicles

    def isEndDepot(i: Int): Boolean = i >= numVars - numVehicles && i < numVars

    def getArrivalTime(vertex: Int, successor: Int): CPIntVar = {
        val dist = travelTimeMatrix(vertex)(successor)
        if (isBeginDepot(vertex)) servingTime(vertex) + dist else servingTime(vertex) + servingDuration(vertex) + dist
    }

    def getCorrespondingRequest(i: Int): Int = if (i < numRequests) i else i - numRequests

    def getCorrespondingPickup(i: Int): Int = i - numRequests

    def getCorrespondingDelivery(i: Int): Int = i + numRequests

    def getCriticalVertex(request: Int): Int = if (timeWindowStarts(request) > 0 || timeWindowEnds(request) < timeHorizon) request else request + numRequests

    def getCorrespondingVertex(i: Int): Int = if (isPickup(i)) getCorrespondingDelivery(i) else getCorrespondingPickup(i)

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
        for (v <- numVehiclesRange) {
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

  case class BranchingChoice(request: Int, cvSucc: Int, ncvSucc: Int, change: Int, vehicle: Int)

  case class Mut[A](var value: A) {}
}

class DarpSol(var succ: Array[Int], var pred: Array[Int], var servingVehicle: Array[Int], var cost: Double, var minServingTime: Array[Double], var maxServingTime: Array[Double])

