package minicp.examples;

import minicp.engine.core.Solver;

import static minicp.cp.Factory.makeSolver;
import minicp.examples.DARPParser.*;
import minicp.examples.DARPDataModel.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class DARPtest {
    static DARPInstance instance;
    int searchTime;
    //    solPath: Option[String],
//    logPath: Option[String],
    boolean firstSolOnly = true;

    public static void main(String[] args) {
        Solver cp = makeSolver();
        instance = DARPParser.parseInstance("C:\\Users\\Utilisateur\\Documents\\UNIF2020\\TFE\\DARP RK\\DARP\\Cordeau\\a2-16.txt");
        System.out.println(instance.nVehicles);
    }
    /*
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


    // finds all possible insertions of request in route v, with critical vertex before cvSucc
    // and compute their cost.
    static void setBestServingTimeFail2(int request, int v, int cvSucc) {
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
                getArrivalTime(cvPred, cv).min());
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
                    getArrivalTime(ncvPred, ncv).min());
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
                getArrivalTime(ncvPred, ncv).min());
        int ncvMaxServingTime, cvMinServingTime;
        if (ncvSucc == cv) { // edge case: we insert drop right after pickup
            cvMinServingTime = Math.max(servingTime[cv].min(),
                    ncvMinServingTime + servingDuration[ncv] + dist[ncv][cv]);
            ncvMaxServingTime = Math.min(servingTime[ncv].max(),
                    cvMaxServingTime - dist[ncv][cv] - servingDuration[ncv]);
        } else {
            cvMinServingTime = Math.max(servingTime[cv].min(),
                    getArrivalTime(cvPred, cv).min());
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
    */
}
