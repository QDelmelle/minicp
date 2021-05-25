package minicp.examples;

import minicp.engine.constraints.*;
import minicp.engine.core.IntVar;
import minicp.engine.core.Solver;
import minicp.search.DFSearch;
import minicp.search.Objective;
import minicp.util.io.InputReader;
import static minicp.cp.Factory.*;
import static minicp.cp.BranchingScheme.*;
import java.lang.management.ManagementFactory;

import java.util.*;
import java.util.stream.IntStream;

public class DialARide {


    public static DialARideSolution solve(int nVehicles, int maxRouteDuration, int vehicleCapacity,
                                          int maxRideTime, ArrayList<RideStop> pickupRideStops, ArrayList<RideStop> dropRideStops,
                                          RideStop depot) {
        // TODO
        // Given a series of dial-a-ride request made by single persons (for request i, pickupRideStops[i] gives the spot
        // where the person wants to be taken, and dropRideStops[i] the spot where (s)he would like to be dropped),
        // minimize the total ride time of all the vehicles.
        // You have nVehicles vehicles, each of them can take at most vehicleCapacity person inside at any time.
        // The maximum time a single person can remain in the vehicle is maxRideTime, and the maximum time a single
        // vehicle can be on the road for a single day is maxRouteDuration.
        // all vehicles start at the depot, and end their day at the depot.
        // Each ride stop must be reached before a given time (window_end) by a vehicle.
        // use distance() to compute the distance between two points.

        // CONSTANTS
        ArrayList<RideStop> allRideStops = pickupRideStops;
        int nRequests = pickupRideStops.size();
        allRideStops.addAll(dropRideStops);
        for(int i=0; i<nVehicles; i++) // this looks stupid but it will be useful later
            allRideStops.add(depot);
        int nNodes = nRequests*2 + nVehicles;
        int[][] disMatrix = makeDistanceMatrix(allRideStops);
        Solver cp = makeSolver(false);

        // nodes 0 -> nRequests-1 = pickups, nodes nRequests -> 2*nRequests-1 = drops
        int[] drop = new int[nRequests];
        for(int i=0; i<nRequests; i++){
            drop[i] = nRequests+i;
        }
        // we duplicate the depot for each vehicle
        // nodes 2*nRequests -> 2*nRequests+nVehicles-1 = depot nodes
        int[] depots = new int[nVehicles];
        for(int i=0; i<nVehicles; i++){
            depots[i] = 2*nRequests+i;
        }

        // VARIABLES

        IntVar[] next = makeIntVarArray(cp, nNodes, nNodes); // next[x] = successor of node x
        IntVar[] prev = makeIntVarArray(cp, nNodes, nNodes); // prev[x] = predecessor of node x
        IntVar[] T = makeIntVarArray(cp, nNodes, maxRouteDuration*nVehicles);
        //T[x] = time at which the node x is reached in the TSP problem
        IntVar[] distPrev = makeIntVarArray(cp, nNodes, 1000); // distance to previous node
        IntVar[] capacity = makeIntVarArray(cp, nNodes, vehicleCapacity);
        //capacity[x] = number of people inside vehicle after x is reached
        IntVar[] chauffeur = makeIntVarArray(cp, nNodes, nVehicles);
        // chauffeur[x] = number of the vehicle who visits node x
        printDomSize(prev, next, T, capacity, chauffeur);

        // OBJECTIVE
        IntVar totalDist = sum(distPrev);
        Objective obj = cp.minimize(totalDist);

        //  CONSTRAINTS
        // circuit for TSP subproblem
        cp.post(new Circuit(next));
        // i = next[prev[i]] = prev[next[i]]
        for(int i=0; i<nNodes; i++) {
            IntVar uselessVar = makeIntVar(cp, i, i);
            cp.post(new Element1DVar(next, prev[i], uselessVar));
            cp.post(new Element1DVar(prev, next[i], uselessVar));
        }

        // capacity and time variation + obj binding
        // vehicle 0 starts at time 0
        T[depots[0]].assign(0);
        for (int i = 0; i < nNodes; i++) {
            if(i == depots[0]){
                cp.post(new Element1DVar(T, prev[i], sum(totalDist, minus(distPrev[i]))));
            }
            else {
                cp.post(new Element1DVar(T, prev[i], sum(T[i], minus(distPrev[i]))));
            }
            cp.post(new Element1DDomainConsistent(disMatrix[i], prev[i], distPrev[i]));
            // capacity[i] = capacity[prev[i]] + type of stop i
            cp.post(new Element1DVar(capacity, prev[i], minus(capacity[i], allRideStops.get(i).type)));
        }
        // maxRouteDuration constraint
        for(int i=0; i<nVehicles-1; i++){
            cp.post(lessOrEqual(T[depots[i+1]], plus(T[depots[i]], maxRouteDuration)));
        }
        cp.post(lessOrEqual(T[depots[0]], plus(T[depots[nVehicles-1]], maxRouteDuration)));
        // can't have a drop after a depot, or a pickup before
        for(int i=0; i<nVehicles; i++){
            for(int j=0; j<nRequests; j++){
                next[depots[i]].remove(drop[j]);
                prev[depots[i]].remove(j);
            }
        }
        // a vehicle returning to depot must be empty
        for(int i=0; i<nVehicles; i++){
            equal(capacity[depots[i]], 0);
        }
        // we decide that vehicle i starts at depots[i]
        for(int i=0; i<nVehicles; i++){
            chauffeur[depots[i]].assign(i);
        }
        cp.fixPoint();

        for(int i=0; i<nRequests; i++){

            // precedency constraint
            //cp.post(lessOrEqual(T[i], T[drop[i]]));

            // maxRideTime constraint
            //cp.post(lessOrEqual(T[drop[i]], plus(T[i], maxRideTime)));

            // chauffeur constraints
            // person i is picked up and dropped by the same vehicle
            cp.post(lessOrEqual(chauffeur[i],chauffeur[drop[i]]));
            cp.post(lessOrEqual(chauffeur[drop[i]],chauffeur[i]));
            // stops on the same route have the same vehicle
            cp.post(new Element1DVar(chauffeur, prev[i], chauffeur[i]));

            // time windows are respected : T[depots[chauffeur[i]]] + window of node i >= T[i]
            IntVar Z = makeIntVar(cp, 0, maxRouteDuration*nVehicles);
            cp.post(new Element1DVar(T,element(depots, chauffeur[i]), Z));
            cp.post(lessOrEqual(minus(T[i], allRideStops.get(i).window_end), Z));
        }

        // maxCapacity constraint -> implicitly taken care of

        printDomSize(prev, next, T, capacity, chauffeur);

        int[] visitNodes = new int[1];
        // LAUNCH SOLVER
        DFSearch dfs = makeDfs(cp, () -> {
            IntVar xs=null;
            int min = Integer.MAX_VALUE;
            int bestNext = -1;
            for (int i = 0; i < nNodes; i++) {
                if (!next[i].isBound()) {
                    int[] tmp = new int[next[i].size()];
                    next[i].fillArray(tmp);
                    int[] nextClosest = minRClosest(i, tmp, disMatrix);
                    int v = disMatrix[i][nextClosest[0]] - disMatrix[i][nextClosest[1]];
                    if (v < min) { // replace next.size to v for closest
                        bestNext = nextClosest[0];
                        xs = next[i];
                        min = v;    //next[i].size();
                    }
                }
            }
            if (xs == null)
                return EMPTY;
            else {
                final int v = bestNext;
                final IntVar xt = xs;
                visitNodes[0]++;
                return branch(() -> equal(xt,v),
                        () -> notEqual(xt, v));
            }
        });

        int[] BestSolution = IntStream.range(0, nNodes).toArray();
        dfs.onSolution(() -> {
            // Update the current best solution
            for (int i = 0; i < nNodes; i++) {
                BestSolution[i] = next[i].min();
            }
            System.out.println("objective:" + totalDist);
        });

        //dfs.solve();

        int nRestarts = 300;
        int failureLimit = 100;
        Random rand = new java.util.Random(0);

        int[] percentage = new int[2];
        percentage[0] = 0; // for simulated annealing % = 50->5
        percentage[1] = 50;

        for (int i = 0; i < nRestarts; i++) {
            int w = i;
            if (i%100==0)
                System.out.println("restart number #"+w);
            // Record the state such that the fragment constraints can be cancelled
            dfs.optimizeSubjectTo(obj,statistics -> statistics.numberOfFailures() >= failureLimit,
                    () -> {
                        // Assign the fragment 10% of the variables randomly chosen
                        for (int j = 0; j < nNodes; j++) {
                            if (rand.nextInt(100) < ((percentage[1]-percentage[0])/((double)(nRestarts))*w+percentage[0])) { // kinda smulated annealing
                                equal(next[j],BestSolution[j]);
                            }
                        }
                    });
        }
        double endTime = ManagementFactory.getThreadMXBean().getCurrentThreadCpuTime()/1000000000.0;
        System.out.println(endTime+": "+visitNodes[0]/endTime);

        // create solution
        DialARideSolution sol = new DialARideSolution(nVehicles, pickupRideStops, dropRideStops, depot,
                vehicleCapacity, maxRideTime, maxRouteDuration);

        for(int i=0; i<nVehicles; i++){
            int node = BestSolution[depots[i]];
            boolean isPickup = allRideStops.get(node).type == 1;
            while(true){
                if(allRideStops.get(node).type == 0)
                    break;
                sol.addStop(i, node%(nRequests+1), isPickup);
                node = next[node].min();
                isPickup = allRideStops.get(node).type == 1;
            }
        }

        return sol;
    }

    /**
     *  a nice function to print the domain sizes of all the variables and see what the propagation did
     */
    public static void printDomSize(IntVar[] prev, IntVar[] next,
                                    IntVar[] T, IntVar[] capacity, IntVar[] chauffeur){
        int prevSize = 0, nextSize = 0, Tsize = 0, capSize = 0, chSize = 0;
        for(int i=0; i<prev.length; i++){
            prevSize += prev[i].size();
            nextSize += next[i].size();
            Tsize += T[i].size();
            capSize += capacity[i].size();
            chSize += chauffeur[i].size();
        }
        int total = chSize + nextSize + prevSize + capSize + Tsize;
        System.out.println("\n Domain sizes: ");
        System.out.println("prev: "+prevSize+", next: "+nextSize+
                ", T: "+Tsize+", capacity: "+capSize+", chauffeur: "+chSize);
        System.out.println("Total: "+total);
    }

    /**
     * Returns the distance between two ride stops
     */
    public static int distance(RideStop a, RideStop b) {
        return (int)(Math.sqrt((a.pos_x - b.pos_x) * (a.pos_x - b.pos_x) + (a.pos_y - b.pos_y) * (a.pos_y - b.pos_y)) * 100);
    }

    /**
     *
     * @param p
     * @param dom = domain of next[p]
     * @param distanceMatrix
     * @return the closest node to p according to distanceMatrix
     */
    private static int closest(int p, int[] dom,int[][] distanceMatrix){
        int closest = -1;
        int min = Integer.MAX_VALUE;
        for(int j=0;j<dom.length;j++){
            if(p==dom[j])
                continue;
            if(distanceMatrix[p][dom[j]]<min){
                min = distanceMatrix[p][dom[j]];
                closest = dom[j];
            }
        }
        return closest;
    }

    /**
     *
     * @param p
     * @param dom = domain of next[p]
     * @param distanceMatrix
     * @return the two closests nodes to p according to distanceMatrix
     */
    private static int[] minRClosest(int p, int[] dom,int[][] distanceMatrix){ // variable needs to be unbound
        int[] closests = new int[2]; // en 0 : plus proche // en 1 suivant
        closests[0] = -1;
        closests[1] = -1;
        int min = Integer.MAX_VALUE;
        int secondMin = Integer.MAX_VALUE;
        for(int j=0;j<dom.length;j++){
            if(p==dom[j])
                continue;
            if(distanceMatrix[p][dom[j]]<min){
                secondMin = min;
                min = distanceMatrix[p][dom[j]];
                closests[1] = closests[0];
                closests[0] = dom[j];
            }
            if(distanceMatrix[p][dom[j]]<secondMin) {
                secondMin = distanceMatrix[p][dom[j]];
                closests[1] = dom[j];
            }
        }
        return closests;
    }

    /**
     *  returns the distance matrix of all the RideStops
     */
    public static int[][] makeDistanceMatrix(List<RideStop> L){
        int[][] t = new int[L.size()][L.size()];
        for(int i=0; i<L.size(); i++){
            for(int j=0; j<L.size(); j++){
                t[i][j] = distance(L.get(i), L.get(j));
            }
        }
        return t;
    }

    /**
     * A solution. To create one, first do new DialARideSolution, then
     * add, for each vehicle, in order, the pickup/drops with addStop(vehicleIdx, rideIdx, isPickup), where
     * vehicleIdx is an integer beginning at 0 and ending at nVehicles - 1, rideIdx is the id of the ride you (partly)
     * fullfill with this stop (from 0 to pickupRideStops.size()-1) and isPickup a boolean indicate if you are beginning
     * or ending the ride. Do not add the last stop to the depot, it is implicit.
     *
     * You can check the validity of your solution with compute(), which returns the total distance, and raises an
     * exception if something is invalid.
     *
     * DO NOT MODIFY THIS CLASS.
     */
    public static class DialARideSolution {
        public ArrayList<Integer>[] stops;
        public ArrayList<RideStop> pickupRideStops;
        public ArrayList<RideStop> dropRideStops;
        public RideStop depot;
        public int capacity;
        public int maxRideTime;
        public int maxRouteDuration;

        public String toString() {
            StringBuilder b = new StringBuilder();
            b.append("Length: ");
            b.append(compute());
            b.append("\n");
            for(int i = 0; i < stops.length; i++) {
                b.append("- ");
                for(int s: stops[i]) {
                    if(s >= pickupRideStops.size()) {
                        b.append(s-pickupRideStops.size());
                        b.append("d, ");
                    }
                    else {
                        b.append(s);
                        b.append("p, ");
                    }
                }
                b.append("\n");
            }
            return b.toString();
        }

        public DialARideSolution(int nVehicles, ArrayList<RideStop> pickupRideStops, ArrayList<RideStop> dropRideStops,
                                 RideStop depot, int vehicleCapacity, int maxRideTime, int maxRouteDuration) {
            stops = new ArrayList[nVehicles];
            for(int i = 0; i < nVehicles; i++)
                stops[i] = new ArrayList<>();

            this.pickupRideStops = pickupRideStops;
            this.dropRideStops = dropRideStops;
            this.depot = depot;
            this.capacity = vehicleCapacity;
            this.maxRideTime = maxRideTime;
            this.maxRouteDuration = maxRouteDuration;
        }

        public void addStop(int vehicleId, int rideId, boolean isPickup) {
            stops[vehicleId].add(rideId + (isPickup ? 0 : pickupRideStops.size()));
        }

        public int compute() {
            int totalLength = 0;
            HashSet<Integer> seenRides = new HashSet<>();

            for(int vehicleId = 0; vehicleId < stops.length; vehicleId ++) {
                HashMap<Integer, Integer> inside = new HashMap<>();
                RideStop current = depot;
                int currentLength = 0;
                for(int next: stops[vehicleId]) {
                    RideStop nextStop;
                    if(next < pickupRideStops.size())
                        nextStop = pickupRideStops.get(next);
                    else
                        nextStop = dropRideStops.get(next - pickupRideStops.size());

                    currentLength += distance(current, nextStop);

                    if(next < pickupRideStops.size()) {
                        if(seenRides.contains(next))
                            throw new RuntimeException("Ride stop visited twice");
                        seenRides.add(next);
                        inside.put(next, currentLength);
                    }
                    else {
                        if(!inside.containsKey(next - pickupRideStops.size()))
                            throw new RuntimeException("Drop before pickup");
                        if(inside.get(next - pickupRideStops.size()) + maxRideTime < currentLength)
                            throw new RuntimeException("Ride time too long");
                        inside.remove(next - pickupRideStops.size());
                    }

                    if(currentLength > nextStop.window_end)
                        throw new RuntimeException("Ride stop visited too late");
                    if(inside.size() > capacity)
                        throw new RuntimeException("Above maximum capacity");

                    current = nextStop;
                }

                currentLength += distance(current, depot);

                if(inside.size() > 0)
                    throw new RuntimeException("Passenger never dropped");
                if(currentLength > maxRouteDuration)
                    throw new RuntimeException("Route too long");

                totalLength += currentLength;
            }

            if(seenRides.size() != pickupRideStops.size())
                throw new RuntimeException("Some rides never fulfilled");

            return totalLength;
        }
    }

    static class RideStop {
        public float pos_x;
        public float pos_y;
        public int type; //0 == depot, 1 == pickup, -1 == drop
        public int window_end;
    }

    public static RideStop readRide(InputReader reader) {
        try {
            RideStop r = new RideStop();
            reader.getInt(); //ignored
            r.pos_x = Float.parseFloat(reader.getString());
            r.pos_y = Float.parseFloat(reader.getString());
            reader.getInt(); //ignored
            r.type = reader.getInt();
            reader.getInt(); //ignored
            r.window_end = reader.getInt() * 100;
            return r;
        } catch (Exception e) {
            return null;
        }
    }

    public static void main(String[] args) {
        // Reading the data

        //TODO change file to test the various instances.
        InputReader reader = new InputReader("data/dialaride/minicustom");

        int nVehicles = reader.getInt();
        reader.getInt(); //ignore
        int maxRouteDuration = reader.getInt()*100;
        int vehicleCapacity = reader.getInt();
        int maxRideTime = reader.getInt()*100;

        RideStop depot = null;
        ArrayList<RideStop> pickupRideStops = new ArrayList<>();
        ArrayList<RideStop> dropRideStops = new ArrayList<>();
        boolean lastWasNotDrop = true;
        while (true) {
            RideStop r = readRide(reader);
            if(r == null)
                break;
            if(r.type == 0) {
                assert depot == null;
                depot = r;
            }
            else if(r.type == 1) {
                assert lastWasNotDrop;
                pickupRideStops.add(r);
            }
            else { //r.type == -1
                lastWasNotDrop = false;
                dropRideStops.add(r);
            }
        }
        assert depot != null;
        assert pickupRideStops.size() == dropRideStops.size();

        System.out.println("nVehicles = "+nVehicles+" nRequests = "+pickupRideStops.size());

        DialARideSolution sol = solve(nVehicles, maxRouteDuration, vehicleCapacity, maxRideTime, pickupRideStops, dropRideStops, depot);
        System.out.println(sol);
    }
}
