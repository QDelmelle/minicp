package minicp.examples;

import java.util.*;

class DARPDataModel {
    static class DARPStop {
        double x;
        double y;
        int service;
        int load;
        int winStart;
        int winEnd;

        public DARPStop(double x, double y, int service, int load, int winStart, int winEnd) {
            this.x = x;
            this.y = y;
            this.service = service;
            this.load = load;
            this.winStart = winStart;
            this.winEnd = winEnd;
        }

        boolean isDepot() {
            return load == 0;
        }

        boolean isStop() {
            return !isDepot();
        }

        boolean isPickup() {
            return load > 0;
        }

        boolean isDrop() {
            return !isPickup();
        }

        public String toString() {
            return "(" + x + ", " + y + ", " + service + ", " + load + ", [" + winStart + ", " + winEnd + "])";
        }
    }

    static class DARPRequest {
        DARPStop pickup;
        DARPStop drop;

        public DARPRequest(DARPStop p, DARPStop d) {
            pickup = p;
            drop = d;
        }

        public String toString() {
            return pickup.toString() + " ---> " + drop.toString();
        }
    }

    static class DARPInstance {
        String name;
        ArrayList<DARPRequest> requests;
        int nVehicles;
        int vCapacity;
        int maxRideTime;
        int maxDuration;
        int startTime;
        int timeHorizon;

        public DARPInstance(String name, int nVehicles, int vCapacity, int maxRideTime, int maxDuration) {
            this.name = name;
            this.nVehicles = nVehicles;
            this.vCapacity = vCapacity;
            this.maxRideTime = maxRideTime;
            this.maxDuration = maxDuration;
            startTime = 0;
            timeHorizon = 0;
            requests = new ArrayList<DARPRequest>();
        }

        public void addRequest(DARPRequest r) {
            requests.add(r);
            timeHorizon = Math.max(timeHorizon, r.drop.winEnd);
        }

        public int nRequests() {
            return requests.size();
        }

        public int nSites() {
            return nRequests() * 2 + nVehicles * 2;
        }

        public DARPStop[] getSites() {
            int n = nSites();
            DARPStop[] sites = new DARPStop[n];
            for (int i = 0; i < nRequests(); i++) {
                sites[i] = requests.get(i).pickup;
                sites[i + nRequests()] = requests.get(i).drop;
            }
            for (int v = 0; v < nVehicles; v++) {
                sites[2 * nRequests() + v] = new DARPStop(0.0, 0.0, 0, 0, 0, timeHorizon);
                sites[n - nVehicles + v] = new DARPStop(0.0, 0.0, 0, 0, 0, timeHorizon);
            }
            return sites;
        }

        public int[][] getDistances(DARPStop[] sites, int SCALING) {
            int n = sites.length;
            int[][] ret = new int[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double dx = (sites[i].x - sites[j].x) * (sites[i].x - sites[j].x);
                    double dy = (sites[i].y - sites[j].y) * (sites[i].y - sites[j].y);
                    ret[i][j] = (int) (Math.round(Math.sqrt(dx + dy) * SCALING));
                }
            }
            return ret;
        }
    }

    /**
     * an instance where vehicles don't necessarily start at the depot, each one has a route in progress.
     */
    static class DynamicDARPInstance extends DARPInstance {
        DARPPath[] paths; // the sites each vehicle has already visited, the last one being its current objective.

        public DynamicDARPInstance(String name, int nVehicles, int vCapacity, int maxRideTime, int maxDuration, int startTime) {
            super(name, nVehicles, vCapacity, maxRideTime, maxDuration);
            this.startTime = startTime;
            paths = new DARPPath[nVehicles];
            for (int i = 0; i < nVehicles; i++) paths[i] = new DARPPath(i);
        }

        void addStep(int vehicle, DARPStep s) {
            paths[vehicle].addStep(s);
        }

        public DynamicDARPInstance copy() {
            DynamicDARPInstance ret = new DynamicDARPInstance(this.name, this.nVehicles, this.vCapacity, this.maxRideTime, this.maxDuration, this.startTime);
            for (DARPRequest r : requests) ret.addRequest(r);
            return ret;
        }
    }

    static class DARPStep {
        int stop;
        int starttime;
        int endtime;

        public DARPStep(int stop, int st, int et) {
            this.stop = stop;
            starttime = st;
            endtime = et;
        }

        boolean isTimeFixed() {
            return starttime == endtime;
        }
    }

    static class DARPPath {
        int vehicle;
        List<DARPStep> steps;
        int len;

        public DARPPath(int v) {
            vehicle = v;
            steps = new ArrayList<DARPStep>();
            len = 0;
        }

        public void addStep(DARPStep s) {
            steps.add(s);
            len++;
        }
    }

    static class DARPSolution {
        public double fails;
        double cost;
        DARPPath[] paths;

        public DARPSolution(DARPPath[] paths, double cost, double fails) {
            this.paths = paths; this.cost = cost; this.fails = fails;
        }

        public boolean isEmpty() {
            return paths[0].len == 0;
        }

        public String toString() {
            String ret = "";
            for (int v = 0; v < paths.length; v++) {
                ret += "vehicle " + v + ": ";
                for (DARPStep s : paths[v].steps) {
                    ret += s.stop + "["+s.starttime+", "+s.endtime+"], ";
                }
                ret += "\n";
            }
            ret += "cost = "+cost;
            return ret;
        }
    }
}

