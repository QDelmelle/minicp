package minicp.examples;

import java.util.*;

class DARPDataModel{
    static class DARPStop{
        int place;
        int request;
        double x;
        double y;
        int service;
        int load;
        int winStart;
        int winEnd;

        public DARPStop(int p, int r, double x, double y, int s, int l, int ws, int we){
            place = p; request = r; this.x = x; this.y = y; service = s; load = l;
            winStart = ws; winEnd = we;
        }

        boolean isDepot() { return load == 0; }

        boolean isStop() { return !isDepot(); }

        boolean isPickup() { return load > 0; }

        boolean isDrop() { return !isPickup(); }
    }

    static class DARPRequest {
        int pickup;
        int drop;

        public DARPRequest(int p, int d){
            pickup = p; drop = d;
        }
    }

    static class DARPVehicle {
        int start;
        int end;
        int capacity;
        public DARPVehicle(int s, int e, int c){
            start = s; end = e; capacity = c;
        }
    }

    static class DARPInstance{
        String name;
        DARPVehicle[] vehicles;
        DARPRequest[] requests;
        DARPStop[] sites;
        int[][] distances;
        int nSites;
        int nVehicles;
        int nRequests;
        int maxRideTime;
        int maxDuration;
        public DARPInstance(String name, DARPVehicle[] vehicles, DARPRequest[] requests,
                            DARPStop[] sites, int[][] distances, int maxRideTime, int maxDuration){
            this.name = name; this.vehicles = vehicles; this.requests = requests;
            this.sites = sites; this.distances = distances;
            nSites = sites.length;
            nVehicles = vehicles.length;
            nRequests = requests.length;
            this.maxRideTime = maxRideTime;
            this.maxDuration = maxDuration;
        }
        int minTravelTime(int request) {
            int pickup = requests[request].pickup;
            int drop = requests[request].drop;
            return sites[pickup].service + distances[pickup][drop] + sites[drop].service;
        }
    }

    static class DARPStep{
        int stop;
        int starttime;
        int endtime;
        public DARPStep(int stop, int st, int et){
            this.stop = stop; starttime = st; endtime = et;
        }
        boolean isTimeFixed(){ return starttime == endtime; }
    }

    static class DARPPath{
        int vehicle;
        List<DARPStep> steps;
        int len;
        public DARPPath(int v){
            vehicle = v; steps = new ArrayList<DARPStep>(); len = 0;
        }
        public void addStep(DARPStep s){
            steps.add(s);len++;
        }
    }

    static class DARPSolution{
        DARPInstance instance;
        DARPPath[] paths;
        public DARPSolution(DARPInstance i, DARPPath[] p){
            instance = i; paths = p;
        }
        /**
         * @return true iif this is a valid solution to the instance and respects all the constraints.
         */
        public boolean isValid(){
            int[] g = new int[instance.nSites];
            int[] vehicle = new int[instance.nSites];
            DARPStep[] stops = new DARPStep[instance.nSites];
            int i = 0;
            for(int v=0;v<paths.length;v++){
                for(DARPStep step : paths[v].steps){
                    stops[step.stop] = step;
                    if (step.starttime < instance.sites[step.stop].winStart || step.endtime > instance.sites[step.stop].winEnd){
                        System.out.println("time window violated for site "+step.stop);
                        return false;
                    }
                    g[step.stop] = i;
                    vehicle[step.stop] = v;
                    i++;
                }
                if (paths[v].steps.get(paths[v].len-1).starttime - paths[v].steps.get(0).endtime > instance.maxDuration){
                    System.out.println("max route duration violated for vehicle "+v);
                    return false;
                }
            }
            if (i != instance.nSites){
                System.out.println("some sites are not visited");
                return false;
            }
            for(int r=0;r<instance.nRequests;r++){
                if (vehicle[r] != vehicle[r + instance.nRequests]){
                    System.out.println("request "+r+" is served by 2 vehicles!");
                    return false;
                }
                if (g[r] > g[r+instance.nRequests]) {
                    System.out.println("request "+r+" is not visited in the right order!");
                    return false;
                }
                if (stops[r+instance.nRequests].starttime - stops[r].endtime - instance.sites[r].service > instance.maxRideTime){
                    System.out.println("max ride time violated for request "+r);
                    return false;
                }
            }
            return true;
        }
    }

    static class DarpSol {
        int[] succ;
        int[] pred;
        int[] servingVehicle;
        double cost;
        double[] minServingTime;
        double[] maxServingTime;

        public DarpSol(int[] succ, int[] pred, int[] servingV, double c, double[] minServingTime, double[] maxServingTime) {
            this.succ = succ; this.pred = pred; servingVehicle = servingV; cost = c;
            this.minServingTime = minServingTime;
            this.maxServingTime = maxServingTime;
        }
    }
}

