package minicp.examples;

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

    class DARPRequest {
        int pickup;
        int drop;
        int maxRideTime;

        public DARPRequest(int p, int d, int mrt){
            pickup = p; drop = d; maxRideTime = mrt;
        }
    }

    class DARPVehicle {
        int start;
        int end;
        int capacity;
        int maxDuration;
        public DARPVehicle(int s, int e, int c, int md){
            start = s; end = e; capacity = c; maxDuration = md;
        }
    }

    class DARPInstance{
        String name;
        DARPVehicle[] vehicles;
        DARPRequest[] requests;
        DARPStop[] sites;
        int[][] distances;
        int nSites = sites.length;
        int nVehicles = vehicles.length;
        int nRequests = requests.length;
        public DARPInstance(String n, DARPVehicle[] v, DARPRequest[] r, DARPStop[] s, int[][] d){
            name = n; vehicles = v; requests = r; sites = s; distances = d;
        }
        int minTravelTime(int request) {
            int pickup = requests[request].pickup;
            int drop = requests[request].drop;
            return sites[pickup].service + distances[pickup][drop] + sites[drop].service;
        }
    }

    class DARPStep{
        int stop;
        int starttime;
        int endtime;
        public DARPStep(int stop, int st, int et){
            this.stop = stop; starttime = st; endtime = et;
        }
        boolean isTimeFixed(){ return starttime == endtime; }
    }

    class DARPPath{
        int vehicle;
        DARPStep[] steps;
        public DARPPath(int v, DARPStep[] s){
            vehicle = v; steps = s;
        }
    }

    class DARPSolution{
        DARPInstance instance;
        DARPPath[] paths;
        public DARPSolution(DARPInstance i, DARPPath[] p){
            instance = i; paths = p;
        }
    }
}

