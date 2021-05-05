package minicp.examples;
import minicp.util.io.InputReader;
import minicp.examples.DARPDataModel.*;

import java.util.ArrayList;

public class DARPParser {
    final static int SCALING = 100;
    static int nRequests = 0;
    public static String getNameFromFilePath(String filePath) {
        String[] tmp = filePath.split("\\\\");
        String name = tmp[tmp.length-1];
        if(name.isEmpty()) return filePath; else return name;
    }

    private static DARPStop lineToStop(String[] line) {
        if (line.length != 7) System.out.println(line.length);
        int place = Integer.parseInt(line[0]);
        int load = Integer.parseInt(line[4]);
        int request = -1;
        if(load < 0) {
            request = place - 1 - nRequests;
        } else if(load > 0){
            request = nRequests;
            nRequests ++;
        }
        return new DARPStop(
                place,
                request,
                Double.parseDouble(line[1]),
                Double.parseDouble(line[2]),
                Integer.parseInt(line[3]) * SCALING,
                load,
                Integer.parseInt(line[5]) * SCALING,
                Integer.parseInt(line[6]) * SCALING
        );
    }

    public static DARPInstance parseInstance(String filePath) {
        InputReader ir = new InputReader(filePath);
        String[] header = ir.getNextLine();
        assert(header.length == 5);

        int nVehicle = Integer.parseInt(header[0]);
        //nRequests = Integer.parseInt(header[1]); // /!\ Not consistent between instances !!!
        int maxRouteDuration = Integer.parseInt(header[2]) * SCALING;
        int vCapacity = Integer.parseInt(header[3]);
        int maxRideTime = Integer.parseInt(header[4]) * SCALING; //Max time in vehicle for request (service excluded)

        DARPStop startDepot = lineToStop(ir.getNextLine()); //start depot

        ArrayList<DARPStop> stoplist = new ArrayList<DARPStop>();
        String[] stopline = ir.getNextLine();
        int nStop = 0;
        while (stopline != null) {
            DARPStop stop = lineToStop(stopline);
            if (stop.load > 0) {
                stoplist.add(stop);
                nStop++;
            }
            stopline = ir.getNextLine();
        }
        DARPStop[] stops = new DARPStop[nStop];
        for (int i = 0; i < nStop; i++) stops[i] = stoplist.get(i);

        DARPStop endDepot = (stoplist.get(nStop-1).load == 0) ?
                stoplist.get(nStop-1) : startDepot;

        // Generating stop for start end depot per vehicle
        int nSite = nVehicle*2 + nStop;
        DARPStop[] sites = new DARPStop[nSite];
        int i=0;
        for (; i<nStop; i++) { sites[i] = stops[i]; }
        for (; i<nStop+nVehicle; i++) { sites[i] = startDepot; }
        for (; i<nSite; i++) { sites[i] = endDepot; }

        int[][] distances = new int[nSite][nSite];
        for (i=0; i<nSite; i++){
            for (int j=0; j<nSite; j++) {
                double dx = (sites[i].x - sites[j].x) * (sites[i].x - sites[j].x);
                double dy = (sites[i].y - sites[j].y) * (sites[i].y - sites[j].y);
                distances[i][j] = (int) Math.round(Math.sqrt(dx + dy) * SCALING);
            }
        }

        DARPVehicle[] vehicles = new DARPVehicle[nVehicle];
        for (int v = 0; v<nVehicle; v++) {
            vehicles[v] = new DARPVehicle(nRequests+v, nRequests+nVehicle+v, vCapacity, maxRouteDuration);
        }

        DARPRequest[] requests = new DARPRequest[nRequests];
        for (int r = 0; r<nRequests; r++) {
            requests[r] = new DARPRequest(r, nRequests+r, maxRideTime);
        }

        return new DARPInstance(getNameFromFilePath(filePath),
                vehicles, requests, sites, distances);
    }
}







