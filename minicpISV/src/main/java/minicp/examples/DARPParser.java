package minicp.examples;
import minicp.util.io.InputReader;
import minicp.examples.DARPDataModel.*;

public class DARPParser {
    final int SCALING = 100;
    int nRequests = 0;
    String getNameFromFilePath(String filePath) {
        String[] tmp = filePath.split("/");
        String name = tmp[tmp.length-1].split("\\.")[0];
        if(name.isEmpty()) return filePath; else return name;
    }

    private DARPStop lineToStop(String line) {
        String[] splitted = line.split(" ");
        assert(splitted.length == 7);
        int place = Integer.parseInt(splitted[0]);
        int load = Integer.parseInt(splitted[4]);
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
                Double.parseDouble(splitted[1]),
                Double.parseDouble(splitted[2]),
                Integer.parseInt(splitted[3]) * SCALING,
                load,
                Integer.parseInt(splitted[5]) * SCALING,
                Integer.parseInt(splitted[6]) * SCALING
        );
    }

    public DARPInstance parseInstance(String filePath) {
        InputReader ir = new InputReader(filePath);
        String[] header = ir.getNextLine();
        assert(header.length == 5);

        int nVehicle = Integer.parseInt(header[0]);
        //nRequests = Integer.parseInt(header[1]); // /!\ Not consistent between instances !!!
        int maxRouteDuration = Integer.parseInt(header[2]) * SCALING;
        int vCapacity = Integer.parseInt(header[3]);
        int maxRideTime = Integer.parseInt(header[4]) * SCALING; //Max time in vehicle for request (service excluded)

        val startDepot: DARPStop = lineToStop(lines.next()) //start depot

        val stopLines: Array[String] = lines.filterNot(_.matches("\\s+")).toArray
        val parsedStops: Array[DARPStop] = stopLines.map(lineToStop) //Next lines are stops /!\ last one might be end depot on some instances

        val stops: Array[DARPStop] = parsedStops.filter(_.load != 0)
        val nStop: Int = stops.length

        val endDepot: DARPStop = if(parsedStops.last.load == 0) parsedStops.last else startDepot //End depot

        // Generating stop for start end depot per vehicle
        val depots: Array[DARPStop] = Array.fill(nVehicle)(startDepot) ++ Array.fill(nVehicle)(endDepot)

                val sites: Array[DARPStop] = stops ++ depots
        val nSite: Int = sites.length

        def dist(i: Int, j: Int): Int = {
                val x = (sites(i).x - sites(j).x) * (sites(i).x - sites(j).x)
                val y = (sites(i).y - sites(j).y) * (sites(i).y - sites(j).y)
                (math.sqrt(x + y) * SCALING).round.toInt
        }

        val distances: Array[Array[Int]] = Array.tabulate(nSite, nSite)((i, j) => {
            dist(i, j)
        })

        DARPInstance(
                getNameFromFilePath(filePath),
                Array.tabulate(nVehicle)(v => DARPVehicle(nRequests+v, nRequests+nVehicle+v, vCapacity, maxRouteDuration)),
        Array.tabulate(nRequests)(r => DARPRequest(r, nRequests+r, maxRideTime)),
        sites,
                distances
        )
    }
}
}






