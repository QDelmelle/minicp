package minicp.examples;
import minicp.util.io.InputReader;
import minicp.examples.DARPDataModel.*;

public class DARPParser {
    final int SCALING = 100;
    String getNameFromFilePath(String filePath) {
        String[] tmp = filePath.split("/");
        String name = tmp[tmp.length-1].split("\\.")[0];
        if(name.isEmpty()) return filePath; else return name;
    }

    public DARPInstance parseInstance(String filePath) {
        InputReader ir = new InputReader(filePath);
        String[] header = ir.getNextLine();
        assert(header.length == 5);

        val nVehicle: Int = header(0).toInt
        //  val nRequests = header(1).toInt // /!\ Not consistent between instances !!!
        val maxRouteDuration: Int = header(2).toInt * SCALING
        val vCapacity: Int = header(3).toInt
        val maxRideTime: Int = header(4).toInt * SCALING //Max time in vehicle for request (service excluded)

        var nRequests = 0

        def lineToStop(line: String): DARPStop = {
                val splitted = line.trim.split("\\s+")
        assert(splitted.length == 7)
        val place = splitted(0).toInt
        val load = splitted(4).toInt
        val request = if(load == 0) -1 else if(load < 0) place - 1 - nRequests else nRequests
        if(load > 0) nRequests +=1
        DARPStop(
                place,
                request,
                splitted(1).toDouble,
                splitted(2).toDouble,
                splitted(3).toInt * SCALING,
                load,
                splitted(5).toInt * SCALING,
                splitted(6).toInt * SCALING
        )}
}





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
