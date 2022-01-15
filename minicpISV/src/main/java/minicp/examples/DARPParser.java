package minicp.examples;

import minicp.util.io.InputReader;
import minicp.examples.DARPDataModel.*;

import java.util.ArrayList;

/**
 * @author Quentin Delmelle qdelmelle@gmail.com
 */

public class DARPParser {
    static int nRequests = 0;

    public static String getNameFromFilePath(String filePath) {
        String[] tmp = filePath.split("\\\\");
        String name = tmp[tmp.length - 1];
        if (name.isEmpty()) return filePath;
        else return name;
    }

    private static DARPStop lineToStop(String[] line) {
        if (line.length != 7) System.out.println("error reading stop: " + line.length);
        int load = Integer.parseInt(line[4]);
        return new DARPStop(
                Double.parseDouble(line[1]),
                Double.parseDouble(line[2]),
                Integer.parseInt(line[3]),
                load,
                Integer.parseInt(line[5]),
                Integer.parseInt(line[6])
        );
    }

    public static DynamicDARPInstance parseInstance(String filePath) {
        InputReader ir = new InputReader(filePath);
        String[] header = ir.getNextLine();
        assert (header.length == 5);

        int nVehicle = Integer.parseInt(header[0]);
        int nRequests = Integer.parseInt(header[1]);
        int maxRouteDuration = Integer.parseInt(header[2]);
        int vCapacity = Integer.parseInt(header[3]);
        int maxRideTime = Integer.parseInt(header[4]); //Max time in vehicle for request (service excluded)

        DynamicDARPInstance ret = new DynamicDARPInstance(getNameFromFilePath(filePath), nVehicle, vCapacity, maxRideTime, maxRouteDuration, 0);
        String[] depot = ir.getNextLine();

        DARPStop[] stops = new DARPStop[2 * nRequests + 2 * nVehicle];
        String[] stopline = ir.getNextLine();
        int i = 0;
        while (stopline != null) {
            stops[i] = lineToStop(stopline);
            stopline = ir.getNextLine();
            i++;
        }

        /*
         */
        for (int r = 0; r < nRequests; r++) {
            ret.addRequest(new DARPRequest(stops[r], stops[r + nRequests]));
        }

        return ret;
    }
}







