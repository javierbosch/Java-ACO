import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime; // Import the LocalDateTime class
import java.time.format.DateTimeFormatter;
import java.util.Scanner;

class Main
{

    private static int nCopies;
    private static int maxIterations;
    private static int[] nAnts;
    private static double[] alphas;
    private static double[] betas;
    private static double[] Qs;
    private static double[] rhos;
    private static int[] sigmas;
    private static String map_file;
    private static String output_file;

    public static void main(String[] args) {

        DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("dd-MM HH:mm:ss");
        System.out.println(LocalDateTime.now().format(dateFormatter) + ": Starting simulation");

        readParameters("tests/" + args[0]);
        WriteHeaders();

        double[][] costMatrix = ReadProblem(map_file);

        int total = nAnts.length * alphas.length * betas.length * sigmas.length * rhos.length * Qs.length * nCopies;
        MultiACO[] colonies = new MultiACO[total];

        for (int n = 0; n < nAnts.length; n++)
        {
            int numAnts = nAnts[n];
            for (int a = 0; a < alphas.length; a++)
            {
                double alpha = alphas[a];
                for (int b = 0; b < betas.length; b++)
                {
                    double beta = betas[b];
                    for (int s = 0; s < sigmas.length; s++)
                    {
                        int sigma = sigmas[s];
                        for (int r = 0; r < rhos.length; r++)
                        {
                            double rho = rhos[r];
                            for (int q = 0; q < Qs.length; q++)
                            {
                                double Q = Qs[q];
                                for (int j = 0; j < nCopies; j++)
                                {
                                    int index = n + a * nAnts.length +
                                            b * nAnts.length * alphas.length +
                                            s * nAnts.length * alphas.length * betas.length +
                                            r * nAnts.length * alphas.length * betas.length * sigmas.length +
                                            q * nAnts.length * alphas.length * betas.length * sigmas.length * rhos.length +
                                            j * nAnts.length * alphas.length * betas.length * sigmas.length * rhos.length * Qs.length;
                                    colonies[index] = new MultiACO(costMatrix, alpha, beta, rho, Q, sigma, numAnts, maxIterations, output_file,total);
                                }
                            }
                        }
                    }
                }
            }
        }
        for (MultiACO colony : colonies) {
            Thread thread = new Thread(colony);
            thread.start();
        }
    }

    static void WriteHeaders(){
        String headers =
                "iteration," +
                        "alpha," +
                        "beta," +
                        "sigma," +
                        "rho," +
                        "Q," +
                        "bestLength," +
                        "bestCurrent," +
                        "nBestPathAnts," +
                        "nBestCurrentPathAnts\n";
        try {
            FileWriter writer = new FileWriter(output_file,false);
            writer.write(headers);
            writer.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    static void readParameters(String filename){
        File file = new File(filename);
        try {
            Scanner scanner = new Scanner(file);
            while (scanner.hasNextLine()){
                String line = scanner.nextLine();
                String[] split = line.split(",");
                String name = split[0];
                double[] values;
                if(split.length > 2){
                    values = Steps(Double.parseDouble(split[1]), Double.parseDouble(split[2]), Integer.parseInt(split[3]));
                }else{
                    try {
                        values = new double[]{Double.parseDouble(split[1])};
                    }catch (Exception e){
                        values = null;
                    }
                }
                switch (name){
                    case "nAnts":
                        nAnts = arrayDoubleToArrayInt(values);
                        break;
                    case "maxIterations":
                        maxIterations = Integer.parseInt(split[1]);
                        break;
                    case "nCopies":
                        nCopies = Integer.parseInt(split[1]);
                        break;
                    case "alphas":
                        alphas = values;
                        break;
                    case "betas":
                        betas = values;
                        break;
                    case "Qs":
                        Qs = values;
                        break;
                    case "rhos":
                        rhos = values;
                        break;
                    case "sigmas":
                        sigmas = arrayDoubleToArrayInt(values);
                        break;
                    case "map_file":
                        map_file = "maps/" + split[1];
                        break;
                    case "output_file":
                        output_file = "results/" + split[1];
                        break;
                }
            }
            scanner.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private static int[] arrayDoubleToArrayInt(double[] values) {
        int[] result = new int[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = (int) values[i];
        }
        return result;
    }

    static double[] Steps(double start, double end, int nSteps) {
        if(nSteps == 1)
            return new double[] {start};
        if(nSteps == 2)
            return new double[] {start, end};
        double[] steps = new double[nSteps];
        double step = (end - start) / (nSteps - 1);
        for (int i = 0; i < nSteps; i++) {
            steps[i] = start + i * step;
        }
        return steps;
    }

    static double[][] ReadProblem(String fileName) {
        File file = new File(fileName);

        String num = "";
        for (int i = 0; i < fileName.length(); i++) {
            char c = fileName.charAt(i);
            if (Character.isDigit(c)) {
                num += c;
            }
        }
        int nCities = Integer.parseInt(num);
        try {
            Scanner scanner = new Scanner(file);
            double[][] coordinates = new double[nCities][];
            int i = 0;
            while (scanner.hasNextLine()){
                String line = scanner.nextLine();
                String x_s = line.split(",")[0];
                String y_s = line.split(",")[1];
                double x = Double.parseDouble(x_s);
                double y = Double.parseDouble(y_s);
                coordinates[i] = new double[] {x, y};
                i++;
            }
            scanner.close();
            return BuildCostMatrix(coordinates, nCities);
        }catch (Exception e){
            e.printStackTrace();
        }
        return null;
    }


    static double[][] BuildCostMatrix(double[][] coordinates, int dim)
    {
        double[][] costMatrix = new double[dim][];
        for (int i = 0; i < dim; i++)
        {
            costMatrix[i] = new double[dim];
            for (int j = 0; j < dim; j++)
            {
                if (i == j)
                {
                    costMatrix[i][j] = 0;
                }
                else
                {
                    double xDiff = coordinates[i][0] - coordinates[j][0];
                    double yDiff = coordinates[i][1] - coordinates[j][1];
                    costMatrix[i][j] = Math.sqrt(xDiff * xDiff + yDiff * yDiff);
                }
            }
        }
        return costMatrix;
    }

}
