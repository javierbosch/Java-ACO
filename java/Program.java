import java.io.File;
import java.io.FileWriter;
import java.time.LocalDateTime; // Import the LocalDateTime class
import java.time.format.DateTimeFormatter;
import java.util.Scanner;
import java.util.stream.IntStream;

class Program
{
    public static void main(String[] args) {
        CreateFile();
        int nCopies = 1;
        int maxIterations = 100;
        DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("dd-MM HH:mm:ss");


        System.out.println(LocalDateTime.now().format(dateFormatter) + ": Starting simulation");
        String filename = "att48_half.txt";
        double[][] costMatrix = ReadProblem(filename, 24);

        int[] nAnts = { 20 };
        double[] alphas = {1};
        double[] betas = {2};
        //double[] Qs = {20};
        double[] rhos = IntStream.range(0, 101).asDoubleStream().map(i -> i * 0.01).toArray();
        double[] Qs = IntStream.range(1, 101).asDoubleStream().map(i -> i).toArray();
        //int[] sigmas = IntStream.range(0, 50).map(i -> 1 + (i * 2)).toArray();
        int[] sigmas = {1};

        ACO[] colonies = new ACO[nAnts.length * alphas.length * betas.length * sigmas.length * rhos.length * Qs.length * nCopies];

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
                                    colonies[index] = new ACO(costMatrix, alpha, beta, rho, Q, sigma, numAnts); ;

                                }
                            }
                        }
                    }
                }
            }
        }

        int iterations = 0;
        while(iterations < maxIterations){
            for (ACO colony : colonies) {
                colony.NextIter();
            }
            iterations++;
            System.out.println(LocalDateTime.now().format(dateFormatter) + String.format(": %.2f", ((double)iterations / (double)maxIterations) * 100));
            if((iterations == 1) | (iterations == 10) | (iterations == 50) | (iterations == 100) | (iterations == 250) | (iterations == 500)){
                Store(colonies, iterations);
            }
        }
        System.out.println(LocalDateTime.now().format(dateFormatter) +  ": All done - good job, Javi");
    }

    static double[][] ReadProblem(String fileName, int nCities) {
        File file = new File(fileName);
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


    static void CreateFile(){
        try {
            FileWriter writer = new FileWriter("../results/output-half2.csv",false);
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
            writer.write(headers);
            writer.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
    static void Store(ACO[] colonies, int iterations) {
        try{
            FileWriter writer = new FileWriter("../results/output-half2.csv",true);
            for (ACO colony : colonies) {
                String line =
                        iterations + "," +
                                colony.alpha + "," +
                                colony.beta + "," +
                                colony.sigma + "," +
                                colony.rho + "," +
                                colony.Q + "," +
                                colony.bestLength + "," +
                                colony.currBestLength + "," +
                                colony.nBestPathAnts + "," +
                                colony.nBestCurrentPathAnts + "\n";
                writer.write(line);
            }
            writer.close();
        }catch (Exception e){
            e.printStackTrace();
        }
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
