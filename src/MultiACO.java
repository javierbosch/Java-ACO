import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

public class MultiACO implements Runnable{

    public double[][] dist;
    public double alpha;
    public double beta;
    public double rho;
    public double Q;
    public int sigma;
    public int numAnts;
    public int maxIterations;

    public double bestLength;
    public int numCities;
    public int[][] ants;
    public double[][] pheromones;
    public int iterations;
    public static FileWriter writer;
    public double currBestLength;
    public static int totalColonies;
    public static int totalDone;

    public MultiACO(double[][] dist, double alpha, double beta, double rho, double Q, int sigma, int numAnts, int maxIterations, String output_file, int totalColonies){
        this.dist = dist;
        this.alpha = alpha;
        this.beta = beta;
        this.rho = rho;
        this.Q = Q;
        this.sigma = sigma;
        this.numAnts = numAnts;
        this.maxIterations = maxIterations;
        try {
            MultiACO.writer = new FileWriter(output_file,true);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        bestLength = Double.MAX_VALUE;
        numCities = dist.length;
        iterations = 0;

        MultiACO.totalColonies = totalColonies;
        MultiACO.totalDone = 0;
    }

    public void run(){
        ants = InitAnts();
        pheromones = InitPheromones();

        for(int i = 0; i<maxIterations; i++){
            UpdateAnts();
            UpdatePheromones();
            currBestLength = CurrBestLength();
            if(currBestLength < bestLength){
                bestLength = currBestLength;
            }
            iterations++;
        }
        totalDone++;
        WriteState();
        // percent done 2 decimal places
        System.out.println(String.format("%.2f", (double)totalDone/totalColonies*100) + "% done");
        if (totalDone == totalColonies) {
            try {
                writer.close();
                System.out.printf("GOOD JOB YOU FINISHED ALL %d COLONIES, MULTITHREADING JAVI!!!\n", totalColonies);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
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
            writer.write(headers);
            writer.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    void WriteState() {
        try {
            String line =
                this.iterations + "," +
                this.alpha + "," +
                this.beta + "," +
                this.sigma + "," +
                this.rho + "," +
                this.Q + "," +
                this.bestLength + "," +
                this.currBestLength + "," + "\n";
            writer.write(line);
            writer.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private double CurrBestLength()
    {
        double tmp = Double.MAX_VALUE;
        for(int i = 0; i<numAnts; i++)
        {
            double len = Length(ants[i]);
            if(tmp > len)
            {
                tmp = len;
            }
        }
        return tmp;
    }


    private int[][] InitAnts()
    {
        int[][] ants = new int[numAnts][];
        for (int k = 0; k <= numAnts - 1; k++)
        {
            ants[k] = RandomTrail();
        }
        return ants;
    }

    private int[] RandomTrail()
    {
        Random rnd = new Random();
        int[] trail = new int[numCities];
        for (int i = 0; i <= numCities - 1; i++)
        {
            trail[i] = i;
        }
        for (int i = 0; i <= numCities - 1; i++)
        {
            int r = rnd.nextInt(numCities);
            int tmp = trail[r];
            trail[r] = trail[i];
            trail[i] = tmp;
        }
        return trail;
    }

    private double[][] InitPheromones()
    {
        double[][] pheromones = new double[numCities][];
        for (int i = 0; i <= numCities - 1; i++)
        {
            pheromones[i] = new double[numCities];
        }
        for (int i = 0; i <= pheromones.length - 1; i++)
        {
            for (int j = 0; j <= pheromones[i].length - 1; j++)
            {
                pheromones[i][j] = 0.01;
            }
        }
        return pheromones;
    }

    private void UpdateAnts()
    {
        Random rnd = new Random();
        for (int k = 0; k <= ants.length - 1; k++)
        {
            int start = rnd.nextInt(numCities);
            int[] newTrail = BuildTrail(start);
            ants[k] = newTrail;
        }
    }

    private int[] BuildTrail(int start)
    {
        int[] trail = new int[numCities];
        boolean[] visited = new boolean[numCities];
        trail[0] = start;
        visited[start] = true;
        for (int i = 0; i < numCities - 1; i++)
        {
            int cityX = trail[i];
            int next = NextCity(cityX, visited);
            trail[i + 1] = next;
            visited[next] = true;
        }
        return trail;
    }

    private int NextCity(int cityX, boolean[] visited)
    {
        Random rnd = new Random();
        double[] moveProb = MoveProb(cityX, visited);

        double[] cumulative = new double[moveProb.length + 1];
        for (int i = 0; i <= moveProb.length - 1; i++)
        {
            cumulative[i + 1] = cumulative[i] + moveProb[i];
        }
        double p = rnd.nextDouble();

        for (int i = 0; i <= cumulative.length - 2; i++)
        {
            if (p >= cumulative[i] && p < cumulative[i + 1])
            {
                return i;
            }
        }
        return NextCity(cityX, visited);
    }


    private double[] MoveProb(int cityX, boolean[] visited)
    {
        int numCities = pheromones.length;
        double[] tau = new double[numCities];
        double sum = 0.0;
        for (int i = 0; i <= tau.length - 1; i++)
        {
            if (i == cityX)
            {
                tau[i] = 0.0;
            }
            else if (visited[i])
            {
                tau[i] = 0.0;
            }
            else
            {
                tau[i] = Math.pow(pheromones[cityX][i], alpha) * Math.pow((1.0 / Distance(cityX, i)), beta);
                if (tau[i] < 0.0001)
                {
                    tau[i] = 0.0001;
                }
                else if (tau[i] > (Double.MAX_VALUE / (numCities * 100)))
                {
                    tau[i] = Double.MAX_VALUE / (numCities * 100);
                }
            }
            sum += tau[i];
        }

        double[] prob = new double[numCities];
        for (int i = 0; i <= prob.length - 1; i++)
        {
            prob[i] = tau[i] / sum;
        }
        return prob;
    }

    private double Distance(int cityX, int cityY)
    {
        return dist[cityX][cityY];
    }



    private void UpdatePheromones()
    {
        for (int i = 0; i <= pheromones.length - 1; i++)
        {
            for (int j = i + 1; j <= pheromones[i].length - 1; j++)
            {
                double sum = 0.0;

                for (int k = 0; k <= ants.length - 1; k++) {
                    double length = Length(ants[k]);
                    if (EdgeInTrail(i, j, ants[k])) {
                        sum +=  Q * (1 / length) ;

                    }
                }
                pheromones[i][j] = rho * pheromones[i][j] + (1 - rho) * sum;

                if (pheromones[i][j] < 0.00001)
                    pheromones[i][j] = 0.00001;
                else if (pheromones[i][j] > 1000000.0)
                    pheromones[i][j] = 1000000.0;

                pheromones[j][i] = pheromones[i][j];
            }
        }
    }


    private boolean EdgeInTrail(int cityX, int cityY, int[] trail)
    {
        // are cityX and cityY adjacent to each other in trail[]?
        int lastIndex = trail.length - 1;
        int idx = IndexOfTarget(trail, cityX);

        if (idx == 0 && trail[1] == cityY)
        {
            return true;
        }
        else if (idx == 0 && trail[lastIndex] == cityY)
        {
            return true;
        }
        else if (idx == 0)
        {
            return false;
        }
        else if (idx == lastIndex && trail[lastIndex - 1] == cityY)
        {
            return true;
        }
        else if (idx == lastIndex && trail[0] == cityY)
        {
            return true;
        }
        else if (idx == lastIndex)
        {
            return false;
        }
        else if (trail[idx - 1] == cityY)
        {
            return true;
        }
        else return trail[idx + 1] == cityY;
    }


    private int IndexOfTarget(int[] trail, int target)
    {
        // helper for RandomTrail
        for (int i = 0; i <= trail.length - 1; i++)
        {
            if (trail[i] == target)
            {
                return i;
            }
        }
        return -1;
    }

    private double Length(int[] trail)
    {
        // total length of a trail
        double result = Distance(trail[0], trail[trail.length-1]);
        for (int i = 0; i < trail.length - 1 ; i++)
        {
            result += Distance(trail[i], trail[i + 1]);
        }
        return result;
    }

}
