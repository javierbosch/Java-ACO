import java.util.Random;

public class ACO {
    public double alpha;
    public double beta;
    public double rho;
    public double Q;
    public double sigma;
    public int numCities;
    public int numAnts;

    public double[][] dists;
    public int[][] ants;
    public double bestLength;
    public double currBestLength;
    public double[][] pheromones;

    public int iterations;

    public int nBestPathAnts;
    public int nBestCurrentPathAnts;

    public int maxIterations;

    public ACO(double[][] dists, double alpha, double beta, double rho, double Q, int sigma, int numAnts)
    {
        this.dists = dists;
        this.alpha = alpha;
        this.beta = beta;
        this.rho = rho;
        this.Q = Q;
        this.sigma = sigma;
        this.numAnts = numAnts;

        bestLength = Double.MAX_VALUE;
        numCities = dists.length;
        ants = InitAnts();
        pheromones = InitPheromones();
        iterations = 0;
    }
    public void NextIter()
    {
        UpdateAnts();
        UpdatePheromones();
        currBestLength = CurrBestLength();
        countBestPathAnts();
        // double pheremonesDiversity = PheromoneDiversity(pheromones);
        // double hammingDistance2 = AverageHammingDistance2(ants);
        // double diversityOfLength = DiversityOfLength(ants);
        iterations++;
    }

    private double CurrBestLength()
    {
        double tmp = Double.MAX_VALUE;
        for(int i = 0; i<numAnts; i++)
        {
            if(tmp > Length(ants[i]))
            {
                tmp = Length(ants[i]);
                nBestCurrentPathAnts = 1;
            }
            else if(tmp == Length(ants[i]))
            {
                nBestCurrentPathAnts++;
            }
        }
        return tmp;
    }
    
    private void countBestPathAnts()
    {
        nBestPathAnts = 0;
        for (int i = 0; i < numAnts; i++)
        {
            if (bestLength == Length(ants[i]))
            {
                nBestPathAnts++;
            }
        }
    }
    
    private int[][] InitAnts()
    {
        int[][] ants = new int[numAnts][];
        for (int k = 0; k <= numAnts - 1; k++)
        {
            Random rnd = new Random();
            int start = rnd.nextInt(numCities);
            ants[k] = RandomTrail(start);
        }
        return ants;
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
                // otherwise first call to UpdateAnts -> BuiuldTrail -> fNode -> MoveProbs => all 0.0 => throws
            }
        }
        return pheromones;
    }

    private int[] RandomTrail(int start)
    {
        Random rnd = new Random();
        // helper for InitAnts
        int[] trail = new int[numCities];

        // sequential
        for (int i = 0; i <= numCities - 1; i++)
        {
            trail[i] = i;
        }

        // Fisher-Yates shuffle
        for (int i = 0; i <= numCities - 1; i++)
        {
            int r = rnd.nextInt(numCities);

            int tmp = trail[r];
            trail[r] = trail[i];
            trail[i] = tmp;
        }

        int idx = IndexOfTarget(trail, start);
        // put start at [0]
        int temp = trail[0];
        trail[0] = trail[idx];
        trail[idx] = temp;

        return trail;
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
        double result = 0.0;
        for (int i = 0; i <= trail.length - 2; i++)
        {
            result += Distance(trail[i], trail[i + 1]);
        }
        return result;
    }

    private void UpdateAnts()
    {
        Random rnd = new Random();
        int numCities = pheromones.length;
        for (int k = 0; k <= ants.length - 1; k++)
        {
            int start = rnd.nextInt(numCities);
            int[] newTrail = BuildTrail(k, start);
            ants[k] = newTrail;
        }
    }

    private int[] BuildTrail(int k, int start)
    {
        int numCities = pheromones.length;
        int[] trail = new int[numCities + 1];
        boolean[] visited = new boolean[numCities];
        trail[0] = start;
        visited[start] = true;
        for (int i = 0; i < numCities - 1; i++)
        {
            int cityX = trail[i];
            int next = NextCity(k, cityX, visited);
            trail[i + 1] = next;
            visited[next] = true;
        }
        trail[numCities] = trail[0];
        return trail;
    }

    private int NextCity(int k, int cityX, boolean[] visited)
    {
        Random rnd = new Random();
        // for ant k (with visited[]), at nodeX, what is next node in trail?
        double[] probs = MoveProbs(k, cityX, visited);

        double[] cumul = new double[probs.length + 1];
        for (int i = 0; i <= probs.length - 1; i++)
        {
            cumul[i + 1] = cumul[i] + probs[i];
            // consider setting cumul[cuml.Length-1] to 1.00
        }
        double p = rnd.nextDouble();

        for (int i = 0; i <= cumul.length - 2; i++)
        {
            if (p >= cumul[i] && p < cumul[i + 1])
            {
                return i;
            }
        }
        return NextCity(k, cityX, visited);
    }

    private double[] MoveProbs(int k, int cityX, boolean[] visited)
    {
        // for ant k, located at nodeX, with visited[], return the prob of moving to each city
        int numCities = pheromones.length;
        double[] taueta = new double[numCities];
        // inclues cityX and visited cities
        double sum = 0.0;
        // sum of all tauetas
        // i is the adjacent city
        for (int i = 0; i <= taueta.length - 1; i++)
        {
            if (i == cityX)
            {
                taueta[i] = 0.0;
                // prob of moving to self is 0
            }
            else if (visited[i])
            {
                taueta[i] = 0.0;
                // prob of moving to a visited city is 0
            }
            else
            {
                taueta[i] = Math.pow((1.0 / Distance(cityX, i)), alpha) * Math.pow(pheromones[cityX][i], beta);
                // could be huge when pheromone[][] is big
                if (taueta[i] < 0.0001)
                {
                    taueta[i] = 0.0001;
                }
                else if (taueta[i] > (Double.MAX_VALUE / (numCities * 100)))
                {
                    taueta[i] = Double.MAX_VALUE / (numCities * 100);
                }
            }
            sum += taueta[i];
        }

        double[] probs = new double[numCities];
        for (int i = 0; i <= probs.length - 1; i++)
        {
            probs[i] = taueta[i] / sum;
            // big trouble if sum = 0.0
        }
        return probs;
    }

    private void UpdatePheromones()
    {
        for (int i = 0; i <= pheromones.length - 1; i++)
        {
            for (int j = i + 1; j <= pheromones[i].length - 1; j++)
            {
                for (int k = 0; k <= ants.length - 1; k++)
                {
                    double length = Length(ants[k]);

                    double decrease = rho * pheromones[i][j];
                    double increase = 0.0;
                    if (EdgeInTrail(i, j, ants[k]))
                    {
                        increase = Q * (1 / length);
                        if (length < bestLength)
                        {
                            bestLength = length;
                            increase = increase * sigma;
                        }

                    }

                    pheromones[i][j] = decrease + increase;

                    if (pheromones[i][j] < 0.0001)
                    {
                        pheromones[i][j] = 0.0001;
                    }
                    else if (pheromones[i][j] > 100000.0)
                    {
                        pheromones[i][j] = 100000.0;
                    }

                    pheromones[j][i] = pheromones[i][j];
                }
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
        else if (trail[idx + 1] == cityY)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    private double Distance(int cityX, int cityY)
    {
        return dists[cityX][cityY];
        }

}
