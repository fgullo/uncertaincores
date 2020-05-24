//implementation with BigDecimals
//probabiities computed with the novel dynamic-programming method proposed in the paper
//incremental initial eta-degrees computation
//TreeSets: experimentally faster than HashSets and without strange time fluctuations

//This is the most reliable implementation

//It differs from the basic KEtaCores_DynamicProgramming_BigDecimals_IncrementalInitialEtaDegrees because it includes the newest Bata-function-based lower-bounding technique 



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
//import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;


public class KEtaCores_DynamicProgramming_BigDecimals_IncrementalInitialEtaDegrees_LB 
{
    //protected static boolean deterministic_mode = true;
    protected static boolean deterministic_mode = false;
    
    //protected static MathContext mc = MathContext.DECIMAL32;
    protected static MathContext mc = MathContext.DECIMAL64;
    //protected static MathContext mc = MathContext.DECIMAL128;
    //protected static MathContext mc = new MathContext(100);
    //protected static MathContext mc = MathContext.UNLIMITED;
    
    protected static MathContext[] mcs = new MathContext[]{MathContext.DECIMAL32,MathContext.DECIMAL64,MathContext.DECIMAL128,new MathContext(70)};
    
    protected static Map<Integer,Map<Integer,Double>> g;//input probabilistic graph
    protected static int n;//vertices in g
    protected static long m;//edges in g
    
    //protected static double epsilon = 0.0001;
    protected static double epsilon = 1.0E-16;
    protected static double epsilon0 = epsilon;
    protected static double epsilon1 = 1.0 - epsilon;
    
    protected static Map<Integer,BigDecimal[]> etaDegProbs;//for each vertex v in g, it contains Pr[deg(v)=k], k=0,...,deg(v)
                                                        //note that the vector is actually filled until k=eta-deg(v) (which could be < deg(v))
    protected static int[] etaDegs;//etaDegs[u] is the (current) eta-degree of vertex u
    
    protected static BigDecimal one;
    protected static BigDecimal zero;
    protected static BigDecimal nepero;
    protected static BigDecimal bigEta;
    
    protected static String matlabOutputPath = "beta_function_values"+File.separator;
    protected static String matlabOutputSuffix = "_beta_matlab.txt_output";
    protected static double[][] betaFunctionValues;
    protected static boolean[] hasExactEtaDegBeenComputed;
    protected static double initialDegreeTime; //in seconds
    
    protected static long initialEtaDegreeThicks = 0;
    
    protected static int runs = 1;
    protected static int startrun = 1;
    
    
    
    public static void main(String[] args) throws IOException
    {
        
        /*
        String dataset = args[0];
        double eta = Double.parseDouble(args[1]);
        run(dataset,eta);
        */
        
        
        //String dataset = "toy.txt";
        //String dataset = "dblp.k100.e1E-3.c2.q0.01.txt";
        String dataset = "biomine_probgraph.txt";
        //String dataset = "biomine_probgraph.txt.graph.part.20_3";
        //String dataset = "Flickr.txt.graph.part.3_0";
        //String dataset = "Flickr.txt";
        //String dataset = "LastFM.txt";
        //String dataset = "DBLP_mu=5.txt";
        //String dataset = "Flickr.txt.graph.part.3_0.graph.part.250_11";
        String suffix = ""+System.currentTimeMillis();
        double[] etas = null;
        
        
        if(args != null)
        {
            if (args.length > 0)
            {
                dataset = args[0];
            }
            
            if (args.length > 1)
            {
                etas = new double[]{Double.parseDouble(args[1])};
            }
            
            if (args.length > 2)
            {
                runs = Integer.parseInt(args[2]);
            }
        }
        
        if(etas == null)
        {
            //etas = new double[]{0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0};
            //etas = new double[]{0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};
            //etas = new double[]{0.9,0.1};
            etas = new double[]{0.1};
            //etas = new double[]{0.2};
        }
        
        runEfficiency_MultiEta_MultiRun(dataset, etas, suffix);
        //runVaryingPrecisionEfficiency_MultiEta_MultiRun(dataset, etas, suffix);
        
        //String output = "output"+File.separator+dataset+"_initialEtaDegrees.txt";
        //BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        for(double eta:etas)
        {
            //run(dataset,eta,suffix);
            //runVaryingPrecision(dataset,eta,suffix);
            //printInitialEtaDegrees(dataset, eta, bw);
        }
        //bw.flush();
        //bw.close();
        
    }
    
    private static void runVaryingPrecision(String dataset, double eta, String suffix) throws IOException
    {
        for(MathContext x:mcs)
        {
            mc = x;
            run(dataset,eta,suffix);
        }
    }
    
    private static void runVaryingPrecisionEfficiency_MultiEta_MultiRun(String dataset, double[] etas, String suffix) throws IOException
    {
        for(MathContext x:mcs)
        {
            mc = x;
            runEfficiency_MultiEta_MultiRun(dataset, etas, suffix);
        }
    }
    
    private static void runEfficiency_MultiEta_MultiRun(String dataset, double[] etas, String suffix) throws IOException
    {
        String head = "dataset"+"\t"+"eta"+"\t"+"runs"+"\t"+"initial eta-deg thicks"+"\t"+"time initial eta-degs (secs)"+"\t"+"time main cycle  (secs)"+"\t"+"time total  (secs)";
 
        new File("output").mkdirs();
        String p = (deterministic_mode)?"":("prec="+mc.getPrecision()+"_");
        String f = "output"+File.separator+"efficiencyTest_LB_runs"+runs+"_"+p+dataset+"_"+suffix;
        BufferedWriter bw = new BufferedWriter(new FileWriter(f));
        bw.write(head);
        bw.newLine();       
        
        for (double eta:etas)
        {
            double time1avg = 0.0;
            double time2avg = 0.0;
            double timeavg = 0.0;
            
            for(int r=startrun; r<=runs; r++)
            {
                System.out.println("******************************************************************");
                System.out.println("RUN "+r+" of "+runs);
                System.out.println("******************************************************************");

                one = new BigDecimal(1.0,mc);
                zero = new BigDecimal(0.0,mc);
                nepero = new BigDecimal(Math.exp(1.0),mc);
                initialEtaDegreeThicks = 0;

                System.out.println("=========================================================");
                String s = (deterministic_mode)?"DETERMINISTIC":("eta="+eta);
                System.out.println(s);
                if(!deterministic_mode)
                {
                    System.out.println("Precision (digits)="+mc.getPrecision());
                }
                loadGraph(dataset);

                if(!deterministic_mode)
                {
                    loadBetaFunctionValuesAndComputeEtaDegLBs(dataset,eta);
                }
                else
                {
                    etaDegs = new int[n];
                    for(int i=0; i<etaDegs.length; i++)
                    {
                        etaDegs[i] = getDegree(i);
                    }
                }

                initialDegreeTime = 0.0;
                bigEta = new BigDecimal(eta,mc);
                //etaDegs = new int[n];
                etaDegProbs = new HashMap<Integer, BigDecimal[]>();
                double start = System.currentTimeMillis();
                //computeMinimumNeighborhoodProbability();
                int c[] = mainCycle(eta);
                double time = System.currentTimeMillis() - start;
                
                
                if (r>0)
                {
                    time /= 1000;
                    double time1 = initialDegreeTime;
                    double time2 = time - time1;

                    time1avg += time1;
                    time2avg += time2;
                    timeavg += time;
                }
                //printOutput(c,eta,dataset,suffix);
                System.out.println("");
                System.out.println("---> Execution time="+(time)+"secs <---");
                System.out.println("---> Initial-degrees-computation time="+initialDegreeTime+"secs <---");
                System.out.println("=========================================================");
                System.out.println("InitialDegreeThicks="+initialEtaDegreeThicks);
                System.out.println("");
                System.out.println("");                 
            }

            time1avg /= runs;
            time2avg /= runs;
            timeavg /= runs;

            String data = dataset+"\t"+eta+"\t"+runs+"\t"+initialEtaDegreeThicks+"\t"+time1avg+"\t"+time2avg+"\t"+timeavg;
            bw.write(data);
            bw.newLine();
            bw.flush();
        }
        bw.flush();
        bw.close();
    }

    private static void run(String dataset, double eta, String suffix) throws IOException
    {
        one = new BigDecimal(1.0,mc);
        zero = new BigDecimal(0.0,mc);
        nepero = new BigDecimal(Math.exp(1.0),mc);
        initialEtaDegreeThicks = 0;
        
        System.out.println("=========================================================");
        String s = (deterministic_mode)?"DETERMINISTIC":("eta="+eta);
        System.out.println(s);
        if(!deterministic_mode)
        {
            System.out.println("Precision (digits)="+mc.getPrecision());
        }
        loadGraph(dataset);
        
        if(!deterministic_mode)
        {
            loadBetaFunctionValuesAndComputeEtaDegLBs(dataset,eta);
        }
        else
        {
            etaDegs = new int[n];
            for(int i=0; i<etaDegs.length; i++)
            {
                etaDegs[i] = getDegree(i);
            }
        }
        
        initialDegreeTime = 0.0;
        bigEta = new BigDecimal(eta,mc);
        //etaDegs = new int[n];
        etaDegProbs = new HashMap<Integer, BigDecimal[]>();
        double start = System.currentTimeMillis();
        //computeMinimumNeighborhoodProbability();
        int c[] = mainCycle(eta);
        double time = System.currentTimeMillis() - start;
        printOutput(c,eta,dataset,suffix);
        System.out.println("");
        System.out.println("---> Execution time="+(time/1000)+"secs <---");
        System.out.println("---> Initial-degrees-computation time="+initialDegreeTime+"secs <---");
        System.out.println("=========================================================");
        System.out.println("InitialDegreeThicks="+initialEtaDegreeThicks);
        System.out.println("");
        System.out.println("");               
    }  
    


    private static void loadGraph(String dataset) throws IOException
    {
        System.out.println("--------------------------------------------------");
        System.out.print("Loading graph ("+dataset+")...");
        long start = System.currentTimeMillis();
        g = new HashMap<Integer, Map<Integer, Double>>();
        n = 0;
        m = 0;
        BufferedReader br = new BufferedReader(new FileReader(dataset));
        String line = br.readLine();
        while(line != null)
        {
            double[] uvp = getUVP(line);
            int u = (int)uvp[0];
            int v = (int)uvp[1];
            double p = uvp[2];
            checkProbability(p);
            if(p>=epsilon0)
            {
                if (u != v)
                {
                    if(!g.containsKey(u))
                    {
                        g.put(u, new HashMap<Integer,Double>());
                    }
                    if(!g.containsKey(v))
                    {
                        g.put(v, new HashMap<Integer,Double>());
                    }
                    g.get(u).put(v,p);
                    g.get(v).put(u,p);

                    m++;
                }
            }
            line = br.readLine();
        }
        n = g.size();
        double time = System.currentTimeMillis() - start;
        System.out.println("DONE! ("+(time/1000)+"secs)");
        System.out.println("--------------------------------------------------");
        //System.out.println("Loading time="+(time/1000)+" secs");
        System.out.println("#vertices (n)="+g.size());
        System.out.println("#edges (m)="+m);
        
        double avgDeg = 0;
        double squareDegSum = 0.0;
        int maxDeg = Integer.MIN_VALUE;
        for(int u:g.keySet())
        {
            int d = g.get(u).size();
            avgDeg += d;
            squareDegSum += d*d;
            if(d > maxDeg)
            {
                maxDeg = d;
            }
        }
        avgDeg /= n;
        System.out.println("max degree (Delta)="+maxDeg);
        System.out.println("avg degree="+avgDeg);
        System.out.println("sum of square degrees (S)="+squareDegSum);
        System.out.println("Delta*m="+maxDeg*m);
        System.out.println("Ratio Delta*m/S="+((double)maxDeg*m)/squareDegSum);
        System.out.println("Ratio S/m="+squareDegSum/m);
        System.out.println("--------------------------------------------------");
    }

    private static void computeExactEtaDegree(int v, double eta) 
    {
        //testBigDecimal2(eta);
        //testBigDecimal(eta);
        //System.out.println("");
        //System.out.println("--------------------------------------------------");
        //System.out.println("Computing initial eta-degrees...");
        long start = System.currentTimeMillis();

        int dv = getDegree(v);
        double[] non1prob = getNon1ProbsDouble(v);
        BigDecimal[] non1probBig = getNon1Probs(non1prob);
        int dvPrime = non1prob.length;//number of non-one probabilities
        BigDecimal[] p = new BigDecimal[dvPrime+1];

        BigDecimal[] X = new BigDecimal[dvPrime+1];
        for(int h=0; h<X.length; h++)
        {
            X[h] = zero;
        }

        BigDecimal prGreaterThanK = one;
        int k=0;
        //while(k<=dv && prGreaterThanK>=eta)
        while(k<=dvPrime && prGreaterThanK.compareTo(bigEta)>0)
        {
            p[k] = computeInitialEtaDegreeProbability(v,k,non1probBig,X);
            BigDecimal bd = p[k];
            String bds = bd.toString();
            prGreaterThanK = prGreaterThanK.subtract(bd,mc);
            //checkProbability(prGreaterThanK);
            String s = prGreaterThanK.toString();
            k++;
        }
        
        initialEtaDegreeThicks += (dvPrime+1)*k;

        etaDegs[v] = k-1+(dv-dvPrime);

        BigDecimal[] realP = new BigDecimal[dv+1];
        for(int i=0; i<dv-dvPrime; i++)
        {
            realP[i] = zero;
        }
        for(int i=dv-dvPrime; i<realP.length; i++)
        {
            realP[i] = p[i-(dv-dvPrime)];
        }

        etaDegProbs.put(v, realP);

        /*
        //if(eta==0.0 && etaDegs[v] != dv)
        if(true)
        //if(false)
        //if(dv<=10 && eta==0.0 && etaDegs[v] != dv)
        {
            System.out.println(v+"\t"+dv+"\t"+etaDegs[v]);
        }
        */
        //System.out.println("");
        //double xxx = Math.pow(maxPeTilde,maxDeg);
        //double yyy = Math.pow(minPeTilde,minDeg);

        double time = System.currentTimeMillis() - start;
        //System.out.println("DONE! ("+(time/1000)+"secs)");
        //System.out.println("--------------------------------------------------");
        
        initialDegreeTime += time/1000;
        
        //System.exit(-1);
        
        //printEtaDegreeStatistics();
    }
    
    private static int getDegree(int v) 
    {
        return g.get(v).size();
    }    


    private static int[] mainCycle(double eta)
    {
        BigDecimal bigEta = new BigDecimal(eta,mc);
        long start = System.currentTimeMillis();
        System.out.println("");
        System.out.println("--------------------------------------------------");
        System.out.println("Main cycle started...");
        
        //initializing data structures
        int[] c = new int[n];
        for(int i=0; i<c.length; i++)
        {
            c[i] = -100;
        }
        int maxEtaDeg = Integer.MIN_VALUE;
        for(Map<Integer,Double> map:g.values())
        {
            int s = map.size();
            if (s > maxEtaDeg)
            {
                maxEtaDeg = s;
            }
        }
            
        /*
        for(int u=0; u<n; u++)
        {
            int d = getDegree(u);
            if (d > maxEtaDeg)
            {
                maxEtaDeg = d;
            }
        }
        */ 
        /*
        for(int etaDeg:etaDegs)
        {
            if (etaDeg > maxEtaDeg)
            {
                maxEtaDeg = etaDeg;
            }
        }
        */
        int maxFilledDPosition = -1;
        Set<Integer>[] D = new Set[maxEtaDeg+1];
        for(int u=0; u<etaDegs.length; u++)
        {
            int etaDeg = etaDegs[u];
            if(etaDeg>maxFilledDPosition)
            {
                maxFilledDPosition = etaDeg;
            }
            if(D[etaDeg]==null)
            {
                //D[etaDeg] = new HashSet<Integer>();
                D[etaDeg] = new TreeSet<Integer>();
            }
            D[etaDeg].add(u);
        }
        
        //main cycle
        int processed = 0;
        int parts = 20;
        int ratio = n/parts;
        for(int k=0; k<D.length && k<=maxFilledDPosition; k++)
        {
            /*
            int scale = Math.min(20,D.length);
            if (k%(D.length/scale)==0)
            {
                double p = ((double)1.0)/scale*100;
                double x = ((double)k)/(D.length/scale)*p;
                System.out.println(x+"%");
            }
            */
            //System.out.print("k="+k+" (of"+D.length+") started...");
            Set<Integer> s = D[k];
            
            /*
            if(s != null && !s.isEmpty())
            {
                System.out.println(k);  
            }
            */
            
            while(s != null && !s.isEmpty())
            {
                //int v = s.iterator().next();
                int v = ((TreeSet<Integer>)s).first();
                
                int newEtaDeg = k;
                
                if(!deterministic_mode && !hasExactEtaDegBeenComputed[v])
                {
                    hasExactEtaDegBeenComputed[v] = true;
                    if(getDegree(v) > k)
                    {
                        computeExactEtaDegree(v,eta);
                        //hasExactEtaDegBeenComputed[v] = true;
                        newEtaDeg = etaDegs[v];
                        if(newEtaDeg>maxFilledDPosition)
                        {
                            maxFilledDPosition = newEtaDeg;
                        }
                        if(etaDegs[v] > k)
                        {
                            s.remove(v);
                            if(D[newEtaDeg]==null)
                            {
                                //D[newEtaDeg] = new HashSet<Integer>();
                                D[newEtaDeg] = new TreeSet<Integer>();
                            }
                            D[newEtaDeg].add(v);
                        }
                    }
                }
                
                if (newEtaDeg <= k)
                {
                    s.remove(v);
                    c[v] = k;
                    processed++;
                    if(processed%ratio==0)
                    {
                        double perc = processed/ratio*(100/parts);
                        System.out.println(perc+"%");
                    }
                    Map<Integer,Double> neighborhood = g.get(v);
                    if(neighborhood != null)
                    {
                        for(int u:neighborhood.keySet())
                        {
                            if(u != v)
                            {
                                int prevEtaDeg = etaDegs[u];
                                if(prevEtaDeg>k)
                                {
                                    //recomputing eta-deg(u)
                                    newEtaDeg = prevEtaDeg-1;
                                    if(!deterministic_mode)
                                    {
                                        if(!hasExactEtaDegBeenComputed[u])
                                        {
                                            newEtaDeg = updateEtaDegLB(u,prevEtaDeg,eta);
                                        }
                                        else
                                        {
                                            newEtaDeg = updateEtaDeg(u,neighborhood.get(u),bigEta);
                                        }
                                    }

                                    if (newEtaDeg != prevEtaDeg)
                                    {
                                        //moving u from D[prevEtaDeg] to D[newEtaDeg]
                                        D[prevEtaDeg].remove(u);
                                        if(D[newEtaDeg]==null)
                                        {
                                            //D[newEtaDeg] = new HashSet<Integer>();
                                            D[newEtaDeg] = new TreeSet<Integer>();
                                        }
                                        D[newEtaDeg].add(u);
                                        etaDegs[u] = newEtaDeg;
                                    }
                                }

                                //removing v from the u's neighbors
                                g.get(u).remove(v);
                            }
                        }
                    }
                }
            }
            //System.out.println("DONE!");
        }
        double time = System.currentTimeMillis() - start;
        System.out.println("DONE! ("+(time/1000)+"secs)");
        System.out.println("--------------------------------------------------");
        
        return c;
    }

    private static int updateEtaDeg(int u, double pe, BigDecimal bigEta) 
    {
        int etaDeg = etaDegs[u];
        
        if(etaDeg==0)
        {
            return 0;
        }
        
        BigDecimal[] p = etaDegProbs.get(u);
        if(pe>=epsilon1)//if the edge removed has probability 1
        {
            for (int i=0; i<etaDeg; i++)
            {
                p[i] = p[i+1];
            }
            return Math.max(0,etaDeg-1);
        }
        
        //update probabilities Pr[deg(u)=i], i=0,...,eta-deg(u)
        //p[0] /= (1.0 - pe);
        BigDecimal bigPe = new BigDecimal(pe,mc);
        BigDecimal bigOneMinusPe = one.subtract(bigPe,mc);
        BigDecimal bigPeDividedByOneMinusPe = divideBigDecimals(bigPe, bigOneMinusPe);
        
        p[0] = divideBigDecimals(p[0],bigOneMinusPe);        
        for(int i=1; i<=etaDeg; i++)
        {
            BigDecimal bd1 = divideBigDecimals(p[i],bigOneMinusPe);
            BigDecimal bd2 = bigPeDividedByOneMinusPe.multiply(p[i-1],mc);
            
            p[i] = bd1.subtract(bd2,mc);
        }
        
        //compute new eta-degree
        BigDecimal prGreaterThanK = one;
        int k = 0;
        int limit = (etaDeg<getDegree(u)-1)?etaDeg:getDegree(u)-1;
        if (limit<0)
        {
            limit = 0;
        }
        while(k<=limit && prGreaterThanK.compareTo(bigEta)>0)
        {
            prGreaterThanK = prGreaterThanK.subtract(p[k],mc);
            k++;
        }
        
        if(k-1<etaDeg-1)//if the new eta-degree decreases of more than one, then it's an error due to numerical stability
        {
            return etaDeg-1;
        }
        
        return k-1;
    }

    private static void printOutput(int[] c, double eta, String d, String suffix) throws IOException
    {
        new File("output").mkdirs();
        String e = (deterministic_mode)?"deterministic":("eta="+eta);
        String p = (deterministic_mode)?"":("prec="+mc.getPrecision()+"_");
        String f = "output"+File.separator+"(k-eta)-coreNumbers_"+e+"_"+p+d+"_"+suffix;
        BufferedWriter bw = new BufferedWriter(new FileWriter(f));
        for(int i=0; i<c.length; i++)
        {
            bw.write(i+"\t"+c[i]);
            bw.newLine();
        }
        bw.flush();
        bw.close();
        
        
        Map<Integer,Integer> map = new HashMap<Integer, Integer>();
        for(int x:c)
        {
            if(!map.containsKey(x))
            {
                map.put(x,0);                        
            }
            int y = map.get(x);
            y++;
            map.put(x, y);
        }
        
        f = "output"+File.separator+"(k-eta)-coreSizes_"+e+"_"+p+d+"_"+suffix;
        bw = new BufferedWriter(new FileWriter(f));
        Integer[] keys = new Integer[map.size()];
        map.keySet().toArray(keys);
        Arrays.sort(keys);
        for(int x:keys)
        {
            bw.write(x+"\t"+map.get(x));
            bw.newLine();
        }
        bw.flush();
        bw.close();
        System.out.println("");
        System.out.println(keys[keys.length-1]+"\t"+map.get(keys[keys.length-1]));
    }

    private static void printEtaDegreeStatistics()
    {
        double avgDeg = 0;
        double squareDegSum = 0.0;
        double squareDegSum2 = 0.0;
        int maxDeg = Integer.MIN_VALUE;
        for(int i=0; i<etaDegs.length; i++)
        {
            int d = etaDegs[i];
            avgDeg += d;
            squareDegSum += d*d;
            squareDegSum2 += d*getDegree(i);
            if(d > maxDeg)
            {
                maxDeg = d;
            }
        }
        avgDeg /= n;
        System.out.println("max eta-degree (eta-Delta)="+maxDeg);
        System.out.println("avg eta-degree="+avgDeg);
        System.out.println("--------------------------------------------------");
        System.out.println("eta-Delta*m="+maxDeg*m);
        System.out.println("Ratio eta-Delta*m/eta-S="+(((double)maxDeg*m)/squareDegSum));
        System.out.println("--------------------------------------------------");
        System.out.println("sum of square eta-degrees (eta-S)="+squareDegSum);
        System.out.println("Ratio eta-S/m="+(squareDegSum/m));
        System.out.println("--------------------------------------------------");
        System.out.println("sum of eta-degree*degree (eta-S_2)="+squareDegSum2);
        System.out.println("Ratio eta-S_2/m="+(squareDegSum2/m));
        System.out.println("--------------------------------------------------");
        
    }

    private static double[] getUVP(String line) 
    {
        StringTokenizer st = new StringTokenizer(line,"\t, ;");
        int u = Integer.parseInt(st.nextToken());
        int v = Integer.parseInt(st.nextToken());
        double p = Double.parseDouble(st.nextToken());
        return new double[]{u,v,p};
    }
    
    public static void checkProbability(double p)
    {
        if(p > 1.0 || p <=0.0)
        {
            throw new RuntimeException("Probability exception: "+p);
        }
    }

    private static void checkLargerThanZero(double p) 
    {
        if(p <0.0)
        {
            throw new RuntimeException("Probability exception: "+p);
        }
    }
    
    public static void computeMinimumNeighborhoodProbability()
    {
        double min = Double.POSITIVE_INFINITY;
        for(int u:g.keySet())
        {
            double p = getMininmumNeighborhoodProbability(u);
            if(p<min)
            {
                min = p;
            }
        }
        System.out.println(min);
    }
    
    public static double getMininmumNeighborhoodProbability(int u)
    {        
        double p = 1.0;
        for(int x:g.get(u).keySet())
        {
            double pr = g.get(u).get(x);
            if (pr<=0.5 && pr > 0.0)
            {
                p *= pr;
            }
            else if (pr < 1.0)
            {
                p *= (1.0-pr);
            }
        }
        
        //checkProbability(p);
        
        if(p==0.0)
        {
            //System.out.println("");
        }
        
        return p;
    }

    private static BigDecimal[] getNon1Probs(double[] v) 
    {        
        BigDecimal[] vec = new BigDecimal[v.length];
        for(int i=0; i<v.length; i++)
        {
            vec[i] = new BigDecimal(v[i],mc);
        }        
        return vec;
    }
    
    private static double[] getNon1ProbsDouble(int v) 
    {
        int number = 0;
        for(double pr:g.get(v).values())
        {
            if (pr <= epsilon1)
            {
                number++;
            }
        }
        
        double[] vec = new double[number];
        int i=0;
        for(double pr:g.get(v).values())
        {
            if (pr <= epsilon1)
            {
                vec[i] = pr;
                i++;
            }
        }
        
        return vec;
    }

    private static BigDecimal divideBigDecimals(BigDecimal dividend, BigDecimal divisor)
    {
        if(mc.equals(MathContext.UNLIMITED))
        {
            return dividend.divide(divisor,BigDecimal.ROUND_HALF_UP);
        }
        
        return dividend.divide(divisor,mc);
    }
    
    
    /*
    private static int computeInitialEtaDegreeProbabilities(int v)
    {
        int dv = getDegree(v);
        double[] non1prob = getNon1ProbsDouble(v);
        BigDecimal[] non1probBig = getNon1Probs(non1prob);
        int dvPrime = non1prob.length;//number of non-one probabilities
                
        BigDecimal[] p = new BigDecimal[dvPrime+1];
        if(p.length == 1)
        {
            p[0] = one;
        }
        else
        {
            p[0] = one.subtract(non1probBig[0],mc);
            p[1] = non1probBig[0];

            for(int i=1; i<dvPrime; i++)
            {
                BigDecimal pe = non1probBig[i];
                BigDecimal oneMinusPe = one.subtract(pe,mc);
                
                BigDecimal[] newP = new BigDecimal[p.length];
                //j=0
                newP[0] = oneMinusPe.multiply(p[0],mc);
                //j \in [1..i]
                for(int j=1; j<=i; j++)
                {                   
                    BigDecimal oldPrJminusOne = p[j-1];
                    BigDecimal oldPrJ = p[j];
                    
                    BigDecimal prod1 = oldPrJminusOne.multiply(pe,mc);
                    BigDecimal prod2 = oldPrJ.multiply(oneMinusPe,mc);
                    
                    newP[j] = prod1.add(prod2,mc);
                }
                //j=i+1
                newP[i+1] = pe.multiply(p[i],mc);
                
                p = newP;
            }
        }
        
        BigDecimal[] realP = new BigDecimal[dv+1];
        for(int i=0; i<dv-dvPrime; i++)
        {
            realP[i] = zero;
        }
        for(int i=dv-dvPrime; i<realP.length; i++)
        {
            realP[i] = p[i-(dv-dvPrime)];
        }
        
        etaDegProbs.put(v, realP);
        
        return dv-dvPrime;
    }   
    */
    
    private static BigDecimal computeInitialEtaDegreeProbability(int v, int k, BigDecimal[] non1prob, BigDecimal[] X)
    {                
        BigDecimal[] Xnew = new BigDecimal[X.length];
        for(int h=0; h<Xnew.length; h++)
        {
            if(h==0)
            {
                if (k==0)
                {
                    Xnew[h] = one;
                }
                else
                {
                    Xnew[h] = zero;
                }
            }
            else if(k>h)
            {
                Xnew[h] = zero;
            }
            else
            {
                BigDecimal peh = non1prob[h-1];
                BigDecimal oneMinusPeh = one.subtract(peh,mc);
                
                BigDecimal bd1 = X[h-1].multiply(peh,mc);
                BigDecimal bd2 = Xnew[h-1].multiply(oneMinusPeh,mc);
                
                Xnew[h] = bd1.add(bd2,mc);
            }
        }
        
        for(int h=0; h<X.length; h++)
        {
            X[h] = Xnew[h];
        }
        
        return X[X.length-1];
    }

    private static void loadBetaFunctionValuesAndComputeEtaDegLBs(String dataset, double eta) throws IOException
    {
        System.out.println("--------------------------------------------------");
        System.out.print("Loading Beta-function values & computing lower bounds on eta-degrees...");
        long start = System.currentTimeMillis();
        
        betaFunctionValues = new double[n][];
        etaDegs = new int[n];
        hasExactEtaDegBeenComputed = new boolean[n];
        for(int i=0; i<hasExactEtaDegBeenComputed.length; i++)
        {
            hasExactEtaDegBeenComputed[i] = false;
        }
        
        String f = matlabOutputPath+dataset+matlabOutputSuffix;
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line = br.readLine();
        while (line != null)
        {
            //compute LB
            StringTokenizer st = new StringTokenizer(line, " \t");
            int v = Integer.parseInt(st.nextToken());
            double pvmin = Double.parseDouble(st.nextToken());
            int dv = Integer.parseInt(st.nextToken());
            
            if (dv == 1)
            {
                double p = g.get(v).values().iterator().next();
                int etaDeg = (p>eta)?1:0;
                etaDegs[v] = etaDeg;
                hasExactEtaDegBeenComputed[v] = true;
            }
            else
            {
                int k=0;
                boolean stop = false;
                while(!stop && k<=dv)
                {
                    double b = Double.parseDouble(st.nextToken());
                    if(b <= eta)
                    {
                        stop = true;
                    }
                    else if(!st.hasMoreTokens())
                    {
                        stop = true;
                        k++;
                    }
                    else
                    {
                        k++;
                    }
                }
                int lbv = k-1;//lb on eta-degree of vertex v
                if(lbv < 0)
                {
                    lbv = 0;
                }
                etaDegs[v] = lbv;
            }
            //load needed Beta values
            line = br.readLine();
            st = new StringTokenizer(line, " \t");
            st.nextToken();
            st.nextToken();
            st.nextToken();
            
            double[] bv = new double[etaDegs[v]+1];
            for(int k=0; k<=etaDegs[v]; k++)
            {
                try
                {
                    bv[k] = Double.parseDouble(st.nextToken());
                }
                catch(java.util.NoSuchElementException e)
                {
                    int x = 0;
                }
            }
            betaFunctionValues[v] = bv;           
            
            line = br.readLine();
        }

        double time = System.currentTimeMillis() - start;
        System.out.println("DONE! ("+(time/1000)+"secs)");
        System.out.println("--------------------------------------------------");
 
    }

    private static int updateEtaDegLB(int u, int prevEtaDeg, double eta) 
    {
        if(prevEtaDeg == 0)
        {
            return 0;
        }
        
        
        /*
        int index = getDegree(u)-prevEtaDeg;
        if(index >= betaFunctionValues[u].length)
        {
            System.out.println("u: "+u);
            System.out.println("dv: "+getDegree(u));
            System.out.println("prev eta-deg: "+prevEtaDeg);
            System.out.println("dv - prev eta-deg: "+index);
            System.out.println("array length: "+betaFunctionValues[u].length);       
        }
        
        double pr = betaFunctionValues[u][getDegree(u)-prevEtaDeg];
        if (pr > eta)
        {
            return prevEtaDeg;
        }
        */
        
        return prevEtaDeg-1;
    }
}
