package lineage;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import org.apache.commons.math3.util.ArithmeticUtils;
import static org.apache.commons.math3.util.ArithmeticUtils.factorial;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class SVIO {
	String[] orProfs;
	public static ArrayList<SVEntry> loadSVFile(String filename, int numSamples, boolean calcVAF) {
		ArrayList<SVEntry> svs = new ArrayList<SVEntry>();
                
                if(!calcVAF) {
                    try {
                        BufferedReader rd = new BufferedReader(new FileReader(filename));
                        String currLine = rd.readLine();
                        currLine = rd.readLine();
                        while (currLine != null) {
                            SVEntry entry = parseSVEntry(currLine, numSamples);
                            
                            svs.add(entry);
                            currLine = rd.readLine();
                        }
                        rd.close();
                        System.out.println("Read " + svs.size() + " SVs from the input file. \n");
                    } catch (IOException e){
                        returnInputFileFormatError("Could not read file: " + filename, null);
                    }

                    //System.out.println(svs.get(2).getVAF(4));
                    
                } else {
                    try {
                    	FileWriter fwV = new FileWriter("testdata/estimatedVAFs.txt");
                        BufferedReader rd = new BufferedReader(new FileReader(filename));
                        String currLine = rd.readLine();
                        currLine = rd.readLine();
                        //currLine = rd.readLine();
                        while(currLine != null) {
                            SVEntry entry = parseSVEntry_vaf(currLine, numSamples, fwV);

                            svs.add(entry);
                            currLine = rd.readLine();
                        }
                        rd.close();
                        fwV.close();
                    } catch (IOException e){
                        returnInputFileFormatError("Could not read file: " + filename, null);
                    }
                }
                return svs;
	}
	
        
	final static int NUM_SV_METADATA_FIELDS = 2;
	public static SVEntry parseSVEntry(String line, int numSamples) {
		String[] entryParts = line.split("\t");
		if(entryParts.length != (numSamples + NUM_SV_METADATA_FIELDS)) {
			returnInputFileFormatError("Expecting the SV description field and " + numSamples + " samples", line);
		} 
		SVEntry entry = new SVEntry(numSamples);
		entry.id = entryParts[0].trim();
		entry.description = entryParts[1].trim();
                
                for(int i = NUM_SV_METADATA_FIELDS; i < entryParts.length; i++) {
			try {
				entry.approxVAF[i - NUM_SV_METADATA_FIELDS] = Double.parseDouble(entryParts[i]);
				//entry.readDepth[i - NUM_SV_METADATA_FIELDS] = Double.parseDouble(entryParts[i + 1]);
			} catch (NumberFormatException e) {
				returnInputFileFormatError("VAF value " + entryParts[i], line);
			} 
                    
                }
                entry.makeProfile(Parameters.MAX_VAF_ABSENT,Parameters.MIN_VAF_PRESENT);
                
                
                /*
		for(int i = NUM_SV_METADATA_FIELDS; i < entryParts.length; i += 2) {
			try {
				entry.readSupport[i - NUM_SV_METADATA_FIELDS] = Integer.parseInt(entryParts[i]);
				entry.readDepth[i - NUM_SV_METADATA_FIELDS] = Double.parseDouble(entryParts[i + 1]);
			} catch (NumberFormatException e) {
				returnInputFileFormatError("VAF value " + entryParts[i + 1], line);
			}
		}
		entry.approximateVAFs();
                */
                
		return entry;
	}
        
	public static SVEntry parseSVEntry_vaf(String line, int numSamples, FileWriter fwV) throws IOException{
		String[] entryParts = line.split("\t");
		if(entryParts.length != (numSamples*2 + 2)) {
			returnInputFileFormatError("Expecting the SV description field and " + numSamples + " samples", line);
		} 
		SVEntry entry = new SVEntry(numSamples);
		entry.id = entryParts[0].trim();
		entry.description = entryParts[1].trim();
                
		int idx = 0;
        for(int i = NUM_SV_METADATA_FIELDS; i < entryParts.length; i++) {
            try {
                i = i+1;
                if(!entryParts[i].equals("0")) {
                    String[] counts = entryParts[i].split(";");
                    System.out.println(entryParts[i]);
                    double vaf = estimateVAF(Integer.parseInt(counts[1]), Integer.parseInt(counts[2]), Integer.parseInt(counts[3]), Integer.parseInt(counts[4]), Integer.parseInt(counts[5]), Integer.parseInt(counts[6]), Integer.parseInt(counts[0]), idx);
                    entry.approxVAF[idx] = vaf;
                } else {
                    entry.approxVAF[idx] = Double.parseDouble("0");
                } 

                    //entry.readDepth[i - NUM_SV_METADATA_FIELDS] = Double.parseDouble(entryParts[i + 1]);
            } catch (NumberFormatException e) {
                returnInputFileFormatError("VAF value " + entryParts[i], line);
            } 
            idx++;   
        }
                //System.exit(0);
                entry.makeProfile(Parameters.MAX_VAF_ABSENT,Parameters.MIN_VAF_PRESENT);
        fwV.write(entry.id+"\t"+Arrays.toString(entry.approxVAF)+"\n");
		return entry;
	}
	
        
        
        
	public static void returnInputFileFormatError(String desc, String entry) {
		System.err.println("[Wrong SV input file format] " + desc);
		if(entry != null) {
			System.err.println("Line: " + entry);
		}
		System.exit(1);
	}
        
        public static float getlambda_n(int t, long g, int read_l, float k) {
            //float estimate = (float) t/g*l_avg*(1-k);
            float estimate = (float) t/g*read_l*(1-k);
            
            return estimate;
        }
        
        public static float getlambda_sr(int t, long g, int read_l, float k) {
            float estimate = (float) t/g*2*read_l*k;
            
            return estimate;
        }
        
        public static float getlambda_crp(int t, long g, int window, float k) {
            float estimate = (float) t/g*window*(1-k);
            
            return estimate;
        }
        
        public static float getlambda_drp(int t, long g, int l_avg, int read_l, float k) {
            float estimate = (float) t/g*(l_avg - 2 * read_l) * k;
            
            return estimate;
        }
        
        public static double getPoissonProbability(float lambda, int obs) {
        	
        	//if(lambda<=0)
        	//	System.out.println("Problem");
            PoissonDistribution pdist = new PoissonDistribution(lambda);
            double prob = pdist.probability(obs);
            //double prob2 = Math.exp(-lambda) * (Math.pow(lambda,obs)/factorial(obs));
          
            return prob;
        }
        
        public static double estimateVAF(int n_a, int n_b, int sr, int drp, int crp, int l_avg, int SVsize, int idx) {
            
            //example values just to test calculations
            //these will be pulled from input files/bams
            //float[] k_arr = {0.1F,0.15F,0.2F,0.25F,0.3F,0.35F,0.4F,0.45F,0.5F,0.55F,0.6F,0.65F,0.7F,0.75F,0.8F,0.85F,0.9F,0.95F};
          //float[] k_arr = {0.1F,0.15F,0.2F,0.25F,0.3F,0.35F,0.4F,0.45F,0.5F,0.55F,0.6F,0.65F,0.7F,0.75F,0.8F,0.85F,0.9F,0.95F};
            //float[] k_arr = {0.01F, 0.02F, 0.03F, 0.04F, 0.05F, 0.06F, 0.07F, 0.08F, 0.09F, 0.1F, 0.11F, 0.12F, 0.13F, 0.14F, 0.15F, 0.16F, 0.17F, 0.18F, 0.19F, 0.2F, 0.21F, 0.22F, 0.23F, 0.24F, 0.25F, 0.26F, 0.27F, 0.28F, 0.29F, 0.3F, 0.31F, 0.32F, 0.33F, 0.34F, 0.35F, 0.36F, 0.37F, 0.38F, 0.39F, 0.4F, 0.41F, 0.42F, 0.43F, 0.44F, 0.45F, 0.46F, 0.47F, 0.48F, 0.49F, 0.5F};
        	float[] k_arr = {0.01F, 0.02F, 0.03F, 0.04F, 0.05F, 0.06F, 0.07F, 0.08F, 0.09F, 0.1F, 0.11F, 0.12F, 0.13F, 0.14F, 0.15F, 
                    0.16F, 0.17F, 0.18F, 0.19F, 0.2F, 0.21F, 0.22F, 0.23F, 0.24F, 0.25F, 0.26F, 0.27F, 0.28F, 0.29F, 0.3F, 0.31F, 0.32F, 
                    0.33F, 0.34F, 0.35F, 0.36F, 0.37F, 0.38F, 0.39F, 0.4F, 0.41F, 0.42F, 0.43F, 0.44F, 0.45F, 0.46F, 0.47F, 0.48F, 0.49F, 
                    0.5F, 0.51F, 0.52F, 0.53F, 0.54F, 0.55F, 0.56F, 0.57F, 0.58F, 0.59F, 0.60F, 0.61F, 0.62F, 0.63F, 0.64F, 0.65F, 0.66F, 
                    0.67F, 0.68F, 0.69F, 0.70F, 0.71F, 0.72F, 0.73F, 0.74F, 0.75F, 0.76F, 0.77F, 0.78F, 0.79F, 0.80F, 0.81F, 0.82F, 0.83F,
                    0.84F, 0.85F, 0.86F, 0.87F, 0.88F, 0.89F, 0.90F, 0.91F, 0.92F, 0.93F, 0.94F, 0.95F, 0.96F, 0.97F, 0.98F, 0.99F};

            
            //int n_a = 17;
            //int n_b =25;
            //int sr = 0;
            //int drp = 13;
            //int crp = 35;
            //int t = 907818900;
            //int t = 1171662930;
            //long g = 3137161264L;
            //int l_avg = 406;
            //int read_l = 151;
            int read_l = 100;
            //long g = 3098825702L;
            //int SVsize = 544;
            long g = 2897293955L;
            //int[] t_val = {907818900, 890042702, 907287041, 902386506, 867859429, 876372688, 867453044};
            //int[] t_val = {1186122611,1109198318,1219667862,1,1};   
            int[] t_val = {1540006722, 1540044314, 1540099738, 1539976823, 1540081622};
            
            idx = idx-1;
            int t = t_val[idx];
            //t_samA = 1186122611;
            //t_samC = 1109198318;
            //t_samD = 1219667862;

            float lambda_n, lambda_sr, lambda_drp, lambda_crp;

            ArrayList<Double> lik_arr = new ArrayList<>();

            for(float i : k_arr) {
                //System.out.println("Likelihood: " + i);
                lambda_n = getlambda_n(t,g,read_l,i);
                lambda_sr = getlambda_sr(t,g,l_avg,i);
                lambda_drp = getlambda_drp(t,g,l_avg, read_l, i);
                lambda_crp = getlambda_crp(t,g,SVsize,i);

                //System.out.println("Lambda_n: " + lambda_n);
                //double lik = getPoissonProbability(lambda_n, n_a);

                double lik = getPoissonProbability(lambda_n, n_a) * getPoissonProbability(lambda_n, n_b) * getPoissonProbability(lambda_sr,sr) * getPoissonProbability(lambda_drp,drp) * getPoissonProbability(lambda_crp,crp); 
                //System.out.println("Likelihood: " + lik);
                lik_arr.add(lik);      
            }

            //System.out.println("Max lik: " + Collections.max(lik_arr));
            int max_idx = lik_arr.indexOf(Collections.max(lik_arr)); 
            float best_k = k_arr[max_idx];
            
            return (double)best_k;
            //System.out.println("Best VAF: " + best_k + "\nLik val: " + Collections.max(lik_arr));
        }
}
