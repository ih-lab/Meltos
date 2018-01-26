package lineage;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import org.apache.commons.math3.util.ArithmeticUtils;
import static org.apache.commons.math3.util.ArithmeticUtils.factorial;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class SSNVIO {
	public static ArrayList<SVEntry> loadSVFile(String filename, int numSamples) {
		ArrayList<SVEntry> svs = new ArrayList<SVEntry>();
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
			System.out.println("Read " + svs.size() + " SSNVs from the input file. \n");
		} catch (IOException e){
			returnInputFileFormatError("Could not read file: " + filename, null);
		}
                
                //System.out.println(svs.get(2).getVAF(4));
		return svs;
	}
	
        
	final static int NUM_SV_METADATA_FIELDS = 4;
	public static SVEntry parseSVEntry(String line, int numSamples) {
		String[] entryParts = line.split("\t");
		if(entryParts.length != (numSamples + NUM_SV_METADATA_FIELDS)) {
			returnInputFileFormatError("Expecting the SV description field and " + numSamples + " samples", line);
		} 
		SVEntry entry = new SVEntry(numSamples);
		entry.id = entryParts[0].trim();
		entry.description = entryParts[1].trim();
                //entry.approxVAF[0]=0;
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
        entry.description="snv";    
		return entry;
	}
	
        
        
        
	public static void returnInputFileFormatError(String desc, String entry) {
		System.err.println("[Wrong SV input file format] " + desc);
		if(entry != null) {
			System.err.println("Line: " + entry);
		}
		System.exit(1);
	}
        
        public static float getlambda_n(int t, long g, int l_avg, float k) {
            float estimate = (float) t/g*l_avg*(1-k);
            
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
            PoissonDistribution pdist = new PoissonDistribution(lambda);
            double prob = pdist.probability(obs);
            //double prob2 = Math.exp(-lambda) * (Math.pow(lambda,obs)/factorial(obs));
          
            return prob;
        }
        
        
}
