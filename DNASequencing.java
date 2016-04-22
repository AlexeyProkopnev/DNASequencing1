import java.util.*;

import static java.lang.Math.*;

//If you submit, comment the following line
import java.io.*;

public class DNASequencing {

    public int initTest(int _testDifficulty) {

        return 0;
    }

    public int preProcessing() {

        return 0;
    }

	public int passReferenceGenome(int chromatidSequenceId, String[] chromatidSequence) {

        return 0;
	}

    String[] getAlignment(int N, double normA, double normS, String[] readName, String[] readSequence) {

        String[] ret = new String[N];

        for (int i = 0; i < N; ++i)
            ret[i] = readName[i] + ",20,1,150,+,0";

        return ret;
    }

    public static void main(String[] args) {

        int testDifficulty = 0;

        //If you submit, comment the following line
        Tester.run(testDifficulty);    
    }
}

//BEGIN TESTER CODE

/*
 * Tester code uses operations with files, so remove it before submitting.
 * 
 * Server message: "Compilation failed because there is no need to make use of 'java.io.File' anywhere in your code, for security purposes."
 *
 **/
class Tester {
    
    /* Constants from the problem statement */

    public static final int MAX_POSITION_DIST = 300;

    public static final double NORM_A_SMALL = -3.392;
    public static final double NORM_A_MEDIUM = -3.962;
    public static final double NORM_A_LARGE = -2.710;
    public static final double MAX_AUC = 0.999999;

    public static final double NORM_S = 0.5;

    /**
     * Read a minisam file and build a map of ground truth
     * @param path	the path of the minisam file storing the ground truth 
     * @return a map[read_name] = read_Position
     */
    static Map<String, Position> parseTruth(String path) throws IOException {

        Map<String, Position> res = new HashMap<String, Position>();

        InputReader reader = new InputReader(new FileInputStream(path));

        while (!reader.isEOF()) {

            String[] tokens = reader.next().split(",");
            
            String qname = tokens[0];
            int chromatid = Integer.parseInt(tokens[1]);
            int from = Integer.parseInt(tokens[2]);
            int to = Integer.parseInt(tokens[3]);
            char strand = tokens[4].charAt(0);
            res.put(qname, new Position(chromatid, from, to, strand));
        }

        return res;
    }

    /**
     * For each string of the results vector, build a read result {confidence, r}
     * @param truth		the map of ground truth position for each read
     * @param results	the vector of results as return by getAlignment
     * @return a vector of ReadResult, that is {confidence, r}
     */
    static List<ReadResult> buildReadResults(Map<String, Position> truth, String[] results) {

        int n = results.length;
        int correct = 0;
        List<ReadResult> readResults = new ArrayList<ReadResult>(n);

        for (String res: results) {

            String[] tokens = res.split(",");
            Position position = truth.get(tokens[0]);
            int r = 1;
            r = (Integer.parseInt(tokens[1]) == position.rname) ? r : 0;
            r = (tokens[4].charAt(0) == position.strand) ? r : 0;
            int start0 = Integer.parseInt(tokens[2]);
            int start1 = position.from;
            r = (abs(start0-start1) < MAX_POSITION_DIST) ? r : 0;
            double confidence = Double.parseDouble(tokens[5]);
            readResults.add(new ReadResult(confidence, r));
            correct += r;
        }

        System.err.printf("Number of correct answers: %d/%d = %.6f\n", correct, n, correct / (double)n);

        return readResults;
    }

    /**
     * Compute the accuracy given the {confidence, r} pairs and the normalization facto
     * @param readResults	a vector of {confidence, r} results
     * @param normA		as described in the problem statement
     * @return	a double, the computed accuracy
     */
    static double computeAccuracy(List<ReadResult> readResults, double normA) {

        int n = readResults.size();

        Collections.sort(readResults);

        // merge results of equal confidence
        List<Integer> cumulSi = new ArrayList<Integer>();
        List<Integer> pos = new ArrayList<Integer>();

        cumulSi.add(readResults.get(0).r);
        pos.add(0);

        for (int i = 1; i < n; ++i) {

            int back = cumulSi.size()-1;
            if (readResults.get(i).confidence == readResults.get(i-1).confidence) {

                cumulSi.set(back, cumulSi.get(back)+readResults.get(i).r); 
                pos.set(back, i);
            }
            else {

                int cumul = cumulSi.get(back) + readResults.get(i).r;
                cumulSi.add(cumul);
                pos.add(i);
            }
        }

        //compute the AuC
        double auc = 0.0;
        double invn = 1.0 / n;
        double invnp1 = 1.0 / (n+1);
        double lfmultiplier = 1.0 / log(n+1);
        int m = cumulSi.size();

        for (int i = 0; i < m; ++i) {

            double fi = 1.0 * (2 + pos.get(i) - cumulSi.get(i)) * invnp1;
            double fi1 = (i == m-1) ? 1.0 : 1.0 * (2 + pos.get(i+1) - cumulSi.get(i+1)) * invnp1;
            double lfi = lfmultiplier * log(fi);
            double lfi1 = lfmultiplier * log(fi1);
            auc += cumulSi.get(i) * (lfi1 - lfi) * invn;
        }

        System.err.printf("auc = %.6f\n", auc);
        double tmp = log(1 - min(auc, MAX_AUC));
        System.err.printf("log(1 - min(auc, MAX_AUC)) = %.6f\n", tmp); 
        System.err.printf("NormA = %.6f\n", normA);
        double accuracy = tmp / normA;
        System.err.printf("accuracy = %.6f\n", accuracy);

        return accuracy;
    }

    /**
     * Perform a single test
     * @param testDifficulty	define the test type (SMALL=0, MEDIUM=1, LARGE=2)
     * @return	alignments in format specified in the problem statement
     */
    static String[] performTest(int testDifficulty, double normA) throws IOException {

        // test data path and description
        String fa1Path = null, fa2Path = null;
        int[] chrIds = null;

        switch (testDifficulty) {
            
            case 0: {
                fa1Path = "../data/small5.fa1";
                fa2Path = "../data/small5.fa2";
                chrIds = new int[] {20}; 
                break;
            }

            case 1: {
                fa1Path = "../data/medium5.fa1";
                fa2Path = "../data/medium5.fa2";
                chrIds = new int[] {1,11,20};
                break;
            }

            case 2: {
                fa1Path = "../data/large5.fa1";
                fa2Path = "../data/large5.fa2";
                chrIds = new int[24];
                for (int i = 1; i <= 24; ++i) 
                    chrIds[i-1] = i;
                break;
            }
        }

        // call the MM DNASequencing methods
        DNASequencing dnaSequencing = new DNASequencing();

        dnaSequencing.initTest(testDifficulty);

        // load chromatid	
        for (int chromatidSequenceId: chrIds) {
            
            List<String> chromatidSequence = new ArrayList<String>();
            String path = "../data/chromatid" + chromatidSequenceId + ".fa";
            InputReader reader = new InputReader(new FileInputStream(path));
            // skip header
            String s = reader.next();
            System.err.println("Skip header: " + s);
            while (!reader.isEOF()) {
                chromatidSequence.add(reader.next());
            }

            dnaSequencing.passReferenceGenome(chromatidSequenceId, chromatidSequence.toArray(new String[chromatidSequence.size()]));		
            chromatidSequence.clear();
        } 

        dnaSequencing.preProcessing();

        // load reads
        List<String> readId = new ArrayList<String>(); 
        List<String> readSeq = new ArrayList<String>();

        InputReader fa1 = new InputReader(new FileInputStream(fa1Path));
        InputReader fa2 = new InputReader(new FileInputStream(fa2Path));

        while (!fa1.isEOF() && !fa2.isEOF()) {
            readId.add(fa1.next().substring(1));
            readId.add(fa2.next().substring(1));
            readSeq.add(fa1.next());
            readSeq.add(fa2.next());
        }

        int nReads = readId.size();
        
        String[] arrayIds = readId.toArray(new String[nReads]);
        readId.clear();

        String[] arraySeq = readSeq.toArray(new String[nReads]);
        readSeq.clear();

        // compute alignments
        String[] results = dnaSequencing.getAlignment(nReads, normA, NORM_S, arrayIds, arraySeq);

        return results;
    }

    public static void run(int testDifficulty) {

        try {

            long start = System.currentTimeMillis();

            String minisamPath = null;
            double normA = 0;

            switch (testDifficulty) {

                case 0: {
                    minisamPath = "../data/small5.minisam";
                    normA = NORM_A_SMALL;
                    break;
                } 
                case 1: {
                    minisamPath = "../data/medium5.minisam";
                    normA = NORM_A_MEDIUM;
                    break;
                } 
                case 2: {
                    minisamPath = "../data/large5.minisam";
                    normA = NORM_A_LARGE;
                    break;
                }
            }

            // perform test
            String[] results = performTest(testDifficulty, normA);
            // load truth
            Map<String, Position> truth = parseTruth(minisamPath);	
            List<ReadResult> readResults = buildReadResults(truth, results);
            // scoring
            double accuracy = computeAccuracy(readResults, normA);

            System.err.printf("\nExecution time: %.3fsec\n", (System.currentTimeMillis() - start) / 1000.0);
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

}

/**
 * Position: describe the position of a read within the genome
 */
class Position {

    public Position(int rname, int from, int to, char strand) {
        this.rname = rname;
        this.from = from;
        this.to = to;
        this.strand = strand;
    }

    int rname;
    int from;
    int to;
    char strand;
}

/**
 * ReadResult: result of a read alignment
 */
class ReadResult implements Comparable<ReadResult> {

    public ReadResult(double confidence, int r) {
        this.confidence = confidence;
        this.r = r;
    }
    
    @Override
    public int compareTo(ReadResult other) {

        return confidence > other.confidence ? -1 : confidence < other.confidence ? 1 : 0;
    }

    double confidence;
    int r;
}

class InputReader {

    private static final int MAX_LINE_LENGTH = 1 << 8; 
    private static final int BUFFER_SIZE = (1 << 14)-1;

    private InputStream stream;

    private int cursor = BUFFER_SIZE;

    private byte buffer[] = new byte[BUFFER_SIZE];


    public InputReader(InputStream stream) throws IOException {

        this.stream = stream;

        readAhead();
    }

    private void readAhead() throws IOException {

        int remaining = BUFFER_SIZE^cursor;
        if (remaining < MAX_LINE_LENGTH) {
            System.arraycopy(buffer, cursor, buffer, 0, remaining);
            remaining += stream.read(buffer, remaining, cursor);
            cursor = 0;
            if (remaining < BUFFER_SIZE)
                Arrays.fill(buffer, remaining, BUFFER_SIZE, (byte) 0);
        }
    }

    private byte read() {
        return buffer[cursor++];
    }

    public String next() throws IOException {

        readAhead();

        int from = cursor;
        byte c = read();

        while (c != '\r' && c != '\n' && c != 0) c = read();

        String ret = new String(buffer, from, cursor-from-1);

        if (buffer[cursor] == '\n') ++cursor;

        return ret;
    }

    public boolean isEOF() {
        return buffer[cursor] == 0;
    }
}

//END TESTER CODE

