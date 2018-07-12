package edu.unc.genomics.nucleosomes;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.apache.log4j.Logger;

import com.beust.jcommander.Parameter;

import edu.unc.genomics.CommandLineToolException;
import edu.unc.genomics.Contig;
import edu.unc.genomics.Interval;
import edu.unc.genomics.ReadablePathValidator;
import edu.unc.genomics.WigMathTool;
import edu.unc.genomics.io.WigFileReader;
import edu.unc.genomics.io.WigFileException;
import edu.unc.utils.FAIREModel;

/**
 * Attempt to predict FAIRE signal from nucleosome occupancy data using a simple
 * probabilistic model
 * 
 * @author timpalpant
 * 
 */
public class PredictFAIRESignal extends WigMathTool {

  private static final Logger log = Logger.getLogger(PredictFAIRESignal.class);

  @Parameter(names = { "-i", "--input" }, description = "Nucleosome occupancy (Wig)", required = true, validateWith = ReadablePathValidator.class)
  public Path inputFile;
  @Parameter(names = { "-s", "--sonication" }, description = "Sonication distribution", required = true, validateWith = ReadablePathValidator.class)
  public Path sonicationFile;
  @Parameter(names = { "-e", "--efficiency" }, description = "FAIRE crosslinking efficiency (0,1)")
  public float efficiency = 0.01f;
  @Parameter(names = { "-x", "--extend" }, description = "Single-end read extension (bp); -1 for paired-end")
  public int extend = -1;

  WigFileReader reader;
  float[] sonication = new float[100];
  double maxValue;

  public static float[] loadSonication(Path s) {
    log.debug("Loading sonication fragment length distribution");
    float[] sonication = new float[100];
    int maxL = 0;
    float total = 0;
    try (BufferedReader reader = Files.newBufferedReader(s, Charset.defaultCharset())) {
      String line;
      while ((line = reader.readLine()) != null) {
        // Parse the line
        String[] entry = line.split("\t");
        if (entry.length != 2) {
          throw new CommandLineToolException("Invalid format for sonication distribution file (length\tpercent)");
        }
        int length = Integer.parseInt(entry[0]);
        float percent = Float.parseFloat(entry[1]);
        // Expand the sonication distribution array if necessary
        if (length >= sonication.length) {
          sonication = Arrays.copyOf(sonication, Math.max(sonication.length + 100, length + 1));
        }
        if (length > maxL) {
          maxL = length;
        }
        sonication[length] = percent;
        total += percent;
      }
    } catch (IOException e) {
      log.fatal("Error loading sonication fragment length distribution");
      e.printStackTrace();
      throw new CommandLineToolException("Error loading sonication fragment length distribution");
    }

    log.debug("Longest fragment length: " + maxL + "bp");
    // Truncate the array to the minimum possible size
    sonication = Arrays.copyOfRange(sonication, 0, maxL + 1);

    // Normalize the sonication distribution so that it has total 1
    for (int i = 0; i < sonication.length; i++) {
      sonication[i] /= total;
    }

    return sonication;
  }

  @Override
  public void setup() {
    try {
      reader = WigFileReader.autodetect(inputFile);
    } catch (IOException e) {
      throw new CommandLineToolException(e);
    }
    addInputFile(reader);
    maxValue = reader.max();
    sonication = loadSonication(sonicationFile);
  }

  @Override
  public float[] compute(Interval chunk) throws IOException, WigFileException {
    int paddedStart = Math.max(chunk.getStart() - sonication.length, reader.getChrStart(chunk.getChr()));
    int paddedStop = Math.min(chunk.getStop() + sonication.length, reader.getChrStop(chunk.getChr()));

    Contig data = reader.query(chunk.getChr(), paddedStart, paddedStop);
    float[] pOcc = data.get(chunk.getStart() - sonication.length, chunk.getStop() + sonication.length);
    // Scale the occupancy by the maximum so that it represents
    // the probability that a nucleosome occupies that base pair
    // You should probably remove outliers (esp. CNVs) first
    for (int i = 0; i < pOcc.length; i++) {
      pOcc[i] /= maxValue;
    }

    FAIREModel model = new FAIREModel(pOcc, sonication, efficiency);
    float[] prediction;
    if (extend > 0) {
      prediction = model.singleEnd(extend);
    } else {
      prediction = model.pairedEnd();
    }

    return Arrays.copyOfRange(prediction, sonication.length, prediction.length - sonication.length);
  }

  /**
   * @param args
   * @throws WigFileException
   * @throws IOException
   */
  public static void main(String[] args) throws IOException, WigFileException {
    new PredictFAIRESignal().instanceMain(args);
  }

}
