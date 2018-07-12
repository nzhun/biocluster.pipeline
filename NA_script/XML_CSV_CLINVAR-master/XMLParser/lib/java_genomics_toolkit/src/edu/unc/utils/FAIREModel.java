package edu.unc.utils;

import org.apache.log4j.Logger;

/**
 * Predict FAIRE from nucleosome maps
 * 
 * @author timpalpant
 *
 */
public class FAIREModel {

  private static final Logger log = Logger.getLogger(FAIREModel.class);

  private final float[] pOcc, sonication;
  private final float efficiency;

  /**
   * Model FAIRE signal from nucleosome maps
   *
   * pOcc must be normalized so that all values represent probabilities in [0,1]
   * sonication must be normalized so that the total == 1
   * 
   */
  public FAIREModel(final float[] pOcc, final float[] sonication, final float efficiency) {
    if (efficiency <= 0 || efficiency >= 1) {
      throw new IllegalArgumentException("Crosslinking efficiency must be in (0,1)");
    }

    this.pOcc = pOcc;
    this.sonication = sonication;
    this.efficiency = efficiency;
  }

  /**
   * Single-end prediction, with artificial uniform extension
   * 
   * @param extend
   *          artificial extension (read length)
   * @return single-end FAIRE signal
   */
  public float[] singleEnd(int extend) {
    if (extend < 1) {
      throw new IllegalArgumentException("Read extension must be >= 1");
    }

    log.debug("Calculating single-end FAIRE prediction");
    float[] watson = new float[pOcc.length];
    float[] crick = new float[pOcc.length];
    // Consider all possible fragment lengths
    for (int l = 1; l < sonication.length; l++) {
      // No need to count if there are no fragments of this length
      if (sonication[l] == 0) {
        continue;
      }

      // Starting at each base pair in the chunk
      for (int x = 0; x < pOcc.length - l; x++) {
        // Calculate the probability that this fragment survives FAIRE
        // and add its probability to the prediction,
        // weighted by the sonication distribution
        float pFAIRE = sonication[l] * (1 - pCrosslinked(x, l));
        watson[x] += pFAIRE;
        crick[x + l - 1] += pFAIRE;
      }
    }

    log.debug("Extending watson fragments");
    float[] prediction = new float[pOcc.length];
    for (int i = 0; i < pOcc.length - extend; i++) {
      for (int j = 0; j <= extend; j++) {
        prediction[i + j] += watson[i];
      }
    }

    log.debug("Extending crick fragments");
    for (int i = pOcc.length - 1; i >= extend; i--) {
      for (int j = 0; j <= extend; j++) {
        prediction[i - j] += crick[i];
      }
    }

    return prediction;
  }

  /**
   * Paired-end, extend to the actual length of each fragment
   * 
   * @return paired-end FAIRE signal for the chunk
   */
  public float[] pairedEnd() {
    log.debug("Calculating paired-end FAIRE prediction");
    float[] prediction = new float[pOcc.length];
    // Consider all possible fragment lengths
    for (int l = 1; l < sonication.length; l++) {
      // No need to count if there are no fragments of this length
      if (sonication[l] == 0) {
        continue;
      }

      // Starting at each base pair in the chunk
      for (int x = 0; x < pOcc.length - l; x++) {
        // Calculate the probability that this fragment survives FAIRE
        // and add its probability to the prediction,
        // weighted by the sonication distribution
        float pFAIRE = sonication[l] * (1 - pCrosslinked(x, l));
        for (int k = 0; k < l; k++) {
          prediction[x + k] += pFAIRE;
        }
      }
    }

    return prediction;
  }

  /**
   * The probability that a fragment of length l, starting at x, is crosslinked
   * at least ncrosslinks times based on the occupancy profile in pOcc
   *
   * For each base b \in [x, x+l): pCrosslinked = U_b ( pOcc(b) and
   * crosslinked(b) ) This model assumes uniform, independent crosslinking at a
   * constant efficiency.
   *
   * @param x
   *          the first base pair of the fragment
   * @param l
   *          the length of the fragment
   * @return single-end FAIRE signal
   */
  public float pCrosslinked(int x, int l) {
    float p = 0;
    // Union of independent events
    for (int b = x; b < x + l; b++) {
      p += efficiency * pOcc[b] * (1 - p);
    }
    return p;
  }

}
