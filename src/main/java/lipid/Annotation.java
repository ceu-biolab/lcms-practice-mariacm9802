package lipid;

import adduct.Adduct;
import adduct.AdductList;

import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IoniationMode ionizationMode;
    private String adduct; // !!TODO The adduct will be detected based on the groupedSignals
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;


    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        // !!TODO This set should be sorted according to help the program to deisotope the signals plus detect the adduct
        this.groupedSignals = new TreeSet<>(Comparator.comparing(Peak::getMz));
        this.groupedSignals.addAll(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
        detectAdductFromPeaks();
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IoniationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }


    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    // !CHECK Take into account that the score should be normalized between -1 and 1
    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    /**
     * @return The normalized score between 0 and 1 that consists on the final number divided into the times that the rule
     * has been applied.
     */
    public double getNormalizedScore() {
        return (double) this.score / this.totalScoresApplied;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }


    // !!TODO Detect the adduct with an algorithm or with drools, up to the user.

    /**
     * Attempts to infer the most likely adduct for this annotation's mz value
     * by comparing grouped peaks and matching their inferred monoisotopic masses
     * within a 10 ppm tolerance. Only non-multimeric positive adducts are considered.
     *
     */
    public void detectAdductFromPeaks() {
        List<Peak> peakList = new ArrayList<>(this.groupedSignals);

        for (Peak p1 : peakList) {
            if (Math.abs(p1.getMz() - this.mz) < 0.01) {
                for (Map.Entry<String, Double> adduct1 : AdductList.MAPMZPOSITIVEADDUCTS.entrySet()) {

                    // Evita aductos multiméricos como [2M+H]+
                    if (Adduct.extractMultimer(adduct1.getKey()) > 1) continue;

                    double m1 = Adduct.getMonoisotopicMassFromMZ(p1.getMz(), adduct1);

                    for (Peak p2 : peakList) {
                        if (p1.equals(p2)) continue;

                        for (Map.Entry<String, Double> adduct2 : AdductList.MAPMZPOSITIVEADDUCTS.entrySet()) {
                            double m2 = Adduct.getMonoisotopicMassFromMZ(p2.getMz(), adduct2);
                            int ppmDiff = Adduct.calculatePPMIncrement(m1, m2);
                            int ppmTolerance = 10;

                            if (ppmDiff <= ppmTolerance) {
                                this.adduct = adduct1.getKey(); // asignar el aducto asociado a this.mz
                                return;
                            }
                        }
                    }
                }
            }
        }
    }
}
