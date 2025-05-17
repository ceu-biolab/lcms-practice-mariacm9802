package adduct;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Adduct {

    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, Map.Entry<String, Double> adduct) {

        Double monoisotopicMass = null;

        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).
        int multimer = extractMultimer(adduct.getKey());
        int charge = extractCharge(adduct.getKey());

        //if charge>1 for an adduct, shift has to be divided by the charge
        Double adductMass = adduct.getValue() / charge;

        // Case 1: Single charge, no multimer (m/z = M +- adductMass).
        if (charge == 1 && multimer == 1) {
            monoisotopicMass = mz + adductMass;
        }
        // Case 2: Multiple charges, no multimer (mz = M/charge +- adductMass).
        else if (charge > 1 && multimer == 1) {
            monoisotopicMass = (mz + adductMass) * charge;
        }
        // Case 3: Multimer, single charge (mz = M * numberOfMultimer +- adductMass).
        else if (charge == 1 && multimer > 1) {
            monoisotopicMass = (mz + adductMass) / multimer;
        }
        // Case 4: Multimer with multiple charges (mz = (M * numberOfMultimer +- adductMass) / charge).
        else {
            monoisotopicMass = ((mz + adductMass) * charge) / multimer;
        }
        return monoisotopicMass;
    }

    /**
     * Extracts the charge of an adduct using a regex pattern.
     *
     * @param adduct
     * @return
     */
    public static int extractCharge(String adduct) {
        // Match the last digit(s) followed by + or − before the closing bracket
        Matcher m = Pattern.compile("([0-9]*)([+-])\\]?$").matcher(adduct);
        if (m.find()) {
            String num = m.group(1);  // May be empty
            return num.isEmpty() ? 1 : Integer.parseInt(num);
        }
        return 1; // Default if no explicit charge found
    }

    /**
     * Extracts the multimer number from the adduct using a regex pattern.
     *
     * @param adduct
     * @return
     */
    public static int extractMultimer(String adduct) {
        // Pattern looks for a number before 'M', e.g. [2M or [M
        Matcher m = Pattern.compile("\\[([0-9]*)M").matcher(adduct);
        if (m.find()) {
            String num = m.group(1); // Group 1 captures the digits before 'M'
            return num.isEmpty() ? 1 : Integer.parseInt(num);
        }
        return 1; // Default: no multimer specified
    }


    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoisotopicMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct) {
        // Obtener la masa del aducto desde AdductList
        Double adductMass = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        if (adductMass == null) {
            throw new IllegalArgumentException("Unknown adduct: " + adduct);
        }
        // Regex para extraer multímero y carga
        Pattern pattern = Pattern.compile("\\[([0-9]*)?M[+-][^\\]]+](\\d*)[+-]?");
        Matcher matcher = pattern.matcher(adduct);

        int multimer = 1;
        int charge = 1;

        if (matcher.matches()) {
            String multimerStr = matcher.group(1);
            String chargeStr = matcher.group(2);

            if (multimerStr != null && !multimerStr.isEmpty()) {
                multimer = Integer.parseInt(multimerStr);
            }

            if (chargeStr != null && !chargeStr.isEmpty()) {
                charge = Integer.parseInt(chargeStr);
            }
        }

        // Aplicar fórmula
        double mass = monoisotopicMass * multimer;
        return (mass + adductMass) / charge;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     * @return the absolute difference un ppm (rounded)
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass) * 1000000 / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the absolute delta (Da) corresponding to the given ppm tolerance for a mass
     *
     * @param experimentalMass Mass measured by MS
     * @param ppm              Tolerance in ppm
     * @return Absolute mass difference (Da) that corresponds to the ppm tolerance
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        deltaPPM =  Math.round(Math.abs((experimentalMass * ppm) / 1000000));
        return deltaPPM;

    }

}
