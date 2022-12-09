package lineageTree.substitutionmodel;


import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


@Description("Allows to calculate transition probabilities for a Typewriter modelled as a Poisson process on the number of edits")
public class TypewriterSubstitutionModelHomogeneous extends SubstitutionModel.Base {
    final public Input<RealParameter> rateInput = new Input<>("rate",
            "Insertion rate for the Poisson process",
            (RealParameter) null);

    /**
     * edit insertion rate  *
     */
    protected double rate;


    @Override
    public void initAndValidate() {
        super.initAndValidate();
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;
        rate = rateInput.get().getValue();


    } // initAndValidate

    protected boolean updateMatrix = true;
    private boolean storedUpdateMatrix = true;


    /**
     * This is to get transition probability between 2 sequences states (with potentially multiple edits having happened)
     */
    public double getSequenceTransitionProbability(List<Integer> start_sequence, List<Integer> end_sequence, double distance) {
        List<Integer> subtracted = new ArrayList<>(end_sequence);

        //substracting start sequence from end sequence: this is what what edited
        start_sequence.forEach(subtracted::remove);

        //upper bound for the poisson process is the number of unedited positions in the start sequence
        List<Integer> zero = Arrays.asList(0);
        start_sequence.removeAll(zero);
        int poisson_up = 5 -  start_sequence.size();

        //initialise the poisson distribution with mean rate * distance
        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(rate*distance);

        //calculate the transition probability for the case where all available positions are edited in:
        // P(max) = 1- sum(P(n))
        if(subtracted.size() == poisson_up ) {
            double sum = 0.0;
            for(int i = 0;  i<=poisson_up; i++) {
                sum += dist.probability(i);
            }

             return (1-sum) * getFrequencyFactor(subtracted);
        }

        //calculate the transition probability for the case where a #edits < avaialable positions
        else{
            double freq = getFrequencyFactor(subtracted);
            double prob = dist.probability(subtracted.size());
            return prob * freq;
        }

    }

    /**
     * Function to obtain the probability factor induced by insert frequencies
     */
    public double getFrequencyFactor(List<Integer> edits) {
        double factor = 1.0;
        for(Integer i : edits){
            factor = factor * frequencies.getFreqs()[i];
        }
        return factor;

    }


    /**
     * CalculationNode implementation follows *
     */
    @Override
    public void store() {
        storedUpdateMatrix = updateMatrix;
//        System.arraycopy(relativeRates, 0, storedRelativeRates, 0, relativeRates.length);

        super.store();
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restore() {

        updateMatrix = storedUpdateMatrix;
// To restore all this stuff just swap the pointers...
//        double[] tmp1 = storedRelativeRates;
//        storedRelativeRates = relativeRates;
//        relativeRates = tmp1;
        super.restore();

    }

    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        updateMatrix = true;
        return true;
    }

    /**
     * this not used in the current implementation because we do not use a matrix format
     */
    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // not implemented

    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return null;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType.getStateCount() != Integer.MAX_VALUE;
    }

} // class GeneralSubstitutionModel
