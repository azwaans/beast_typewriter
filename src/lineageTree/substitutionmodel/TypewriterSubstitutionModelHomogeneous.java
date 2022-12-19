package lineageTree.substitutionmodel;


import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
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
    final public Input<RealParameter> frequenciesInput = new Input<>("editfrequencies",
            "Edit frequencies for the typewriter process",
            (RealParameter) null);

    /**
     * edit insertion rate  *
     */
    protected RealParameter insertFrequencies;



    @Override
    public void initAndValidate() {
        super.initAndValidate();
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;
        insertFrequencies = frequenciesInput.get();


    } // initAndValidate

    protected boolean updateMatrix = true;
    private boolean storedUpdateMatrix = true;


    /**
     * This is to get transition probability between 2 sequences states (with potentially multiple edits having happened)
     */
    public double getSequenceTransitionProbability(final List<Integer> start_sequence, final List<Integer> end_sequence, double distance) {
        //Log.info.println("start sequence"+ start_sequence);
        //Log.info.println("end sequence"+ end_sequence);

        List<Integer> startstate = new ArrayList(start_sequence);
        List<Integer> endstate = new ArrayList(end_sequence);

        List<Integer> zero = Arrays.asList(0);

        //removing all unedited sites from each sequence
        startstate.removeAll(zero);
        endstate.removeAll(zero);

        //subtracting start sequence from end sequence: edits introduced
        startstate.forEach(endstate::remove);

        //if endstate is less edited than the start state, violates ordering
        if(startstate.size() > endstate.size() ){
            return 0.0;
        }

        //available positions are 5 - number of edited positions
        int poisson_up = 5 -  startstate.size();

        //initialise the poisson distribution with mean rate * distance
        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(distance);

        //calculate the transition probability for the case where all available positions are edited in:
        // P(max) = 1- sum(P(n))
        if(endstate.size() == poisson_up ) {
            double sum = 0.0;
            for(int i = 0;  i< poisson_up; i++) {
                sum += dist.probability(i);
            }
             return (1-sum) * getFrequencyFactor(endstate);
        }

        //calculate the transition probability for the case where a #edits < avaialable positions
        else{
            double freq = getFrequencyFactor(endstate);
            double prob = dist.probability(endstate.size());
            return prob * freq;
        }

    }

    /**
     * Function to obtain the probability factor induced by insert frequencies
     */
    public double getFrequencyFactor(List<Integer> edits) {
        double factor = 1.0;
        double[] insertFrequenciesValue = insertFrequencies.getDoubleValues();

        for(Integer i : edits){
            //Log.info.println("edit for which we are trying to find the freq" + i);
            factor = factor * insertFrequenciesValue[i-1];
        }
        return factor;

    }

    /**
     * Function to obtain the probability factor induced by insert frequencies
     */
    public double[] getInsertionProbs() {
        double[] insertFrequenciesValue = insertFrequencies.getDoubleValues();

        return insertFrequenciesValue;

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
