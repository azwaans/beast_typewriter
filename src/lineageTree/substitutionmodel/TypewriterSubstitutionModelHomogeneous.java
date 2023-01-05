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

import static java.lang.Math.E;
import static java.lang.Math.log;


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
        //TODO nrOfStates +1 because of unedited
        nrOfStates = frequencies.getFreqs().length;
        insertFrequencies = frequenciesInput.get();


    } // initAndValidate

    protected boolean updateMatrix = true;
    private boolean storedUpdateMatrix = true;


    /**
     * This is to get transition probability between 2 sequences states (with potentially multiple edits having happened)
     * start sequence is sequence of parent node
     * end sequence is sequence of a child node
     */
    public double getSequenceTransitionProbability(final List<Integer> start_sequence, final List<Integer> end_sequence, double distance) {

        List<Integer> startstate = new ArrayList(start_sequence);
        List<Integer> endstate = new ArrayList(end_sequence);

        //create an unedited state to subtract from sequences to get only edited sites
        List<Integer> zero = Arrays.asList(0);

        //removing all unedited sites from each sequence
        startstate.removeAll(zero);
        endstate.removeAll(zero);

        //if endstate is less edited than the start state, violates ordering
        if(startstate.size() > endstate.size() ){
            return 0.0;
        }
        //subtracting start sequence from end sequence: edits introduced
        // if start state has identical elements to end state remove
        startstate.forEach(endstate::remove);
        //TODO foor loop over positions, make explicit checks, compare if that has the same result as above

        //available positions are 5 - number of edited positions
        // TODO rename nrOfPossibleInserts
        int poisson_up = 5 - startstate.size();

        //initialise the poisson distribution with mean rate * distance
        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(distance);

        //calculate the transition probability for the case where all available positions are edited in:
        // P(max) = 1- sum(P(n))
        //TODO rename endstate to newInserts
        if(endstate.size() == poisson_up ) {

            double sum = 0.0;
            for(int i = 0;  i<poisson_up; i++) {
                sum += dist.probability(i);
            }
            return  (1 - sum)*getFrequencyFactor(endstate);
        }

        //calculate the transition probability for the case where a #edits < available positions
        else{
            // Or keep in structure as in the if statement before
            // rename combine probability
            double freq = getFrequencyFactor(endstate);
            // todo rename poisson probability
            double prob = dist.probability(endstate.size());
            return prob * freq;
        }

    }

    //TODO function getNewInserts

    //TODO function calculate truncated poisson probability
// double sum = 0.0;
//            for(int i = 0;  i<poisson_up; i++) {
//                sum += dist.probability(i);
//            }
    /**
     * Function to obtain the probability factor induced by insert frequencies
     * combineInsertProbabilities
     */
    public double getFrequencyFactor(List<Integer> edits) {
        //TODO rename combinedProbability
        double factor = 1.0;
        double[] insertFrequenciesValue = insertFrequencies.getDoubleValues();

        for(Integer i : edits){
            //inserts are in {1, ..., nInserts}; insertProbabilities are in {0, ..., nInserts - 1}
            factor = factor * insertFrequenciesValue[i-1];
        }
        return factor;

    }

    /**
     * Function to obtain the insert probabilities
     * TODO : Check in initAndValidate: a) value in [0,1]; b) sum=1;
     */
    public double[] getInsertProbabilities() {

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
     * TODO: Ask Tim: How to deal with this class in the inheritance scheme and what to do with these
     * methods that we do not need.
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
