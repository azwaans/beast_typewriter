package lineageTree.substitutionmodel;


import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



@Description("Calculates transition probabilities sequence transitions from an ancestral node to a child node. " +
        "Assume that the sequences proposed represent valid transition pairs.")

public class TypewriterSubstitutionModel extends SubstitutionModel.Base {
    final public Input<RealParameter> editProbabilitiesInput = new Input<>("editProbabilities",
            "Edit probabilities for the typewriter process", Input.Validate.REQUIRED);

    /**
     * edit insertion rate  *
     */
    protected RealParameter editProbabilities;
    double[] editProbs;
     int targetBClength;


    @Override
    public void initAndValidate() {

        super.initAndValidate();
        //Here, nrOfStates represents the number of possible edits at each position of the TargetBCs + the unedited state
        nrOfStates = frequencies.getFreqs().length;

        // Get edit frequencies and check correct input
        editProbabilities = editProbabilitiesInput.get();
        editProbs = editProbabilities.getDoubleValues();

        if (nrOfStates != (editProbs.length + 1)){
            throw new RuntimeException(String.format("The edit probability vector has to be one element " +
                    "shorter than the vector of the possible state frequencies!"));
        }

        double insertProbabilitiesSum = Arrays.stream(editProbs).sum();
        if (Math.abs(insertProbabilitiesSum - 1.0) > 1e-6) {
            throw new IllegalArgumentException(String.format(
                    "sum of insert probabilities is not 1"));

        }

        for( double insertProbability : editProbs) {
            if(insertProbability>1 || insertProbability<0) {
                throw new IllegalArgumentException(String.format(
                        "insert probabilities is not between 0 and 1"));
            }
        }



    }


    /**
     * This function calculates the probability of transitioning between 2 sequences states in given evolutionary time (distance)
     * (with potentially multiple edits having happened)
     *
     * @param startSequence  is a sequence state at a parent node
     * @param endSequence is a sequence state at a child node
     */
    public double getSequenceTransitionProbability(final List<Integer> startSequence, final List<Integer> endSequence, double distance) {

        List<Integer> startState = new ArrayList(startSequence);
        List<Integer> endState = new ArrayList(endSequence);

        //create an unedited state to subtract from sequences to get only edited sites
        List<Integer> zero = Arrays.asList(0);

        //removing all unedited sites from each sequence
        startState.removeAll(zero);
        endState.removeAll(zero);

        //if endState is less edited than the start state, violates ordering
        if(startState.size() > endState.size() ){
            return 0.0;
        }

        //subtracting start sequence from end sequence: edits introduced
        // if start state has identical elements to end state remove
        startState.forEach(endState::remove);
        List<Integer> newInserts = endState;

        //available positions are targetBClength length - number of edited positions
        int nrOfPossibleInserts = targetBClength - startState.size();

        //initialise the poisson distribution with mean rate * distance
        org.apache.commons.math.distribution.PoissonDistribution poissonDistribution = new PoissonDistributionImpl(distance);

        //calculate the transition probability for the case where all available positions are edited in
        // This is the absorbing state in the poisson process
        // P(max) = 1- sum(P(n)) * probability of this insert combination
        if(newInserts.size() == nrOfPossibleInserts ) {

            return calculateAbsorbingStateProbability(poissonDistribution, nrOfPossibleInserts) * combinedInsertProbabilities(newInserts);
        }
        //calculate the transition probability for the case where a #edits < available positions
        //this is a regular draw from the poisson process * probability of this insert combination
        else if (newInserts.size() < nrOfPossibleInserts){

            return poissonDistribution.probability(newInserts.size()) * combinedInsertProbabilities(newInserts);

        } else{

            throw new RuntimeException("Error! Number of new inserts is larger than nr of possible inserts!");
        }
    }



    /**
     * This function calculates the probability of reaching/editing the last unedited position in the barcode with nbrOfPossibleInserts
     * available positions
     *
     * @param nbrOfPossibleInserts  is the number of available positions until the absorbing state is reached
     * @param dist is the poisson
     */
    public double calculateAbsorbingStateProbability(org.apache.commons.math.distribution.PoissonDistribution dist,int nbrOfPossibleInserts) {
        double absorbingStateProbability = 1.0;
        //TODO use cumulative distribution until nrOfPossible inserts
        for(int i = 0;  i < nbrOfPossibleInserts ; i++) {
            absorbingStateProbability -= dist.probability(i);
        }
        return absorbingStateProbability;

    }


    /**
     * Function to obtain the probability factor induced by insert frequencies
     * combineInsertProbabilities
     */
    public double combinedInsertProbabilities(List<Integer> inserts) {

        double factor = 1.0;
        double[] editProbabilitiesValue = editProbs;

        for(Integer i : inserts){
            //inserts are in {1, ..., nInserts}; insertProbabilities are in {0, ..., nInserts - 1}
            factor = factor * editProbabilitiesValue[i-1];
        }
        return factor;

    }

    /**
     * Function to obtain the array of insert probabilities
     *
     */
    public double[] getInsertProbabilities() {
        return editProbs;
    }

    /**
     * Function to set the array length for a particular alignment.
     *
     */
    public void setTargetBClength(int length) {
        targetBClength = length;
    }


    /**
     * CalculationNode implementation follows *
     */
    @Override
    public void store() {
        editProbs = editProbabilities.getDoubleValues();
        super.store();
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restore() {
        editProbs = editProbabilities.getDoubleValues();
        super.restore();

    }

    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        editProbs = editProbabilities.getDoubleValues();
        return true;
    }

    /**
     * TODO: Ask Tim: How to deal with this class in the inheritance scheme and what to do with these
     * methods that we do not need.
     * this not used in the current implementation because we do not use a matrix format
     */
    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return null;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType.getStateCount() != Integer.MAX_VALUE;
    }

}
