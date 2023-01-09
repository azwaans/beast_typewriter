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

        double[] insertProbabilities = insertFrequencies.getDoubleValues();
        if( Arrays.stream(insertProbabilities).sum() != 1 ) {

            throw new IllegalArgumentException(String.format(
                    "sum of insert probabilities is not 1"));

        }

        for( double insertProbability : insertProbabilities) {
            if(insertProbability>1 || insertProbability<0) {
                throw new IllegalArgumentException(String.format(
                        "insert probabilities is not between 0 and 1"));
            }
        }



    }

    protected boolean updateMatrix = true;
    private boolean storedUpdateMatrix = true;


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
        //todo here, we check if we have any sequence transitions that violate ordering/do weird things
        //specifically, if there is already editing in startState, all edits in startState should be found in end state.
//        else if (startState.size() != 0){
//            //if there is already editing in startState, all edits in startState should be found in end state.
//            boolean sharedEdits = true;
//            for(int i=0; i<startState.size(); i++){
//                if(endState.get(i) != startState.get(i)) {
//                    sharedEdits = false;
//                }
//            }
//            if(sharedEdits == false) {
//                Log.info.println("there are no shared edits where there should be");
//                Log.info.println("startseq : " + startSequence);
//                Log.info.println("endseq : " + endSequence);
//            }

//        }
        //subtracting start sequence from end sequence: edits introduced
        // if start state has identical elements to end state remove
        startState.forEach(endState::remove);
        List<Integer> newInserts = endState;
        //TODO foor loop over positions, make explicit checks, compare if that has the same result as above



        //available positions are 5 - number of edited positions
        int nrOfPossibleInserts = 5 - startState.size();

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

    //TODO function getNewInserts

    /**
     * This function calculates the probability of reaching/editing the last unedited position in the barcode with nbrOfPossibleInserts
     * available positions
     *
     * @param nbrOfPossibleInserts  is the number of available positions until the absorbing state is reached
     * @param dist is the poisson
     */
    public double calculateAbsorbingStateProbability(org.apache.commons.math.distribution.PoissonDistribution dist,int nbrOfPossibleInserts) {
        double absorbingStateProbability = 1.0;
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
        double[] insertFrequenciesValue = insertFrequencies.getDoubleValues();

        for(Integer i : inserts){
            //inserts are in {1, ..., nInserts}; insertProbabilities are in {0, ..., nInserts - 1}
            factor = factor * insertFrequenciesValue[i-1];
        }
        return factor;

    }

    /**
     * Function to obtain the array of insert probabilities
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
//      System.arraycopy(relativeRates, 0, storedRelativeRates, 0, relativeRates.length);

        super.store();
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restore() {

        updateMatrix = storedUpdateMatrix;
//      To restore all this stuff just swap the pointers...
//      double[] tmp1 = storedRelativeRates;
//      storedRelativeRates = relativeRates;
//      relativeRates = tmp1;
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
