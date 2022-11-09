package lineageTree.substitutionmodel;


import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.util.BEASTClassLoader;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.DoubleStream;


@Description("Specifies transition probability vector for a Typewriter model")
public class TypewriterSubstitutionModel extends SubstitutionModel.Base {
    final public Input<RealParameter> ratesInput = new Input<>("rates",
            "Rates at which each target is cut in the barcode",
            (RealParameter) null);
    /**
     * an m_nStates vector current rates  *
     */
    protected Double[] rateVector;

    /**
     * Used for precalculations, sum of rates
     */
    private double sumOfRates;
    private boolean updateIntermediates = true;

    public void calculateIntermediates() {
        this.sumOfRates = arraySum(rateVector);
    }

    /**
     * The sum of the values in the array
     */
    public static double arraySum(Double[] xs) {
        double tmp = 0.0;
        for (double x : xs) {
            tmp += x;
        }
        return(tmp);
    }



    @Override
    public void initAndValidate() {
        super.initAndValidate();
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;
        rateVector = ratesInput.get().getValues();
        if (rateVector.length != nrOfStates ) {
            throw new IllegalArgumentException("Dimension of input 'rates' is " + rateVector.length + " but a " +
                    "rate vector of dimension " + nrOfStates + " was " +
                    "expected");
        }

    } // initAndValidate



    protected boolean updateMatrix = true;
    private boolean storedUpdateMatrix = true;

    //This is to get transition probabilities for a single editing event
    public double getEditTransitionProbability(int edit, double rate, double startTime, double endTime) {
        if (updateMatrix) {
            calculateIntermediates();
            updateMatrix = false;
        }
        double distance = (startTime - endTime) * rate;

        double pb = (rateVector[edit] - rateVector[edit] * Math.exp(-distance * sumOfRates)) / sumOfRates;
        return pb;
    }

    //This is to get transition probability between sequences (with potentially multiple edits)
    public double getSequenceTransitionProbability(List<Integer> start_sequence, List<Integer> end_sequence, double distance) {
        List<Integer> subtracted = new ArrayList<>(end_sequence);

        //substracting start sequence from end sequence
        start_sequence.forEach(subtracted::remove);
        double transition_prob = 1;

        //the sequence is unedited, only the probability of staying unedited
        if(subtracted.size() == 0 ) {
            return getTransitionProbability(0, distance);
        }

        //the subtraction is non empty, editing happened, we multiply edit probabilities
        //todo here would be where things would change for model variant 2
        for(Integer edit: subtracted) {
            transition_prob = transition_prob * getTransitionProbability(edit, distance);
        }

        return transition_prob;
    }

    //This is to get transition probability for a single position
    //todo remove redundancy with geteditprobability
    public double getTransitionProbability(Integer edit, double distance) {
//        if (updateMatrix) {
//            calculateIntermediates();
//            updateMatrix = false;
//        }
        //edit0 is staying in the unedited state, return the 1st element
        if(edit == 0) {
            return Math.exp(-distance *sumOfRates);
        }

        //there is an edit, return the corresponding transition probability
        else {
        double pb = (rateVector[edit-1] - rateVector[edit-1] * Math.exp(-distance * sumOfRates)) / sumOfRates;
        return pb;}
    }


    /**
     * access to (copy of) rate vector *
     */
    public Double[] getrateVector() {
        return rateVector.clone();
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


    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        double[] full_transition_probabilities = new double[rateVector.length];
        double branchtime = (endTime - startTime) * rate;
        for(int i=1 ; i<=rateVector.length; i++){
            full_transition_probabilities[i] = getTransitionProbability(i,branchtime);
        }
        matrix = full_transition_probabilities;

    }

    public double [] getInsertionProbs(){
        calculateIntermediates();

        double [] insertionProbs = new double [rateVector.length];

        for (int i=0; i<insertionProbs.length; i++){
            insertionProbs[i] = rateVector[i] / sumOfRates;
        }

        return insertionProbs;
    }


    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return null;
    }

    /**
     * This function returns the Eigen vectors.
     *
     * @return the array
     */


    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType.getStateCount() != Integer.MAX_VALUE;
    }

} // class GeneralSubstitutionModel
